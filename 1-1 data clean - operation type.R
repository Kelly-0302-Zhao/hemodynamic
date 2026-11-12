## =============================================================================
## 1-1 data clean - operation type.R -- baseline tables for the three hospitals: unify formats, fill in departments, merge, keep single-drug patients, assign surgery type
## =============================================================================

# Packages
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(flextable)
library(officer)
library(tidyverse)
library(tableone)
library(kableExtra)
library(data.table)
library(readr)
library(ggplot2)
library(stringi)
library(forestplot)
library(stringdist)
library(scales)
library(MuMIn)
library(geepack)
library(lubridate)
library(cobalt)
library(MatchIt)
library(broom)
library(grid)
library(purrr)
library(survey)
library(Hmisc)
library(openxlsx)
setwd("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827/JCA_run") #Set the storage location
#===============================================================================
#Data cleaning, baseline description
#===============================================================================
data1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/master.csv", 
                  locale = locale(encoding = "GBK")) # 24649obs baseline patient data from Nanjing First Hospital
data2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/master.csv", 
                  locale = locale(encoding = "GBK")) # 50572obs baseline patient data from AnYi hospital
data3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/master.csv", 
                  locale = locale(encoding = "UTF-8")) # 37452obs baseline patient data from ZhongDa hospital
#=====First fill operating dept with dept stayed; some surgery names cannot identify the surgery category, so departments are used to fill them
#Define the standardization function: set up a function that specifies the format of each variable and use it to unify the data formats of the three hospitals
standardize_df <- function(df) {
  # Explicitly set the Shanghai time zone to avoid an 8-hour offset
  my_tz <- "Asia/Shanghai"
  df %>%
    mutate(
      # --- ID variables: add trimming ---
      # Some centers have invisible spaces after the ID, which can cause merge failures
      across(c(PATIENT_ID, VISIT_ID, OPER_ID, INP_NO), ~ trimws(as.character(.x))),
      # --- Department and classification columns ---
      OPERATING_DEPT = as.character(OPERATING_DEPT),
      DEPT_STAYED    = as.character(DEPT_STAYED),
      ASA_GRADE      = as.character(ASA_GRADE),
      EMERGENCY_INDICATOR = as.numeric(EMERGENCY_INDICATOR),
      # --- Time column handling ---
      across(
        c(IN_HOSPITAL_TIME, IN_DATE_TIME, ANES_START_TIME, START_DATE_TIME, 
          END_DATE_TIME, ANES_END_TIME, OUT_DATE_TIME, IN_ICU_TIME, 
          OUT_ICU_TIME, OUT_HOSPITAL_TIME),
        ~ {
          # Convert to character first, then try several parsing orders
          tmp_time <- parse_date_time(as.character(.x), 
                                      orders = c("ymd HMS", "ymd HM", "mdy HMS", "dmy HMS"), 
                                      tz = my_tz)
          # Give a warning if a value that was not NA becomes NA after conversion
          return(tmp_time) } ),
      # --- Demographics and drugs: guard against non-numeric values ---
      across(c(AGE, HEIGHT, WEIGHT, BMI, PROPOFOL, CIPROFOL), as.numeric),
      # --- Character columns: standardize NA ---
      across(c(DIAG_BEFORE_OPERATION, OPERATION_NAME, ANESTHESIA_METHOD, 
               OPERATION_HISTORY, OPERATION_NAME_OF_HISTORY, 
               RECENTLY_TAKING_MEDICATION, PREOPERATIVE_BLOOD_PRESSURE), as.character) )}
# Perform the standardization
## Run standardize_df once for each of the three hospitals to unify formats
data1 <- standardize_df(data1)
data2 <- standardize_df(data2)
data3 <- standardize_df(data3)
attr(data1$ANES_START_TIME, "tzone")
attr(data2$ANES_START_TIME, "tzone")
attr(data3$ANES_START_TIME, "tzone")

colSums(is.na(data1))#patient 8 missing
colSums(is.na(data2))#patient 48 missing
colSums(is.na(data3))

#Add the hospital label and drop rows with missing key fields (a missing patientid causes errors in later matching)
data1 <- data1 %>% filter(!is.na(PATIENT_ID)) %>% mutate(Hospital = "NanYi") # 24641obs
data2 <- data2 %>% filter(!is.na(PATIENT_ID)) %>% mutate(Hospital = "AnYi") # 50524obs
data3 <- data3 %>% filter(!is.na(PATIENT_ID)) %>% mutate(Hospital = "ZhongDa") # 37452obs

#Numbers are found in dept; convert them to NA first
## This function treats department names that are purely numeric as dirty data and sets them to NA
clean_dept <- function(df) {
  df %>%
    mutate(
      OPERATING_DEPT = as.character(OPERATING_DEPT),
      
      # Purely numeric -> NA
      OPERATING_DEPT = ifelse(
        grepl("^\\d+$", trimws(OPERATING_DEPT)),
        NA,
        OPERATING_DEPT
      )
    )
}

## Run clean_dept once for each of the three hospitals (sets purely numeric departments to NA)
data1 <- clean_dept(data1)
data2 <- clean_dept(data2)
data3 <- clean_dept(data3)

## ============================================================================
## Step 1: fill missing operating department (OPERATING_DEPT) with the patient's department (DEPT_STAYED)
## ============================================================================
fill_dept <- function(df) {
  df %>%
    mutate(
      OPERATING_DEPT = coalesce(OPERATING_DEPT, DEPT_STAYED),
      OPERATING_DEPT = as.character(OPERATING_DEPT),
      OPERATION_NAME = as.character(OPERATION_NAME)
    )
}

## Run fill_dept once for each of the three hospitals (fills operating department with DEPT_STAYED)
data1 <- fill_dept(data1)
data2 <- fill_dept(data2)
data3 <- fill_dept(data3)


## ============================================================================
## Step 2: build a "surgery name -> department" mapping dictionary (using the mode, not distinct values)
## ============================================================================
mapping_dict <- data1 %>%
  filter(!is.na(OPERATION_NAME) & !is.na(OPERATING_DEPT)) %>%
  count(OPERATION_NAME, OPERATING_DEPT) %>%
  group_by(OPERATION_NAME) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(OPERATING_DEPT = as.character(OPERATING_DEPT))


## ============================================================================
## Step 3: apply exact matching with the "surgery name -> department" dictionary, then fill another batch of missing departments
## ============================================================================
safe_fill_dept <- function(target_df, mapping_dict) {
  
  target_df %>%
    left_join(mapping_dict, by = "OPERATION_NAME", suffix = c("", "_mapped")) %>%
    mutate(
      OPERATING_DEPT = coalesce(OPERATING_DEPT, OPERATING_DEPT_mapped)
    ) %>%
    select(-OPERATING_DEPT_mapped)
}

## Apply exact dictionary matching to fill departments again for AnYi / ZhongDa
data2 <- safe_fill_dept(data2, mapping_dict)
data3 <- safe_fill_dept(data3, mapping_dict)

## ============================================================================
## Step 4: for departments still missing after exact dictionary matching, fill them again with fuzzy matching (character similarity, maximum distance 1)
## ============================================================================
fuzzy_fill_dept <- function(target_df, mapping_dict) {
  
  mapping_dict <- mapping_dict %>%
    mutate(clean_name = str_replace_all(OPERATION_NAME, "[[:punct:]\\s]", ""))
  
  target_df <- target_df %>%
    mutate(
      clean_name = str_replace_all(OPERATION_NAME, "[[:punct:]\\s]", ""),
      OPERATING_DEPT = as.character(OPERATING_DEPT)
    )
  
  matched <- target_df %>% filter(!is.na(OPERATING_DEPT))
  to_fix  <- target_df %>% filter(is.na(OPERATING_DEPT))
  
  if(nrow(to_fix) > 0) {
    
    idx <- amatch(
      to_fix$clean_name,
      mapping_dict$clean_name,
      maxDist = 1   # Tighten the matching condition
    )
    
    # Fill only the successfully matched values
    to_fix$OPERATING_DEPT <- ifelse(
      is.na(idx),
      NA,
      mapping_dict$OPERATING_DEPT[idx]
    )
  }
  
  bind_rows(matched, to_fix) %>%
    select(-clean_name)
}

## For departments still missing, fill again with fuzzy matching (only successful matches are filled, to avoid wrong fills)
data2 <- fuzzy_fill_dept(data2, mapping_dict)
data3 <- fuzzy_fill_dept(data3, mapping_dict)

## ============================================================================
## Step 5: merge the three hospitals into one table and unify the types of the department and surgery name columns (112617 rows)
## ============================================================================
data1 <- bind_rows(data1, data2, data3) %>%
  mutate(
    OPERATING_DEPT = as.character(OPERATING_DEPT),
    OPERATION_NAME = as.character(OPERATION_NAME)
  )


## ============================================================================
## Step 6: keep only patients who received ciprofol alone or propofol alone, and assign the Group label (105574 rows)
## ============================================================================
data1 <- data1 %>%
  filter((CIPROFOL == 1 & PROPOFOL == 0) |
           (CIPROFOL == 0 & PROPOFOL == 1)) %>%
  mutate(Group = ifelse(CIPROFOL == 1, "Ciprofol_only", "Propofol_only"))


## ============================================================================
## Step 7: assign each surgery a surgery type label by "department name + surgery name" (orthopedics/general surgery/urology...)
## ============================================================================
data1 <- data1 %>%
  mutate(operation_type = case_when(
    
    # Orthopedics
    str_detect(OPERATING_DEPT, "骨科|关节|脊柱|创伤|手外|运动医学") 
    ~ "Orthopaedics",
    
    # General surgery
    str_detect(OPERATING_DEPT, "普外|肝胆|胃肠|肛肠|乳腺|甲状腺|烧伤|整形|肿瘤外科") 
    ~ "General Surgery",
    
    # Urology
    str_detect(OPERATING_DEPT, "泌尿") 
    ~ "Urology",
    
    # Cardiothoracic surgery
    str_detect(OPERATING_DEPT, "胸外|肺外科|心血管外科") 
    ~ "Thoracic Surgery",
    
    # Obstetrics and gynecology
    str_detect(OPERATING_DEPT, "妇科|产科|爱婴") 
    ~ "Gynecology",
    
    # Neurosurgery
    str_detect(OPERATING_DEPT, "神经外科") 
    ~ "Neurosurgery",
    
    # Second fallback (using surgery name)
    TRUE ~ case_when(
      str_detect(OPERATION_NAME, "骨折|关节置换|脊柱|椎间盘|韧带") ~ "Orthopaedics",
      str_detect(OPERATION_NAME, "胆囊|胃癌|结肠|直肠|乳腺|甲状腺") ~ "General Surgery",
      str_detect(OPERATION_NAME, "前列腺|肾脏|输尿管|膀胱") ~ "Urology",
      str_detect(OPERATION_NAME, "肺|胸腔|食管") ~ "Thoracic Surgery",
      str_detect(OPERATION_NAME, "子宫|卵巢|剖宫产|妇科") ~ "Gynecology",
      str_detect(OPERATION_NAME, "颅内|脑肿瘤|脑出血|神经外科") ~ "Neurosurgery",
      TRUE ~ "Others"
    )
  ))


## ============================================================================
## Step 8: set surgery type as a factor with a fixed level order
## ============================================================================
op_levels <- c(
  "General Surgery", "Gynecology",
  "Neurosurgery", "Orthopaedics", "Thoracic Surgery",
  "Urology", "Others"
)

data1$operation_type <- factor(data1$operation_type, levels = op_levels)

table(data1$operation_type)
# General Surgery：28317；Gynecology：3456；Neurosurgery：3336；Orthopaedics：24230
# Thoracic Surgery：9482；Urology：21913；Others：14840
