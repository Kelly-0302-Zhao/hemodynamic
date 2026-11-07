## =============================================================================
## 1-1 data clean - operation type.R —— 三家医院基础表：统一格式、填补科室、合并、筛单药、打手术类型
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
setwd("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827/JCA_run") #设置存储位置
#===============================================================================
#数据清洗，基线描述
#===============================================================================
data1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/master.csv", 
                  locale = locale(encoding = "GBK")) # 24649obs 南京市第一医院患者基础数据
data2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/master.csv", 
                  locale = locale(encoding = "GBK")) # 50572obs 安徽附一患者基础数据
data3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/master.csv", 
                  locale = locale(encoding = "UTF-8")) # 37452obs 中大患者基础数据
#=====首先用dept stayed 填补operating dept 由于部分手术名称无法分辨手术分类，因此用科室填补
#定义标准化函数：设置一个函数，规定每个变量的格式，用此函数去统一三个医院的数据格式
standardize_df <- function(df) {
  # 显式定义上海时区，避免 8 小时偏差
  my_tz <- "Asia/Shanghai"
  df %>%
    mutate(
      # --- ID 变量：增加去空格处理 ---
      # 有些中心数据 ID 后面可能带不可见空格，导致 merge 失败
      across(c(PATIENT_ID, VISIT_ID, OPER_ID, INP_NO), ~ trimws(as.character(.x))),
      # --- 科室与分类列 ---
      OPERATING_DEPT = as.character(OPERATING_DEPT),
      DEPT_STAYED    = as.character(DEPT_STAYED),
      ASA_GRADE      = as.character(ASA_GRADE),
      EMERGENCY_INDICATOR = as.numeric(EMERGENCY_INDICATOR),
      # --- 时间列处理 ---
      across(
        c(IN_HOSPITAL_TIME, IN_DATE_TIME, ANES_START_TIME, START_DATE_TIME, 
          END_DATE_TIME, ANES_END_TIME, OUT_DATE_TIME, IN_ICU_TIME, 
          OUT_ICU_TIME, OUT_HOSPITAL_TIME),
        ~ {
          # 先转字符，再尝试多种顺序解析
          tmp_time <- parse_date_time(as.character(.x), 
                                      orders = c("ymd HMS", "ymd HM", "mdy HMS", "dmy HMS"), 
                                      tz = my_tz)
          # 如果原本不是 NA 但转完变 NA 了，给出警告
          return(tmp_time) } ),
      # --- 人口学与药物：防止非数值干扰 ---
      across(c(AGE, HEIGHT, WEIGHT, BMI, PROPOFOL, CIPROFOL), as.numeric),
      # --- 字符列：确保 NA 统一 ---
      across(c(DIAG_BEFORE_OPERATION, OPERATION_NAME, ANESTHESIA_METHOD, 
               OPERATION_HISTORY, OPERATION_NAME_OF_HISTORY, 
               RECENTLY_TAKING_MEDICATION, PREOPERATIVE_BLOOD_PRESSURE), as.character) )}
# 执行转换
## 三家医院各自跑一遍 standardize_df，统一格式
data1 <- standardize_df(data1)
data2 <- standardize_df(data2)
data3 <- standardize_df(data3)
attr(data1$ANES_START_TIME, "tzone")
attr(data2$ANES_START_TIME, "tzone")
attr(data3$ANES_START_TIME, "tzone")

colSums(is.na(data1))#patient 8缺失
colSums(is.na(data2))#patient 48缺失
colSums(is.na(data3))

#增加医院标签并剔除关键 key 缺失的行（patientid缺失会导致后续匹配出现error）
data1 <- data1 %>% filter(!is.na(PATIENT_ID)) %>% mutate(Hospital = "NanYi") # 24641obs
data2 <- data2 %>% filter(!is.na(PATIENT_ID)) %>% mutate(Hospital = "AnYi") # 50524obs
data3 <- data3 %>% filter(!is.na(PATIENT_ID)) %>% mutate(Hospital = "ZhongDa") # 37452obs

#发现dept里有数字，先将数字变成NA
## 这个函数：把"纯数字"的科室名当成脏数据，清成 NA
clean_dept <- function(df) {
  df %>%
    mutate(
      OPERATING_DEPT = as.character(OPERATING_DEPT),
      
      # 纯数字 → NA
      OPERATING_DEPT = ifelse(
        grepl("^\\d+$", trimws(OPERATING_DEPT)),
        NA,
        OPERATING_DEPT
      )
    )
}

## 三家医院各自跑一遍 clean_dept（把"纯数字"科室清成 NA）
data1 <- clean_dept(data1)
data2 <- clean_dept(data2)
data3 <- clean_dept(data3)

## ============================================================================
## 第 1 步：用"病人所在科室 DEPT_STAYED"填补缺失的手术科室（OPERATING_DEPT）
## ============================================================================
fill_dept <- function(df) {
  df %>%
    mutate(
      OPERATING_DEPT = coalesce(OPERATING_DEPT, DEPT_STAYED),
      OPERATING_DEPT = as.character(OPERATING_DEPT),
      OPERATION_NAME = as.character(OPERATION_NAME)
    )
}

## 三家医院各自跑一遍 fill_dept（用 DEPT_STAYED 填补手术科室）
data1 <- fill_dept(data1)
data2 <- fill_dept(data2)
data3 <- fill_dept(data3)


## ============================================================================
## 第 2 步：构建"手术名 → 科室"映射字典（取众数，而不是 distinct）
## ============================================================================
mapping_dict <- data1 %>%
  filter(!is.na(OPERATION_NAME) & !is.na(OPERATING_DEPT)) %>%
  count(OPERATION_NAME, OPERATING_DEPT) %>%
  group_by(OPERATION_NAME) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(OPERATING_DEPT = as.character(OPERATING_DEPT))


## ============================================================================
## 第 3 步：用"手术名 → 科室"字典做精确匹配，再补一批缺失科室
## ============================================================================
safe_fill_dept <- function(target_df, mapping_dict) {
  
  target_df %>%
    left_join(mapping_dict, by = "OPERATION_NAME", suffix = c("", "_mapped")) %>%
    mutate(
      OPERATING_DEPT = coalesce(OPERATING_DEPT, OPERATING_DEPT_mapped)
    ) %>%
    select(-OPERATING_DEPT_mapped)
}

## 用字典精确匹配，给 AnYi / ZhongDa 再补一次科室
data2 <- safe_fill_dept(data2, mapping_dict)
data3 <- safe_fill_dept(data3, mapping_dict)

## ============================================================================
## 第 4 步：字典精确匹配后仍缺失的科室，用模糊匹配（字符相似度，最大距离 1）再补一次
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
      maxDist = 1   # 收紧匹配条件
    )
    
    # 只填充成功匹配的
    to_fix$OPERATING_DEPT <- ifelse(
      is.na(idx),
      NA,
      mapping_dict$OPERATING_DEPT[idx]
    )
  }
  
  bind_rows(matched, to_fix) %>%
    select(-clean_name)
}

## 仍缺科室的，再用模糊匹配补一次（只填匹配成功的，避免误填）
data2 <- fuzzy_fill_dept(data2, mapping_dict)
data3 <- fuzzy_fill_dept(data3, mapping_dict)

## ============================================================================
## 第 5 步：三家医院合并成一张大表，并统一科室/手术名两列的类型（112617 行）
## ============================================================================
data1 <- bind_rows(data1, data2, data3) %>%
  mutate(
    OPERATING_DEPT = as.character(OPERATING_DEPT),
    OPERATION_NAME = as.character(OPERATION_NAME)
  )


## ============================================================================
## 第 6 步：只保留"单用环泊酚 / 单用丙泊酚"的患者，并打上 Group 标签（105574 行）
## ============================================================================
data1 <- data1 %>%
  filter((CIPROFOL == 1 & PROPOFOL == 0) |
           (CIPROFOL == 0 & PROPOFOL == 1)) %>%
  mutate(Group = ifelse(CIPROFOL == 1, "Ciprofol_only", "Propofol_only"))


## ============================================================================
## 第 7 步：按"科室名 + 手术名"给每台手术打上"手术类型"标签（骨科/普外/泌尿…）
## ============================================================================
data1 <- data1 %>%
  mutate(operation_type = case_when(
    
    # 骨科
    str_detect(OPERATING_DEPT, "骨科|关节|脊柱|创伤|手外|运动医学") 
    ~ "Orthopaedics",
    
    # 普外
    str_detect(OPERATING_DEPT, "普外|肝胆|胃肠|肛肠|乳腺|甲状腺|烧伤|整形|肿瘤外科") 
    ~ "General Surgery",
    
    # 泌尿
    str_detect(OPERATING_DEPT, "泌尿") 
    ~ "Urology",
    
    # 心胸
    str_detect(OPERATING_DEPT, "胸外|肺外科|心血管外科") 
    ~ "Thoracic Surgery",
    
    # 妇产
    str_detect(OPERATING_DEPT, "妇科|产科|爱婴") 
    ~ "Gynecology",
    
    # 神外
    str_detect(OPERATING_DEPT, "神经外科") 
    ~ "Neurosurgery",
    
    # 二次兜底（用手术名）
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
## 第 8 步：把手术类型设成"有固定顺序"的因子
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
