## =============================================================================
## 1-4 data clean use the diagnosis table and preoperative medications to build the four comorbidity flags
## =============================================================================

#Organize and clean the preoperative baseline diseases
diag1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/diagnosis.csv", 
                  locale = locale(encoding = "GBK")) # 217040obs
diag2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/diagnosis.csv", 
                  locale = locale(encoding = "UTF-8")) # 216651obs
diag3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/diagnosis.csv", 
                  locale = locale(encoding = "UTF-8")) # 981477obs
#Unify the formats of the three tables
standardize_diag <- function(df, hospital_name) {
  df %>%
    mutate(
      # Force all to trimmed strings --- resolves the conflict between dbl and chr and the inability to merge lgl NA
      across(any_of(c("PATIENT_ID", "INP_NO", "VISIT_ID", "OPER_ID")), 
             ~ trimws(as.character(.x))),
      # --- Diagnosis codes and names: convert all to character ---
      across(any_of(c("P_DIAGNOSIS_CODE", "P_DIAGNOSIS_NAME", "P_DIAGNOSIS_DESC", "P_DIAGNOSIS_TYPE")), 
             as.character),
      # --- Sequence-number columns: convert all to numeric ---
      across(any_of(c("P_DIAGNOSIS_NO", "diagnosticId")), 
             ~ as.numeric(as.character(.x))),
      # Add the hospital label
      Hospital = hospital_name
    ) %>% # Select only the needed core columns and drop hospital-specific extra columns (such as diagDateFormat)
    select(any_of(c("PATIENT_ID", "INP_NO", "VISIT_ID", "OPER_ID", 
                    "P_DIAGNOSIS_TYPE", "P_DIAGNOSIS_NAME", "P_DIAGNOSIS_CODE", 
                    "P_DIAGNOSIS_DESC", "Hospital")))}
# Perform the standardization
diag1 <- standardize_diag(diag1, "NanYi")
diag2 <- standardize_diag(diag2, "AnYi") #1884 patients were found to have multiple surgery records
diag3 <- standardize_diag(diag3, "ZhongDa")
# Merge (at this point diag1 holds the data of the three hospitals)
diag1 <- bind_rows(diag1, diag2, diag3) # 1415168 obs

## ============================================================================
## Step 2: keep only preoperative diagnoses such as "preoperative/admission/initial", then screen the four chronic disease keywords and rename them
## ============================================================================
#Filter data recorded before surgery
diag1 <- diag1 %>% 
  # Keep only preoperative/admission/initial diagnoses
  filter(P_DIAGNOSIS_TYPE %in% c("术前诊断", "入院诊断", "初步诊断","入院西医诊断","门诊诊断","确诊诊断","住院诊断","门诊诊断Outpatient","确诊诊断Confirmed","术前诊断Preoperative","入院诊断Admitting"))
keywords <- c("高血压", "糖尿病", "冠状动脉", "冠心病", 
              "心力衰竭", "心衰", "心功能", "心功能不全")
pattern <- paste(keywords, collapse = "|")  # Build the regular expression
diag1 <- diag1 %>% #Filter data containing the key diseases
  filter(grepl(pattern, P_DIAGNOSIS_NAME))  # 118550 obs
diag1 <- diag1 %>% #Data cleaning
  mutate(
    P_DIAGNOSIS_NAME = case_when(
      grepl("高血压", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Hypertension",
      grepl("糖尿病", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Diabetes",
      grepl("冠状动脉|冠心病", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Coronary artery disease",
      grepl("心力衰竭|心衰|心功能", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Heart failure",
      TRUE ~ P_DIAGNOSIS_NAME  )) 
## ============================================================================
## Step 3: Convert "longitudinal diagnoses" into four binary comorbidity flags (0/1) with one row per patient, join back to data1, and replace NA with 0
## ============================================================================
diag_flags <- diag1 %>% #Add one column per disease
  mutate(
    Hypertension = if_else(grepl("高血压|Hypertension", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0),
    Diabetes = if_else(grepl("糖尿病|Diabetes", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0),
    Coronary_artery_disease = if_else(grepl("冠状动脉|Coronary artery disease", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0),
    Heart_failure = if_else(grepl("心力衰竭|Heart failure", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0))
diag_flags_summary <- diag_flags %>% #Set the binary variables
  group_by(PATIENT_ID) %>%
  summarise(
    Hypertension = as.integer(max(Hypertension)),
    Diabetes = as.integer(max(Diabetes)),
    Coronary_artery_disease = as.integer(max(Coronary_artery_disease)),
    Heart_failure = as.integer(max(Heart_failure)),
    .groups = "drop" )
data1 <- data1 %>% #Join the result back into data1
  left_join(diag_flags_summary, by = c("PATIENT_ID"))
disease_cols <- c("Hypertension", "Diabetes", "Coronary_artery_disease", "Heart_failure")
data1 <- data1 %>% # Replace NA with 0
  mutate(across(all_of(disease_cols), ~replace_na(.x, 0)))
table(data1$Hypertension) # 18957
table(data1$Diabetes) # 6859
table(data1$Coronary_artery_disease) # 727
table(data1$Heart_failure) # 245

## ============================================================================
## Step 4: Use the category terms in the "preoperative medication" field to recover missed comorbidity patients, and add one more round (marked as 1)
## ============================================================================
table(data1$RECENTLY_TAKING_MEDICATION, useNA = "ifany") # Print the distribution of the field values
# Recovery keyword list: category terms + all common generic names (the stems "dihydropyridine/sartan/olol/pril/glitide/gliflozin/glutide" each cover a whole drug class)
keywords_hypertension <- "降压药|地平|氨氯地平|硝苯地平|非洛地平|尼群地平|贝尼地平|左旋氨氯地平|沙坦|缬沙坦|厄贝沙坦|氯沙坦|替米沙坦|坎地沙坦|奥美沙坦|阿利沙坦|普利|倍他乐克|洛尔|美托洛尔|比索洛尔|阿替洛尔|拉贝洛尔|卡维地洛|氢氯噻嗪|吲达帕胺|利血平|复方利血平|乌拉地尔|特拉唑嗪|多沙唑嗪"
keywords_diabetes     <- "糖尿病药|胰岛素|二甲双胍|格列|阿卡波糖|伏格列波糖|米格列醇|吡格列酮|罗格列酮|列净|恩格列净|达格列净|卡格列净|艾托格列净|利拉鲁肽|司美格鲁肽|度拉糖肽|艾塞那肽|鲁肽"
keywords_cad          <- "硝酸甘油|硝酸异山梨酯|单硝酸异山梨酯|抗血小板药|阿司匹林|波立维|氯吡格雷|替格瑞洛|曲美他嗪|尼可地尔|复方丹参滴丸|速效救心丸|麝香保心丸|地奥心血康"
keywords_hf           <- "利尿剂|呋塞米|速尿|托拉塞米|布美他尼|地高辛|螺内酯|安体舒通|依普利酮|诺欣妥|沙库巴曲|伊伐布雷定|芪苈强心"

# 1. Run the recovery logic (argument name: ignore_case)
data1 <- data1 %>%
  mutate(
    # Hypertension
    Hypertension = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_hypertension, ignore_case = TRUE)) %in% TRUE,
      1, Hypertension ),
    # Diabetes
    Diabetes = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_diabetes, ignore_case = TRUE)) %in% TRUE,
      1, Diabetes ),
    # Coronary artery disease
    Coronary_artery_disease = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_cad, ignore_case = TRUE)) %in% TRUE,
      1, Coronary_artery_disease ),
    # Heart failure
    Heart_failure = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_hf, ignore_case = TRUE)) %in% TRUE,
      1, Heart_failure ) )
# Check the recovery effect: see whether the proportion marked as 1 has increased
table(data1$Hypertension) # Diagnosis 18957 -> after recovery 28922 (+9965, all from "antihypertensive drugs")
table(data1$Diabetes) # 6859 -> 10047 (+3188, all from "antidiabetic drugs")
table(data1$Coronary_artery_disease) # 727 -> 727 (no matching category in the field, recovered as 0)
table(data1$Heart_failure) # 245 -> 245 (no matching category in the field, recovered as 0)
