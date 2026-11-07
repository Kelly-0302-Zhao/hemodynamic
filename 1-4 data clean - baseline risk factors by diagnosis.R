## =============================================================================
## 1-4 data clean 用诊断表 + 术前用药，造四大合并症标志
## =============================================================================

#术前基础疾病整理清洗
diag1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/diagnosis.csv", 
                  locale = locale(encoding = "GBK")) # 217040obs
diag2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/diagnosis.csv", 
                  locale = locale(encoding = "UTF-8")) # 216651obs
diag3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/diagnosis.csv", 
                  locale = locale(encoding = "UTF-8")) # 981477obs
#统一三个表格的格式
standardize_diag <- function(df, hospital_name) {
  df %>%
    mutate(
      # 全部强制转为去空格的字符串 --- 解决 dbl 与 chr 冲突，以及 lgl 型 NA 无法合并的问题
      across(any_of(c("PATIENT_ID", "INP_NO", "VISIT_ID", "OPER_ID")), 
             ~ trimws(as.character(.x))),
      # --- 诊断编码与名称：统一转字符 ---
      across(any_of(c("P_DIAGNOSIS_CODE", "P_DIAGNOSIS_NAME", "P_DIAGNOSIS_DESC", "P_DIAGNOSIS_TYPE")), 
             as.character),
      # --- 序号类：统一转数值 ---
      across(any_of(c("P_DIAGNOSIS_NO", "diagnosticId")), 
             ~ as.numeric(as.character(.x))),
      # 增加医院标签
      Hospital = hospital_name
    ) %>% # 只选择需要的核心列，丢弃各医院特有的杂乱列（如 diagDateFormat）
    select(any_of(c("PATIENT_ID", "INP_NO", "VISIT_ID", "OPER_ID", 
                    "P_DIAGNOSIS_TYPE", "P_DIAGNOSIS_NAME", "P_DIAGNOSIS_CODE", 
                    "P_DIAGNOSIS_DESC", "Hospital")))}
# 执行转换
diag1 <- standardize_diag(diag1, "NanYi")
diag2 <- standardize_diag(diag2, "AnYi") #发现有1884名患者是出现多次手术的数据
diag3 <- standardize_diag(diag3, "ZhongDa")
# 合并（此时diag1为三家医院的数据）
diag1 <- bind_rows(diag1, diag2, diag3) # 1415168 obs

## ============================================================================
## 第 2 步：只保留"术前/入院/初步"等术前诊断，再筛出四大慢病关键词并重命名
## ============================================================================
#筛选手术之前的数据
diag1 <- diag1 %>% 
  # 只保留术前/入院/初步诊断
  filter(P_DIAGNOSIS_TYPE %in% c("术前诊断", "入院诊断", "初步诊断","入院西医诊断","门诊诊断","确诊诊断","住院诊断","门诊诊断Outpatient","确诊诊断Confirmed","术前诊断Preoperative","入院诊断Admitting"))
keywords <- c("高血压", "糖尿病", "冠状动脉", "冠心病", 
              "心力衰竭", "心衰", "心功能", "心功能不全")
pattern <- paste(keywords, collapse = "|")  # 构建正则表达式
diag1 <- diag1 %>% #筛选含有关键疾病的数据
  filter(grepl(pattern, P_DIAGNOSIS_NAME))  # 118550 obs
diag1 <- diag1 %>% #数据清洗
  mutate(
    P_DIAGNOSIS_NAME = case_when(
      grepl("高血压", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Hypertension",
      grepl("糖尿病", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Diabetes",
      grepl("冠状动脉|冠心病", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Coronary artery disease",
      grepl("心力衰竭|心衰|心功能", P_DIAGNOSIS_NAME, ignore.case = TRUE) ~ "Heart failure",
      TRUE ~ P_DIAGNOSIS_NAME  )) 
## ============================================================================
## 第 3 步：把"纵向诊断"转成"每个病人一行"的四大合并症标志(0/1)，join 回 data1，NA 补 0
## ============================================================================
diag_flags <- diag1 %>% #给每个疾病放一列
  mutate(
    Hypertension = if_else(grepl("高血压|Hypertension", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0),
    Diabetes = if_else(grepl("糖尿病|Diabetes", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0),
    Coronary_artery_disease = if_else(grepl("冠状动脉|Coronary artery disease", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0),
    Heart_failure = if_else(grepl("心力衰竭|Heart failure", P_DIAGNOSIS_NAME, ignore.case = TRUE), 1, 0))
diag_flags_summary <- diag_flags %>% #设置二分类变量
  group_by(PATIENT_ID) %>%
  summarise(
    Hypertension = as.integer(max(Hypertension)),
    Diabetes = as.integer(max(Diabetes)),
    Coronary_artery_disease = as.integer(max(Coronary_artery_disease)),
    Heart_failure = as.integer(max(Heart_failure)),
    .groups = "drop" )
data1 <- data1 %>% #将结果join到data1里
  left_join(diag_flags_summary, by = c("PATIENT_ID"))
disease_cols <- c("Hypertension", "Diabetes", "Coronary_artery_disease", "Heart_failure")
data1 <- data1 %>% # 将 NA 替换为 0
  mutate(across(all_of(disease_cols), ~replace_na(.x, 0)))
table(data1$Hypertension) # 18957
table(data1$Diabetes) # 6859
table(data1$Coronary_artery_disease) # 727
table(data1$Heart_failure) # 245

## ============================================================================
## 第 4 步：用"术前用药"字段里的类目词，打捞漏标的合并症患者，再补一轮（标成 1）
## ============================================================================
table(data1$RECENTLY_TAKING_MEDICATION, useNA = "ifany") # 打印字段取值分布
# 打捞词库：类目词 + 全部常见通用名（词根"地平/沙坦/洛尔/普利/格列/列净/鲁肽"覆盖一整类药）
keywords_hypertension <- "降压药|地平|氨氯地平|硝苯地平|非洛地平|尼群地平|贝尼地平|左旋氨氯地平|沙坦|缬沙坦|厄贝沙坦|氯沙坦|替米沙坦|坎地沙坦|奥美沙坦|阿利沙坦|普利|倍他乐克|洛尔|美托洛尔|比索洛尔|阿替洛尔|拉贝洛尔|卡维地洛|氢氯噻嗪|吲达帕胺|利血平|复方利血平|乌拉地尔|特拉唑嗪|多沙唑嗪"
keywords_diabetes     <- "糖尿病药|胰岛素|二甲双胍|格列|阿卡波糖|伏格列波糖|米格列醇|吡格列酮|罗格列酮|列净|恩格列净|达格列净|卡格列净|艾托格列净|利拉鲁肽|司美格鲁肽|度拉糖肽|艾塞那肽|鲁肽"
keywords_cad          <- "硝酸甘油|硝酸异山梨酯|单硝酸异山梨酯|抗血小板药|阿司匹林|波立维|氯吡格雷|替格瑞洛|曲美他嗪|尼可地尔|复方丹参滴丸|速效救心丸|麝香保心丸|地奥心血康"
keywords_hf           <- "利尿剂|呋塞米|速尿|托拉塞米|布美他尼|地高辛|螺内酯|安体舒通|依普利酮|诺欣妥|沙库巴曲|伊伐布雷定|芪苈强心"

# 1. 执行打捞逻辑 (参数名为 ignore_case)
data1 <- data1 %>%
  mutate(
    # 高血压
    Hypertension = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_hypertension, ignore_case = TRUE)) %in% TRUE,
      1, Hypertension ),
    # 糖尿病
    Diabetes = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_diabetes, ignore_case = TRUE)) %in% TRUE,
      1, Diabetes ),
    # 冠心病
    Coronary_artery_disease = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_cad, ignore_case = TRUE)) %in% TRUE,
      1, Coronary_artery_disease ),
    # 心力衰竭
    Heart_failure = if_else(
      str_detect(RECENTLY_TAKING_MEDICATION, regex(keywords_hf, ignore_case = TRUE)) %in% TRUE,
      1, Heart_failure ) )
# 检查打捞效果：看看标记为 1 的比例是否有增加
table(data1$Hypertension) # 诊断 18957 -> 打捞后 28922（+9965，全靠"降压药"）
table(data1$Diabetes) # 6859 -> 10047（+3188，全靠"糖尿病药"）
table(data1$Coronary_artery_disease) # 727 -> 727（字段里没有相关类目，打捞为 0）
table(data1$Heart_failure) # 245 -> 245（字段里没有相关类目，打捞为 0）
