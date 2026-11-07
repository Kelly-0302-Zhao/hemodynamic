## =============================================================================
## 1-2 data clean - ASA&GENERAL&DURATION.R —— ASA 统一分级、只保留全麻、补齐麻醉起止时间
## =============================================================================

#规范ASA分级
## ============================================================================
## 第 1 步：把各家医院五花八门的 ASA 写法统一成 Ⅰ/Ⅱ/Ⅲ（识别不出记 NA）
## ============================================================================
data1 <- data1 %>%
  mutate(ASA_GRADE_raw = ASA_GRADE)
data1 <- data1 %>%
  mutate(
    ASA_GRADE = case_when(
      
      str_detect(
        ASA_GRADE_raw,
        regex(
          "^\\s*(?:ASA\\s*)?(?:3|Ⅲ|III)\\s*(?:\\+?\\s*E)?\\s*(?:级|类|A)?\\s*$",
          ignore_case = TRUE
        )
      ) ~ "Ⅲ",
      
      str_detect(
        ASA_GRADE_raw,
        regex(
          "^\\s*(?:ASA\\s*)?(?:2|Ⅱ|II)\\s*(?:\\+?\\s*E)?\\s*(?:级|类|A)?\\s*$",
          ignore_case = TRUE
        )
      ) ~ "Ⅱ",
      
      str_detect(
        ASA_GRADE_raw,
        regex(
          "^\\s*(?:ASA\\s*)?(?:1|Ⅰ|I)\\s*(?:\\+?\\s*E)?\\s*(?:级|类|A)?\\s*$",
          ignore_case = TRUE
        )
      ) ~ "Ⅰ",
      
      TRUE ~ NA_character_
    )
  )
table(data1$ASA_GRADE, useNA = "ifany") 
# Ⅰ     Ⅱ     Ⅲ 
# 1959 69735 33880 
## ============================================================================
## 第 2 步：只保留"全身麻醉"的手术（剔除非全麻 / 不能确认全麻者， 2682 人）
## ============================================================================
general <- read_csv(
  "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/General Anesthesia.csv",
  locale = locale(encoding = "GBK")) %>% 
  distinct() #导入临床协助整理的全麻列表
# 1. 预处理字典表（确保文字不含空格，以免影响匹配）
general <- general %>%
  mutate(ANESTHESIA_METHOD = trimws(as.character(ANESTHESIA_METHOD))) %>%
  distinct()
# 2. 对合并后的 data1 进行筛选（筛选了2682名非全麻或不能确认全麻的患者，剩余102892个患者）
data1 <- data1 %>%
  # 关联全麻字典
  left_join(general, by = "ANESTHESIA_METHOD") %>%
  filter(
    GENERAL == 1 | 
      grepl("全麻|全身麻醉|静吸复合|插管|喉罩", ANESTHESIA_METHOD)) %>%
  # 移除辅助列
  select(-GENERAL)

## ============================================================================
## 第 3 步：补齐麻醉开始/结束时间（缺失时用入室/出室时间兜底），为算时长做准备
## ============================================================================
#计算手术时长
data1 <- data1 %>% #手术时长清洗缺失值，这里定义的手术时长是按照麻醉时间来定义的
  mutate( # 补充麻醉开始时间
    ANES_START_TIME = coalesce(
      ANES_START_TIME,
      START_DATE_TIME,
      IN_DATE_TIME
    ),#补充麻醉结束时间
    ANES_END_TIME = coalesce(
      ANES_END_TIME,
      END_DATE_TIME,
      OUT_DATE_TIME ))