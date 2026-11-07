## =============================================================================
## 1-6 data clean —— 再筛 data1（监测/事件对得上）、清年龄性别、算急诊/ICU、写出 clean 表
## =============================================================================
# 提取在术中时间窗内有监测点的 ID 清单
valid_analysis_ids <- monitor1 %>%
  inner_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID, START_DATE_TIME, END_DATE_TIME), 
             by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")) %>%
  # 核心：必须有术中监测点
  filter(TIME_POINT >= START_DATE_TIME & TIME_POINT <= END_DATE_TIME) %>%
  distinct(PATIENT_ID, VISIT_ID, OPER_ID) # 85445obs
# 匹配有监测点的人群，才能进一步分析 
## ============================================================================
## 第 1 步：只保留"术中时间窗内有监测点"的手术
## ============================================================================
data1 <- data1 %>%
  semi_join(valid_analysis_ids, by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")) # 85445obs

event_id <- event %>%
  mutate(
    PATIENT_ID = as.character(PATIENT_ID),
    VISIT_ID   = as.character(VISIT_ID),
    OPER_ID    = as.character(OPER_ID)
  ) %>%
  distinct(PATIENT_ID, VISIT_ID, OPER_ID) # 112618obs
## ============================================================================
## 第 2 步：再只保留能在 event 表对上 ID 的手术（85445 行）
## ============================================================================
data1 <- data1 %>%
  semi_join(
    event_id,
    by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")
  ) # 85445 obs

## ============================================================================
## 第 3 步：删掉 AGE 缺失、SEX 非男/女的行（85397 行）
## ============================================================================
data1 <- data1 %>% filter(!is.na(AGE)  & SEX %in% c("男","女")) # 85397obs
#判断索引值是否有NA
sum(is.na(data1$PATIENT_ID))
sum(is.na(data1$VISIT_ID))
sum(is.na(data1$OPER_ID)) # 都没有

## ============================================================================
## 第 4 步：算"是否进 ICU / 是否急诊"标志，并生成每台手术唯一 cluster_id
## ============================================================================
data1 <- data1 %>%
  mutate(
    # 是否进入 ICU
    IN_ICU_INDICATOR = if_else(
      !is.na(IN_ICU_TIME),
      1L,
      0L
    ),
    
    # 急诊指标只依据原始 EMERGENCY_INDICATOR
    EMERGENCY_INDICATOR = replace_na(
      as.numeric(EMERGENCY_INDICATOR),
      0
    ),
    
    EMERGENCY_INDICATOR = if_else(
      EMERGENCY_INDICATOR > 0,
      1L,
      0L
    )
  )


data1 <- data1 %>%
  mutate(
    PATIENT_ID   = as.character(PATIENT_ID),
    VISIT_ID     = as.character(VISIT_ID),
    OPER_ID      = as.character(OPER_ID),
    cluster_id   = paste(PATIENT_ID, VISIT_ID, OPER_ID, sep = "_")  # 每个手术唯一ID
  )
length(unique(data1$cluster_id)) # 85397 obs

## ============================================================================
## 第 5 步：把清洗好的四张表写出 CSV（落在 20260827\JCA_run\）
## ============================================================================
fwrite(data1, "master_clean_85397obs.csv")
fwrite(monitor1, "monitor_clean_72460818obs.csv")
fwrite(event, "event_clean_3681170obs.csv")
fwrite(diag1, "术前诊断_118550obs.csv")

