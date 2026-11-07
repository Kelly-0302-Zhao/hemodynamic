## =============================================================================
## 1-3 data clean 读术中 event 表、用"有事件记录"筛 data1、补麻醉结束时间、算手术时长
## =============================================================================

#用event数据匹配master，筛选出数据全面的人群
event1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/event.csv", 
                   locale = locale(encoding = "GBK")) # 891543obs
event2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/event.csv", 
                   locale = locale(encoding = "GBK")) # 1564897obs
event3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/event.csv", 
                   locale = locale(encoding = "UTF-8")) # 1224730obs

##  event 表整理成统一格式（ID 转字符、时间统一时区、只留核心列）
standardize_event <- function(df, hospital_name) {df %>%
    mutate(
      # 1. 统一 ID 类型为字符型，确保匹配一致性
      PATIENT_ID = as.character(PATIENT_ID),
      VISIT_ID   = as.character(VISIT_ID),
      OPER_ID    = as.character(OPER_ID),
      # 2. 统一时间格式与时区
      # 针对 event2 为字符的情况，以及 event1/3 为 dttm 的情况
      START_TIME = as.POSIXct(START_TIME, format = "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai"),
      END_TIME   = as.POSIXct(END_TIME,   format = "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai"),
      # 3. 统一数据类型
      DOSAGE = as.numeric(DOSAGE),
      ITEM_NAME = as.character(ITEM_NAME),
      # 4. 增加中心标记
      Hospital = hospital_name) %>%
    # 5. 仅保留分析核心列
    select(PATIENT_ID, VISIT_ID, OPER_ID, ITEM_NAME, DOSAGE, DOSAGE_UNITS, START_TIME, END_TIME, Hospital)}
# 执行清洗并合并 3681170 obs
event <- bind_rows(
  standardize_event(event1, "NanYi"),
  standardize_event(event2, "AnYi"),
  standardize_event(event3, "ZhongDa"))
# 清理原始表以释放内存
rm(event1, event2, event3); gc()
event_id <- event %>% #提取event有的ID
  mutate(
    PATIENT_ID = as.character(PATIENT_ID),
    VISIT_ID   = as.character(VISIT_ID),
    OPER_ID    = as.character(OPER_ID)
  ) %>%
  distinct(PATIENT_ID, VISIT_ID, OPER_ID) # 112618obs
## ============================================================================
## 第 2 步：只保留"能在 event 表里对上 ID"的手术（保证每人都有术中事件记录）102892obs
## ============================================================================
data1 <- data1 %>%
  semi_join(
    event_id,
    by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")
  ) 

## ============================================================================
## 第 3 步：用每台手术"最后一次 event 的时间"去补缺失的麻醉结束时间
## ============================================================================
event_last_time <- event %>%
  group_by(PATIENT_ID, VISIT_ID, OPER_ID) %>%
  summarise(# 先把 END_TIME 和 START_TIME 合并成一个长向量，取其中最大的非空值
    last_event_time = {
      all_times <- c(START_TIME, END_TIME)
      all_times <- all_times[!is.na(all_times)] # 剔除 NA
      if (length(all_times) > 0) max(all_times) else as.POSIXct(NA) },
    .groups = "drop")
event_last_time <- event_last_time %>%
  mutate(
    last_event_time = as.POSIXct(
      last_event_time ))
data1 <- data1 %>%
  left_join(event_last_time, by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")) %>%
  mutate(# 只有当原始结束时间缺失，且补救时间有效时才补
    ANES_END_TIME = coalesce(ANES_END_TIME, last_event_time) )
sum(is.na(data1$ANES_START_TIME)) # 1
sum(is.na(data1$ANES_END_TIME)) # 7

# 检查是否存在异常的时间逻辑
error_times <- data1 %>% 
  filter(ANES_END_TIME <= ANES_START_TIME)
nrow(error_times)#存在5个异常，判断为录入错误（结束时间早于开始时间，由于前期已经标正时间和时差，故不存在时间偏差）
## ============================================================================
## 第 4 步：剔除时间逻辑错误（结束早于开始），正式算手术时长（分钟）
## ============================================================================
data1 <- data1 %>%
  # 1. 剔除时间逻辑错误 (5例) 
  filter(ANES_END_TIME > ANES_START_TIME) %>%
  # 2. 正式计算手术时长（分钟）
  mutate(anes_duration_mins = as.numeric(difftime(ANES_END_TIME, ANES_START_TIME, units = "mins"))
  )
# 再次检查：确认没有负数或 NA（无NA，但最大值为4685mins，存在极端值）
summary(data1$anes_duration_mins)
# 看看超过 12 小时（720分钟）的手术有多少例（38例）
nrow(data1 %>% filter(anes_duration_mins > 720))
# 根据临床意义，剔除或截断超长手术（超过 24 小时）（患者人数从102879降至102872，剔除了7位极端时间患者）
## 剔除超过 24 小时（1440 分钟）的极端手术时长（视为记录错误）
data1 <- data1 %>%
  filter(anes_duration_mins < 1440) # 剔除超过 24 小时的样本，通常认为这些是记录错误

## ============================================================================
## 第 5 步：按手术时长分箱（<60 / 60-120 / 120-180 / >180 分钟），方便描述
## ============================================================================
data1 <- data1 %>% #分箱法查看不同时间手术时长的人数
  mutate(
    anesthesia_duration_bin = case_when(
      anes_duration_mins < 60 ~ "<60 min",
      anes_duration_mins >= 60 & anes_duration_mins < 120 ~ "60-120 min",
      anes_duration_mins >= 120 & anes_duration_mins < 180 ~ "120-180 min",
      anes_duration_mins >= 180 ~ ">180 min",
      TRUE ~ NA_character_ ))
