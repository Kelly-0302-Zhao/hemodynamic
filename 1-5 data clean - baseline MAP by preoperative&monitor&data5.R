## =============================================================================
## 1-5 data clean -清洗术中监测大表 monitor1 + 生成"术前已有 MAP<65"手术名单 bp.exclude1
## =============================================================================

## ============================================================================
## 第 1 步：读三家术中监测 monitor 表，标化合并成一张大表
## ============================================================================
monitor1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/monitor.csv",
                     locale = locale(encoding = "UTF-8")) # 16154084obs
monitor2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/monitor.csv",
                     locale = locale(encoding = "UTF-8")) # 39061431obs
monitor3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/monitor.csv",
                     locale = locale(encoding = "UTF-8")) # 39457555obs
# --- 1. 定义与 data1 一致的时区 ---
MY_TZ <- "Asia/Shanghai"
convert_monitor_df <- function(df, hospital_name) {
  message("正在标化: ", hospital_name)
  # 1. 转换为 data.table 
  setDT(df)
  # 2. 统一列名大写，确保合并时不产生重复列
  setnames(df, toupper(names(df)))
  # 3. 核心：统一数据类型与时区
  df[, `:=`(
    # ID 类全部强制转为去空格的 Character
    PATIENT_ID = trimws(as.character(PATIENT_ID)),
    VISIT_ID   = trimws(as.character(VISIT_ID)),
    OPER_ID    = trimws(as.character(OPER_ID)),
    Hospital   = hospital_name,
    # 使用 lubridate::parse_date_time 兼容不同中心的文本格式
    TIME_POINT = parse_date_time(as.character(TIME_POINT),
                                 orders = c("ymd HMS", "ymd HM", "mdy HMS"),
                                 tz = MY_TZ),
    # 数值类：强制转 numeric
    PSBP = as.numeric(PSBP),
    PDBP = as.numeric(PDBP),
    MAP  = as.numeric(MAP),
    HR   = as.numeric(HR))]
  # 4. MAP 自动打捞 逻辑：只要 SBP/DBP 都在，就算出 MAP
  df[is.na(MAP) & !is.na(PSBP) & !is.na(PDBP),
     MAP := round((1/3) * PSBP + (2/3) * PDBP, 1)]
  # 5. 过滤：只剔除三项血压指标全空的真正无效行
  res <- df[!(is.na(PSBP) & is.na(PDBP) & is.na(MAP))]
  # 6. 精简列，释放不必要的字段
  keep_cols <- c("PATIENT_ID", "VISIT_ID", "OPER_ID", "Hospital", "TIME_POINT", "PSBP", "PDBP", "MAP", "HR")
  return(res[, ..keep_cols])}
m1 <- convert_monitor_df(monitor1, "nanyi"); rm(monitor1); gc() # 9138215
m2 <- convert_monitor_df(monitor2, "anyi"); rm(monitor2); gc() # 28975255
m3 <- convert_monitor_df(monitor3, "zhongda"); rm(monitor3); gc() # 34573998
# --- 最终合并 ---
monitor1 <- rbindlist(list(m1, m2, m3), use.names = TRUE)
# 再次释放中间变量
rm(m1, m2, m3); gc()

## ============================================================================
## 第 2 步：清洗 monitor1（修正 MAP 异常、补算 MAP、去掉重复时间戳、算测点间隔）
## 间隔 >10 分钟不再删点，监测中断只用来断开“事件"，真正要删的只有"超出麻醉开始—结束时间窗"的测点
## ============================================================================
# 清洗：去掉SBP和DBP为空的行
monitor1 <- monitor1 %>%
  mutate(
    SBP = PSBP,
    DBP = PDBP,
    MAP = if_else(MAP < 20 | MAP > 225, NA_real_, MAP),
    # MAP 缺失时用公式补充
    MAP = if_else(
      is.na(MAP) & !is.na(SBP) & !is.na(DBP),
      SBP / 3 + DBP * 2 / 3,
      MAP   )  ) %>%
  # 删除计算后MAP仍为NA的行
  filter(!is.na(MAP)) %>%
  select(
    PATIENT_ID, VISIT_ID, OPER_ID,
    TIME_POINT,
    SBP, DBP, MAP )
# 同一时间戳的多条记录全部保留，不去重：监测仪在一分钟之内可能测好几次血压，这是正常的真实测量，
monitor1 <- monitor1 %>%
  arrange(PATIENT_ID, VISIT_ID, OPER_ID, TIME_POINT)
# 计算相邻节点间隔 dt —— 只用来核对监测中断的比例，不据此删任何行
monitor1 <- monitor1 %>%
  group_by(PATIENT_ID, VISIT_ID, OPER_ID) %>%
  mutate(
    dt = as.numeric(difftime(lead(TIME_POINT), TIME_POINT, units = "mins"))  ) %>%
  ungroup()
nrow(monitor1) # 不去重的测点数：约 72,460,818

## ============================================================================
## 第 3 步：标出"手术开始前就存在 MAP<65"的手术（易感人群名单）
## ============================================================================
bp.exclude1 <- monitor1 %>%
  left_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID, ANES_START_TIME, START_DATE_TIME, END_DATE_TIME)) %>%
  semi_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID)) %>%
  filter(TIME_POINT <= START_DATE_TIME) %>%
  filter(MAP < 65) %>%
  mutate(IOH_EXCLUDE = 1) %>%
  distinct(PATIENT_ID, VISIT_ID, OPER_ID) # 手术开始前存在MAP<65
total_dupes <- sum(duplicated(data1)) # 判断是否有全部重复的行，结果为0
