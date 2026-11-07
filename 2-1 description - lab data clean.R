## =============================================================================
## 2-1 description - lab data clean.R —— 读三家检验 lab 表、标化合并、写出，并跑数据质量摘要
## =============================================================================
demographics <- data1
operation <- event
lab1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/lab.csv", 
                 locale = locale(encoding = "GBK")) # 3489306obs
lab2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/lab.csv", 
                 locale = locale(encoding = "UTF-8")) # 9683536obs
lab3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/lab.csv", 
                 locale = locale(encoding = "UTF-8")) # 4575869obs
vital <- monitor1 # 72460818obs（同一分钟多条记录不去重）
setDT(vital)

# 定义统一时区
MY_TZ <- "Asia/Shanghai"

standardize_lab_data <- function(df, hospital_name) {
  message("正在标准化 Lab 数据: ", hospital_name)
  # 转为 data.table 原地处理
  setDT(df)
  # 统一列名大写，防止合并冲突
  setnames(df, toupper(names(df)))
  # 强制主键和 ID 为字符型
  id_cols <- intersect(c("PATIENT_ID", "VISIT_ID", "OPER_ID", "INP_NO"), names(df))
  df[, (id_cols) := lapply(.SD, function(x) trimws(as.character(x))), .SDcols = id_cols]
  # 统一时间格式与时区
  if ("TIME_POINT" %in% names(df)) {
    df[, TIME_POINT := with_tz(as.POSIXct(TIME_POINT), tzone = MY_TZ)]  }
  # 3. 统一结果和项目名称
  df[, `:=`(
    REPORT_ITEM_NAME = as.character(REPORT_ITEM_NAME),
    RESULT = as.character(RESULT),
    UNITS = as.character(UNITS),
    Hospital = hospital_name  )]
  # 4. 只选择核心列，丢弃杂质列
  keep_cols <- c("PATIENT_ID", "VISIT_ID", "OPER_ID", "Hospital", 
                 "TIME_POINT", "REPORT_ITEM_NAME", "RESULT", "UNITS")
  # 仅保留存在的列
  res <- df[, .SD, .SDcols = intersect(keep_cols, names(df))]
  return(res)}
# 清理 lab
lab1_clean <- standardize_lab_data(lab1, "NanYi")
rm(lab1); gc()
lab2_clean <- standardize_lab_data(lab2, "AnYi")
rm(lab2); gc()
lab3_clean <- standardize_lab_data(lab3, "ZhongDa")
rm(lab3); gc()
## 合并三家 lab 表，写出 clean 大表（约 1775 万行），随后释放中间变量
labs <- bind_rows(lab1_clean, lab2_clean, lab3_clean)
#输出
fwrite(labs, "lab_clean_17748711obs.csv")
gc()
rm(lab1_clean,lab2_clean,lab3_clean)

## 对一张表算"行数 / 唯一首列数 / 唯一前三列组合 / 各列有效数与缺失率"
get_summary <- function(df) {
  data <- get(df, envir = .GlobalEnv)
  unique.1 <- length(unique(data[[1]]))
  unique.2 <- nrow(unique(data[, 1:3]))
  valid_count <- nrow(data) - colSums(is.na(data))
  fill_rate <- round(100 - (colSums(is.na(data)) / nrow(data) * 100),2)
  column_summary <- data.frame(
    column = names(data),
    valid_count = valid_count,
    fill_rate = fill_rate)
  result <- data.frame(
    num_rows = nrow(data),
    unique_first_col = unique.1,
    unique_first_three_cols = unique.2,
    column_summary = column_summary)
  return(result)}
## 对四张核心表（人口学 / 术中事件 / 检验 / 监测）跑数据质量摘要
s1 <- get_summary("demographics")
s2 <- get_summary("operation")
s3 <- get_summary("labs")
s5 <- get_summary("vital")