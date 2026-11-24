## =============================================================================
## 2-1 description - lab data clean.R - read the lab tables from the three hospitals, standardize and merge them, write them out, and run a data quality summary
## =============================================================================
demographics <- data1
operation <- event
lab1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/lab.csv", 
                 locale = locale(encoding = "GBK")) # 3489306obs
lab2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/lab.csv", 
                 locale = locale(encoding = "UTF-8")) # 9683536obs
lab3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/lab.csv", 
                 locale = locale(encoding = "UTF-8")) # 4575869obs
vital <- monitor1 # 72460818obs (multiple records within the same minute are not deduplicated)
setDT(vital)

# Define a unified time zone
MY_TZ <- "Asia/Shanghai"

standardize_lab_data <- function(df, hospital_name) {
  message("Standardizing Lab data: ", hospital_name)
  # Convert to data.table for in-place processing
  setDT(df)
  # Uppercase column names to avoid merge conflicts
  setnames(df, toupper(names(df)))
  # Coerce the primary key and ID columns to character
  id_cols <- intersect(c("PATIENT_ID", "VISIT_ID", "OPER_ID", "INP_NO"), names(df))
  df[, (id_cols) := lapply(.SD, function(x) trimws(as.character(x))), .SDcols = id_cols]
  # Standardize the time format and time zone
  if ("TIME_POINT" %in% names(df)) {
    df[, TIME_POINT := with_tz(as.POSIXct(TIME_POINT), tzone = MY_TZ)]  }
  # 3. Standardize the result and item names
  df[, `:=`(
    REPORT_ITEM_NAME = as.character(REPORT_ITEM_NAME),
    RESULT = as.character(RESULT),
    UNITS = as.character(UNITS),
    Hospital = hospital_name  )]
  # 4. Keep only the core columns and drop the extraneous ones
  keep_cols <- c("PATIENT_ID", "VISIT_ID", "OPER_ID", "Hospital", 
                 "TIME_POINT", "REPORT_ITEM_NAME", "RESULT", "UNITS")
  # Keep only the columns that exist
  res <- df[, .SD, .SDcols = intersect(keep_cols, names(df))]
  return(res)}
# Clean the lab data
lab1_clean <- standardize_lab_data(lab1, "NanYi")
rm(lab1); gc()
lab2_clean <- standardize_lab_data(lab2, "AnYi")
rm(lab2); gc()
lab3_clean <- standardize_lab_data(lab3, "ZhongDa")
rm(lab3); gc()
## Merge the lab tables from the three hospitals, write out the clean master table (about 17.75 million rows), then release the intermediate variables
labs <- bind_rows(lab1_clean, lab2_clean, lab3_clean)
#Output
fwrite(labs, "lab_clean_17748711obs.csv")
gc()
rm(lab1_clean,lab2_clean,lab3_clean)

## Compute "row count / number of unique values in the first column / number of unique combinations of the first three columns / valid count and missing rate for each column" for one table
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
## Run a data quality summary on the four core tables (demographics / intraoperative events / lab / monitoring)
s1 <- get_summary("demographics")
s2 <- get_summary("operation")
s3 <- get_summary("labs")
s5 <- get_summary("vital")