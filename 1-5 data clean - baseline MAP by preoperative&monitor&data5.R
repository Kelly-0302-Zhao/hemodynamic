## =============================================================================
## 1-5 data clean - clean the large intraoperative monitoring table monitor1 + generate the list of procedures with "preoperative MAP<65" (bp.exclude1)
## =============================================================================

## ============================================================================
## Step 1: Read the intraoperative monitoring monitor tables of the three hospitals, standardize them and combine them into one large table
## ============================================================================
monitor1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/monitor.csv",
                     locale = locale(encoding = "UTF-8")) # 16154084obs
monitor2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/monitor.csv",
                     locale = locale(encoding = "UTF-8")) # 39061431obs
monitor3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/monitor.csv",
                     locale = locale(encoding = "UTF-8")) # 39457555obs
# --- 1. Define the time zone consistent with data1 ---
MY_TZ <- "Asia/Shanghai"
convert_monitor_df <- function(df, hospital_name) {
  message("Standardizing: ", hospital_name)
  # 1. Convert to data.table 
  setDT(df)
  # 2. Uppercase the column names to avoid duplicate columns when combining
  setnames(df, toupper(names(df)))
  # 3. Core: unify the data types and the time zone
  df[, `:=`(
    # Convert all ID columns to trimmed Character
    PATIENT_ID = trimws(as.character(PATIENT_ID)),
    VISIT_ID   = trimws(as.character(VISIT_ID)),
    OPER_ID    = trimws(as.character(OPER_ID)),
    Hospital   = hospital_name,
    # Use lubridate::parse_date_time to handle the text formats of different centers
    TIME_POINT = parse_date_time(as.character(TIME_POINT),
                                 orders = c("ymd HMS", "ymd HM", "mdy HMS"),
                                 tz = MY_TZ),
    # Numeric columns: force conversion to numeric
    PSBP = as.numeric(PSBP),
    PDBP = as.numeric(PDBP),
    MAP  = as.numeric(MAP),
    HR   = as.numeric(HR))]
  # 4. Automatic MAP recovery logic: compute MAP whenever both SBP and DBP are present
  df[is.na(MAP) & !is.na(PSBP) & !is.na(PDBP),
     MAP := round((1/3) * PSBP + (2/3) * PDBP, 1)]
  # 5. Filtering: remove only the truly invalid rows where all three blood pressure measures are missing
  res <- df[!(is.na(PSBP) & is.na(PDBP) & is.na(MAP))]
  # 6. Trim the columns and drop the unnecessary fields
  keep_cols <- c("PATIENT_ID", "VISIT_ID", "OPER_ID", "Hospital", "TIME_POINT", "PSBP", "PDBP", "MAP", "HR")
  return(res[, ..keep_cols])}
m1 <- convert_monitor_df(monitor1, "nanyi"); rm(monitor1); gc() # 9138215
m2 <- convert_monitor_df(monitor2, "anyi"); rm(monitor2); gc() # 28975255
m3 <- convert_monitor_df(monitor3, "zhongda"); rm(monitor3); gc() # 34573998
# --- Final combination ---
monitor1 <- rbindlist(list(m1, m2, m3), use.names = TRUE)
# Free the intermediate variables again
rm(m1, m2, m3); gc()

## ============================================================================
## Step 2: Clean monitor1 (correct abnormal MAP, fill in MAP, remove duplicate timestamps, compute the interval between measurements)
## Intervals > 10 minutes no longer delete points; monitoring interruptions are only used to break "episodes", and the only points to delete are those "outside the anesthesia start-to-end time window"
## ============================================================================
# Cleaning: remove the rows with missing SBP or DBP
monitor1 <- monitor1 %>%
  mutate(
    SBP = PSBP,
    DBP = PDBP,
    MAP = if_else(MAP < 20 | MAP > 225, NA_real_, MAP),
    # Fill in MAP by formula when it is missing
    MAP = if_else(
      is.na(MAP) & !is.na(SBP) & !is.na(DBP),
      SBP / 3 + DBP * 2 / 3,
      MAP   )  ) %>%
  # Remove the rows where MAP is still NA after the computation
  filter(!is.na(MAP)) %>%
  select(
    PATIENT_ID, VISIT_ID, OPER_ID,
    TIME_POINT,
    SBP, DBP, MAP )
# Keep all records with the same timestamp without deduplication: the monitor may measure blood pressure several times within one minute, which is a normal real measurement,
monitor1 <- monitor1 %>%
  arrange(PATIENT_ID, VISIT_ID, OPER_ID, TIME_POINT)
# Compute the interval dt between adjacent points - used only to check the proportion of monitoring interruptions, and no rows are deleted based on it
monitor1 <- monitor1 %>%
  group_by(PATIENT_ID, VISIT_ID, OPER_ID) %>%
  mutate(
    dt = as.numeric(difftime(lead(TIME_POINT), TIME_POINT, units = "mins"))  ) %>%
  ungroup()
nrow(monitor1) # Number of measurement points without deduplication: about 72,460,818

## ============================================================================
## Step 3: Flag the procedures with "MAP<65 already present before the start of surgery" (susceptible population list)
## ============================================================================
bp.exclude1 <- monitor1 %>%
  left_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID, ANES_START_TIME, START_DATE_TIME, END_DATE_TIME)) %>%
  semi_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID)) %>%
  filter(TIME_POINT <= START_DATE_TIME) %>%
  filter(MAP < 65) %>%
  mutate(IOH_EXCLUDE = 1) %>%
  distinct(PATIENT_ID, VISIT_ID, OPER_ID) # MAP<65 present before the start of surgery
total_dupes <- sum(duplicated(data1)) # Check whether there are fully duplicated rows; the result is 0
