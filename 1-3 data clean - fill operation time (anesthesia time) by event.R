## =============================================================================
## 1-3 data clean read the intraoperative event table, filter data1 by having event records, fill in anesthesia end times, compute surgery duration
## =============================================================================

#Match event data with master and keep patients with complete data
event1 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/nanyi/event.csv", 
                   locale = locale(encoding = "GBK")) # 891543obs
event2 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/anyi/event.csv", 
                   locale = locale(encoding = "GBK")) # 1564897obs
event3 <- read_csv("D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/zhongda/event.csv", 
                   locale = locale(encoding = "UTF-8")) # 1224730obs

##  Standardize the event table into a unified format (ID to character, unified time zone, core columns only)
standardize_event <- function(df, hospital_name) {df %>%
    mutate(
      # 1. Unify ID type as character to ensure consistent matching
      PATIENT_ID = as.character(PATIENT_ID),
      VISIT_ID   = as.character(VISIT_ID),
      OPER_ID    = as.character(OPER_ID),
      # 2. Unify time format and time zone
      # Handle the case where event2 is character and where event1/3 are dttm
      START_TIME = as.POSIXct(START_TIME, format = "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai"),
      END_TIME   = as.POSIXct(END_TIME,   format = "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai"),
      # 3. Unify data types
      DOSAGE = as.numeric(DOSAGE),
      ITEM_NAME = as.character(ITEM_NAME),
      # 4. Add the center label
      Hospital = hospital_name) %>%
    # 5. Keep only the core analysis columns
    select(PATIENT_ID, VISIT_ID, OPER_ID, ITEM_NAME, DOSAGE, DOSAGE_UNITS, START_TIME, END_TIME, Hospital)}
# Run the cleaning and merge, 3681170 obs
event <- bind_rows(
  standardize_event(event1, "NanYi"),
  standardize_event(event2, "AnYi"),
  standardize_event(event3, "ZhongDa"))
# Clear the original tables to free memory
rm(event1, event2, event3); gc()
event_id <- event %>% #Extract the IDs present in event
  mutate(
    PATIENT_ID = as.character(PATIENT_ID),
    VISIT_ID   = as.character(VISIT_ID),
    OPER_ID    = as.character(OPER_ID)
  ) %>%
  distinct(PATIENT_ID, VISIT_ID, OPER_ID) # 112618obs
## ============================================================================
## Step 2: keep only surgeries whose ID matches in the event table (so every patient has intraoperative event records), 102892obs
## ============================================================================
data1 <- data1 %>%
  semi_join(
    event_id,
    by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")
  ) 

## ============================================================================
## Step 3: use the time of the last event of each surgery to fill in missing anesthesia end times
## ============================================================================
event_last_time <- event %>%
  group_by(PATIENT_ID, VISIT_ID, OPER_ID) %>%
  summarise(# First combine END_TIME and START_TIME into one long vector and take the largest non-missing value
    last_event_time = {
      all_times <- c(START_TIME, END_TIME)
      all_times <- all_times[!is.na(all_times)] # Remove NA
      if (length(all_times) > 0) max(all_times) else as.POSIXct(NA) },
    .groups = "drop")
event_last_time <- event_last_time %>%
  mutate(
    last_event_time = as.POSIXct(
      last_event_time ))
data1 <- data1 %>%
  left_join(event_last_time, by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")) %>%
  mutate(# Fill only when the original end time is missing and the backfilled time is valid
    ANES_END_TIME = coalesce(ANES_END_TIME, last_event_time) )
sum(is.na(data1$ANES_START_TIME)) # 1
sum(is.na(data1$ANES_END_TIME)) # 7

# Check for abnormal time logic
error_times <- data1 %>% 
  filter(ANES_END_TIME <= ANES_START_TIME)
nrow(error_times)#5 abnormalities were found and judged to be data entry errors (end time earlier than start time; since times and time differences were already corrected earlier, there is no time offset)
## ============================================================================
## Step 4: remove time logic errors (end earlier than start) and formally compute surgery duration (minutes)
## ============================================================================
data1 <- data1 %>%
  # 1. Remove time logic errors (5 cases) 
  filter(ANES_END_TIME > ANES_START_TIME) %>%
  # 2. Formally compute surgery duration (minutes)
  mutate(anes_duration_mins = as.numeric(difftime(ANES_END_TIME, ANES_START_TIME, units = "mins"))
  )
# Check again: confirm there are no negatives or NA (no NA, but the maximum is 4685mins, so extreme values exist)
summary(data1$anes_duration_mins)
# See how many surgeries exceed 12 hours (720 minutes) (38 cases)
nrow(data1 %>% filter(anes_duration_mins > 720))
# Based on clinical meaning, remove or truncate extremely long surgeries (over 24 hours) (patients decreased from 102879 to 102872, removing 7 patients with extreme durations)
## Remove extreme surgery durations over 24 hours (1440 minutes) (treated as recording errors)
data1 <- data1 %>%
  filter(anes_duration_mins < 1440) # Remove samples over 24 hours, which are generally regarded as recording errors

## ============================================================================
## Step 5: bin the surgery duration (<60 / 60-120 / 120-180 / >180 minutes) for description
## ============================================================================
data1 <- data1 %>% #Use binning to view patient counts by surgery duration
  mutate(
    anesthesia_duration_bin = case_when(
      anes_duration_mins < 60 ~ "<60 min",
      anes_duration_mins >= 60 & anes_duration_mins < 120 ~ "60-120 min",
      anes_duration_mins >= 120 & anes_duration_mins < 180 ~ "120-180 min",
      anes_duration_mins >= 180 ~ ">180 min",
      TRUE ~ NA_character_ ))
