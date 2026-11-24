## =============================================================================
## 1-6 data clean - further filter data1 (monitoring/events match), clean age and sex, compute the emergency/ICU flags, and write the clean tables
## =============================================================================
# Extract the list of IDs with monitoring points within the intraoperative time window
valid_analysis_ids <- monitor1 %>%
  inner_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID, START_DATE_TIME, END_DATE_TIME), 
             by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")) %>%
  # Core: must have intraoperative monitoring points
  filter(TIME_POINT >= START_DATE_TIME & TIME_POINT <= END_DATE_TIME) %>%
  distinct(PATIENT_ID, VISIT_ID, OPER_ID) # 85445obs
# Match the population with monitoring points before further analysis 
## ============================================================================
## Step 1: Keep only the procedures with "monitoring points within the intraoperative time window"
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
## Step 2: Then keep only the procedures whose IDs match the event table (85445 rows)
## ============================================================================
data1 <- data1 %>%
  semi_join(
    event_id,
    by = c("PATIENT_ID", "VISIT_ID", "OPER_ID")
  ) # 85445 obs

## ============================================================================
## Step 3: Remove the rows with missing AGE or with SEX other than male/female (85397 rows)
## ============================================================================
data1 <- data1 %>% filter(!is.na(AGE)  & SEX %in% c("男","女")) # 85397obs
#Check whether the index values contain NA
sum(is.na(data1$PATIENT_ID))
sum(is.na(data1$VISIT_ID))
sum(is.na(data1$OPER_ID)) # None of them

## ============================================================================
## Step 4: Compute the "ICU admission / emergency" flags and generate a unique cluster_id for each procedure
## ============================================================================
data1 <- data1 %>%
  mutate(
    # Whether the patient entered ICU
    IN_ICU_INDICATOR = if_else(
      !is.na(IN_ICU_TIME),
      1L,
      0L
    ),
    
    # The emergency indicator relies only on the original EMERGENCY_INDICATOR
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
    cluster_id   = paste(PATIENT_ID, VISIT_ID, OPER_ID, sep = "_")  # Unique ID for each procedure
  )
length(unique(data1$cluster_id)) # 85397 obs

## ============================================================================
## Step 5: Write the four cleaned tables to CSV (into 20260827\JCA_run\)
## ============================================================================
fwrite(data1, "master_clean_85397obs.csv")
fwrite(monitor1, "monitor_clean_72460818obs.csv")
fwrite(event, "event_clean_3681170obs.csv")
fwrite(diag1, "术前诊断_118550obs.csv")

