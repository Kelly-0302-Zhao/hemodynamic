## =============================================================================
## 1-2 data clean - ASA&GENERAL&DURATION.R -- unify ASA grading, keep only general anesthesia, complete anesthesia start and end times
## =============================================================================

#Standardize ASA grading
## ============================================================================
## Step 1: unify the various ASA notations across hospitals into I/II/III (record NA when unrecognized)
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
## Step 2: keep only "general anesthesia" surgeries (excluding those that are not general anesthesia or cannot be confirmed, 2682 patients)
## ============================================================================
general <- read_csv(
  "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/General Anesthesia.csv",
  locale = locale(encoding = "GBK")) %>% 
  distinct() #Import the general anesthesia list compiled with clinical assistance
# 1. Preprocess the dictionary table (make sure the text contains no spaces, so that matching is not affected)
general <- general %>%
  mutate(ANESTHESIA_METHOD = trimws(as.character(ANESTHESIA_METHOD))) %>%
  distinct()
# 2. Filter the merged data1 (2682 patients who were not general anesthesia or could not be confirmed were excluded, leaving 102892 patients)
data1 <- data1 %>%
  # Join the general anesthesia dictionary
  left_join(general, by = "ANESTHESIA_METHOD") %>%
  filter(
    GENERAL == 1 | 
      grepl("全麻|全身麻醉|静吸复合|插管|喉罩", ANESTHESIA_METHOD)) %>%
  # Remove helper columns
  select(-GENERAL)

## ============================================================================
## Step 3: complete anesthesia start/end times (fall back to room-in/room-out times when missing) in preparation for computing duration
## ============================================================================
#Compute the surgery duration
data1 <- data1 %>% #Clean missing values in surgery duration; surgery duration is defined here by anesthesia time
  mutate( # Fill in anesthesia start times
    ANES_START_TIME = coalesce(
      ANES_START_TIME,
      START_DATE_TIME,
      IN_DATE_TIME
    ),#Fill in anesthesia end times
    ANES_END_TIME = coalesce(
      ANES_END_TIME,
      END_DATE_TIME,
      OUT_DATE_TIME ))