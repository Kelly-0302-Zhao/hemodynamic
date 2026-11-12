## =============================================================================
## 1-7 data clean - group of susceptible.R - build the "susceptible population" subgroup (MAP<65 before drug administration)
## =============================================================================
#Build the population in which hypotension occurred before drug administration as the susceptible population, for subsequent subgroup analysis
monitor_susceptible <- monitor1 %>%
  semi_join(bp.exclude1, by = c("PATIENT_ID", "VISIT_ID", "OPER_ID"))

## ============================================================================
## Step 1: Link the susceptible monitoring points to the "procedure timeline + drug group" and assign cluster_id
## ============================================================================
monitor_susceptible_final <- monitor_susceptible %>%
  left_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID, 
                             Group, Hospital, 
                             START_DATE_TIME, END_DATE_TIME), 
            by = c("PATIENT_ID", "VISIT_ID", "OPER_ID"))
# Count the data volume of this subset ---
message("Total number of measurement points in the susceptible population: ", nrow(monitor_susceptible_final)) # 25338320 obs
monitor_susceptible_final <- monitor_susceptible_final %>%
  mutate(cluster_id   = paste(PATIENT_ID, VISIT_ID, OPER_ID, sep = "_") )
# Save as a separate subgroup analysis file ---
# Output CSV
fwrite(monitor_susceptible_final, "monitor_below65_before_start.csv")
#Obtain the baseline data table of the susceptible population
## ============================================================================
## Step 2: Take the unique procedure list of the susceptible population, build the baseline table from the first record of each procedure, and write two CSV files
## ============================================================================
susceptible_patients_list <- monitor_susceptible_final %>%
  distinct(cluster_id)
message("Number of unique procedures in the susceptible population: ", nrow(susceptible_patients_list)) # 27594 obs
# 2. Extract the baseline data of these patients from data1 and deduplicate; for patients undergoing multiple surgeries, the unit is the procedure
baseline_data_susceptible <- data1 %>%
  semi_join(susceptible_patients_list, by = "cluster_id") %>%
  arrange(cluster_id, START_DATE_TIME) %>%
  group_by(cluster_id) %>%
  slice(1) %>%
  ungroup()
# 3. Report the patient distribution of this subgroup
message("Susceptible population baseline table (N=", nrow(baseline_data_susceptible), ") built.") # n=27449
print(table(baseline_data_susceptible$Group)) # Ciprofol_only 1051 : Propofol_only 26398 
# Save for subsequent subgroup analysis
fwrite(baseline_data_susceptible, "monitor_below65_before_start_baseline.csv")

## ============================================================================
## Step 3: Selection bias check - what percentage of each group is the "susceptible population" that was already hypotensive before surgery
## ============================================================================
## 3.1 Flag each procedure with the "susceptible" status and the year
sb <- data1 %>%
  mutate(
    susceptible = as.integer(cluster_id %in% baseline_data_susceptible$cluster_id),
    year        = as.integer(substr(as.character(START_DATE_TIME), 1, 4))
  )

## 3.2 Two functions: one computes n and the percentage for each group, the other computes the between-group RR (ciprofol vs propofol)
sb_count <- function(dt, level, hosp = NA_character_, yr = NA_integer_) {
  dt %>%
    filter(Group %in% c("Ciprofol_only", "Propofol_only")) %>%
    group_by(Group) %>%
    summarise(n_surgery       = n(),
              n_susceptible   = sum(susceptible),
              pct_susceptible = round(100 * mean(susceptible), 2),
              .groups = "drop") %>%
    mutate(level = level, Hospital = hosp, year = yr) %>%
    select(level, Hospital, year, Group, n_surgery, n_susceptible, pct_susceptible)
}
sb_rr <- function(dt, level, hosp = NA_character_, yr = NA_integer_) {
  dt <- dt %>% filter(Group %in% c("Ciprofol_only", "Propofol_only"))
  n1 <- sum(dt$Group == "Ciprofol_only"); a <- sum(dt$Group == "Ciprofol_only" & dt$susceptible == 1L)
  n0 <- sum(dt$Group == "Propofol_only");  b <- sum(dt$Group == "Propofol_only"  & dt$susceptible == 1L)
  if (n1 == 0 || n0 == 0) return(NULL)                     # This stratum contains only one drug, so RR cannot be computed
  p1 <- a / n1; p0 <- b / n0
  if (a == 0 || b == 0) {                                  # When the numerator is 0, use the +0.5 correction and do not give CI
    rr <- ((a + 0.5) / (n1 + 1)) / ((b + 0.5) / (n0 + 1)); se <- NA_real_
  } else {
    rr <- p1 / p0; se <- sqrt(1/a - 1/n1 + 1/b - 1/n0)
  }
  data.frame(level = level, Hospital = hosp, year = yr,
             n_ciprofol   = n1, pct_ciprofol = round(100 * p1, 2),
             n_propofol   = n0, pct_propofol = round(100 * p0, 2),
             RR = round(rr, 3),
             CI_low  = round(exp(log(rr) - 1.96 * se), 3),
             CI_high = round(exp(log(rr) + 1.96 * se), 3))
}

## 3.3 Assemble three tables: all three hospitals / stratified by hospital / ZhongDa by year
sb_year_list <- sort(unique(sb$year[sb$Hospital == "ZhongDa" & !is.na(sb$year)]))
sb_share <- bind_rows(
  sb_count(sb, "全部三院"),
  sb_count(sb %>% filter(Hospital == "ZhongDa"), "仅 ZhongDa（两药都有）", "ZhongDa"),
  sb_count(sb %>% filter(Hospital == "NanYi"),   "仅 NanYi（环泊酚仅 2024 年）", "NanYi"),
  sb_count(sb %>% filter(Hospital == "AnYi"),    "仅 AnYi（零环泊酚）", "AnYi"),
  bind_rows(lapply(sb_year_list, function(y)
    sb_count(sb %>% filter(Hospital == "ZhongDa", year == y), "ZhongDa 按年份", "ZhongDa", y)))
)
sb_rr_tab <- bind_rows(
  sb_rr(sb, "全部三院"),
  sb_rr(sb %>% filter(Hospital == "ZhongDa"), "仅 ZhongDa（两药都有）", "ZhongDa"),
  sb_rr(sb %>% filter(Hospital == "ZhongDa", year %in% 2022:2024), "ZhongDa 2022-2024（同期）", "ZhongDa"),
  bind_rows(lapply(sb_year_list, function(y)
    sb_rr(sb %>% filter(Hospital == "ZhongDa", year == y), "ZhongDa 按年份", "ZhongDa", y)))
)

## 3.4 Print
cat("\n[Selection bias check] Susceptible population (preoperative MAP<65) as a percentage of procedures in each group\n")
print(as.data.frame(sb_share), row.names = FALSE)
cat("\n[Between-group risk ratio] RR = susceptible share in the Ciprofol group / susceptible share in the Propofol group (RR>1 means more severe patients in the Ciprofol group)\n")
print(as.data.frame(sb_rr_tab), row.names = FALSE) #Significant in all three hospitals or in ZhongDa only, with RR < 1 in both, indicating that more severely affected patients are in the propofol group
fwrite(sb_share,  "Susceptible_Share_By_Group.csv")
fwrite(sb_rr_tab, "Susceptible_Share_RR.csv")
