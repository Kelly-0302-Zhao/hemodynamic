## =============================================================================
## 3-1 GEE analysis - build the surgery-level analysis table
## =============================================================================
## ============================================================================
## Step 1: preparation (load packages, set paths, hard-code the constants)
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # read and process tables
  library(geepack)      # Used later for the GEE models
})

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
DATA <- file.path(D, "JCA_run")      # Cleaned intermediate tables
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## Three ID columns. The ID column types may differ across the three tables, so convert them all to character before use
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## Parse all timestamps as UTC (raw data carry a Z suffix)
MY_TZ <- "UTC"

## ============================================================================
## Step 2: read the surgery base table, parse the anesthesia start/end time, and arrange it as "one row = one surgery"
## ============================================================================
parse_iso <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  as.POSIXct(gsub("T", " ", sub("Z$", "", as.character(x))),
             tz = MY_TZ, format = "%Y-%m-%d %H:%M:%S")
}
to_chr_ids <- function(DT) {
  DT[, (ID_VARS) := lapply(.SD, as.character), .SDcols = ID_VARS]
  DT
}

master <- fread(file.path(DATA, "master_clean_85397obs.csv"), encoding = "UTF-8")
master <- to_chr_ids(master)
master[, anes_start := parse_iso(ANES_START_TIME)]
master[, anes_end   := parse_iso(ANES_END_TIME)]
master[, year       := as.integer(format(as.Date(anes_start), "%Y"))]
## Group factor: propofol as the reference group (ciprofol vs propofol)
master[, Group := relevel(factor(as.character(Group)), ref = "Propofol_only")]
## Convert the following categorical variables to factors for the subsequent models
for (v in c("ASA_GRADE", "operation_type", "Hospital", "SEX")) {
  if (v %in% names(master)) master[[v]] <- factor(master[[v]])
}

## ============================================================================
## Step 3: read the surgery-level hypotension outcome table and join it to the base table by ID
## ============================================================================
surg <- setDT(readRDS(file.path(OUT, "Final_Analysis_Surgery.rds")))

## Select the needed columns from the base table
keep <- c(ID_VARS, "Group", "Hospital", "year", "anes_start", "anes_end",
          "anes_duration_mins", "AGE", "SEX", "BMI", "ASA_GRADE",
          "EMERGENCY_INDICATOR", "operation_type",
          "Hypertension", "Diabetes", "Coronary_artery_disease",
          "Heart_failure", "IN_ICU_INDICATOR")

d <- master[, ..keep][surg[, !c("anes_start", "anes_end", "anes_duration_mins"),
                           with = FALSE], on = ID_VARS]



## ============================================================================
## Step 4: add the derived columns needed for the analysis to the master table
## ============================================================================
d[, lg_dur := log(anes_duration_mins)]
d[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
d[, year_f  := factor(year)]
d[, any_hypo := as.integer(n_epi >= 1L)]

cat("Derived columns added: lg_dur / sid_int / year_f / any_hypo\n")


## ============================================================================
## Step 5: unit check (print "procedures vs patients" and look at group × center / year)
## ============================================================================
cat("=============== Unit check (very important)===============\n")
cat("Analysis table rows (= procedures):", nrow(d), "\n")
cat("Of which, patients                :", uniqueN(d$PATIENT_ID), "\n")
cat("Patients with multiple procedures :",
    d[, .N, by = PATIENT_ID][N > 1, .N], "patients (a total of",
    d[, .N, by = PATIENT_ID][N > 1, sum(N)], "procedures)\n")
cat("Procedures with analyzable outcome:", d[!is.na(prolonged), .N], "\n")
cat("Patients with analyzable outcome  :", uniqueN(d[!is.na(prolonged)]$PATIENT_ID), "\n\n")

cat("---- Group × center (procedures)----\n")
print(dcast(d[!is.na(prolonged)], Hospital ~ Group, value.var = "prolonged", fun.aggregate = length))
cat("\n---- Group × year (procedures)----\n")
print(dcast(d[!is.na(prolonged)], year ~ Group, value.var = "prolonged", fun.aggregate = length))


## ============================================================================
## Step 6: save the master surgery-level analysis table
## ============================================================================
saveRDS(d, file.path(OUT, "Final_Analysis_Long_Table.rds"))
fwrite(d[, .(PATIENT_ID, VISIT_ID, OPER_ID, Hospital, year, Group, n_epi, n_long,
             prolonged, prolonged_cum, multi_epi, short_only, any_hypo,
             MAP65_TIME, MAP65_AUC, n_map, rate_map, anes_duration_mins,
             AGE, SEX, BMI, ASA_GRADE, Hypertension, Diabetes,
             Coronary_artery_disease, Heart_failure)],
       file.path(OUT, "Final_Analysis_Long_Table.csv"))
