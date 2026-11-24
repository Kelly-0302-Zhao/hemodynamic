## =============================================================================
## 2-3 description - descriptive statistics of hypotension (Table 2)
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)   # read and process tables
  library(geepack)      # run GEE models
})

## Root directory. All inputs and outputs are placed under this directory
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
DATA <- file.path(D, "JCA_run")   # Home of the cleaned intermediate tables (master / monitor / event / lab)
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")     # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

## Three ID columns. The ID column types are inconsistent across the three tables, so convert them all to character before merging (otherwise the merge will misalign)
ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
MY_TZ    <- "UTC"      # Raw timestamps carry a Z suffix; parse everything as UTC throughout to keep the internal representation consistent
GAP_LIMIT <- 10        # Minutes: a gap of more than 10 minutes between consecutive measurements is treated as a "monitoring gap"; hypotension duration is not accumulated across a gap

## : convert ISO time strings such as "2024-01-01T08:30:00Z" to POSIXct
parse_iso <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  as.POSIXct(gsub("T", " ", sub("Z$", "", as.character(x))),
             tz = MY_TZ, format = "%Y-%m-%d %H:%M:%S")
}
## Convert the three ID columns to character
to_chr_ids <- function(DT) {
  DT[, (ID_VARS) := lapply(.SD, as.character), .SDcols = ID_VARS]
  DT
}

cat("[Common setup] ready | Output directory:", OUT, "\n")


## ============================================================================
## Step 1: read the surgery table / event table and add group, center, and year
## ============================================================================

d   <- setDT(readRDS(file.path(OUT, "Final_Analysis_Surgery.rds")))
epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## Read the surgery base table (one row = one surgery) and take group / center / year to split the two groups for the descriptive summary
master <- fread(file.path(DATA, "master_clean_85397obs.csv"), encoding = "UTF-8")
master <- to_chr_ids(master)
master[, anes_start := parse_iso(ANES_START_TIME)]
master[, anes_end   := parse_iso(ANES_END_TIME)]
master[, year       := as.integer(format(as.Date(anes_start), "%Y"))]
master[, Group      := relevel(factor(as.character(Group)), ref = "Propofol_only")]
for (v in c("ASA_GRADE", "operation_type", "Hospital", "SEX")) {
  if (v %in% names(master)) master[[v]] <- factor(master[[v]])
}
info <- master[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, Hospital, year)]

## Attach group / center / year to the surgery table, then delete the surgeries with missing "prolonged hypotension" information
d   <- info[d, on = ID_VARS]
d   <- d[!is.na(prolonged)]
d[, any_hypo := as.integer(n_epi >= 1L)]

cat("Surgery table:", nrow(d), "procedures (prolonged missing excluded)\n")
cat("Event table:", nrow(epi), "hypotensive episodes\n\n")


## ============================================================================
## Step 2: compute the surgery-level descriptive table
## ============================================================================

by_surg <- d[, .(
  `Surgical procedures, n`                = .N,
  `Patients, n`                    = uniqueN(PATIENT_ID),
  `Any IOH (MAP <65 mm Hg), n (%)`        = sprintf("%d (%.1f%%)", sum(any_hypo),
                                          100 * mean(any_hypo)),
  `Prolonged IOH (>10 min), n (%)`     = sprintf("%d (%.1f%%)", sum(prolonged),
                                          100 * mean(prolonged)),
  `Cumulative IOH (>10 min), n (%)`           = sprintf("%d (%.1f%%)", sum(prolonged_cum),
                                          100 * mean(prolonged_cum)),
  `>=2 episodes, n (%)`             = sprintf("%d (%.1f%%)", sum(multi_epi),
                                          100 * mean(multi_epi)),
  `Short episodes only, n (%)`           = sprintf("%d (%.1f%%)", sum(short_only),
                                          100 * mean(short_only)),
  `Episodes per procedure, mean`          = sprintf("%.2f", mean(n_epi)),
  `Cumulative IOH duration, min, median (IQR)` = sprintf("%.1f (%.1f, %.1f)",
                                          median(MAP65_TIME),
                                          quantile(MAP65_TIME, .25),
                                          quantile(MAP65_TIME, .75)),
  `Cumulative IOH burden (AUC), median (IQR)`     = sprintf("%.1f (%.1f, %.1f)",
                                          median(MAP65_AUC),
                                          quantile(MAP65_AUC, .25),
                                          quantile(MAP65_AUC, .75)),
  `Intraoperative MAP readings, median`       = sprintf("%.0f", median(n_map)),
  `Anesthesia duration, min, median`       = sprintf("%.0f", median(anes_duration_mins)),
  `Monitoring density, readings/min, median`    = sprintf("%.2f", median(rate_map))
), by = Group]

cat("---- Table 2A: procedure-level metrics (denominator = procedures)----\n")
print(by_surg)
cat(sprintf("\nUnit note: %d procedures | %d patients\n\n", nrow(d), uniqueN(d$PATIENT_ID)))


## ============================================================================
## Step 3: compute the event-level descriptive table (Table 2B, denominator = one hypotension episode)
## ============================================================================

de <- merge(epi, d[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group)], by = ID_VARS)

by_epi <- de[, .(
  `Hypotensive episodes, n`            = .N,
  `Procedures involved, n`            = uniqueN(paste(PATIENT_ID, VISIT_ID, OPER_ID)),
  `Duration per episode, min, median (IQR)`  = sprintf("%.1f (%.1f, %.1f)", median(duration),
                                          quantile(duration, .25), quantile(duration, .75)),
  `Episodes >10 min, n (%)`        = sprintf("%d (%.1f%%)", sum(duration > 10),
                                          100 * mean(duration > 10)),
  `Burden per episode (AUC), median (IQR)`      = sprintf("%.1f (%.1f, %.1f)", median(event_auc),
                                          quantile(event_auc, .25), quantile(event_auc, .75)),
  `Maximum depth, mm Hg, median`      = sprintf("%.1f", median(max_map_depth))
), by = Group]

cat("---- Table 2B: episode-level metrics (denominator = one hypotensive episode)----\n")
print(by_epi)


## ============================================================================
## Step 4: save the two tables (csv + Word each)
## ============================================================================

fwrite(by_surg, file.path(OUT, "Table2_SurgeryLevel.csv"))
fwrite(by_epi,  file.path(OUT, "Table2_EpisodeLevel.csv"))

## Save the Word version
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft1 <- flextable(as.data.frame(by_surg)) |> theme_booktabs() |> autofit() |>
         set_caption(caption = "Table 2A. Intraoperative hypotension by treatment group. Unit of analysis: surgical procedure (one anaesthesia episode).")
  print(read_docx() |> body_add_flextable(ft1),
        target = file.path(OUT, "Table2_SurgeryLevel.docx"))
  ft2 <- flextable(as.data.frame(by_epi)) |> theme_booktabs() |> autofit() |>
         set_caption(caption = "Table 2B. Episode-level characteristics of intraoperative hypotension. Unit of analysis: one hypotensive episode.")
  print(read_docx() |> body_add_flextable(ft2),
        target = file.path(OUT, "Table2_EpisodeLevel.docx"))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

