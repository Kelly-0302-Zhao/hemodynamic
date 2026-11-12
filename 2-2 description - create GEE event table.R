## =============================================================================
## 2-2 description - create GEE event table.R - rebuild the intraoperative hypotension event table (event-level + surgery-level)
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

## Three ID columns. The ID column types are inconsistent across the three tables, so convert them all to character before merging
ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
MY_TZ    <- "UTC"      # Raw timestamps carry a Z suffix; parse everything as UTC throughout to keep the internal representation consistent
GAP_LIMIT <- 10        # Minutes: a gap of more than 10 minutes between consecutive measurements is treated as a "monitoring gap"; hypotension duration is not accumulated across a gap


## Convert ISO time strings such as "2024-01-01T08:30:00Z" to POSIXct
parse_iso <- function(x) {
  ## Return directly if the value is already POSIXct, without passing through as.character():
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
## Step 1: read the two base tables (surgery table + intraoperative monitoring table)
## ============================================================================

## Surgery table: one row = one surgery, already containing group, covariates, anesthesia start/end time, and anesthesia duration
master <- fread(file.path(DATA, "master_clean_85397obs.csv"), encoding = "UTF-8") 
master <- to_chr_ids(master)
master[, anes_start := parse_iso(ANES_START_TIME)]
master[, anes_end   := parse_iso(ANES_END_TIME)]
master[, year       := as.integer(format(as.Date(anes_start), "%Y"))]
master[, Group      := relevel(factor(as.character(Group)), ref = "Propofol_only")]
for (v in c("ASA_GRADE", "operation_type", "Hospital", "SEX")) {
  if (v %in% names(master)) master[[v]] <- factor(master[[v]])
}
cat("Surgery table:", nrow(master), "procedures,", uniqueN(master$PATIENT_ID), "patients\n")

## Intraoperative monitoring table: about 72.46 million raw records (multiple records within the same minute are all kept, not deduplicated). Read only the needed columns (three ID columns + time + MAP)
bp <- fread(file.path(DATA, "monitor_clean_72460818obs.csv"), encoding = "UTF-8",
            select = c(ID_VARS, "TIME_POINT", "MAP"))
bp <- to_chr_ids(bp)
setnames(bp, "TIME_POINT", "t_raw")
bp[, t := parse_iso(t_raw)]
bp[, t_raw := NULL]
cat("Monitor table:", nrow(bp), "measurement points\n\n")


## ============================================================================
## Step 2: run a time benchmark self-check first, then apply the intraoperative time window filter
## ============================================================================

## Attach each surgery's anesthesia start/end time to every measurement
bp <- master[, .(PATIENT_ID, VISIT_ID, OPER_ID, anes_start, anes_end)][bp, on = ID_VARS]
n_all <- nrow(bp)
bp[, in_window    := t >= anes_start & t <= anes_end]
bp[, pos_min      := as.numeric(difftime(t, anes_start, units = "mins"))]
bp[, before_start := t < anes_start]
bp[, after_end    := t > anes_end]

cat("---- Time reference self-check (position of measurement points relative to anesthesia start)----\n")
cat(sprintf("Missing ANES_START_TIME in master: %d procedures; missing ANES_END_TIME: %d\n",
            master[is.na(anes_start), .N], master[is.na(anes_end), .N]))
cat(sprintf("Measurement position: before anesthesia start %.2f%% | within the intraoperative window %.2f%% | after anesthesia end %.2f%%\n",
            100 * mean(bp$before_start, na.rm = TRUE),
            100 * mean(bp$in_window,    na.rm = TRUE),
            100 * mean(bp$after_end,    na.rm = TRUE)))
cat("Quantiles of measurement point position relative to anesthesia start (minutes; negative = before anesthesia):\n")
print(round(quantile(bp$pos_min, c(.001, .01, .05, .25, .50, .75, .95, .99, .999),
                     na.rm = TRUE), 1))

## Check surgery by surgery whether there are any measurements inside the intraoperative window
ov <- bp[!is.na(anes_start) & !is.na(anes_end),
         .(n_point = .N, n_in = sum(in_window, na.rm = TRUE)), by = ID_VARS]
pct_cov <- 100 * ov[n_in > 0, .N] / nrow(ov)
cat(sprintf("Among procedures with anesthesia start and end times, proportion with at least 1 measurement point in the intraoperative window: %.2f%%\n", pct_cov))

if (!is.finite(pct_cov) || pct_cov < 50) {
  stop("Almost no measurement points in the intraoperative window - the time reference of monitor and master is probably inconsistent (time zone or unit difference).\n",
       "       Please check the time column definitions of both tables first; do not apply the intraoperative window filter directly.")
}

bp[, bin_lo := pmin(pmax(floor(pos_min / 10) * 10, -180), 480)]
chk <- bp[!is.na(bin_lo), .(n_point = .N), by = bin_lo][order(bin_lo)]
chk[, pct := round(100 * n_point / sum(n_point), 3)]
fwrite(chk, file.path(OUT, "Check_Window_Position.csv"))

cat(sprintf("\nTime window filter: %d → %d measurement points (%.1f%% removed)\n",
            n_all, sum(bp$in_window, na.rm = TRUE),
            100 * (1 - sum(bp$in_window, na.rm = TRUE) / n_all)))
bp <- bp[in_window == TRUE]
bp[, `:=`(in_window = NULL, pos_min = NULL, before_start = NULL,
          after_end = NULL, bin_lo = NULL)]


## ============================================================================
## Step 3: multiple records within the same minute -> collapse into "one measurement per minute" (take the lowest MAP within that minute)
## ============================================================================
bp[, MAP := {
  v <- MAP[!is.na(MAP)]
  if (length(v)) min(v) else NA_real_
}, by = c(ID_VARS, "t")]
bp[, n_raw := .N, by = c(ID_VARS, "t")]

n_before_agg <- nrow(bp)
bp <- unique(bp, by = c(ID_VARS, "t"))
cat(sprintf("Collapsing multiple records within the same minute: %d → %d measurement points (%.2f raw records per point on average); %d with missing MAP\n\n",
            n_before_agg, nrow(bp), n_before_agg / nrow(bp), sum(is.na(bp$MAP))))


## ============================================================================
## Step 4: compute the sampling interval and assign a "coverage duration" to each measurement
## ============================================================================

setorder(bp, PATIENT_ID, VISIT_ID, OPER_ID, t)
bp[, idx   := seq_len(.N), by = ID_VARS]
bp[, next_t := shift(t, type = "lead"), by = ID_VARS]
bp[, delta  := as.numeric(difftime(next_t, t, units = "mins"))]

## The "typical sampling interval" of that surgery (count only the points whose interval is within the normal range 0 < interval <= 10 minutes)
bp[, med_delta := as.numeric(median(delta[delta > 0 & delta <= GAP_LIMIT], na.rm = TRUE)),
   by = ID_VARS]
bp[!is.finite(med_delta) | med_delta <= 0, med_delta := 5]

## Single-point coverage duration: use the true interval for normal intervals; for missing or above-threshold intervals, use the typical interval only
bp[, cover := ifelse(is.na(delta) | delta <= 0 | delta > GAP_LIMIT, med_delta, delta)]
bp[cover > GAP_LIMIT, cover := GAP_LIMIT]      # Fallback cap: any single point represents at most 10 minutes

cat(sprintf("Median typical sampling interval: %.1f minutes\n", median(bp$med_delta, na.rm = TRUE)))
cat(sprintf("Proportion of measurement intervals >%d minutes (classified as monitoring interruption, which breaks episodes): %.3f%%; %d measurement points with missing MAP\n",
            GAP_LIMIT, 100 * mean(bp$delta > GAP_LIMIT, na.rm = TRUE), sum(is.na(bp$MAP))))

## ============================================================================
## Step 5: split into "continuous segments" and cut out each hypotension episode (event-level table)
##   Event definition: consecutive measurements with MAP<65 count as one episode (a single point also counts as one episode); a new episode starts only when MAP recovers to >=65 or when the interval to the adjacent measurement is >10 minutes (monitoring interruption).
## ============================================================================

bp[, is_low := as.integer(MAP < 65)]
bp[, depth  := pmax(0, 65 - MAP)]
bp[, step_auc := depth * cover]

## Conditions for starting a new segment: the first point of this surgery / a change in the hypotension state / a true monitoring interruption (interval > 10 minutes)
bp[, brk := as.integer(idx == 1L | is_low != shift(is_low) | (!is.na(delta) & delta > GAP_LIMIT))]
bp[, seg := cumsum(brk), by = ID_VARS]

## Event-level: one row = one hypotension episode
epi_episodes <- bp[is_low == 1L, .(
  start_time = min(t),
  end_time   = max(t),
  duration   = sum(cover, na.rm = TRUE),        # Sum the per-point coverage durations (not across interruptions)
  event_auc  = sum(step_auc, na.rm = TRUE),
  n_point    = .N,
  max_map_depth = max(depth, na.rm = TRUE)
), by = c(ID_VARS, "seg")]
epi_episodes[, seg := NULL]

cat("\nEpisode-level table (one row = one hypotensive episode):", nrow(epi_episodes), "episodes\n")
cat(sprintf("  Episode duration: median %.1f minutes, share >10 minutes: %.1f%%\n",
            median(epi_episodes$duration), 100 * mean(epi_episodes$duration > 10)))
cat(sprintf("  Single-point episodes with only 1 measurement point: %d (%.1f%%) - each still counts as 1 episode\n",
            sum(epi_episodes$n_point == 1L), 100 * mean(epi_episodes$n_point == 1L)))


## ============================================================================
## Step 6: aggregate the surgery-level table
## ============================================================================

surg <- bp[, .(
  n_epi        = uniqueN(seg[is_low %in% 1L]),                 # ★ Number of episodes (outcome for M4)
  MAP65_TIME   = sum(cover * is_low, na.rm = TRUE),            # Minutes accumulated below 65
  MAP65_AUC    = sum(step_auc, na.rm = TRUE),                  # Cumulative hypotension area (burden)
  n_map        = .N,                                           # Number of valid intraoperative measurements
  first_low_t  = suppressWarnings(min(t[is_low == 1L]))
), by = ID_VARS]

## Add back the surgeries with "no hypotension at all" (episode count = 0)
surg <- merge(master[, .(PATIENT_ID, VISIT_ID, OPER_ID, anes_start, anes_end,
                         anes_duration_mins)], surg, by = ID_VARS, all.x = TRUE)

surg[, no_monitor := is.na(n_map)]

surg[is.na(n_epi) & !is.na(anes_start), `:=`(n_epi = 0L, MAP65_TIME = 0,
                                             MAP65_AUC = 0, n_map = 0L)]
surg[, rate_map := as.numeric(NA)]
surg[!no_monitor & !is.na(anes_duration_mins),
     rate_map := n_map / anes_duration_mins]

## Outcome variables (all defined at the level of "one surgery")
surg[, n_long     := NA_integer_]
surg[, prolonged  := NA_integer_]     # Whether that surgery had hypotension lasting >10 minutes
surg[, prolonged_cum := NA_integer_]  # Whether the cumulative time below 65 in that surgery was >10 minutes
surg[, multi_epi  := NA_integer_]     # Whether that surgery had >=2 hypotensive episodes
surg[, short_only := NA_integer_]     # That surgery had only short episodes and no long episode

agg <- epi_episodes[, .(n_long = sum(duration > 10, na.rm = TRUE)), by = ID_VARS]
surg <- merge(surg, agg, by = ID_VARS, all.x = TRUE, suffixes = c("", ".epi"))
surg[, n_long := as.integer(fcoalesce(n_long.epi, 0L))][, n_long.epi := NULL]

surg[!is.na(n_epi), `:=`(
  prolonged     = as.integer(n_long >= 1L),
  prolonged_cum = as.integer(MAP65_TIME > 10),
  multi_epi     = as.integer(n_epi >= 2L),
  short_only    = as.integer(n_epi >= 1L & n_long == 0L))]

## Surgeries with a missing anesthesia time window or with no measurements inside the intraoperative window: the outcome cannot be defined, set to NA (all subsequent models automatically exclude these surgeries and report the counts faithfully)
bad <- surg[is.na(anes_start) | is.na(anes_end) | no_monitor]
if (nrow(bad) > 0) {
  cat(sprintf("\n[Note] %d procedures cannot have an intraoperative outcome defined (missing anesthesia start/end times: %d + no measurement points intraoperatively: %d):\n",
              nrow(bad), surg[is.na(anes_start) | is.na(anes_end), .N], surg[no_monitor == TRUE, .N]))
  print(master[bad[, .(PATIENT_ID, VISIT_ID, OPER_ID)], on = ID_VARS,
               nomatch = NULL][, .N, by = .(Hospital, Group)][order(-N)])
}
surg[is.na(anes_start) | is.na(anes_end) | no_monitor, `:=`(
  n_epi = NA_integer_, n_long = NA_integer_, prolonged = NA_integer_,
  prolonged_cum = NA_integer_, multi_epi = NA_integer_, short_only = NA_integer_)]


## ============================================================================
## Step 7: cross-check at the surgery level against the old pipeline (data_ioh.csv) to confirm the numbers are reasonable
## ============================================================================

old_path <- file.path(D, "data_ioh.csv")
if (file.exists(old_path)) {
  old <- fread(old_path, encoding = "UTF-8")
  old <- to_chr_ids(old)
  cmp <- merge(surg[, .(PATIENT_ID, VISIT_ID, OPER_ID, new_time = MAP65_TIME,
                        new_auc = MAP65_AUC, new_n = n_epi)],
               old[, .(PATIENT_ID, VISIT_ID, OPER_ID, old_time = MAP65_TIME,
                       old_auc = MAP65_AUC, old_ioh = MAP65_IOH)],
               by = ID_VARS, all.x = TRUE)
  cat("\n---- Old vs new specification comparison (procedure level, medians)----\n")
  cat(sprintf("  Cumulative duration below 65: old %.1f minutes → new %.1f minutes\n",
              median(cmp$old_time, na.rm = TRUE), median(cmp$new_time, na.rm = TRUE)))
  cat(sprintf("  Cumulative hypotension AUC  : old %.1f → new %.1f\n",
              median(cmp$old_auc, na.rm = TRUE), median(cmp$new_auc, na.rm = TRUE)))
  cat(sprintf("  Procedures with IOH (old specification): %d (%.1f%%)\n",
              sum(cmp$old_ioh == 1, na.rm = TRUE), 100 * mean(cmp$old_ioh == 1, na.rm = TRUE)))
  cat(sprintf("  Procedures with IOH (new specification): %d (%.1f%%)\n",
              sum(cmp$new_n > 0, na.rm = TRUE), 100 * mean(cmp$new_n > 0, na.rm = TRUE)))
  fwrite(cmp, file.path(OUT, "Check_EventTable_Old_vs_New.csv"))
} else {
  cat("\n[Note] Old pipeline comparison table not found:", old_path, "\n")
  cat("       The most recent old-pipeline data_ioh.csv is under 20260623/; its cohort snapshot differs from this run,\n")
  cat("       so no forced comparison is made here and the old vs new specification comparison is skipped (this does not affect any main result).\n")
  cat("       For this QC, copy data_ioh.csv from the matching snapshot to ", D, " and rerun.\n")
}

cat("\n=========== Procedure-level outcome distribution (one row = one procedure) ===========\n")
print(surg[!is.na(prolonged), .(n_surgery = .N,
                                IOH率 = sprintf("%.1f%%", 100 * mean(n_epi > 0)),
                                `Prolonged IOH rate >10min` = sprintf("%.1f%%", 100 * mean(prolonged == 1)),
                                `Rate of ≥2 episodes` = sprintf("%.1f%%", 100 * mean(multi_epi == 1)),
                                `Episodes/procedure` = sprintf("%.2f", mean(n_epi)),
                                `Median cumulative duration` = sprintf("%.1f", median(MAP65_TIME)))])
cat("\nPatient counts also given:", nrow(surg[!is.na(prolonged)]), "procedures,",
    uniqueN(surg[!is.na(prolonged)]$PATIENT_ID), "patients\n")


## ============================================================================
## Step 8: save the two main tables
## ============================================================================

saveRDS(epi_episodes, file.path(OUT, "Final_Analysis_Events.rds"))
saveRDS(surg,         file.path(OUT, "Final_Analysis_Surgery.rds"))

