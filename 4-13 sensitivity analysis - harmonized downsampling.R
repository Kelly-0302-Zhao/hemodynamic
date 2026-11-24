## =============================================================================
## 4-13 sensitivity analysis - harmonized downsampling of invasive monitoring
## =============================================================================
## Purpose: downsample invasive-layer MAP readings to 5-min intervals (simulating
##          non-invasive sampling frequency), rebuild the event table, and refit
##          GEE to test whether the invasive-layer effect is driven by detection bias.
## Method: for each invasive-layer procedure, round timestamps to 5-min windows,
##         take the lowest MAP within each 5-min window, and reproduce the event
##         segmentation logic from script 2-2 (GAP_LIMIT = 10 min unchanged).
## =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(geepack)
})

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
DATA <- file.path(D, "JCA_run")
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")
MY_TZ    <- "UTC"
GAP_LIMIT <- 10
DOWNSAMPLE_MIN <- 5

to_chr_ids <- function(DT) {
  DT[, (ID_VARS) := lapply(.SD, as.character), .SDcols = ID_VARS]
  DT
}
parse_iso <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  as.POSIXct(gsub("T", " ", sub("Z$", "", as.character(x))),
             tz = MY_TZ, format = "%Y-%m-%d %H:%M:%S")
}

## ============================================================================
## Step 1: load invasive-layer procedure list + weights + event table (for GEE)
## ============================================================================
cat("[1] Loading invasive-layer procedure list...\n")
layer <- fread(file.path(OUT, "Modality_Layer_For_Subgroup_v4.csv"), encoding = "UTF-8")
layer <- to_chr_ids(layer)
invasive_ops <- layer[layer_final == "invasive", .(PATIENT_ID, VISIT_ID, OPER_ID)]
cat(sprintf("  Invasive-layer procedures: %d\n", nrow(invasive_ops)))

## Surgery table + weights
dw <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
dw <- to_chr_ids(dw)
w  <- setDT(readRDS(file.path(OUT, "IPTW_Weights.rds")))
w  <- to_chr_ids(w)
dw <- dw[w[, .(PATIENT_ID, VISIT_ID, OPER_ID, sw, in_overlap)], on = ID_VARS, nomatch = 0]
dw <- dw[!is.na(prolonged)]

## Invasive-layer surgery table (with weights)
dw_inv <- dw[invasive_ops, on = ID_VARS, nomatch = 0]
cat(sprintf("  Invasive-layer weighted procedures: %d (ciprofol %d / propofol %d)\n",
            nrow(dw_inv),
            sum(dw_inv$Group == "Ciprofol_only"),
            sum(dw_inv$Group == "Propofol_only")))

## Original event table (for comparison)
epi_orig <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))
epi_orig <- to_chr_ids(epi_orig)

## ============================================================================
## Step 2: load monitor_clean, keep invasive-layer procedures and intraoperative window only
## ============================================================================
cat("[2] Loading monitor_clean (72M rows, reading 6 columns only)...\n")
bp <- fread(file.path(DATA, "monitor_clean_72460818obs.csv"), encoding = "UTF-8",
            select = c(ID_VARS, "TIME_POINT", "MAP"))
bp <- to_chr_ids(bp)
setnames(bp, "TIME_POINT", "t_raw")
bp[, t := parse_iso(t_raw)]
bp[, t_raw := NULL]

## Keep invasive-layer procedures only
bp <- bp[invasive_ops, on = ID_VARS, nomatch = 0]
cat(sprintf("  Invasive-layer raw readings: %d (%.1f ten-thousandths)\n", nrow(bp), nrow(bp) / 1e4))

## Time-window filtering (using master anes_start/anes_end)
master <- fread(file.path(DATA, "master_clean_85397obs.csv"), encoding = "UTF-8")
master <- to_chr_ids(master)
master[, anes_start := parse_iso(ANES_START_TIME)]
master[, anes_end   := parse_iso(ANES_END_TIME)]
anes <- master[invasive_ops, .(PATIENT_ID, VISIT_ID, OPER_ID, anes_start, anes_end), on = ID_VARS]

bp <- merge(bp, anes, by = ID_VARS, all.x = TRUE)
bp[, in_window := !is.na(t) & !is.na(anes_start) & t >= anes_start & t <= anes_end]
bp <- bp[in_window == TRUE]
bp[, `:=`(anes_start = NULL, anes_end = NULL, in_window = NULL)]
cat(sprintf("  After intraoperative window filter: %d readings\n", nrow(bp)))

## ============================================================================
## Step 3: downsample -- take lowest MAP per 5-min window
## ============================================================================
cat("[3] Downsampling to 5-min intervals...\n")

## Round timestamp to 5-min windows
bp[, t_5min := as.POSIXct(floor(as.numeric(t) / (DOWNSAMPLE_MIN * 60)) * (DOWNSAMPLE_MIN * 60),
                          origin = "1970-01-01", tz = MY_TZ)]

## Take lowest MAP within each 5-min window (same logic as 2-2 "lowest per minute")
bp[, MAP := {
  v <- MAP[!is.na(MAP)]
  if (length(v)) min(v) else NA_real_
}, by = c(ID_VARS, "t_5min")]

## Deduplicate (one row per procedure per 5-min window)
bp <- unique(bp, by = c(ID_VARS, "t_5min"))
setnames(bp, "t_5min", "t")
cat(sprintf("  After downsampling: %d readings (avg %.1f readings per procedure)\n",
            nrow(bp), nrow(bp) / uniqueN(bp, by = ID_VARS)))

## ============================================================================
## Step 4: reproduce 2-2 event segmentation logic (interval + cover + episode split)
## ============================================================================
cat("[4] Rebuilding event table...\n")

setorder(bp, PATIENT_ID, VISIT_ID, OPER_ID, t)
bp[, idx   := seq_len(.N), by = ID_VARS]
bp[, next_t := shift(t, type = "lead"), by = ID_VARS]
bp[, delta  := as.numeric(difftime(next_t, t, units = "mins"))]

## Typical sampling interval -- should be 5 min after downsampling
bp[, med_delta := as.numeric(median(delta[delta > 0 & delta <= GAP_LIMIT], na.rm = TRUE)),
   by = ID_VARS]
bp[!is.finite(med_delta) | med_delta <= 0, med_delta := DOWNSAMPLE_MIN]

## Single-point coverage duration
bp[, cover := ifelse(is.na(delta) | delta <= 0 | delta > GAP_LIMIT, med_delta, delta)]
bp[cover > GAP_LIMIT, cover := GAP_LIMIT]

cat(sprintf("  Median typical sampling interval after downsampling: %.1f min\n", median(bp$med_delta, na.rm = TRUE)))

## Split episodes
bp[, is_low := as.integer(MAP < 65)]
bp[, depth  := pmax(0, 65 - MAP)]
bp[, step_auc := depth * cover]

bp[, brk := as.integer(idx == 1L | is_low != shift(is_low) | (!is.na(delta) & delta > GAP_LIMIT))]
bp[, seg := cumsum(brk), by = ID_VARS]

## Episode-level table
epi_ds <- bp[is_low == 1L, .(
  start_time = min(t),
  end_time   = max(t),
  duration   = sum(cover, na.rm = TRUE),
  event_auc  = sum(step_auc, na.rm = TRUE),
  n_point    = .N,
  max_map_depth = max(depth, na.rm = TRUE)
), by = c(ID_VARS, "seg")]
epi_ds[, seg := NULL]

cat(sprintf("  Downsampled event table: %d episodes (original invasive-layer %d)\n",
            nrow(epi_ds),
            nrow(epi_orig[invasive_ops, on = ID_VARS, nomatch = 0])))
cat(sprintf("  Episode duration median %.1f min, >10min proportion %.1f%%\n",
            median(epi_ds$duration), 100 * mean(epi_ds$duration > 10)))

## ============================================================================
## Step 5: run GEE on downsampled invasive layer (M1-M3 episode-level)
## ============================================================================
cat("[5] Downsampled invasive-layer GEE...\n")

## Attach weights to event table
ep_ds <- epi_ds[dw_inv[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                            Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
                on = ID_VARS, nomatch = NULL]
ep_ds <- ep_ds[!is.na(Group) & is.finite(w) & w > 0]
ep_ds[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

## M1: prolonged episode OR
ep_ds[, long_epi := as.integer(duration > 10)]
d1 <- ep_ds[complete.cases(ep_ds[, .(long_epi, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
m1 <- geeglm(long_epi ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d1, id = sid_int, weights = d1$w,
             family = binomial(), corstr = "exchangeable")
s1 <- summary(m1)$coefficients
e1 <- s1["GroupCiprofol_only", "Estimate"]; se1 <- s1["GroupCiprofol_only", "Std.err"]

## M2: AUC
ep_ds[, log_eauc := log(event_auc + 0.1)]
d2 <- ep_ds[complete.cases(ep_ds[, .(log_eauc, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
m2 <- geeglm(log_eauc ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d2, id = sid_int, weights = d2$w,
             family = gaussian(), corstr = "exchangeable")
s2 <- summary(m2)$coefficients
e2 <- s2["GroupCiprofol_only", "Estimate"]; se2 <- s2["GroupCiprofol_only", "Std.err"]

## M3: duration
ep_ds[, log_dur := log(duration + 0.1)]
d3 <- ep_ds[complete.cases(ep_ds[, .(log_dur, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
m3 <- geeglm(log_dur ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d3, id = sid_int, weights = d3$w,
             family = gaussian(), corstr = "exchangeable")
s3 <- summary(m3)$coefficients
e3 <- s3["GroupCiprofol_only", "Estimate"]; se3 <- s3["GroupCiprofol_only", "Std.err"]

## M4: episode rate (procedure-level)
adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  vv <- dw_inv[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]

## Downsampled procedure-level episode count
epi_per_surg <- ep_ds[, .(n_epi_ds = .N), by = ID_VARS]
dw_ds <- merge(dw_inv, epi_per_surg, by = ID_VARS, all.x = TRUE)
dw_ds[is.na(n_epi_ds), n_epi_ds := 0L]
dw_ds[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
dw_ds[, lg_dur := log(anes_duration_mins)]
d4 <- dw_ds[complete.cases(dw_ds[, c("n_epi_ds", "Group", adj_ok, "lg_dur", "sw"), with = FALSE])]
m4 <- geeglm(as.formula(paste("n_epi_ds ~ Group +", paste(adj_ok, collapse = " + "), "+ offset(lg_dur)")),
             data = d4, id = sid_int, weights = d4$sw,
             family = poisson(), corstr = "exchangeable")
s4 <- summary(m4)$coefficients
e4 <- s4["GroupCiprofol_only", "Estimate"]; se4 <- s4["GroupCiprofol_only", "Std.err"]

## ============================================================================
## Step 6: summary comparison -- original invasive vs downsampled invasive vs non-invasive
## ============================================================================
cat("\n[6] Summary comparison...\n")

## Original invasive + non-invasive results (from Subgroup_Monitoring_Modality_Results.csv)
sub_grp <- fread(file.path(OUT, "Subgroup_Monitoring_Modality_Results.csv"), encoding = "UTF-8")
get_eff <- function(dt, pop, m_key) {
  r <- dt[population == pop & grepl(m_key, model)]
  if (nrow(r) == 0) return(NA_character_)
  r$effect[1]
}

orig_inv_m1   <- get_eff(sub_grp, "invasive",     "M1")
orig_inv_m2   <- get_eff(sub_grp, "invasive",     "M2")
orig_inv_m3   <- get_eff(sub_grp, "invasive",     "M3")
orig_inv_m4   <- get_eff(sub_grp, "invasive",     "M4")
orig_non_m1   <- get_eff(sub_grp, "non_invasive", "M1")
orig_non_m2   <- get_eff(sub_grp, "non_invasive", "M2")
orig_non_m3   <- get_eff(sub_grp, "non_invasive", "M3")
orig_non_m4   <- get_eff(sub_grp, "non_invasive", "M4")
orig_inv_nepi <- sub_grp[population == "invasive" & grepl("M1", model), n_unit][1]
orig_non_nepi <- sub_grp[population == "non_invasive" & grepl("M1", model), n_unit][1]

## Downsampled results
ds_tab <- data.table(
  population = "Invasive downsampled (5-min)",
  M1_OR  = sprintf("%.2f (%.2f-%.2f)", exp(e1), exp(e1-1.96*se1), exp(e1+1.96*se1)),
  M1_p   = 2 * pnorm(-abs(e1/se1)),
  M2_pct = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                   100*(exp(e2)-1), 100*(exp(e2-1.96*se2)-1), 100*(exp(e2+1.96*se2)-1)),
  M2_p   = 2 * pnorm(-abs(e2/se2)),
  M3_pct = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                   100*(exp(e3)-1), 100*(exp(e3-1.96*se3)-1), 100*(exp(e3+1.96*se3)-1)),
  M3_p   = 2 * pnorm(-abs(e3/se3)),
  M4_RR  = sprintf("%.2f (%.2f-%.2f)", exp(e4), exp(e4-1.96*se4), exp(e4+1.96*se4)),
  M4_p   = 2 * pnorm(-abs(e4/se4)),
  n_episodes = nrow(ep_ds),
  n_procedures = nrow(dw_ds))

cat("\n==== Downsampled invasive-layer results ====\n")
print(ds_tab)

## Comparison table
cmp <- data.table(
  Metric = c("M1 OR", "M2 AUC %", "M3 Duration %", "M4 RR", "n_episodes"),
  `Original Invasive (1-min)` = c(orig_inv_m1, orig_inv_m2, orig_inv_m3,
                                  orig_inv_m4, as.character(orig_inv_nepi)),
  `Invasive Downsampled (5-min)` = c(ds_tab$M1_OR, ds_tab$M2_pct, ds_tab$M3_pct,
                                      ds_tab$M4_RR, as.character(ds_tab$n_episodes)),
  `Original Non-invasive (5-min)` = c(orig_non_m1, orig_non_m2, orig_non_m3,
                                      orig_non_m4, as.character(orig_non_nepi)))

cat("\n==== Three-layer comparison ====\n")
print(cmp)

## Write outputs
fwrite(ds_tab, file.path(OUT, "Harmonized_Downsampling_Invasive_5min.csv"))
fwrite(cmp,   file.path(OUT, "Harmonized_Downsampling_Comparison.csv"))

cat("\n==== Verdict ====\n")
m1_orig <- as.numeric(sub(" .*", "", orig_inv_m1))
m1_ds   <- exp(e1)
m1_non  <- as.numeric(sub(" .*", "", orig_non_m1))
cat(sprintf("  M1 OR: original invasive(1min) %.2f -> downsampled(5min) %.2f -> non-invasive(5min) %.2f\n", m1_orig, m1_ds, m1_non))
if (abs(m1_ds - m1_non) < abs(m1_ds - m1_orig)) {
  cat("  -> Downsampled result approaches non-invasive layer: detection bias is likely a driver\n")
} else {
  cat("  -> Downsampled result remains close to original invasive layer: not driven by detection bias\n")
}
