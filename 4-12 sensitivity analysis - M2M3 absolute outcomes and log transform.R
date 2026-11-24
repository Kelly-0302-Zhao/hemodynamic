## =============================================================================
## 4-12 sensitivity analysis - M2/M3 absolute outcomes + log transform without +0.1
## =============================================================================
##   1. Report weighted geometric mean AUC, arithmetic mean AUC, median duration,
##      mean duration + bootstrap CI for both groups
##      (reviewer request: report absolute standardized outcomes)
##   2. Compare log(AUC+0.1) vs log(AUC) for M2/M3 (verify whether +0.1 affects
##      conclusions; AUC min = 0.20, no zeros, so log(AUC) is safe)
## =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(geepack)
})

## Root directory
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## ---- Load event table + weights (same as 3-6 steps 2-3) ----
d   <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
w   <- setDT(readRDS(file.path(OUT, "IPTW_Weights.rds")))
epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

for (v in ID_VARS) {
  set(d,   j = v, value = as.character(d[[v]]))
  set(w,   j = v, value = as.character(w[[v]]))
  set(epi, j = v, value = as.character(epi[[v]]))
}

dw <- d[w[, .(PATIENT_ID, VISIT_ID, OPER_ID, ps_hat, sw, in_overlap)],
        on = ID_VARS, nomatch = 0]
dw <- dw[!is.na(prolonged)]

## Attach weights to event table
ep <- epi[dw[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                 Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
          on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]
ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

cat("Analysis population:", nrow(ep), "episodes /", uniqueN(ep$PATIENT_ID), "patients\n\n")

## ============================================================================
## Part 1: weighted absolute outcomes for both groups + bootstrap CI
## ============================================================================

## bootstrap by procedure: resample procedures (same as 3-6 step 9)
boot_abs <- function(epi_dt, metric, B = 2000, seed = 20260928) {
  ## metric: "geomean_auc", "arithmean_auc", "median_dur", "meandur"
  sid <- epi_dt$sid_int
  np  <- max(sid)
  isc <- as.integer(epi_dt$Group == "Ciprofol_only")
  wt  <- epi_dt$w
  auc <- epi_dt$event_auc
  dur <- epi_dt$duration

  ## Aggregate by procedure (per-procedure total weight + weighted statistic)
  dt <- data.table(sid, isc, wt, auc, dur)

  calc_stat <- function(sub, grp) {
    s <- sub[isc == grp]
    w <- s$wt
    switch(metric,
      "geomean_auc"   = exp(weighted.mean(log(s$auc), w)),
      "arithmean_auc" = weighted.mean(s$auc, w),
      "median_dur"    = weightedMedian(s$dur, w),
      "meandur"       = weighted.mean(s$dur, w),
      NA_real_)
  }

  ## True sample
  cip_true <- calc_stat(dt, 1)
  pro_true <- calc_stat(dt, 0)

  ## bootstrap
  set.seed(seed)
  cip_reps <- numeric(B)
  pro_reps <- numeric(B)
  for (b in seq_len(B)) {
    mult <- tabulate(sample.int(np, np, replace = TRUE), np)
    idx  <- which(mult[sid] > 0)
    sub  <- dt[idx]
    ## Multiply each procedure's weight by its resample count
    sub[, wmult := wt * mult[sid]]
    ## Temporarily replace weight
    orig_w <- sub$wt
    sub$wt <- sub$wmult
    cip_reps[b] <- tryCatch(calc_stat(sub, 1), error = function(e) NA_real_)
    pro_reps[b] <- tryCatch(calc_stat(sub, 0), error = function(e) NA_real_)
    sub$wt <- orig_w
  }
  cip_reps <- cip_reps[is.finite(cip_reps)]
  pro_reps <- pro_reps[is.finite(pro_reps)]

  list(
    cip_est = cip_true,
    cip_lo  = as.numeric(quantile(cip_reps, .025)),
    cip_hi  = as.numeric(quantile(cip_reps, .975)),
    pro_est = pro_true,
    pro_lo  = as.numeric(quantile(pro_reps, .025)),
    pro_hi  = as.numeric(quantile(pro_reps, .975))
  )
}

## weightedMedian (dependency-free implementation)
weightedMedian <- function(x, w) {
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  idx <- which(cw >= 0.5)[1]
  if (idx == 1) return(x[1])
  ## Linear interpolation
  w_lo <- cw[idx - 1]
  w_hi <- cw[idx]
  if (w_lo == w_hi) return(x[idx])
  x[idx - 1] + (0.5 - w_lo) / (w_hi - w_lo) * (x[idx] - x[idx - 1])
}

cat("==== Part 1: weighted absolute outcomes (bootstrap by procedure, B=2000) ====\n\n")

metrics <- list(
  list(name = "Geometric mean AUC",          key = "geomean_auc",   fmt = "%.2f"),
  list(name = "Arithmetic mean AUC",         key = "arithmean_auc", fmt = "%.2f"),
  list(name = "Median duration (min)",       key = "median_dur",    fmt = "%.1f"),
  list(name = "Mean duration (min)",         key = "meandur",       fmt = "%.2f"))

abs_out <- list()
for (m in metrics) {
  r <- boot_abs(ep, m$key)
  abs_out[[length(abs_out) + 1]] <- data.table(
    Metric = m$name,
    Ciprofol = sprintf(m$fmt, r$cip_est),
    Ciprofol_95CI = sprintf("(%.2f to %.2f)", r$cip_lo, r$cip_hi),
    Propofol = sprintf(m$fmt, r$pro_est),
    Propofol_95CI = sprintf("(%.2f to %.2f)", r$pro_lo, r$pro_hi),
    Ratio_or_diff = if (m$key %in% c("geomean_auc", "arithmean_auc")) {
      sprintf("%.2f%%", 100 * (r$cip_est / r$pro_est - 1))
    } else {
      sprintf("%+.2f min", r$cip_est - r$pro_est)
    })
  cat(sprintf("  %s: Ciprofol %s %s | Propofol %s %s | %s\n",
              m$name,
              sprintf(m$fmt, r$cip_est), sprintf("(%.2f to %.2f)", r$cip_lo, r$cip_hi),
              sprintf(m$fmt, r$pro_est), sprintf("(%.2f to %.2f)", r$pro_lo, r$pro_hi),
              if (m$key %in% c("geomean_auc", "arithmean_auc"))
                sprintf("%.2f%%", 100 * (r$cip_est / r$pro_est - 1))
              else
                sprintf("%+.2f min", r$cip_est - r$pro_est)))
}

abs_tab <- rbindlist(abs_out)
fwrite(abs_tab, file.path(OUT, "M2M3_Absolute_Outcomes.csv"))

## ============================================================================
## Part 2: compare log(AUC+0.1) vs log(AUC) for M2/M3
## ============================================================================

cat("\n==== Part 2: log(AUC+0.1) vs log(AUC) comparison ====\n\n")

## log(AUC+0.1) -- original transform (same as 3-6)
ep[, log_eauc_01 := log(event_auc + 0.1)]
d2a <- ep[complete.cases(ep[, .(log_eauc_01, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
m2a <- geeglm(log_eauc_01 ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
              data = d2a, id = sid_int, weights = d2a$w,
              family = gaussian(), corstr = "exchangeable")
s2a <- summary(m2a)$coefficients
e2a <- s2a["GroupCiprofol_only", "Estimate"]
se2a <- s2a["GroupCiprofol_only", "Std.err"]

## log(AUC) -- without +0.1
ep[, log_eauc_raw := log(event_auc)]
d2b <- ep[complete.cases(ep[, .(log_eauc_raw, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
m2b <- geeglm(log_eauc_raw ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
              data = d2b, id = sid_int, weights = d2b$w,
              family = gaussian(), corstr = "exchangeable")
s2b <- summary(m2b)$coefficients
e2b <- s2b["GroupCiprofol_only", "Estimate"]
se2b <- s2b["GroupCiprofol_only", "Std.err"]

## log(duration+0.1) -- original transform
ep[, log_dur_01 := log(duration + 0.1)]
d3a <- ep[complete.cases(ep[, .(log_dur_01, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
m3a <- geeglm(log_dur_01 ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
              data = d3a, id = sid_int, weights = d3a$w,
              family = gaussian(), corstr = "exchangeable")
s3a <- summary(m3a)$coefficients
e3a <- s3a["GroupCiprofol_only", "Estimate"]
se3a <- s3a["GroupCiprofol_only", "Std.err"]

## log(duration) -- without +0.1 (duration min = 1.00, safe)
ep[, log_dur_raw := log(duration)]
d3b <- ep[complete.cases(ep[, .(log_dur_raw, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
m3b <- geeglm(log_dur_raw ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
              data = d3b, id = sid_int, weights = d3b$w,
              family = gaussian(), corstr = "exchangeable")
s3b <- summary(m3b)$coefficients
e3b <- s3b["GroupCiprofol_only", "Estimate"]
se3b <- s3b["GroupCiprofol_only", "Std.err"]

## Comparison table
cmp_tab <- data.table(
  Model = c("M2 AUC (log(x+0.1))", "M2 AUC (log(x))",
            "M3 Duration (log(x+0.1))", "M3 Duration (log(x))"),
  Transform = c("log(AUC+0.1)", "log(AUC)", "log(dur+0.1)", "log(dur)"),
  Estimate = c(e2a, e2b, e3a, e3b),
  SE = c(se2a, se2b, se3a, se3b),
  Pct_change = c(100*(exp(e2a)-1), 100*(exp(e2b)-1), 100*(exp(e3a)-1), 100*(exp(e3b)-1)),
  Pct_lo = c(100*(exp(e2a-1.96*se2a)-1), 100*(exp(e2b-1.96*se2b)-1),
             100*(exp(e3a-1.96*se3a)-1), 100*(exp(e3b-1.96*se3b)-1)),
  Pct_hi = c(100*(exp(e2a+1.96*se2a)-1), 100*(exp(e2b+1.96*se2b)-1),
             100*(exp(e3a+1.96*se3a)-1), 100*(exp(e3b+1.96*se3b)-1)))
cmp_tab[, Effect := sprintf("%+.1f%% (%+.1f%% to %+.1f%%)", Pct_change, Pct_lo, Pct_hi)]
cmp_tab[, P := 2 * pnorm(-abs(Estimate / SE))]

print(cmp_tab[, .(Model, Transform, Effect, P = sprintf("%.3g", P))])
fwrite(cmp_tab, file.path(OUT, "M2M3_LogTransform_Comparison.csv"))

cat("\nComparison summary:\n")
cat(sprintf("  M2 AUC:  log(AUC+0.1) = %+.1f%%  vs  log(AUC) = %+.1f%%  -> difference %.2f percentage points\n",
            100*(exp(e2a)-1), 100*(exp(e2b)-1), abs(100*(exp(e2a)-1) - 100*(exp(e2b)-1))))
cat(sprintf("  M3 Dur:  log(dur+0.1) = %+.1f%%  vs  log(dur) = %+.1f%%  -> difference %.2f percentage points\n",
            100*(exp(e3a)-1), 100*(exp(e3b)-1), abs(100*(exp(e3a)-1) - 100*(exp(e3b)-1))))
