## =============================================================================
## 4-5 subgroup analysis - age 65-75.R
## -----------------------------------------------------------------------------
## 4-5 Subgroup analysis: procedures in patients aged 65-74 only
## =============================================================================


## ============================================================================
## Step 1: preparation
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # read and process tables
  library(geepack)      # run GEE models
})

LEVEL_EN <- c("女" = "Female", "男" = "Male",
              "Ⅰ" = "I", "Ⅱ" = "II", "Ⅲ" = "III")
lvl_en <- function(x) {
  out <- as.character(x)
  hit <- LEVEL_EN[out]
  out[!is.na(hit)] <- hit[!is.na(hit)]
  out
}


## root directory, output directory
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")                 # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

## Three ID columns. The ID column types differ across the three tables, so convert them all to character before use
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## The 10 baseline variables used for the propensity score (same as the whole-cohort main analysis; PREOPERATIVE_MAP excluded since 2026-09-30)
PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR",
             "operation_type", "Hypertension", "Diabetes", "Coronary_artery_disease",
             "Heart_failure")

ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## Label of this subgroup
LABEL <- "AGE 65-75"

## Read two tables
tab     <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))  # Surgery table
epi_all <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))       # Event table

## Convert ID columns to character
for (v in ID_VARS) {
  set(tab,     j = v, value = as.character(tab[[v]]))
  set(epi_all, j = v, value = as.character(epi_all[[v]]))
}

cat("Loaded: surgery table", nrow(tab), "procedures | event table", nrow(epi_all), "episodes\n")


## ============================================================================
## Step 2: extract the procedures in this subgroup, run an unweighted (crude) comparison, and check whether the sample size is sufficient
## ============================================================================

## Keep only procedures in patients aged 65-74
x <- tab[AGE >= 65 & AGE < 75]
cat("Subgroup included: surgery", nrow(x), "procedures |", uniqueN(x$PATIENT_ID), "patients\n")

if (nrow(x) < 200 || uniqueN(x$Group) < 2) {
  cat("  Too few procedures or only one arm; skipped.\n")
} else {

  cat(sprintf("  Group: Ciprofol %d procedures / Propofol %d procedures\n",
              sum(x$Group == "Ciprofol_only"), sum(x$Group == "Propofol_only")))

  ## Unweighted comparison: raw differences between the two groups on a few summary metrics (no weighting, for reference only)
  crude <- x[, .(n_surgery    = .N,
                 n_patient    = uniqueN(PATIENT_ID),
                 prolonged_pct = 100 * mean(prolonged, na.rm = TRUE),
                 MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
                 n_epi_mean   = mean(n_epi, na.rm = TRUE)),
              by = Group]
  cat("\n  ---- Step 2: unweighted (crude) comparison ----\n")
  print(crude)


  ## ==========================================================================
  ## Step 3: re-fit the propensity score within this subgroup and compute the stabilized weights
  ## ==========================================================================
  xw <- copy(x)
  xw[, treat := as.integer(Group == "Ciprofol_only")]

  use <- character(0)
  for (v in PS_VARS) {
    if (!v %in% names(xw)) next
    vv <- xw[[v]]
    if (is.numeric(vv)) {
      if (is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0) use <- c(use, v)
    } else {
      if (uniqueN(as.character(vv[!is.na(vv)])) > 1) use <- c(use, v)
    }
  }
  dropped <- setdiff(PS_VARS, use)
  if (length(dropped))
    cat("   [PS] constant within this subgroup, dropped automatically:", paste(dropped, collapse = ", "), "\n")

  ## Fit only among the patients with all 11 variables complete
  cc  <- complete.cases(xw[, ..use])
  xp  <- xw[cc]
  if (nrow(xp) < 100 || min(table(xp$treat)) < 10)
    stop("Too few procedures in this subgroup, or one arm is too small; propensity score cannot be fitted")

  ## Logistic regression for the propensity score, then compute the stabilized weights (treatment group = treatment rate/PS, control group = (1-treatment rate)/(1-PS))
  ps <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
            data = xp, family = binomial())
  xp[, ps_hat := predict(ps, type = "response")]
  prev <- mean(xp$treat)
  xp[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

  xp[, `:=`(ql = as.numeric(quantile(sw_raw, 0.01)),
            qh = as.numeric(quantile(sw_raw, 0.99))), by = Group]
  xp[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
  xp[, c("ql", "qh") := NULL]

  lo_overlap <- max(min(xp$ps_hat[xp$treat == 1]), min(xp$ps_hat[xp$treat == 0]))
  hi_overlap <- min(max(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 0]))
  xp[, in_overlap := ps_hat >= lo_overlap & ps_hat <= hi_overlap]

  ## Covariate balance before and after weighting (SMD = between-group mean difference / pooled SD; <0.1 counts as balanced)
  wsmd <- function(v, w) {
    v <- as.numeric(v)
    m1 <- weighted.mean(v[xp$treat == 1], w[xp$treat == 1])
    m0 <- weighted.mean(v[xp$treat == 0], w[xp$treat == 0])
    v1 <- sum(w[xp$treat == 1] * (v[xp$treat == 1] - m1)^2) / sum(w[xp$treat == 1])
    v0 <- sum(w[xp$treat == 0] * (v[xp$treat == 0] - m0)^2) / sum(w[xp$treat == 0])
    sp <- sqrt((v1 + v0) / 2)
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    abs(m1 - m0) / sp
  }
  bal <- rbindlist(lapply(use, function(v) {
    vv <- xp[[v]]
    if (is.numeric(vv)) {
      data.table(variable   = v,
                 smd_before = wsmd(vv, rep(1, nrow(xp))),
                 smd_after  = wsmd(vv, xp$sw))
    } else {
      rbindlist(lapply(setdiff(unique(as.character(vv)), NA), function(l) {
        ind <- as.integer(as.character(vv) == l)
        data.table(variable   = paste0(v, ": ", lvl_en(l)),
                   smd_before = wsmd(ind, rep(1, nrow(xp))),
                   smd_after  = wsmd(ind, xp$sw))
      }))
    }
  }))

  ## Effective sample size (Kish): when the weights are imbalanced, judge how much real information is left behind the inflated weights
  ess <- function(w) sum(w)^2 / sum(w^2)
  ess_tab <- xp[, .(Group       = as.character(first(Group)),
                    n_procedures = .N,
                    n_patients   = uniqueN(PATIENT_ID),
                    weight_median = median(sw),
                    weight_max = max(sw),
                    ESS      = ess(sw)), by = Group][, Group := NULL]

  cat(sprintf("   [PS] Procedures included in weighting %d (patients %d); PS range Ciprofol %.3f-%.3f, Propofol %.3f-%.3f\n",
              nrow(xp), uniqueN(xp$PATIENT_ID),
              min(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 1]),
              min(xp$ps_hat[xp$treat == 0]), max(xp$ps_hat[xp$treat == 0])))
  cat(sprintf("   [Balance] Max SMD before weighting %.3f → max SMD after weighting %.3f\n",
              max(bal$smd_before, na.rm = TRUE), max(bal$smd_after, na.rm = TRUE)))
  e_c <- ess_tab[Group == "Ciprofol_only", ESS]
  e_p <- ess_tab[Group == "Propofol_only", ESS]
  cat(sprintf("   [Weight] Ciprofol ESS %.0f / %d procedures, Propofol ESS %.0f / %d procedures\n",
              e_c, sum(xp$treat == 1), e_p, sum(xp$treat == 0)))
  cat(sprintf("   [Overlap] PS overlap region between groups = [%.3f, %.3f]; outside the region %d procedures (%.1f%%)\n",
              lo_overlap, hi_overlap, sum(!xp$in_overlap), 100 * mean(!xp$in_overlap)))


  ## ==========================================================================
  ## Step 4: prepare event-level data + model M1
  ## ==========================================================================

  ## Attach the group, weight, and adjustment variables from the surgery table to each event → ep: one row = one episode
  ep <- epi_all[xp[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                       Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
                on = ID_VARS, nomatch = NULL]
  ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

  ## GEE requires the cluster id to be an integer; convert the surgery id (patient + visit + surgery) into an integer id here
  ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

  ## Adjustment variables: drop those that are constant within this subgroup
  adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
    vv <- xp[[v]]
    if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
    else uniqueN(as.character(vv[!is.na(vv)])) > 1
  }, logical(1))]
  ADJ <- if (length(adj_ok)) paste(adj_ok, collapse = " + ") else "1"

  ## Outcome: 1 if this episode lasted longer than 10 minutes, 0 otherwise
  ep[, long_epi := as.integer(duration > 10)]
  keep1 <- c("long_epi", "Group", adj_ok, "sid_int")   # data.table's with=FALSE only accepts a character vector
  d1 <- ep[complete.cases(ep[, keep1, with = FALSE])]

  m1 <- geeglm(as.formula(paste("long_epi ~ Group +", ADJ)),
               data = d1, id = sid_int, weights = d1$w,
               family = binomial(), corstr = "exchangeable")

  ## Extract the "ciprofol vs propofol" coefficient (control = propofol)
  s1   <- summary(m1)$coefficients
  est1 <- s1["GroupCiprofol_only", "Estimate"]    # log-scale coefficient
  se1  <- s1["GroupCiprofol_only", "Std.err"]     # standard error
  p1   <- 2 * pnorm(-abs(est1 / se1))            # P value

  ## OR and confidence interval = exponentiated coefficient
  r1 <- data.table(
    analysis   = LABEL,
    model      = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d1),
    n_patient  = uniqueN(d1$PATIENT_ID),
    estimate   = exp(est1),
    lo         = exp(est1 - 1.96 * se1),
    hi         = exp(est1 + 1.96 * se1),
    p          = p1,
    effect_type = "OR",
    effect     = sprintf("%.2f (%.2f-%.2f)", exp(est1), exp(est1 - 1.96 * se1), exp(est1 + 1.96 * se1)))
  result <- r1

  cat("M1 proportion of episodes >10 min  OR =", r1$effect, "  P =", sprintf("%.3g", p1), "\n")
  cat("   （", nrow(d1), "episodes, involving", uniqueN(d1$PATIENT_ID), "patients)\n")


  ## ==========================================================================
  ## Step 5: model M2 - severity score (AUC) of each hypotension episode, comparison between the two groups
  ## ==========================================================================

  ep[, log_eauc := log(event_auc + 0.1)]   # +0.1 avoids log(0)
  keep2 <- c("log_eauc", "Group", adj_ok, "sid_int")
  d2 <- ep[complete.cases(ep[, keep2, with = FALSE])]

  m2 <- geeglm(as.formula(paste("log_eauc ~ Group +", ADJ)),
               data = d2, id = sid_int, weights = d2$w,
               family = gaussian(), corstr = "exchangeable")

  s2   <- summary(m2)$coefficients
  est2 <- s2["GroupCiprofol_only", "Estimate"]
  se2  <- s2["GroupCiprofol_only", "Std.err"]
  p2   <- 2 * pnorm(-abs(est2 / se2))

  ## Convert the log coefficient to a percentage change: 100 * (exp(coefficient) - 1)
  r2 <- data.table(
    analysis   = LABEL,
    model      = "M2 Per-episode hypotension burden (AUC)",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d2),
    n_patient  = uniqueN(d2$PATIENT_ID),
    estimate   = 100 * (exp(est2) - 1),
    lo         = 100 * (exp(est2 - 1.96 * se2) - 1),
    hi         = 100 * (exp(est2 + 1.96 * se2) - 1),
    p          = p2,
    effect_type = "Pct",
    effect     = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                         100 * (exp(est2) - 1),
                         100 * (exp(est2 - 1.96 * se2) - 1),
                         100 * (exp(est2 + 1.96 * se2) - 1)))
  result <- rbind(result, r2)

  cat("M2 per-episode AUC change              ", r2$effect, "  P =", sprintf("%.3g", p2), "\n")


  ## ==========================================================================
  ## Step 6: model M3 - duration of each hypotension episode in minutes, comparison between the two groups
  ## ==========================================================================

  ep[, log_duration := log(duration + 0.1)]
  keep3 <- c("log_duration", "Group", adj_ok, "sid_int")
  d3 <- ep[complete.cases(ep[, keep3, with = FALSE])]

  m3 <- geeglm(as.formula(paste("log_duration ~ Group +", ADJ)),
               data = d3, id = sid_int, weights = d3$w,
               family = gaussian(), corstr = "exchangeable")

  s3   <- summary(m3)$coefficients
  est3 <- s3["GroupCiprofol_only", "Estimate"]
  se3  <- s3["GroupCiprofol_only", "Std.err"]
  p3   <- 2 * pnorm(-abs(est3 / se3))

  r3 <- data.table(
    analysis   = LABEL,
    model      = "M3 Duration of each hypotensive episode",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d3),
    n_patient  = uniqueN(d3$PATIENT_ID),
    estimate   = 100 * (exp(est3) - 1),
    lo         = 100 * (exp(est3 - 1.96 * se3) - 1),
    hi         = 100 * (exp(est3 + 1.96 * se3) - 1),
    p          = p3,
    effect_type = "Pct",
    effect     = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                         100 * (exp(est3) - 1),
                         100 * (exp(est3 - 1.96 * se3) - 1),
                         100 * (exp(est3 + 1.96 * se3) - 1)))
  result <- rbind(result, r3)

  cat("M3 per-episode duration change         ", r3$effect, "  P =", sprintf("%.3g", p3), "\n")


  ## ==========================================================================
  ## Step 7: model M4 - number of hypotension episodes per procedure, comparison between the two groups
  ## ==========================================================================

  xp[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
  xp[, lg_dur  := log(anes_duration_mins)]   # Log-transform anesthesia duration and use it as an offset

  keep4 <- c("n_epi", "Group", adj_ok, "lg_dur", "sid_int")
  d4 <- xp[complete.cases(xp[, keep4, with = FALSE])]

  m4 <- geeglm(as.formula(paste("n_epi ~ Group +", ADJ, "+ offset(lg_dur)")),
               data = d4, id = sid_int, weights = d4$sw,
               family = poisson(), corstr = "exchangeable")

  s4   <- summary(m4)$coefficients
  est4 <- s4["GroupCiprofol_only", "Estimate"]
  se4  <- s4["GroupCiprofol_only", "Std.err"]
  p4   <- 2 * pnorm(-abs(est4 / se4))

  r4 <- data.table(
    analysis   = LABEL,
    model      = "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)",
    unit       = "Procedure (one anesthesia)",
    n_unit     = nrow(d4),
    n_patient  = uniqueN(d4$PATIENT_ID),
    estimate   = exp(est4),
    lo         = exp(est4 - 1.96 * se4),
    hi         = exp(est4 + 1.96 * se4),
    p          = p4,
    effect_type = "RR",
    effect     = sprintf("%.2f (%.2f-%.2f)", exp(est4), exp(est4 - 1.96 * se4), exp(est4 + 1.96 * se4)))
  result <- rbind(result, r4)

  cat("M4 episode rate (per hour)         RR =", r4$effect, "  P =", sprintf("%.3g", p4), "\n")
  cat("\n  ---- Results after IPTW weighting within the subgroup ----\n")
  print(result[, .(model, unit, effect, effect_type, p)])


  ## ==========================================================================
  ## Step 8: weighted absolute metrics (weighted incidence and other summaries for each group)
  ## ==========================================================================

  abs_tab <- xp[, .(n_surgery    = .N,
                    n_patient    = uniqueN(PATIENT_ID),
                    prolonged_pct = 100 * weighted.mean(prolonged, sw, na.rm = TRUE),
                    MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
                    n_epi_mean   = weighted.mean(n_epi, sw, na.rm = TRUE)),
                 by = Group]
  cat("\n  ---- Absolute metrics after weighting ----\n")
  print(abs_tab)


  ## ==========================================================================
  ## Step 9: save the results to disk
  ## ==========================================================================

  fwrite(crude,   file.path(OUT, paste0("Sensitivity_Analysis_Crude_",     LABEL, ".csv")))
  fwrite(abs_tab, file.path(OUT, paste0("Sensitivity_Analysis_Absolute_", LABEL, ".csv")))
  fwrite(bal,     file.path(OUT, paste0("Table1_IPTW_",                    LABEL, ".csv")))
  fwrite(ess_tab, file.path(OUT, paste0("Sensitivity_Analysis_ESS_",       LABEL, ".csv")))
  fwrite(result,  file.path(OUT, paste0("Sensitivity_Analysis_Results_",  LABEL, ".csv")))

  tryCatch({
    suppressPackageStartupMessages({ library(flextable); library(officer) })
    ft <- flextable(as.data.frame(result[, .(
             Model = model, Unit = unit, Sample_size = n_unit, n_patients = n_patient,
             Estimate = effect, Type = effect_type,
             P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])) |>
         theme_booktabs() |> autofit() |>
         set_caption(caption = sprintf(paste(
           "Sensitivity analysis - %s. Weighted GEE models after re-fitting IPTW within this subgroup.",
           "M1-M3: one hypotensive episode per row; M4: one surgical procedure per row;",
           "all clustered by surgical procedure."), LABEL))
    print(read_docx() |> body_add_flextable(ft),
          target = file.path(OUT, paste0("Sensitivity_Analysis_Results_", LABEL, ".docx")))

    ft2 <- flextable(as.data.frame(bal)) |> theme_booktabs() |> autofit() |>
           set_caption(caption = sprintf("Table 1. Covariate balance before and after IPTW - %s.", LABEL))
    print(read_docx() |> body_add_flextable(ft2),
          target = file.path(OUT, paste0("Table1_IPTW_", LABEL, ".docx")))
  }, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

}

