## =============================================================================
## 4-11 subgroup analysis - sex male.R
## -----------------------------------------------------------------------------
## 4-11 Subgroup analysis: male only
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


## Root directory. All inputs and outputs are placed under this directory
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## Three ID columns. The ID column types differ across the three tables, so convert them all to character before use
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## File label: determines all output file names
lab <- "male"

PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR",
             "operation_type", "Hypertension", "Diabetes",
             "Coronary_artery_disease", "Heart_failure")

ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

cat("[Setup] ready | Output directory:", OUT, "\n")


## ============================================================================
## Step 2: read the procedure-level analysis table, subset by "male only", and report the cohort size
## ============================================================================

## Surgery table: one row = one procedure
tab <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))

## Convert the ID columns of all three tables to character
for (v in ID_VARS) set(tab, j = v, value = as.character(tab[[v]]))


x <- tab[as.character(SEX) == "男"]   

cat("Sensitivity cohort included (male only):",
    sprintf("Procedures %d | patients %d", nrow(x), uniqueN(x$PATIENT_ID)), "\n")
cat(sprintf("  Group: Ciprofol %d procedures / Propofol %d procedures\n",
            sum(x$Group == "Ciprofol_only"), sum(x$Group == "Propofol_only")))


## ============================================================================
## Step 3: unweighted (crude) comparison 
## ============================================================================
crude <- x[, .(n_surgery = .N,
               n_patient = uniqueN(PATIENT_ID),
               prolonged_pct    = 100 * mean(prolonged, na.rm = TRUE),
               MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
               n_epi_mean = mean(n_epi, na.rm = TRUE)),
             by = Group]
print(crude)


## ============================================================================
## Step 4: re-fit the propensity score within this subgroup and recompute the weights (IPTW)
## ============================================================================

## Flag each surgery with whether it received ciprofol
x[, treat := as.integer(Group == "Ciprofol_only")]

use <- character(0)
for (v in PS_VARS) {
  if (!v %in% names(x)) next
  vv <- x[[v]]
  if (is.numeric(vv)) {
    if (is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0) use <- c(use, v)
  } else {
    if (uniqueN(as.character(vv[!is.na(vv)])) > 1) use <- c(use, v)
  }
}
dropped <- setdiff(PS_VARS, use)
if (length(dropped))
  cat("   [PS] constant within this subgroup, dropped automatically:", paste(dropped, collapse = ", "), "\n")

cc <- complete.cases(x[, ..use])
xp <- x[cc]
if (nrow(xp) < 100 || min(table(xp$treat)) < 10)
  stop("Too few procedures in this subgroup, or one arm is too small; propensity score cannot be fitted")

ps_fit <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
              data = xp, family = binomial())
xp[, ps_hat := predict(ps_fit, type = "response")]

## Stabilized weights: treatment group = treatment rate/PS, control group = (1-treatment rate)/(1-PS)
prev <- mean(xp$treat)
xp[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

xp[, `:=`(ql = as.numeric(quantile(sw_raw, 0.01)),
          qh = as.numeric(quantile(sw_raw, 0.99))), by = Group]
xp[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
xp[, c("ql", "qh") := NULL]

## Common support (overlap region) flag: mark procedures whose PS does not overlap between the two groups
lo_overlap <- max(min(xp$ps_hat[xp$treat == 1]), min(xp$ps_hat[xp$treat == 0]))
hi_overlap <- min(max(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 0]))
xp[, in_overlap := ps_hat >= lo_overlap & ps_hat <= hi_overlap]

## Balance before and after weighting (SMD = between-group mean difference / pooled SD; <0.1 counts as balanced)
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
bal_rows <- list()
for (v in use) {
  vv <- xp[[v]]
  if (is.numeric(vv)) {
    bal_rows[[length(bal_rows) + 1]] <- data.table(
      variable   = v,
      smd_before = wsmd(vv, rep(1, nrow(xp))),
      smd_after  = wsmd(vv, xp$sw))
  } else {
    for (lv in setdiff(unique(as.character(vv)), NA)) {
      ind <- as.integer(as.character(vv) == lv)
      bal_rows[[length(bal_rows) + 1]] <- data.table(
        variable   = paste0(v, ": ", lvl_en(lv)),
        smd_before = wsmd(ind, rep(1, nrow(xp))),
        smd_after  = wsmd(ind, xp$sw))
    }
  }
}
balance <- rbindlist(bal_rows)

## Effective sample size (ESS) and weight distribution, computed by group
ess_tab <- xp[, .(Group = as.character(first(Group)),
                   n_procedures = .N,
                   n_patients   = uniqueN(PATIENT_ID),
                   weight_median = median(sw),
                   weight_max = max(sw),
                   ESS      = sum(sw)^2 / sum(sw^2)),
              by = Group][, Group := NULL]

cat(sprintf("   [PS] Procedures included in weighting %d (patients %d); PS range Ciprofol %.3f-%.3f, Propofol %.3f-%.3f\n",
            nrow(xp), uniqueN(xp$PATIENT_ID),
            min(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 1]),
            min(xp$ps_hat[xp$treat == 0]), max(xp$ps_hat[xp$treat == 0])))
cat(sprintf("   [Balance] Max SMD before weighting %.3f (%d variables >0.1) → max SMD after weighting %.3f (%d variables)\n",
            max(balance$smd_before, na.rm = TRUE), sum(balance$smd_before > .1, na.rm = TRUE),
            max(balance$smd_after,  na.rm = TRUE), sum(balance$smd_after  > .1, na.rm = TRUE)))
e_c <- ess_tab[Group == "Ciprofol_only", ESS]
e_p <- ess_tab[Group == "Propofol_only", ESS]
cat(sprintf("   [Weight] Within-group truncation quantiles 1%%/99%%; Ciprofol ESS %.0f / %d procedures, Propofol ESS %.0f / %d procedures\n",
            e_c, sum(xp$treat == 1), e_p, sum(xp$treat == 0)))
cat(sprintf("   [Overlap] PS overlap region between groups = [%.3f, %.3f]; procedures outside the region %d (%.1f%%)\n",
            lo_overlap, hi_overlap, sum(!xp$in_overlap), 100 * mean(!xp$in_overlap)))
if (sum(!xp$in_overlap) > 0 && "Hospital" %in% names(xp)) {
  cat("          Center distribution of procedures outside the region:\n")
  print(xp[in_overlap == FALSE, .N, by = .(Hospital, Group)][order(-N)])
}

## Weighted analysis population: xp (with ps_hat, sw)
xw <- xp


## ============================================================================
## Step 5: read the event-level table, keep only the events of this cohort, and attach the group/weight/adjustment variables
## ============================================================================

epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))
for (v in ID_VARS) set(epi, j = v, value = as.character(epi[[v]]))

## Keep only the events belonging to surgeries in this cohort
keep_ids <- xw[, .(PATIENT_ID, VISIT_ID, OPER_ID)]
ep <- epi[keep_ids, on = ID_VARS, nomatch = NULL]

## Attach the group, weight, and adjustment variables of the surgery each event belongs to → ep: one row = one episode
ep <- ep[xw[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
         on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

cat("Event table included:", nrow(ep), "episodes, involving", uniqueN(ep$PATIENT_ID), "patients\n")


## ============================================================================
## Step 6: model M1 - proportion of hypotension episodes lasting longer than 10 minutes, comparison between the two groups
## ============================================================================

adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  if (!v %in% names(xw)) return(FALSE)
  vv <- xw[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]
ADJ <- if (length(adj_ok)) paste(adj_ok, collapse = " + ") else "1"

## The event table also needs an integer cluster id (the id in geeglm must be an integer)
ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

## Outcome: 1 if this episode lasted longer than 10 minutes, 0 otherwise
ep[, long_epi := as.integer(duration > 10)]

## Drop events with missing values for variables used in the formula
d1 <- ep[complete.cases(ep[, c("long_epi", "Group", adj_ok, "sid_int"), with = FALSE])]
result_sens <- data.table(
  analysis = character(), model = character(), unit = character(),
  n_unit = integer(), n_patient = integer(),
  estimate = numeric(), lo = numeric(), hi = numeric(),
  p = numeric(), effect_type = character(), effect = character())
if (nrow(d1) > 50 && uniqueN(d1$Group) > 1) {
  f1  <- as.formula(paste("long_epi ~ Group +", ADJ))
  m1  <- geeglm(f1, data = d1, id = sid_int, weights = d1$w,
                family = binomial(), corstr = "exchangeable")
  s1  <- summary(m1)$coefficients
  est1 <- s1["GroupCiprofol_only", "Estimate"]
  se1  <- s1["GroupCiprofol_only", "Std.err"]
  p1   <- 2 * pnorm(-abs(est1 / se1))
  result_sens <- rbind(result_sens, data.table(
    analysis   = lab,
    model      = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d1), n_patient = uniqueN(d1$PATIENT_ID),
    estimate   = exp(est1),
    lo  = exp(est1 - 1.96 * se1), hi = exp(est1 + 1.96 * se1),
    p          = p1,
    effect_type = "OR",
    effect     = sprintf("%.2f (%.2f-%.2f)", exp(est1), exp(est1 - 1.96 * se1), exp(est1 + 1.96 * se1))))
  cat("M1 proportion of episodes >10 min  OR =", result_sens[model %like% "M1", effect],
      "  P =", sprintf("%.3g", p1), "\n")
} else {
  cat("M1 too few episodes, skipped\n")
}


## ============================================================================
## Step 7: model M2 - severity score (AUC) of each hypotension episode, comparison between the two groups
## ============================================================================
ep[, log_eauc := log(event_auc + 0.1)]   # +0.1 avoids log(0)
d2 <- ep[complete.cases(ep[, c("log_eauc", "Group", adj_ok, "sid_int"), with = FALSE])]
if (nrow(d2) > 50 && uniqueN(d2$Group) > 1) {
  f2  <- as.formula(paste("log_eauc ~ Group +", ADJ))
  m2  <- geeglm(f2, data = d2, id = sid_int, weights = d2$w,
                family = gaussian(), corstr = "exchangeable")
  s2  <- summary(m2)$coefficients
  est2 <- s2["GroupCiprofol_only", "Estimate"]
  se2  <- s2["GroupCiprofol_only", "Std.err"]
  p2   <- 2 * pnorm(-abs(est2 / se2))
  result_sens <- rbind(result_sens, data.table(
    analysis   = lab,
    model      = "M2 Per-episode hypotension burden (AUC)",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d2), n_patient = uniqueN(d2$PATIENT_ID),
    estimate   = 100 * (exp(est2) - 1),
    lo  = 100 * (exp(est2 - 1.96 * se2) - 1),
    hi  = 100 * (exp(est2 + 1.96 * se2) - 1),
    p          = p2,
    effect_type = "Pct",
    effect     = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                         100 * (exp(est2) - 1),
                         100 * (exp(est2 - 1.96 * se2) - 1),
                         100 * (exp(est2 + 1.96 * se2) - 1))))
  cat("M2 per-episode AUC change              ", result_sens[model %like% "M2", effect],
      "  P =", sprintf("%.3g", p2), "\n")
} else {
  cat("M2 too few episodes, skipped\n")
}


## ============================================================================
## Step 8: model M3 - duration of each hypotension episode in minutes, comparison between the two groups
## ============================================================================
ep[, log_duration := log(duration + 0.1)]
d3 <- ep[complete.cases(ep[, c("log_duration", "Group", adj_ok, "sid_int"), with = FALSE])]
if (nrow(d3) > 50 && uniqueN(d3$Group) > 1) {
  f3  <- as.formula(paste("log_duration ~ Group +", ADJ))
  m3  <- geeglm(f3, data = d3, id = sid_int, weights = d3$w,
                family = gaussian(), corstr = "exchangeable")
  s3  <- summary(m3)$coefficients
  est3 <- s3["GroupCiprofol_only", "Estimate"]
  se3  <- s3["GroupCiprofol_only", "Std.err"]
  p3   <- 2 * pnorm(-abs(est3 / se3))
  result_sens <- rbind(result_sens, data.table(
    analysis   = lab,
    model      = "M3 Duration of each hypotensive episode",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d3), n_patient = uniqueN(d3$PATIENT_ID),
    estimate   = 100 * (exp(est3) - 1),
    lo  = 100 * (exp(est3 - 1.96 * se3) - 1),
    hi  = 100 * (exp(est3 + 1.96 * se3) - 1),
    p          = p3,
    effect_type = "Pct",
    effect     = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                         100 * (exp(est3) - 1),
                         100 * (exp(est3 - 1.96 * se3) - 1),
                         100 * (exp(est3 + 1.96 * se3) - 1))))
  cat("M3 per-episode duration change         ", result_sens[model %like% "M3", effect],
      "  P =", sprintf("%.3g", p3), "\n")
} else {
  cat("M3 too few episodes, skipped\n")
}


## ============================================================================
## Step 9: model M4 - number of hypotension episodes per procedure, comparison between the two groups
## ============================================================================
xw[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
xw[, lg_dur  := log(anes_duration_mins)]   # Log-transform anesthesia duration and use it as an offset

d4 <- xw[complete.cases(xw[, c("n_epi", "Group", adj_ok, "lg_dur", "sid_int"), with = FALSE])]
m4  <- geeglm(as.formula(paste("n_epi ~ Group +", ADJ, "+ offset(lg_dur)")),
              data = d4, id = sid_int, weights = d4$sw,
              family = poisson(), corstr = "exchangeable")
s4  <- summary(m4)$coefficients
est4 <- s4["GroupCiprofol_only", "Estimate"]
se4  <- s4["GroupCiprofol_only", "Std.err"]
p4   <- 2 * pnorm(-abs(est4 / se4))
result_sens <- rbind(result_sens, data.table(
  analysis   = lab,
  model      = "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)",
  unit       = "Procedure (one anesthesia)",
  n_unit     = nrow(d4), n_patient = uniqueN(d4$PATIENT_ID),
  estimate   = exp(est4),
  lo  = exp(est4 - 1.96 * se4), hi = exp(est4 + 1.96 * se4),
  p          = p4,
  effect_type = "RR",
  effect     = sprintf("%.2f (%.2f-%.2f)", exp(est4), exp(est4 - 1.96 * se4), exp(est4 + 1.96 * se4))))
cat("M4 episode rate (per hour)         RR =", result_sens[model %like% "M4", effect],
    "  P =", sprintf("%.3g", p4), "\n")

cat("\n==================== Four-model results after weighting (male only) ====================\n")
print(result_sens[, .(model, unit, n_unit, n_patient, type = effect_type, effect,
                      p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])


## ============================================================================
## Step 10: weighted absolute metrics - weighted incidence and mean number of episodes for each group
## ============================================================================
abs_tab <- xw[, .(n_surgery = .N,
                  n_patient = uniqueN(PATIENT_ID),
                  prolonged_pct     = 100 * weighted.mean(prolonged, sw, na.rm = TRUE),
                  MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
                  n_epi_mean = weighted.mean(n_epi, sw, na.rm = TRUE)),
                by = Group]
print(abs_tab)


## ============================================================================
## Step 11: save the results to disk
## ============================================================================
fwrite(crude,       file.path(OUT, paste0("Sensitivity_Analysis_Crude_",    lab, ".csv")))
fwrite(abs_tab,     file.path(OUT, paste0("Sensitivity_Analysis_Absolute_", lab, ".csv")))
fwrite(balance,     file.path(OUT, paste0("Table1_IPTW_",                   lab, ".csv")))
fwrite(ess_tab,     file.path(OUT, paste0("Sensitivity_Analysis_ESS_",      lab, ".csv")))
fwrite(result_sens, file.path(OUT, paste0("Sensitivity_Analysis_Results_",  lab, ".csv")))

## docx: results table
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(result_sens[, .(
           Model = model, Unit = unit, Sample_size = n_unit, n_patients = n_patient,
           Estimate = effect, Type = effect_type,
           P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])) |>
       theme_booktabs() |> autofit() |>
       set_caption(caption = paste(
         "Sensitivity analysis - male. Weighted GEE models after re-fitting IPTW within this subgroup.",
         "M1-M3: one hypotensive episode per row; M4: one surgical procedure per row; all clustered by surgical procedure."))
  print(read_docx() |> body_add_flextable(ft),
        target = file.path(OUT, paste0("Sensitivity_Analysis_Results_", lab, ".docx")))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

## docx: balance table
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft2 <- flextable(as.data.frame(balance)) |> theme_booktabs() |> autofit() |>
        set_caption(caption = paste0("Table 1. Covariate balance before and after IPTW - ", lab, "."))
  print(read_docx() |> body_add_flextable(ft2),
        target = file.path(OUT, paste0("Table1_IPTW_", lab, ".docx")))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

