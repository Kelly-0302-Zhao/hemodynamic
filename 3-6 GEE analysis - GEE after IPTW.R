## =============================================================================
## 3-6 GEE analysis -  the four main weighted models (main results of the paper)
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
DATA <- file.path(D, "JCA_run")      # The cleaned intermediate tables are placed here (output of the 1-x series of scripts)
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## Three ID columns.
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 4 baseline variables adjusted for in the GEE model
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## Parse all timestamps as UTC
MY_TZ <- "UTC"

## ============================================================================
## Step 2: read the three tables, assemble the "weighted analysis population" and report how many procedures and patients were included
## ============================================================================

## Surgery table: one row = one procedure
d <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))

## Weight table: the propensity score (ps_hat) and stabilized weight (sw) computed in 3-5
w <- setDT(readRDS(file.path(OUT, "IPTW_Weights.rds")))

## Episode table: one row = one hypotensive episode
epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## Convert the ID columns of all three tables to character (merging misaligns when the types differ - a pitfall already hit)
for (v in ID_VARS) {
  set(d,   j = v, value = as.character(d[[v]]))
  set(w,   j = v, value = as.character(w[[v]]))
  set(epi, j = v, value = as.character(epi[[v]]))
}

## Attach the weights to the procedure table; then drop procedures with missing "prolonged hypotension" data
dw <- d[w[, .(PATIENT_ID, VISIT_ID, OPER_ID, ps_hat, sw, in_overlap)],
        on = ID_VARS, nomatch = 0]
dw <- dw[!is.na(prolonged)]

cat("Weighted analysis population:", nrow(dw), "procedures |", uniqueN(dw$PATIENT_ID), "patients\n")
cat(sprintf("  (total %d procedures in the surgery table; the difference is procedures with missing propensity score covariates)\n\n",
            d[!is.na(prolonged), .N]))


## ============================================================================
## Step 3: model M1 - proportion of hypotension episodes lasting longer than 10 minutes, comparison between the two groups
##   Data: episode table (one row = one hypotensive episode), each episode inherits the weight of its parent procedure
##   Model: logistic GEE, reporting the OR (ciprofol vs propofol)
## ============================================================================

## Attach the group, weight and adjustment variables of the procedure table to each episode -> ep: one row = one episode
ep <- epi[dw[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                 Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
          on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

## GEE requires the cluster id to be an integer; convert the surgery id (patient + visit + surgery) into an integer id here
ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

## Outcome: 1 if this episode lasted longer than 10 minutes, 0 otherwise
ep[, long_epi := as.integer(duration > 10)]

## Drop events with missing values for variables used in the formula
d1 <- ep[complete.cases(ep[, .(long_epi, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]

## Run the model
m1 <- geeglm(long_epi ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d1, id = sid_int, weights = d1$w,
             family = binomial(), corstr = "exchangeable")

## Extract the "ciprofol vs propofol" coefficient from the model (control = propofol)
s1 <- summary(m1)$coefficients
est1 <- s1["GroupCiprofol_only", "Estimate"]   # log-scale coefficient
se1  <- s1["GroupCiprofol_only", "Std.err"]    # standard error
p1   <- 2 * pnorm(-abs(est1 / se1))            # P value

## OR and confidence interval = exponentiated coefficient
result_main <- data.table(
  analysis  = "IPTW-weighted (full cohort)",
  model     = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
  unit      = "Episode (one hypotensive episode)",
  n_unit    = nrow(d1),                       # how many episodes
  n_patient = uniqueN(d1$PATIENT_ID),         # how many patients are involved
  estimate  = exp(est1), lo = exp(est1 - 1.96 * se1), hi = exp(est1 + 1.96 * se1),
  p         = p1,
  effect_type = "OR",
  effect    = sprintf("%.2f (%.2f-%.2f)", exp(est1), exp(est1 - 1.96 * se1), exp(est1 + 1.96 * se1)))

cat("M1 proportion of episodes >10 min  OR =", result_main$effect,
    "  P =", sprintf("%.3g", p1), "\n")
cat("   （", nrow(d1), "episodes, involving", uniqueN(d1$PATIENT_ID), "patients)\n")


## ============================================================================
## Step 4: model M2 - severity score (AUC) of each hypotension episode, comparison between the two groups
##   AUC = total area of "depth x duration" for this episode, measuring the overall severity of this episode
##   Because the AUC is strongly skewed, take the log before comparison; results are converted back to "percentage change"
## ============================================================================

ep[, log_eauc := log(event_auc + 0.1)]   # +0.1 avoids log(0)
d2 <- ep[complete.cases(ep[, .(log_eauc, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]

m2 <- geeglm(log_eauc ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d2, id = sid_int, weights = d2$w,
             family = gaussian(), corstr = "exchangeable")

s2   <- summary(m2)$coefficients
est2 <- s2["GroupCiprofol_only", "Estimate"]
se2  <- s2["GroupCiprofol_only", "Std.err"]
p2   <- 2 * pnorm(-abs(est2 / se2))

## Convert the log coefficient to a percentage change: 100 * (exp(coefficient) - 1)
result_main <- rbind(result_main, data.table(
  analysis  = "IPTW-weighted (full cohort)",
  model     = "M2 Per-episode hypotension burden (AUC)",
  unit      = "Episode (one hypotensive episode)",
  n_unit    = nrow(d2), n_patient = uniqueN(d2$PATIENT_ID),
  estimate  = 100 * (exp(est2) - 1),
  lo        = 100 * (exp(est2 - 1.96 * se2) - 1),
  hi        = 100 * (exp(est2 + 1.96 * se2) - 1),
  p         = p2,
  effect_type = "Pct",
  effect    = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                      100 * (exp(est2) - 1),
                      100 * (exp(est2 - 1.96 * se2) - 1),
                      100 * (exp(est2 + 1.96 * se2) - 1))))

cat("M2 per-episode AUC change              ", result_main[model %like% "M2", effect],
    "  P =", sprintf("%.3g", p2), "\n")


## ============================================================================
## Step 5: model M3 - duration of each hypotension episode in minutes, comparison between the two groups
##   Same approach as M2, only with episode duration as the dependent variable
## ============================================================================

ep[, log_duration := log(duration + 0.1)]
d3 <- ep[complete.cases(ep[, .(log_duration, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]

m3 <- geeglm(log_duration ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d3, id = sid_int, weights = d3$w,
             family = gaussian(), corstr = "exchangeable")

s3   <- summary(m3)$coefficients
est3 <- s3["GroupCiprofol_only", "Estimate"]
se3  <- s3["GroupCiprofol_only", "Std.err"]
p3   <- 2 * pnorm(-abs(est3 / se3))

result_main <- rbind(result_main, data.table(
  analysis  = "IPTW-weighted (full cohort)",
  model     = "M3 Duration of each hypotensive episode",
  unit      = "Episode (one hypotensive episode)",
  n_unit    = nrow(d3), n_patient = uniqueN(d3$PATIENT_ID),
  estimate  = 100 * (exp(est3) - 1),
  lo        = 100 * (exp(est3 - 1.96 * se3) - 1),
  hi        = 100 * (exp(est3 + 1.96 * se3) - 1),
  p         = p3,
  effect_type = "Pct",
  effect    = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                      100 * (exp(est3) - 1),
                      100 * (exp(est3 - 1.96 * se3) - 1),
                      100 * (exp(est3 + 1.96 * se3) - 1))))

cat("M3 per-episode duration change         ", result_main[model %like% "M3", effect],
    "  P =", sprintf("%.3g", p3), "\n")


## ============================================================================
## Step 6: model M4 - number of hypotension episodes per procedure, comparison between the two groups
##   Data: procedure table (one row = one surgery), dependent variable = n_epi (number of episodes in this surgery)
##   Surgeries differ in length and longer ones naturally accumulate more episodes, so add offset(log anesthesia duration),
##   What the model actually compares is "how many per hour", i.e. the "rate" of episode counts
## ============================================================================
adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  vv <- dw[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]

dw[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
dw[, lg_dur  := log(anes_duration_mins)]   # Log-transform anesthesia duration and use it as an offset

d4 <- dw[complete.cases(dw[, c("n_epi", "Group", adj_ok, "lg_dur", "sw"), with = FALSE])]

f4 <- as.formula(paste("n_epi ~ Group +", paste(adj_ok, collapse = " + "), "+ offset(lg_dur)"))
m4 <- geeglm(f4, data = d4, id = sid_int, weights = d4$sw,
             family = poisson(), corstr = "exchangeable")

s4   <- summary(m4)$coefficients
est4 <- s4["GroupCiprofol_only", "Estimate"]
se4  <- s4["GroupCiprofol_only", "Std.err"]
p4   <- 2 * pnorm(-abs(est4 / se4))

result_main <- rbind(result_main, data.table(
  analysis  = "IPTW-weighted (full cohort)",
  model     = "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)",
  unit      = "Procedure (one anesthesia)",
  n_unit    = nrow(d4), n_patient = uniqueN(d4$PATIENT_ID),
  estimate  = exp(est4), lo = exp(est4 - 1.96 * se4), hi = exp(est4 + 1.96 * se4),
  p         = p4,
  effect_type = "RR",
  effect    = sprintf("%.2f (%.2f-%.2f)", exp(est4), exp(est4 - 1.96 * se4), exp(est4 + 1.96 * se4))))

cat("M4 episode rate (per hour)         RR =", result_main[model %like% "M4", effect],
    "  P =", sprintf("%.3g", p4), "\n")

## ============================================================================
## Step 7: save the results of the four models (one csv + one Word file each)
##   Expected output: out\Analysis_Results_main.csv / .docx
## ============================================================================

fwrite(result_main, file.path(OUT, "Analysis_Results_main.csv"))

## Save a Word version
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft1 <- flextable(as.data.frame(result_main[, .(
          Model = model, Unit = unit, Sample_size = n_unit, n_patients = n_patient,
          Estimate = effect, Type = effect_type,
          P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])) |>
       theme_booktabs() |> autofit() |>
       set_caption(caption = paste(
         "Table 5. Outcomes after IPTW weighting (M1-M3: one hypotensive episode per row;",
         "M4: one surgical procedure per row; all clustered by surgical procedure)"))
  doc <- read_docx() |> body_add_flextable(ft1)
  print(doc, target = file.path(OUT, "Analysis_Results_main.docx"))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

cat("\n==================== Main results (after IPTW weighting)====================\n")
print(result_main[, .(model, unit, n_unit, n_patient, type = effect_type,
                      effect, p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])

## ============================================================================
## Step 8: control analysis - restricted to the range where the two groups' characteristics overlap, re-weight and redo the four models
##   AnYi hospital has never used ciprofol; its propofol patients are too far from the ciprofol patients on the "propensity score" and fall in the region where the two groups' characteristics do not overlap
## ============================================================================

if (dw[in_overlap == FALSE, .N] > 0) {

  cat("\n==================== Control analysis: PS overlap region only (re-weighted with IPTW)====================\n")
  cat(sprintf("Excluded %d non-overlapping procedures (%.1f%%), of which %d were at AnYi hospital\n",
              dw[in_overlap == FALSE, .N], 100 * mean(!dw$in_overlap),
              dw[in_overlap == FALSE & Hospital == "AnYi", .N]))
  dw2 <- dw[in_overlap == TRUE]
  print(dw2[, .N, by = .(Hospital, Group)][order(Hospital, Group)])

  ## ---- 8a. refit the propensity score and recompute the weights within the overlap region----
  dw2[, treat := as.integer(Group == "Ciprofol_only")]

  ## 10 baseline variables used for the propensity score
  PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR",
               "operation_type", "Hypertension", "Diabetes",
               "Coronary_artery_disease", "Heart_failure")
  use <- character(0)
  for (v in PS_VARS) {
    if (!v %in% names(dw2)) next
    vv <- dw2[[v]]
    if (is.numeric(vv)) {
      if (is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0) use <- c(use, v)
    } else {
      if (uniqueN(as.character(vv[!is.na(vv)])) > 1) use <- c(use, v)
    }
  }
  ps_fit <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
                data = dw2[complete.cases(dw2[, ..use])], family = binomial())
  dw2 <- dw2[complete.cases(dw2[, ..use])]
  dw2[, ps_hat := predict(ps_fit, type = "response")]

  ## Stabilized weights: treatment group = treatment rate/PS, control group = (1-treatment rate)/(1-PS)
  prev <- mean(dw2$treat)
  dw2[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

  ## Weight truncation: within each group, cut extreme weights at the 1% and 99% quantiles
  dw2[, ql := as.numeric(quantile(sw_raw, 0.01)), by = Group]
  dw2[, qh := as.numeric(quantile(sw_raw, 0.99)), by = Group]
  dw2[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
  dw2[, c("ql", "qh") := NULL]

  ## Covariate balance before and after weighting (SMD = between-group mean difference / pooled SD; <0.1 counts as balanced)
  wsmd <- function(v, tr, w) {
    v <- as.numeric(v)
    m1_ <- weighted.mean(v[tr == 1], w[tr == 1]); m0_ <- weighted.mean(v[tr == 0], w[tr == 0])
    v1_ <- sum(w[tr == 1] * (v[tr == 1] - m1_)^2) / sum(w[tr == 1])
    v0_ <- sum(w[tr == 0] * (v[tr == 0] - m0_)^2) / sum(w[tr == 0])
    sp_ <- sqrt((v1_ + v0_) / 2)
    if (!is.finite(sp_) || sp_ == 0) return(NA_real_)
    abs(m1_ - m0_) / sp_
  }
  bal_rows <- list()
  for (v in use) {
    vv <- dw2[[v]]
    if (is.numeric(vv)) {
      bal_rows[[length(bal_rows) + 1]] <- data.table(
        variable = v, smd_before = wsmd(vv, dw2$treat, rep(1, nrow(dw2))),
        smd_after = wsmd(vv, dw2$treat, dw2$sw))
    } else {
      for (lv in setdiff(unique(as.character(vv)), NA)) {
        ind <- as.integer(as.character(vv) == lv)
        bal_rows[[length(bal_rows) + 1]] <- data.table(
          variable = paste0(v, ": ", lvl_en(lv)),
          smd_before = wsmd(ind, dw2$treat, rep(1, nrow(dw2))),
          smd_after = wsmd(ind, dw2$treat, dw2$sw))
      }
    }
  }
  balance_overlap <- rbindlist(bal_rows)
  fwrite(balance_overlap, file.path(OUT, "Table1_IPTW_overlap_Balance_SMD.csv"))
  cat(sprintf("Max SMD before weighting %.3f → max SMD after weighting %.3f\n",
              max(balance_overlap$smd_before, na.rm = TRUE),
              max(balance_overlap$smd_after, na.rm = TRUE)))

  ## ---- 8b. rerun the four models within the overlap region
  ep2 <- epi[dw2[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                     Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
             on = ID_VARS, nomatch = NULL]
  ep2 <- ep2[!is.na(Group) & is.finite(w) & w > 0]
  ep2[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

  result_ov <- NULL

  ## M1
  ep2[, long_epi := as.integer(duration > 10)]
  d1o <- ep2[complete.cases(ep2[, .(long_epi, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
  m1o <- geeglm(long_epi ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
                data = d1o, id = sid_int, weights = d1o$w,
                family = binomial(), corstr = "exchangeable")
  so <- summary(m1o)$coefficients
  e1 <- so["GroupCiprofol_only", "Estimate"]; s1o <- so["GroupCiprofol_only", "Std.err"]
  result_ov <- rbind(result_ov, data.table(
    analysis = "PS overlap region only", model = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
    unit = "Episode (one hypotensive episode)", n_unit = nrow(d1o), n_patient = uniqueN(d1o$PATIENT_ID),
    estimate = exp(e1), lo = exp(e1 - 1.96 * s1o), hi = exp(e1 + 1.96 * s1o),
    p = 2 * pnorm(-abs(e1 / s1o)), effect_type = "OR",
    effect = sprintf("%.2f (%.2f-%.2f)", exp(e1), exp(e1 - 1.96 * s1o), exp(e1 + 1.96 * s1o))))

  ## M2
  ep2[, log_eauc := log(event_auc + 0.1)]
  d2o <- ep2[complete.cases(ep2[, .(log_eauc, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
  m2o <- geeglm(log_eauc ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
                data = d2o, id = sid_int, weights = d2o$w,
                family = gaussian(), corstr = "exchangeable")
  so <- summary(m2o)$coefficients
  e2 <- so["GroupCiprofol_only", "Estimate"]; s2o <- so["GroupCiprofol_only", "Std.err"]
  result_ov <- rbind(result_ov, data.table(
    analysis = "PS overlap region only", model = "M2 Per-episode hypotension burden (AUC)",
    unit = "Episode (one hypotensive episode)", n_unit = nrow(d2o), n_patient = uniqueN(d2o$PATIENT_ID),
    estimate = 100 * (exp(e2) - 1), lo = 100 * (exp(e2 - 1.96 * s2o) - 1), hi = 100 * (exp(e2 + 1.96 * s2o) - 1),
    p = 2 * pnorm(-abs(e2 / s2o)), effect_type = "Pct",
    effect = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                     100 * (exp(e2) - 1), 100 * (exp(e2 - 1.96 * s2o) - 1), 100 * (exp(e2 + 1.96 * s2o) - 1))))

  ## M3
  ep2[, log_duration := log(duration + 0.1)]
  d3o <- ep2[complete.cases(ep2[, .(log_duration, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]
  m3o <- geeglm(log_duration ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
                data = d3o, id = sid_int, weights = d3o$w,
                family = gaussian(), corstr = "exchangeable")
  so <- summary(m3o)$coefficients
  e3 <- so["GroupCiprofol_only", "Estimate"]; s3o <- so["GroupCiprofol_only", "Std.err"]
  result_ov <- rbind(result_ov, data.table(
    analysis = "PS overlap region only", model = "M3 Duration of each hypotensive episode",
    unit = "Episode (one hypotensive episode)", n_unit = nrow(d3o), n_patient = uniqueN(d3o$PATIENT_ID),
    estimate = 100 * (exp(e3) - 1), lo = 100 * (exp(e3 - 1.96 * s3o) - 1), hi = 100 * (exp(e3 + 1.96 * s3o) - 1),
    p = 2 * pnorm(-abs(e3 / s3o)), effect_type = "Pct",
    effect = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                     100 * (exp(e3) - 1), 100 * (exp(e3 - 1.96 * s3o) - 1), 100 * (exp(e3 + 1.96 * s3o) - 1))))

  ## M4
  ## Recheck in the overlap-region population: whether the 4 adjustment variables still vary across individuals
  adj_ok2 <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
    vv <- dw2[[v]]
    if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
    else uniqueN(as.character(vv[!is.na(vv)])) > 1
  }, logical(1))]
  dw2[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
  dw2[, lg_dur := log(anes_duration_mins)]
  d4o <- dw2[complete.cases(dw2[, c("n_epi", "Group", adj_ok2, "lg_dur", "sw"), with = FALSE])]
  m4o <- geeglm(as.formula(paste("n_epi ~ Group +", paste(adj_ok2, collapse = " + "), "+ offset(lg_dur)")),
                data = d4o, id = sid_int, weights = d4o$sw,
                family = poisson(), corstr = "exchangeable")
  so <- summary(m4o)$coefficients
  e4 <- so["GroupCiprofol_only", "Estimate"]; s4o <- so["GroupCiprofol_only", "Std.err"]
  result_ov <- rbind(result_ov, data.table(
    analysis = "PS overlap region only", model = "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)",
    unit = "Procedure (one anesthesia)", n_unit = nrow(d4o), n_patient = uniqueN(d4o$PATIENT_ID),
    estimate = exp(e4), lo = exp(e4 - 1.96 * s4o), hi = exp(e4 + 1.96 * s4o),
    p = 2 * pnorm(-abs(e4 / s4o)), effect_type = "RR",
    effect = sprintf("%.2f (%.2f-%.2f)", exp(e4), exp(e4 - 1.96 * s4o), exp(e4 + 1.96 * s4o))))

  print(result_ov[, .(model, unit, effect, p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])
  fwrite(result_ov, file.path(OUT, "Analysis_Results_overlap.csv"))

  tryCatch({
    suppressPackageStartupMessages({ library(flextable); library(officer) })
    ft <- flextable(as.data.frame(result_ov[, .(Model = model, Unit = unit, Estimate = effect,
                                                P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])) |>
         theme_booktabs() |> autofit() |>
         set_caption(caption = paste("Outcomes after IPTW within the overlap region only",
                                     "(PS re-fitted in this subset). All clustered by surgical procedure."))
    print(read_docx() |> body_add_flextable(ft),
          target = file.path(OUT, "Analysis_Results_overlap.docx"))
    ft2 <- flextable(as.data.frame(balance_overlap)) |> theme_booktabs() |> autofit() |>
           set_caption(caption = "Covariate balance before and after IPTW within the overlap region only.")
    print(read_docx() |> body_add_flextable(ft2),
          target = file.path(OUT, "Table1_IPTW_overlap_Balance_SMD.docx"))
  }, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))
}

## ============================================================================
## Step 9: weighted incidence in each group and the incidence difference (absolute effect size)
##   Three quantities per outcome - weighted ciprofol incidence, weighted propofol incidence, and their difference
##   The 95% CI of the difference is computed by bootstrap: resample "patients" (not surgeries) randomly with replacement 2000 times
## ============================================================================

cat("\n==================== Absolute effect estimates (weighted risks + bootstrap)====================\n")
ocs <- c("prolonged", "prolonged_cum", "multi_epi", "short_only", "any_hypo")

## Weights are also attached to the episode table to compute the proportion of episodes >10 min (same specification as M1)
epw_abs <- epi[dw[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, sw)],
               on = ID_VARS, nomatch = 0]
epw_abs[, long_epi := as.integer(duration > 10)]

boot_rd <- function(x, oc, wvar, B = 2000, seed = 20260928) {
  w  <- x[[wvar]]
  ok <- !is.na(x[[oc]]) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]

  y   <- as.numeric(x[[oc]])
  isc <- as.integer(x$Group == "Ciprofol_only")
  sid <- as.integer(factor(paste(x$PATIENT_ID, x$VISIT_ID, x$OPER_ID))); np <- max(sid)

  ## First sum the total weight and the weighted outcome within each surgery (one surgery is the analysis unit)
  per <- data.table(sid, isc, w, y)[, .(W = sum(w), Y = sum(w * y)), by = .(sid, isc)]
  CW <- numeric(np); CY <- numeric(np); PW <- numeric(np); PY <- numeric(np)
  CW[per$sid[per$isc == 1]] <- per$W[per$isc == 1]; CY[per$sid[per$isc == 1]] <- per$Y[per$isc == 1]
  PW[per$sid[per$isc == 0]] <- per$W[per$isc == 0]; PY[per$sid[per$isc == 0]] <- per$Y[per$isc == 0]

  ## Incidence difference in the observed sample
  est <- sum(CY) / sum(CW) - sum(PY) / sum(PW)

  ## bootstrap: resample patients 2000 times, recomputing the incidence difference each time
  set.seed(seed)
  reps <- numeric(B)
  for (b in seq_len(B)) {
    mult <- tabulate(sample.int(np, np, replace = TRUE), np)
    sc <- sum(mult * CW); sp <- sum(mult * PW)
    reps[b] <- if (sc > 0 && sp > 0) sum(mult * CY) / sc - sum(mult * PY) / sp else NA_real_
  }
  reps <- reps[is.finite(reps)]
  c(est = est, lo = as.numeric(quantile(reps, .025)), hi = as.numeric(quantile(reps, .975)))
}

## Procedure-level outcome (one row = one surgery)
abs_rows <- list()
for (oc in ocs) {
  a  <- boot_rd(dw, oc, wvar = "sw")
  cw <- dw$Group == "Ciprofol_only"; pw <- dw$Group == "Propofol_only"
  abs_rows[[length(abs_rows) + 1]] <- data.table(
    Outcome = oc, Unit = "Procedure (one anesthesia)",
    Ciprofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(dw[[oc]][cw], dw$sw[cw])),
    Propofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(dw[[oc]][pw], dw$sw[pw])),
    Risk_difference      = sprintf("%+.1f pp", 100 * a["est"]),
    Risk_difference_95CI = sprintf("(%.1f to %.1f)", 100 * a["lo"], 100 * a["hi"]),
    Ciprofol_n = sum(cw), Propofol_n = sum(pw))
}
abs_tab <- rbindlist(abs_rows)

## Episode-level outcome: proportion of episodes >10 min (denominator = hypotensive episodes that occurred; this is not the same as the procedure-level incidence, and the two cannot be compared directly)
cw2 <- epw_abs$Group == "Ciprofol_only"; pw2 <- epw_abs$Group == "Propofol_only"
a_e <- boot_rd(epw_abs, "long_epi", wvar = "sw")
abs_epi <- data.table(
  Outcome = "long_epi", Unit = "Episode (one hypotensive episode)",
  Ciprofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(epw_abs$long_epi[cw2], epw_abs$sw[cw2])),
  Propofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(epw_abs$long_epi[pw2], epw_abs$sw[pw2])),
  Risk_difference      = sprintf("%+.1f pp", 100 * a_e["est"]),
  Risk_difference_95CI = sprintf("(%.1f to %.1f)", 100 * a_e["lo"], 100 * a_e["hi"]),
  Ciprofol_n = sum(cw2), Propofol_n = sum(pw2))
abs_tab <- rbind(abs_epi, abs_tab)   # Put the episode-level row on top

print(abs_tab)
fwrite(abs_tab, file.path(OUT, "Analysis_Results_absolute_main.csv"))

tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(abs_tab)) |> theme_booktabs() |> autofit() |>
        set_caption(caption = paste("Weighted absolute risks and risk differences, full cohort.",
          "First row is episode-level (M1 definition); remaining rows are surgical-procedure level.",
          "Cluster bootstrap by surgical procedure."))
  print(read_docx() |> body_add_flextable(ft),
        target = file.path(OUT, "Analysis_Results_absolute_main.docx"))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))


## ============================================================================
## Step 10: same as Step 9, but done only within the overlap region
## ============================================================================

if (exists("dw2") && !is.null(dw2$sw)) {
  cat("\n==================== Absolute effect estimates (PS overlap region only)====================\n")
  abs_rows_b <- list()
  for (oc in ocs) {
    a  <- boot_rd(dw2, oc, wvar = "sw")
    cw <- dw2$Group == "Ciprofol_only"; pw <- dw2$Group == "Propofol_only"
    abs_rows_b[[length(abs_rows_b) + 1]] <- data.table(
      Outcome = oc, Unit = "Procedure (one anesthesia)",
      Ciprofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(dw2[[oc]][cw], dw2$sw[cw])),
      Propofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(dw2[[oc]][pw], dw2$sw[pw])),
      Risk_difference      = sprintf("%+.1f pp", 100 * a["est"]),
      Risk_difference_95CI = sprintf("(%.1f to %.1f)", 100 * a["lo"], 100 * a["hi"]),
      Ciprofol_n = sum(cw), Propofol_n = sum(pw))
  }
  abs_tab_b <- rbindlist(abs_rows_b)

  epb <- epi[dw2[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, sw)],
             on = ID_VARS, nomatch = 0]
  epb[, long_epi := as.integer(duration > 10)]
  cw3 <- epb$Group == "Ciprofol_only"; pw3 <- epb$Group == "Propofol_only"
  a_eb <- boot_rd(epb, "long_epi", wvar = "sw")
  abs_epi_b <- data.table(
    Outcome = "long_epi", Unit = "Episode (one hypotensive episode)",
    Ciprofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(epb$long_epi[cw3], epb$sw[cw3])),
    Propofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(epb$long_epi[pw3], epb$sw[pw3])),
    Risk_difference      = sprintf("%+.1f pp", 100 * a_eb["est"]),
    Risk_difference_95CI = sprintf("(%.1f to %.1f)", 100 * a_eb["lo"], 100 * a_eb["hi"]),
    Ciprofol_n = sum(cw3), Propofol_n = sum(pw3))
  abs_tab_b <- rbind(abs_epi_b, abs_tab_b)

  print(abs_tab_b)
  fwrite(abs_tab_b, file.path(OUT, "Analysis_Results_absolute_overlap.csv"))

  tryCatch({
    suppressPackageStartupMessages({ library(flextable); library(officer) })
    ft <- flextable(as.data.frame(abs_tab_b)) |> theme_booktabs() |> autofit() |>
          set_caption(caption = paste("Weighted absolute risks and risk differences within the",
            "overlap region only. First row is episode-level (M1 definition);",
            "remaining rows are surgical-procedure level. Cluster bootstrap by surgical procedure."))
    print(read_docx() |> body_add_flextable(ft),
          target = file.path(OUT, "Analysis_Results_absolute_overlap.docx"))
  }, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))
}


## ============================================================================
## Step 11: E-value - how strong an "unmeasured confounder" would be needed to overturn our conclusion
##   The larger the value, the more robust the conclusion. The algorithm is pure arithmetic, two lines are enough
## ============================================================================

cat("\n==================== E-value (for the M1 primary outcome)====================\n")
## E-value formula: v + sqrt(v * (v - 1)); if v<1, take the reciprocal first
evalue_of <- function(v) { if (v < 1) v <- 1 / v; v + sqrt(v * (v - 1)) }

m1_row <- result_main[model %like% "M1"]
or_pt <- m1_row$estimate                        # OR point estimate
or_ci <- if (or_pt < 1) m1_row$hi else m1_row$lo  # The end of the confidence interval closest to 1

cat(sprintf("  M1 point estimate: OR = %.2f → RR ≈ %.3f → E-value = %.2f\n",
            or_pt, sqrt(or_pt), evalue_of(sqrt(or_pt))))
cat(sprintf("  M1 CI limit      : OR = %.2f → RR ≈ %.3f → E-value = %.2f\n",
            or_ci, sqrt(or_ci), evalue_of(sqrt(or_ci))))
cat(sprintf("  (under the loose assumption of OR≈RR, the two E-values are %.2f / %.2f)\n",
            evalue_of(or_pt), evalue_of(or_ci)))


## ============================================================================
## Step 12: combine the main results table and the absolute risk table into one Word document
## ============================================================================

tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft1 <- flextable(as.data.frame(result_main[, .(
          Model = model, Unit = unit, Sample_size = n_unit, n_patients = n_patient,
          Estimate = effect, Type = effect_type,
          P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])) |>
       theme_booktabs() |> autofit() |>
       set_caption(caption = paste(
         "Table 5. Outcomes after IPTW weighting (M1-M3: one hypotensive episode per row;",
         "M4: one surgical procedure per row; all clustered by surgical procedure)"))
  ft2 <- flextable(as.data.frame(abs_tab)) |> theme_booktabs() |> autofit() |>
       set_caption(caption = paste(
         "Table 6. Weighted absolute risks and risk differences (cluster bootstrap by surgical procedure)"))
  doc <- read_docx() |> body_add_flextable(ft1) |> body_add_par("") |> body_add_flextable(ft2)
  print(doc, target = file.path(OUT, "Analysis_Results_main.docx"))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

