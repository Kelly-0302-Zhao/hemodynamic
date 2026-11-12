## =============================================================================
## 4-2 subgroup analysis - operation duration.R
## -----------------------------------------------------------------------------
## Stratification by procedure (anesthesia) duration + interaction analysis
## =============================================================================

suppressPackageStartupMessages({ library(data.table); library(geepack) })

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
DATA <- file.path(D, "JCA_run")
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")
PS_VARS  <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR", "operation_type",
              "Hypertension", "Diabetes", "Coronary_artery_disease", "Heart_failure")

## ============================================================================
## Reusable function: refit PS + stabilized weights + truncation within the subset
## ============================================================================
fit_iptw <- function(dt) {
  dt <- copy(dt)
  dt[, treat := as.integer(Group == "Ciprofol_only")]
  use <- character(0)
  for (v in PS_VARS) {
    if (!v %in% names(dt)) next
    vv <- dt[[v]]
    if (is.numeric(vv)) { if (is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0) use <- c(use, v) }
    else { if (uniqueN(as.character(vv[!is.na(vv)])) > 1) use <- c(use, v) }
  }
  dt <- dt[complete.cases(dt[, ..use])]
  ps <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))), data = dt, family = binomial())
  dt[, ps_hat := predict(ps, type = "response")]
  prev <- mean(dt$treat)
  dt[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]
  dt[, ql := as.numeric(quantile(sw_raw, 0.01)), by = Group]
  dt[, qh := as.numeric(quantile(sw_raw, 0.99)), by = Group]
  dt[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
  dt[, c("ql", "qh", "sw_raw") := NULL]
  dt
}

epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## ============================================================================
## Reusable function: run M1 (episode-level OR) on a set of already weighted procedures
## ============================================================================
run_m1 <- function(dw_sub, tag, n_stratum = NA_integer_) {
  ep <- epi[dw_sub[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
            on = ID_VARS, nomatch = NULL]
  ep <- ep[!is.na(Group) & is.finite(w) & w > 0]
  if (nrow(ep) < 30) return(NULL)
  ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
  ep[, long_epi := as.integer(duration > 10)]
  adj <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
    vv <- ep[[v]]
    if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
    else uniqueN(as.character(vv[!is.na(vv)])) > 1
  }, logical(1))]
  d1 <- ep[complete.cases(ep[, c("long_epi", "Group", adj, "w"), with = FALSE])]
  m <- geeglm(as.formula(paste("long_epi ~ Group +", paste(adj, collapse = " + "))),
              data = d1, id = d1$sid_int, weights = d1$w,
              family = binomial(), corstr = "exchangeable")
  cm <- as.matrix(summary(m)$coefficients)
  est <- cm["GroupCiprofol_only", "Estimate"]
  se  <- cm["GroupCiprofol_only", "Std.err"]
  data.table(stratum = tag, n_surgery = n_stratum, n_event = nrow(d1),
             OR = round(exp(est), 3),
             lo = round(exp(est - 1.96 * se), 3),
             hi = round(exp(est + 1.96 * se), 3),
             p = round(2 * pnorm(-abs(est / se)), 4))
}

## Joint Wald test: given a geeglm model and a set of "interaction" coefficient names, return W, df, p
wald_joint <- function(fit, terms) {
  V <- vcov(fit); b <- coef(fit)
  terms <- intersect(terms, names(b))
  if (length(terms) == 0) return(list(W = NA_real_, df = 0, p = NA_real_))
  bb <- b[terms]; VV <- V[terms, terms]
  W <- as.numeric(t(bb) %*% solve(VV) %*% bb)
  list(W = W, df = length(terms), p = pchisq(W, df = length(terms), lower.tail = FALSE))
}

## ============================================================================
## Read the table and build duration groups
## ============================================================================
d <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
for (v in ID_VARS) set(d, j = v, value = as.character(d[[v]]))
d <- d[!is.na(prolonged) & is.finite(anes_duration_mins) & anes_duration_mins > 0]

## Tertile cut points
q1 <- unname(quantile(d$anes_duration_mins, 1/3)); q2 <- unname(quantile(d$anes_duration_mins, 2/3))
cat(sprintf("\nAnaesthesia duration tertile cut-points: T1 <%s | T2 %s-%s | T3 >%s minutes (median %s)\n",
            round(q1), round(q1), round(q2), round(q2), round(median(d$anes_duration_mins))))
d[, dur_q := fcase(anes_duration_mins <= q1, "T1_short",
                   anes_duration_mins <= q2, "T2_mid",
                   default = "T3_long")]
d[, dur_q := factor(dur_q, levels = c("T1_short", "T2_mid", "T3_long"))]

## ============================================================================
## Part 1: stratified estimates (weights refit within each stratum)
## ============================================================================
strata_rows <- list()
for (lv in levels(d$dur_q)) {
  sub <- fit_iptw(d[dur_q == lv])
  r <- run_m1(sub, lv, n_stratum = nrow(d[dur_q == lv]))
  if (!is.null(r)) strata_rows[[length(strata_rows) + 1]] <- r
}
## Also add an overall row as a reference
strata_rows[[length(strata_rows) + 1]] <- run_m1(fit_iptw(d), "Overall", nrow(d))
strata_tab <- rbindlist(strata_rows)
cat("\n=== Stratification by duration tertile (M1 OR, re-weighted within stratum) ===\n")
print(strata_tab)
fwrite(strata_tab, file.path(OUT, "Subgroup_Duration_Strata.csv"))

## ============================================================================
## Part 2: interaction test (whole cohort, whole-cohort weights)
##   2a. Interaction term between categorical duration stratum and Group, joint Wald test
##   2b. Interaction between continuous log duration and Group (1 degree of freedom)
## ============================================================================
dwf <- fit_iptw(d)                       # Whole-cohort weights (for the interaction models)
epf <- epi[dwf[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, Hypertension, AGE, SEX, ASA_GRADE,
                   dur_q, anes_duration_mins, w = sw)], on = ID_VARS, nomatch = NULL]
epf <- epf[!is.na(Group) & is.finite(w) & w > 0]
epf[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
epf[, long_epi := as.integer(duration > 10)]
adjf <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  vv <- epf[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]
epf <- epf[complete.cases(epf[, c("long_epi", "Group", adjf, "dur_q", "w"), with = FALSE])]

## 2a Categorical interaction
m_int <- geeglm(as.formula(paste("long_epi ~ Group * dur_q +", paste(adjf, collapse = " + "))),
                data = epf, id = epf$sid_int, weights = epf$w,
                family = binomial(), corstr = "exchangeable")
inter_terms <- grep(":", names(coef(m_int)), value = TRUE)
w_int <- wald_joint(m_int, inter_terms)

## 2b Continuous interaction (log duration x Group; anes_duration_mins has already been joined into epf)
epf[, lg_dur := log(anes_duration_mins)]
epf <- epf[is.finite(lg_dur)]
m_cont <- geeglm(as.formula(paste("long_epi ~ Group * lg_dur +", paste(adjf, collapse = " + "))),
                 data = epf, id = epf$sid_int, weights = epf$w,
                 family = binomial(), corstr = "exchangeable")
w_cont <- wald_joint(m_cont, grep(":", names(coef(m_cont)), value = TRUE))

inter_tab <- rbind(
  data.table(test = "分组 × 时长三分位（分类交互）", terms = paste(inter_terms, collapse = "; "),
             Wald_W = round(w_int$W, 3), df = w_int$df, p = signif(w_int$p, 3),
             verdict = fifelse(w_int$p < 0.05, "交互显著：药效随时长变化，需按时长解读",
                               "交互不显著：没有证据表明药效随时长改变")),
  data.table(test = "分组 × log(时长)（连续交互）", terms = paste(grep(":", names(coef(m_cont)), value = TRUE), collapse = "; "),
             Wald_W = round(w_cont$W, 3), df = w_cont$df, p = signif(w_cont$p, 3),
             verdict = fifelse(w_cont$p < 0.05, "交互显著", "交互不显著：药效不随时长单调变化")))
cat("\n=== Interaction test (full cohort) ===\n")
print(inter_tab)
fwrite(inter_tab, file.path(OUT, "Subgroup_Duration_Interaction.csv"))

## ============================================================================
## Part 3: sensitivity analysis with a different set of clinical cut points (<120 / 120-240 / >=240 minutes)
## ============================================================================
d[, dur_clin := fcase(anes_duration_mins < 120, "A_<120",
                      anes_duration_mins < 240, "B_120-240",
                      default = "C_>=240")]
d[, dur_clin := factor(dur_clin, levels = c("A_<120", "B_120-240", "C_>=240"))]
clin_rows <- list()
for (lv in levels(d$dur_clin)) {
  sub <- fit_iptw(d[dur_clin == lv])
  r <- run_m1(sub, lv, n_stratum = nrow(d[dur_clin == lv]))
  if (!is.null(r)) clin_rows[[length(clin_rows) + 1]] <- r
}
clin_tab <- rbindlist(clin_rows)
cat("\n=== Clinical cut-point sensitivity (<120 / 120-240 / >=240 minutes) ===\n")
print(clin_tab)
fwrite(clin_tab, file.path(OUT, "Subgroup_Duration_Strata_ClinicalCut.csv"))
