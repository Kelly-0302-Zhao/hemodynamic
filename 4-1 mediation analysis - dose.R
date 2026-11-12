## =============================================================================
## 4-1 mediation analysis - dose.R
## -----------------------------------------------------------------------------
## Mediation analysis of the induction dose (g-formula / g-computation, risk difference scale, parametric bootstrap 2000 times)
## =============================================================================

suppressPackageStartupMessages({ library(data.table) })

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## Baseline confounders shared by the mediator model and the outcome model (excluding postoperative/post-exposure variables other than the mediator itself)
CONF <- c("AGE", "SEX", "BMI", "ASA_GRADE", "Hypertension", "Diabetes",
          "Coronary_artery_disease", "Heart_failure", "EMERGENCY_INDICATOR",
          "operation_type", "anes_duration_mins")
## 11 baseline variables used for the propensity score
PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR", "operation_type",
             "Hypertension", "Diabetes", "Coronary_artery_disease", "Heart_failure",
             "anes_duration_mins")

## =============================================================================
## Re-fit the propensity score on the given data + stabilized weights + 1%/99% truncation
## =============================================================================
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
  ps <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
            data = dt, family = binomial())
  dt[, ps_hat := predict(ps, type = "response")]
  prev <- mean(dt$treat)
  dt[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]
  dt[, ql := as.numeric(quantile(sw_raw, 0.01)), by = Group]
  dt[, qh := as.numeric(quantile(sw_raw, 0.99)), by = Group]
  dt[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
  dt[, c("ql", "qh", "sw_raw") := NULL]
  dt
}

## =============================================================================
## g-formula mediation analysis (VanderWeele's g-computation)
##   med_col specifies the column name of the mediator (log1p first, then standardize)
##   Mediator model (exposure -> mediator): linear weighting
##   Outcome model (with mediator): weighted logistic
##   Four counterfactual means (risk difference scale):
##     TE  = ψ(1,m1) − ψ(0,m0)   total effect
##     NIE = ψ(1,m1) − ψ(1,m0)   natural indirect effect (via the mediator)
##     NDE = ψ(1,m0) − ψ(0,m0)   natural direct effect (not via the mediator)
##     IIE = ψ(1,m1) − ψ(0,m1)   interaction-induced effect
##   95%CI: parametric bootstrap - draw one set of coefficients from the coefficient covariance matrix of each of the two models,
##   Recompute the counterfactual means, repeat B times, and take the 2.5%/97.5% percentiles.
## =============================================================================
mvrnorm1 <- function(mu, V) {
  V <- (V + t(V)) / 2
  ev <- eigen(V, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = length(ev$values))
  out <- mu + as.numeric(L %*% rnorm(length(mu)))
  names(out) <- names(mu)
  out
}

gf_mediation <- function(dat, med_col, fe_col = NULL, B = 2000,
                         seed = 20261106, label = "", mediator_name = "") {
  dt <- copy(dat)
  conf_use <- intersect(CONF, names(dt))
  dt[, m_raw := log1p(get(med_col))]
  sd_m <- sd(dt$m_raw, na.rm = TRUE)
  dt[, m_z := (m_raw - mean(m_raw, na.rm = TRUE)) / sd_m]
  keepvars <- c("prolonged", "Group", "m_z", conf_use)
  dt <- dt[complete.cases(dt[, ..keepvars])]
  dt[, treat := as.integer(Group == "Ciprofol_only")]

  fconf <- paste(conf_use, collapse = " + ")
  feT   <- if (!is.null(fe_col)) paste0(" + factor(", fe_col, ")") else ""
  f_med <- as.formula(paste("m_z ~ treat +", fconf, feT))
  f_out <- as.formula(paste("prolonged ~ treat + m_z +", fconf, feT))

  mm <- lm(f_med, data = dt, weights = dt$sw)
  mo <- glm(f_out, data = dt, family = binomial(), weights = dt$sw)

  ## ---- Mediator model: split out X_med (drop the treat column); replace treat with a for the counterfactual ----
  df_med <- model.matrix(f_med, data = dt); cn_med <- colnames(df_med)
  i_tm <- which(cn_med == "treat")
  X_med <- df_med[, -i_tm, drop = FALSE]
  b_med <- coef(mm)
  a_treat <- unname(b_med[i_tm]); b_med_rest <- b_med[-i_tm]

  ## ---- Outcome model: X_out drops treat / m_z; the remaining columns are aligned with b_rest ----
  df_out <- model.matrix(f_out, data = dt); cn_out <- colnames(df_out)
  i_treat_out <- which(cn_out == "treat")
  i_mz_out    <- which(cn_out == "m_z")
  drop_cols  <- c(i_treat_out, i_mz_out)
  X_out      <- df_out[, -drop_cols, drop = FALSE]
  b_out      <- coef(mo)
  b_treat    <- unname(b_out[i_treat_out])
  b_mz       <- unname(b_out[i_mz_out])
  b_rest     <- b_out[-drop_cols]

  lo <- min(dt$m_z); hi <- max(dt$m_z)

  ## ---- Four counterfactual means ----
  psi_set <- function(bmr, at, br, bt, bm) {
    m1 <- pmin(pmax(as.numeric(X_med %*% bmr) + at, lo), hi)
    m0 <- pmin(pmax(as.numeric(X_med %*% bmr) + 0,  lo), hi)
    lin <- function(a, m) plogis(as.numeric(X_out %*% br) + bt * a + bm * m)
    c(TE  = mean(lin(1, m1)) - mean(lin(0, m0)),
      NIE = mean(lin(1, m1)) - mean(lin(1, m0)),
      NDE = mean(lin(1, m0)) - mean(lin(0, m0)),
      IIE = mean(lin(1, m1)) - mean(lin(0, m1)),
      psi11 = mean(lin(1, m1)), psi00 = mean(lin(0, m0)))
  }

  pt <- psi_set(b_med_rest, a_treat, b_rest, b_treat, b_mz)

  ## ---- Parametric bootstrap ----
  set.seed(seed)
  bm_n <- names(b_med); bo_n <- names(b_out)
  med_rest_names <- setdiff(bm_n, "treat")
  out_rest_names <- setdiff(bo_n, c("treat", "m_z"))
  Vm <- vcov(mm); Vo <- vcov(mo)
  bt_mat <- matrix(NA_real_, nrow = B, ncol = 5,
                   dimnames = list(NULL, c("TE","NIE","NDE","IIE","prop")))
  for (b in seq_len(B)) {
    bms <- mvrnorm1(b_med, Vm); bos <- mvrnorm1(b_out, Vo)
    v <- psi_set(bms[med_rest_names], unname(bms["treat"]), bos[out_rest_names],
                 unname(bos["treat"]), unname(bos["m_z"]))
    bt_mat[b, 1:4] <- v[1:4]
    bt_mat[b, 5]   <- if (abs(v["TE"]) > 1e-10) v["NIE"] / v["TE"] else NA_real_
  }

  q <- function(x) unname(quantile(x, c(.025, .975), na.rm = TRUE))
  prop_pt  <- if (abs(pt[["TE"]]) > 1e-10) pt[["NIE"]] / pt[["TE"]] else NA_real_
  prop_lo  <- q(bt_mat[, "prop"])[1]; prop_hi <- q(bt_mat[, "prop"])[2]
  prop_ciw <- prop_hi - prop_lo

  data.table(
    mediator     = mediator_name,
    spec         = label,
    n_surgery    = nrow(dt),
    n_cipro      = sum(dt$treat),
    n_prop       = sum(dt$treat == 0),
    a_exposure_to_mediator = round(unname(coef(mm)["treat"]), 4),
    a_p           = unname(summary(mm)$coefficients["treat", 4]),
    b_mediator_to_outcome  = round(b_mz, 4),
    b_p           = unname(summary(mo)$coefficients["m_z", 4]),
    TE_RD         = round(pt[["TE"]],  5),
    TE_lo         = round(q(bt_mat[, "TE"])[1],  5), TE_hi = round(q(bt_mat[, "TE"])[2],  5),
    NIE_RD        = round(pt[["NIE"]], 5),
    NIE_lo        = round(q(bt_mat[, "NIE"])[1],  5), NIE_hi = round(q(bt_mat[, "NIE"])[2], 5),
    NDE_RD        = round(pt[["NDE"]], 5),
    NDE_lo        = round(q(bt_mat[, "NDE"])[1],  5), NDE_hi = round(q(bt_mat[, "NDE"])[2], 5),
    IIE_RD        = round(pt[["IIE"]], 5),
    IIE_lo        = round(q(bt_mat[, "IIE"])[1],  5), IIE_hi = round(q(bt_mat[, "IIE"])[2], 5),
    or_total_model = round(exp(unname(coef(mo)["treat"])), 3),
    or_total_lo    = round(exp(unname(confint(mo)["treat", 1])), 3),
    or_total_hi    = round(exp(unname(confint(mo)["treat", 2])), 3),
    prop_mediated = round(prop_pt, 3), prop_lo = round(prop_lo, 3), prop_hi = round(prop_hi, 3),
    interpretable = is.finite(prop_lo) && is.finite(prop_hi) && prop_ciw <= 2.0,
    verdict = if (!is.finite(prop_lo) || !is.finite(prop_hi))
      "中介比例不可解释: 总效应趋零导致比值发散（请以 NIE/NDE 的绝对值为准）"
    else if (prop_ciw > 2.0)
      "中介比例不可解释: 95%CI 跨度过宽（>2.0），间接效应与总效应无法分离"
    else if (prop_lo <= 0 && prop_hi >= 0)
      "中介比例 95%CI 包含 0"
    else if (abs(prop_pt) < 0.05)
      sprintf("中介比例实质为零(%.1f%%): 幅度 <5%%", 100 * prop_pt)
    else "存在显著中介"
  )
}

## =============================================================================
## Main workflow
## =============================================================================
cat("====================================================================\n")
cat("4-1 Mediation analysis of the induction dose (propofol equivalent ind_equiv)\n")
cat("(g-formula, risk difference scale, parametric bootstrap with 2000 replications)\n")
cat("====================================================================\n")

dw <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
for (v in ID_VARS) set(dw, j = v, value = as.character(dw[[v]]))
dw <- dw[!is.na(prolonged)]

## ---- dose (procedure level, specification v12)----
## The dose table is organized by (surgery, drug) => first pick "the drug corresponding to the group that this surgery belongs to", then LEFT JOIN it back
dose_file <- "_dose_metrics_v12.csv"
dose_path <- file.path(OUT, dose_file)
if (!file.exists(dose_path)) {
  cand <- list.files(OUT, pattern = "^_dose_metrics_v[0-9]+\\.csv$")
  if (!length(cand)) stop("Dose table _dose_metrics_v*.csv not found")
  dose_file <- cand[order(as.integer(sub("^_dose_metrics_v([0-9]+)\\.csv$", "\\1", cand)), decreasing = TRUE)][1]
  dose_path <- file.path(OUT, dose_file)
  cat(sprintf("!! %s does not exist, falling back to %s\n", "_dose_metrics_v12.csv", dose_file))
}
cat(sprintf("Dose table: %s\n", dose_file))
dose <- fread(dose_path)
if (nrow(dose) < 80000) stop(sprintf("Unexpected number of rows in the dose table (%d); the wrong file may have been selected", nrow(dose)))
for (v in ID_VARS) set(dose, j = v, value = as.character(dose[[v]]))
cat(sprintf("Dose table %d rows; drug distribution:", nrow(dose)))
print(dose[, .N, by = drug])

dw[, want_drug := fifelse(Group == "Ciprofol_only", "ciprofol", "propofol")]
dose_pick <- merge(dw[, c(..ID_VARS, "want_drug"), with = FALSE], dose,
                   by = ID_VARS, all.x = TRUE, sort = FALSE)
dose_pick <- dose_pick[drug == want_drug]
dose_pick <- unique(dose_pick, by = ID_VARS)
dw <- merge(dw, dose_pick[, c(..ID_VARS, "ind_mgkg", "inf_mgkg_h", "ind_mgkg_adm"), with = FALSE],
            by = ID_VARS, all.x = TRUE, sort = FALSE)
dw <- unique(dw, by = ID_VARS)

## ---- potency ratio -> propofol equivalent---
RATIO <- median(dw[Group == "Propofol_only" & !is.na(ind_mgkg), ind_mgkg]) /
         median(dw[Group == "Ciprofol_only"  & !is.na(ind_mgkg), ind_mgkg])
dw[, ind_equiv := ind_mgkg * fifelse(Group == "Ciprofol_only", RATIO, 1)]
cat(sprintf("\nPotency ratio (median ratio of induction dose mg/kg, Propofol:Ciprofol) = %.4f (literature 4.1:1)\n", RATIO))
cat("Between-group distribution of the induction dose (raw value vs propofol equivalent):\n")
print(rbind(
  dw[!is.na(ind_mgkg), .(口径 = "raw value mg/kg", 台数 = .N,
                         中位 = round(median(ind_mgkg), 3),
                         q25 = round(quantile(ind_mgkg, .25), 3),
                         q75 = round(quantile(ind_mgkg, .75), 3)), by = Group],
  dw[!is.na(ind_equiv), .(口径 = "propofol equivalent mg/kg", 台数 = .N,
                          中位 = round(median(ind_equiv), 3),
                          q25 = round(quantile(ind_equiv, .25), 3),
                          q75 = round(quantile(ind_equiv, .75), 3)), by = Group]))

dw[, op_year := as.numeric(as.character(year))]
dw[, cell := paste(Hospital, op_year, sep = "_")]

cat(sprintf("\nFull cohort: %d procedures\n", nrow(dw)))
cat(sprintf("  Induction dose recorded in %d procedures (%.1f%%)  Ciprofol %d (%.1f%%) / Propofol %d (%.1f%%)\n",
            sum(!is.na(dw$ind_mgkg)), 100 * mean(!is.na(dw$ind_mgkg)),
            sum(!is.na(dw$ind_mgkg) & dw$Group == "Ciprofol_only"),
            100 * mean(!is.na(dw$ind_mgkg) & dw$Group == "Ciprofol_only"),
            sum(!is.na(dw$ind_mgkg) & dw$Group == "Propofol_only"),
            100 * mean(!is.na(dw$ind_mgkg) & dw$Group == "Propofol_only")))

## =============================================================================
## Main analysis: induction dose (propofol equivalent) mediation
## =============================================================================
cat("\n\n############ Induction dose (propofol-equivalent) mediation ############\n")
dose_sub <- dw[!is.na(ind_equiv) & !is.na(prolonged)]
cat(sprintf("Dose subset: %d procedures (ciprofol %d / propofol %d), covering %d centre x year cells\n",
            nrow(dose_sub), sum(dose_sub$Group == "Ciprofol_only"),
            sum(dose_sub$Group == "Propofol_only"), uniqueN(dose_sub$cell)))
cat("Centre x year cells where ciprofol appears in this subset:",
    paste(sort(unique(dose_sub[Group == "Ciprofol_only"]$cell)), collapse = ", "), "\n")

## ---- Cells with both arms present: center × year cells with >=20 procedures each ----
find_both_cells <- function(sub) {
  cc <- sub[, .N, by = .(cell, Group)]
  piv <- dcast(cc, cell ~ Group, value.var = "N", fill.value = 0)
  for (cn in c("Ciprofol_only", "Propofol_only"))
    if (!cn %in% names(piv)) piv[, (cn) := 0L]
  piv[, Ciprofol_only := as.integer(as.character(Ciprofol_only))]
  piv[, Propofol_only  := as.integer(as.character(Propofol_only))]
  piv[is.na(Ciprofol_only), Ciprofol_only := 0L]
  piv[is.na(Propofol_only),  Propofol_only  := 0L]
  piv[, within := pmin(Ciprofol_only, Propofol_only)]
  piv[within >= 20]$cell
}

both_dose <- find_both_cells(dose_sub)
cat(sprintf("\nCentre x year cells with both arms present and >=20 procedures each (dose): %s\n", paste(both_dose, collapse = ", ")))

ds    <- fit_iptw(dose_sub)
ds_co <- fit_iptw(dose_sub[cell %in% both_dose])

res_ds <- rbindlist(list(
  gf_mediation(ds, med_col = "ind_equiv", fe_col = NULL, B = 2000, seed = 20261110,
               label = "① 粗口径", mediator_name = "诱导剂量(丙泊酚当量)"),
  gf_mediation(ds, med_col = "ind_equiv", fe_col = "cell", B = 2000, seed = 20261111,
               label = "② ①+中心×年FE", mediator_name = "诱导剂量(丙泊酚当量)"),
  gf_mediation(ds_co, med_col = "ind_equiv", fe_col = NULL, B = 2000, seed = 20261112,
               label = "③ 两臂共存格子", mediator_name = "诱导剂量(丙泊酚当量)"),
  gf_mediation(ds_co, med_col = "ind_equiv", fe_col = "cell", B = 2000, seed = 20261113,
               label = "④ ③+中心×年FE", mediator_name = "诱导剂量(丙泊酚当量)")
), fill = TRUE)
cat("\n---- Induction dose - four specifications ----\n")
print(res_ds[, .(mediator, spec, n_surgery, n_cipro, n_prop,
                 a_exposure_to_mediator, b_mediator_to_outcome,
                 TE_RD, TE_lo, TE_hi, NIE_RD, NIE_lo, NIE_hi,
                 NDE_RD, NDE_lo, NDE_hi,
                 prop_mediated, prop_lo, prop_hi, or_total_model, verdict)])

## =============================================================================
## Auxiliary diagnostics
## =============================================================================
cat("\n---- Auxiliary diagnostic 1: exposure x mediator interaction test (is the no-interaction assumption valid) ----\n")
for (spec_lab in c("粗口径", "①+中心×年FE")) {
  dat <- copy(ds)
  fc <- if (spec_lab == "粗口径") NULL else "cell"
  conf_use <- intersect(CONF, names(dat))
  dat[, m_z := as.numeric(scale(log1p(ind_equiv)))]
  kv <- c("prolonged", "Group", "m_z", conf_use)
  dat <- dat[complete.cases(dat[, ..kv])]
  dat[, treat := as.integer(Group == "Ciprofol_only")]
  fconf <- paste(conf_use, collapse = " + ")
  feT <- if (!is.null(fc)) paste0(" + factor(", fc, ")") else ""
  mo_i <- glm(as.formula(paste("prolonged ~ treat * m_z +", fconf, feT)),
              data = dat, family = binomial(), weights = dat$sw)
  cat(sprintf("   [%s] treat:m_z = %.4f (SE %.4f, p = %.3g)\n", spec_lab,
              unname(coef(mo_i)["treat:m_z"]),
              unname(summary(mo_i)$coefficients["treat:m_z", 2]),
              unname(summary(mo_i)$coefficients["treat:m_z", 4])))
}

cat("\n---- Auxiliary diagnostic 2: attenuation of the a path (exposure -> equivalent dose) ----\n")
da <- copy(dose_sub)
da[, treat := as.integer(Group == "Ciprofol_only")]
da[, m_z := as.numeric(scale(log1p(ind_equiv)))]
d_crude <- lm(m_z ~ treat, data = da)
d_fe    <- lm(m_z ~ treat + factor(cell), data = da)
a_crude <- unname(coef(d_crude)["treat"]); a_fe <- unname(coef(d_fe)["treat"])
cat(sprintf("   Crude a = %.4f (SE %.4f, p=%.3g)\n", a_crude,
            summary(d_crude)$coefficients["treat", 2], summary(d_crude)$coefficients["treat", 4]))
cat(sprintf("   With FE a = %.4f (SE %.4f, p=%.3g)  ← attenuation %.0f%%\n", a_fe,
            summary(d_fe)$coefficients["treat", 2], summary(d_fe)$coefficients["treat", 4],
            100 * (1 - a_fe / a_crude)))
cat("   => After conversion to propofol-equivalent units the two dose distributions almost overlap => the a path is inherently tiny,\n")
cat("     Dose mediation identifies the residual pharmacological difference AT AN EQUIVALENT DOSE, not a difference in dose quantity.\n")

cat("\n---- Auxiliary diagnostic 3: scale artefact of the raw mg/kg specification (control) ----\n")
d_raw <- fit_iptw(dw[!is.na(ind_mgkg) & !is.na(prolonged)])
rr <- gf_mediation(d_raw, med_col = "ind_mgkg", fe_col = NULL, B = 2000, seed = 20261114,
                   label = "对照", mediator_name = "诱导剂量(原始数值mg/kg·尺度假象)")
print(rr[, .(mediator, spec, a_exposure_to_mediator, b_mediator_to_outcome,
             TE_RD, NIE_RD, NIE_lo, NIE_hi, prop_mediated, prop_lo, prop_hi)])
cat("   WARNING: an a path of magnitude -1.6 in that row means drug identity was treated as a dose channel - a scale artefact,\n")
cat("     it must not be cited as a cross-drug conclusion; the only reportable specification is the propofol-equivalent one.\n")

## ---- Delivery table (4 rows = four specifications)----
write_atomic <- function(tab, fname) {
  tmp <- file.path(OUT, paste0("_stag_", fname))
  fwrite(tab, tmp)
  dst <- file.path(OUT, fname)
  ok <- suppressWarnings(file.copy(tmp, dst, overwrite = TRUE))
  if (ok) { unlink(tmp); cat("\nWritten to:", dst, "\n") }
  else {
    cat("\n!! Target file is locked (may be open in Excel/WPS) and could not be overwritten:", dst, "\n")
    cat("!! Results saved to a temporary file:", tmp, "\n")
  }
}
write_atomic(res_ds, "Mediation_dose.csv")
cat(sprintf("\nTable rows: %d (1 mediator x 4 specifications)\n", nrow(res_ds)))
