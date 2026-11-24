## =============================================================================
## 3-5 GEE analysis -  compute the propensity score (PS) and IPTW weights
## =============================================================================

## ============================================================================
## Step 1: preparation (load packages, set paths, hard-code a few constants shared by the whole file)
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # read and process tables
  library(geepack)      # This script only builds the weights; geeglm is used in 3-6, so load it here in advance
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
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")   # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

## Three ID columns (patient / visit / surgery).
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## Parse all timestamps as UTC (raw data carry a Z suffix)
MY_TZ <- "UTC"

## The 10 baseline variables used for the propensity score.
## The main analysis does not adjust for center (Hospital) or year (year): structural differences across centers/years are checked in the sensitivity analysis.
PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR",
             "operation_type", "Hypertension", "Diabetes", "Coronary_artery_disease",
             "Heart_failure")

## Adjustment variables shared across GEE models
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## Weight truncation quantiles: within each group, cut extreme weights at the 1% / 99% quantiles (see the pitfall note in Step 4)
TRUNC <- c(.01, .99)

## ============================================================================
## Step 2: read the procedure-level analysis table and report how many procedures and patients the full cohort contributed
## ============================================================================

## Read the procedure-level analysis table produced by 3-1 (one row = one surgery)
p_tab <- file.path(OUT, "Final_Analysis_Long_Table.rds")
if (!file.exists(p_tab)) stop("Final_Analysis_Long_Table.rds not found; please run 3-1 first.")
d <- setDT(readRDS(p_tab))

## Drop procedures with missing "prolonged hypotension" data
n_before <- nrow(d)
d <- d[!is.na(prolonged)]

cat("====================================================================\n")
cat("3-5 Propensity score and IPTW weights\n")
cat(sprintf("   Read %d rows from the surgery table → drop %d rows with missing prolonged → %d rows enter the analysis\n",
            n_before, n_before - nrow(d), nrow(d)))
cat(sprintf("   Full-cohort analysis population: %d procedures | %d patients\n",
            nrow(d), uniqueN(d$PATIENT_ID)))

## ============================================================================
## Step 3: flag each patient as treatment group or not, and select the variables eligible for the propensity score
## ============================================================================

## Treatment group = 1 (ciprofol), control group = 0 (propofol)
d[, treat := as.integer(Group == "Ciprofol_only")]

## Check each variable: for numeric variables, whether the standard deviation is > 0; for categorical variables, whether there is more than one level
use <- character(0)
for (v in PS_VARS) {
  if (!v %in% names(d)) next
  vv <- d[[v]]
  if (is.numeric(vv)) {
    if (is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0) use <- c(use, v)
  } else {
    if (uniqueN(as.character(vv[!is.na(vv)])) > 1) use <- c(use, v)
  }
}
dropped <- setdiff(PS_VARS, use)

cat(sprintf("   Variables included in the propensity score (%d): %s\n", length(use), paste(use, collapse = ", ")))
if (length(dropped))
  cat(sprintf("   Constant within this population and dropped automatically: %s\n", paste(dropped, collapse = ", ")))
cat(sprintf("   Treatment group (Ciprofol) %d procedures / Control group (Propofol) %d procedures\n",
            sum(d$treat == 1), sum(d$treat == 0)))

## ============================================================================
## Step 4: fit the propensity score, compute stabilized weights, truncate within group, flag the overlap region
## ============================================================================

## (1) complete-case subset
cc <- complete.cases(d[, ..use])
dw <- d[cc]
if (nrow(dw) < 100 || min(table(dw$treat)) < 10)
  stop("Too few procedures or one group too small; the propensity score cannot be fitted")

## (2) logistic regression to fit the PS
ps_fit <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
              data = dw, family = binomial())
dw[, ps_hat := predict(ps_fit, type = "response")]

## (3) stabilized weight: the true proportion of ciprofol use in the current population
prev <- mean(dw$treat)
dw[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

## (4) truncation is done [within group] (pitfall: truncating by overall quantiles flattens the treatment group)
dw[, `:=`(ql = as.numeric(quantile(sw_raw, TRUNC[1])),
          qh = as.numeric(quantile(sw_raw, TRUNC[2]))), by = Group]
dw[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
dw[, c("ql", "qh") := NULL]

## (5) common support (overlap region) flag
## AnYi hospital has never used ciprofol; the PS of its propofol cases falls close to 0, and IPTW can only give them weights of about 1, so they cannot be down-weighted and removed - this is an inherent limitation of IPTW and is only flagged here.
lo_overlap <- max(min(dw$ps_hat[dw$treat == 1]), min(dw$ps_hat[dw$treat == 0]))
hi_overlap <- min(max(dw$ps_hat[dw$treat == 1]), max(dw$ps_hat[dw$treat == 0]))
dw[, in_overlap := ps_hat >= lo_overlap & ps_hat <= hi_overlap]

cat(sprintf("   Fitting complete: %d procedures included in the weighted analysis (%d patients)\n",
            nrow(dw), uniqueN(dw$PATIENT_ID)))
cat(sprintf("   PS range  Propofol %.3f-%.3f | Ciprofol %.3f-%.3f\n",
            min(dw$ps_hat[dw$treat == 0]), max(dw$ps_hat[dw$treat == 0]),
            min(dw$ps_hat[dw$treat == 1]), max(dw$ps_hat[dw$treat == 1])))
cat(sprintf("   Stabilized weights (within-group truncation %.0f%%/%.0f%%) range: %.3f - %.1f\n",
            100 * TRUNC[1], 100 * TRUNC[2], min(dw$sw), max(dw$sw)))


## ============================================================================
## Step 5: compute the balance table before and after weighting (standardized mean difference SMD for each covariate)
## ============================================================================

## Weighted SMD of a single variable. Note that v must be numeric; categorical variables are first turned into 0/1 indicator variables before computation
wsmd <- function(v, w) {
  v <- as.numeric(v)
  m1 <- weighted.mean(v[dw$treat == 1], w[dw$treat == 1])
  m0 <- weighted.mean(v[dw$treat == 0], w[dw$treat == 0])
  v1 <- sum(w[dw$treat == 1] * (v[dw$treat == 1] - m1)^2) / sum(w[dw$treat == 1])
  v0 <- sum(w[dw$treat == 0] * (v[dw$treat == 0] - m0)^2) / sum(w[dw$treat == 0])
  sp <- sqrt((v1 + v0) / 2)
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  abs(m1 - m0) / sp
}

## For the selected variables, compute numeric ones directly; for categorical ones, compute once per level
bal_rows <- list()
for (v in use) {
  vv <- dw[[v]]
  if (is.numeric(vv)) {
    bal_rows[[length(bal_rows) + 1]] <- data.table(
      variable  = v,
      smd_before = wsmd(vv, rep(1, nrow(dw))),
      smd_after  = wsmd(vv, dw$sw))
  } else {
    for (lv in setdiff(unique(as.character(vv)), NA)) {
      ind <- as.integer(as.character(vv) == lv)
      bal_rows[[length(bal_rows) + 1]] <- data.table(
        variable  = paste0(v, ": ", lvl_en(lv)),
        smd_before = wsmd(ind, rep(1, nrow(dw))),
        smd_after  = wsmd(ind, dw$sw))
    }
  }
}
bal <- rbindlist(bal_rows)

cat("\n---- Covariate balance table (SMD: before weighting → after weighting)----\n")
print(bal)
cat(sprintf("   Max SMD before weighting %.3f (%d variables >0.1) → max SMD after weighting %.3f (%d variables)\n",
            max(bal$smd_before, na.rm = TRUE), sum(bal$smd_before > .1, na.rm = TRUE),
            max(bal$smd_after,  na.rm = TRUE), sum(bal$smd_after  > .1, na.rm = TRUE)))

fwrite(bal, file.path(OUT, "Table1_IPTW_Balance_SMD.csv"))


## ============================================================================
## Step 6: compute the effective sample size (ESS) table and print all diagnostic information
## ============================================================================

## Kish effective sample size
ess <- function(w) sum(w)^2 / sum(w^2)

## ESS table: by group, procedures / patients / median weight / max weight / ESS
ess_tab <- dw[, .(Group = as.character(first(Group)),
                   n_procedures = .N,
                   n_patients = uniqueN(PATIENT_ID),
                   weight_median = median(sw),
                   weight_max = max(sw),
                   ESS = ess(sw)), by = Group][, Group := NULL]

## The row order of ess_tab depends on the order in which groups appear, so values must be taken by group name, not by position
e_c <- ess_tab[Group == "Ciprofol_only", ESS]
e_p <- ess_tab[Group == "Propofol_only", ESS]

cat(sprintf("   [PS]   Procedures included in weighting %d (patients %d); PS range Ciprofol %.3f-%.3f, Propofol %.3f-%.3f\n",
            nrow(dw), uniqueN(dw$PATIENT_ID),
            min(dw$ps_hat[dw$treat == 1]), max(dw$ps_hat[dw$treat == 1]),
            min(dw$ps_hat[dw$treat == 0]), max(dw$ps_hat[dw$treat == 0])))
cat(sprintf("   [Balance] Max SMD before weighting %.3f (%d variables >0.1) → max SMD after weighting %.3f (%d variables)\n",
            max(bal$smd_before, na.rm = TRUE), sum(bal$smd_before > .1, na.rm = TRUE),
            max(bal$smd_after,  na.rm = TRUE), sum(bal$smd_after  > .1, na.rm = TRUE)))
cat(sprintf("   [Weight] Within-group truncation quantiles %.0f%%/%.0f%%; Ciprofol ESS %.0f / %d procedures, Propofol ESS %.0f / %d procedures\n",
            100 * TRUNC[1], 100 * TRUNC[2],
            e_c, sum(dw$treat == 1), e_p, sum(dw$treat == 0)))
cat(sprintf("   [Overlap] PS overlap region between groups = [%.3f, %.3f]; procedures outside the region %d (%.1f%%)\n",
            lo_overlap, hi_overlap, sum(!dw$in_overlap),
            100 * mean(!dw$in_overlap)))
if (sum(!dw$in_overlap) > 0 && "Hospital" %in% names(dw)) {
  cat("          Center distribution of procedures outside the region:\n")
  print(dw[in_overlap == FALSE, .N, by = .(Hospital, Group)][order(-N)])
}

cat("\n---- Weight description (by center)----\n")
print(dw[, .(n_procedures = .N, n_patients = uniqueN(PATIENT_ID),
             weight_median = round(median(sw), 2), weight_max = round(max(sw), 1),
             ESS = round(ess(sw))), by = .(Hospital, Group)][order(Hospital, Group)])

cat("\n---- Contribution of each center to the weighting (ESS/procedures)----\n")
print(dw[, .(n_procedures = .N, ESS = round(ess(sw)),
             ESS占比 = sprintf("%.1f%%", 100 * ess(sw) / .N)), by = Hospital][order(-n_procedures)])

cat("\nNote: AnYi hospital has only the Propofol group, so the ESS after weighting is very limited; the main results are driven mainly by the ZhongDa hospital and the NanYi hospital; this should be stated clearly in the main text.\n")

cat("\n---- Effective sample size (ESS) table ----\n")
print(ess_tab)
fwrite(ess_tab, file.path(OUT, "Table1_IPTW_ESS.csv"))


## ============================================================================
## Step 7: generate Table 1 (before / after weighting, comparison of baseline characteristics between the two groups)
## ============================================================================

## Lists of continuous and categorical variables (consistent with Table 1 of the paper)
CONT <- c("AGE", "BMI", "anes_duration_mins", "year")
CATS <- c("SEX", "ASA_GRADE", "operation_type", "EMERGENCY_INDICATOR",
          "Hypertension", "Diabetes", "Coronary_artery_disease", "Heart_failure",
          "Hospital")
LABELS <- c(AGE = "Age (yr)", BMI = "Body mass index (kg/m2)",
            anes_duration_mins = "Anesthesia duration (min)", year = "Year of surgery",
            SEX = "Sex", ASA_GRADE = "ASA physical status", operation_type = "Surgery type",
            EMERGENCY_INDICATOR = "Emergency surgery", Hypertension = "Hypertension",
            Diabetes = "Diabetes", Coronary_artery_disease = "Coronary artery disease",
            Heart_failure = "Heart failure", Hospital = "Center")

## Post-weighting mean(SD) strings (used internally in this step)
wmean_sd <- function(v, w) {
  m <- weighted.mean(v, w, na.rm = TRUE)
  s <- sqrt(sum(w * (v - m)^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
  sprintf("%.1f (%.1f)", m, s)
}

## Get the SMD 
pick_smd <- function(keys, col) {
  v <- bal[variable %in% keys, get(col)]
  if (length(v)) as.numeric(v[1]) else NA_real_
}
rows <- list()
for (v in CONT) {
  if (!v %in% names(dw)) next
  a <- dw[Group == "Ciprofol_only"]; b <- dw[Group == "Propofol_only"]
  rows[[length(rows) + 1]] <- data.table(
    Variable = LABELS[[v]], Level = "",
    Ciprofol_unweighted = sprintf("%.1f (%.1f)", mean(a[[v]], na.rm = TRUE), sd(a[[v]], na.rm = TRUE)),
    Propofol_unweighted = sprintf("%.1f (%.1f)", mean(b[[v]], na.rm = TRUE), sd(b[[v]], na.rm = TRUE)),
    Ciprofol_weighted   = wmean_sd(a[[v]], a$sw),
    Propofol_weighted   = wmean_sd(b[[v]], b$sw),
    SMD_before = round(pick_smd(v, "smd_before"), 3),
    SMD_after = round(pick_smd(v, "smd_after"), 3))
}
for (v in CATS) {
  if (!v %in% names(dw)) next
  lv <- levels(factor(dw[[v]]))
  for (l in lv) {
    a <- dw[Group == "Ciprofol_only"]; b <- dw[Group == "Propofol_only"]
    ind_a <- as.integer(as.character(a[[v]]) == l)
    ind_b <- as.integer(as.character(b[[v]]) == l)
    key <- paste0(v, ": ", lvl_en(l))
    rows[[length(rows) + 1]] <- data.table(
      Variable = LABELS[[v]], Level = lvl_en(l),
      Ciprofol_unweighted = sprintf("%d (%.1f%%)", sum(ind_a), 100 * mean(ind_a)),
      Propofol_unweighted = sprintf("%d (%.1f%%)", sum(ind_b), 100 * mean(ind_b)),
      ## After weighting report percentages only: after IPTW the "weighted procedures" is a fraction and displaying it would be misread as a real count of procedures
      Ciprofol_weighted   = sprintf("%.1f%%", 100 * weighted.mean(ind_a, a$sw)),
      Propofol_weighted   = sprintf("%.1f%%", 100 * weighted.mean(ind_b, b$sw)),
      ## First look up "variable: level" (true categorical variables); if not found, look up by variable name (0/1 comorbidities)
      SMD_before = round(pick_smd(c(key, v), "smd_before"), 3),
      SMD_after = round(pick_smd(c(key, v), "smd_after"), 3))
  }
}
tab1 <- rbindlist(rows, fill = TRUE)

cat(sprintf("\nTable header: Ciprofol %d procedures / %d patients; Propofol %d procedures / %d patients\n",
            dw[Group == "Ciprofol_only", .N], uniqueN(dw[Group == "Ciprofol_only"]$PATIENT_ID),
            dw[Group == "Propofol_only", .N], uniqueN(dw[Group == "Propofol_only"]$PATIENT_ID)))
cat("\n---- Table 1: baseline characteristics before / after weighting ---- \n")
print(tab1)
fwrite(tab1, file.path(OUT, "Table1_IPTW.csv"))


## ============================================================================
## Step 8: save Table 1 also as a Word version
## ============================================================================
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(tab1)) |> theme_booktabs() |> autofit() |>
    set_header_labels(Variable = "Variable", Level = "Level",
                      Ciprofol_unweighted = "Ciprofol, unweighted",
                      Propofol_unweighted = "Propofol, unweighted",
                      Ciprofol_weighted = "Ciprofol, IPTW-weighted",
                      Propofol_weighted = "Propofol, IPTW-weighted",
                      SMD_before = "SMD before", SMD_after = "SMD after") |>
    set_caption(caption = sprintf(
      "Table 1. Baseline characteristics before and after IPTW. Unit of analysis: surgical procedure (one anaesthesia episode). Ciprofol %d procedures (%d patients); Propofol %d procedures (%d patients).",
      dw[Group == "Ciprofol_only", .N], uniqueN(dw[Group == "Ciprofol_only"]$PATIENT_ID),
      dw[Group == "Propofol_only", .N], uniqueN(dw[Group == "Propofol_only"]$PATIENT_ID)))
  doc <- read_docx() |> body_add_flextable(ft)
  print(doc, target = file.path(OUT, "Table1_IPTW.docx"))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))


## ============================================================================
## Step 9: draw PS overlap diagnostic plots (overlaid PS distributions of the two groups + boxplot of the weight distribution)
## ============================================================================
png(file.path(OUT, "Figure_PS_Overlap.png"), width = 1400, height = 700, res = 150)
par(mfrow = c(1, 2))
hist(dw[Group == "Propofol_only", ps_hat], breaks = 40, col = rgb(.2, .4, .8, .5),
     border = "white", main = "倾向评分重叠（叠加权重前后）", xlab = "PS",
     ylim = c(0, max(table(cut(dw$ps_hat, 40)))))
hist(dw[Group == "Ciprofol_only", ps_hat], breaks = 40, col = rgb(.8, .3, .2, .5),
     border = "white", add = TRUE)
legend("topright", c("丙泊酚", "环泊酚"), fill = c(rgb(.2,.4,.8,.5), rgb(.8,.3,.2,.5)), bty = "n")
boxplot(list(Propofol = dw[Group == "Propofol_only", sw], Ciprofol = dw[Group == "Ciprofol_only", sw]),
        col = c(rgb(.2,.4,.8,.5), rgb(.8,.3,.2,.5)), main = "IPTW 权重分布", ylab = "权重")
dev.off()

## ============================================================================
## Step 10: bin the PS distribution (quantify the overlap region)
## ============================================================================

## Bin breaks
br <- c(0, 0.0338, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001)
lab <- c("<0.034", "0.034-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4",
         "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
dw[, ps_bin := cut(ps_hat, breaks = br, right = FALSE, labels = lab)]

## For each PS bin, count the procedures in each group and the percentage within each group
bin_tab <- function(x) {
  t <- dcast(x[, .N, by = .(ps_bin, Group)], ps_bin ~ Group, value.var = "N", fill = 0)
  setnames(t, c("Propofol_only", "Ciprofol_only"), c("Propofol_n", "Ciprofol_n"))
  t[, `:=`(Propofol_pct = round(100 * Propofol_n / sum(Propofol_n), 1),
           Ciprofol_pct = round(100 * Ciprofol_n / sum(Ciprofol_n), 1))]
  t[]
}

cat("\n---- PS bins for the full cohort (procedures + within-group proportion %)----\n")
tab_all <- bin_tab(dw)
print(tab_all)
fwrite(tab_all, file.path(OUT, "Check_PS_Distribution.csv"))
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(tab_all)) |> theme_booktabs() |> autofit() |>
    set_caption(caption = "Propensity score distribution by treatment group (full cohort).")
  print(read_docx() |> body_add_flextable(ft),
        target = file.path(OUT, "Check_PS_Distribution.docx"))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

cat("\n---- PS bins within the overlap region only (in_overlap == TRUE) ----\n")
tab_ov <- bin_tab(dw[in_overlap == TRUE])
print(tab_ov)
fwrite(tab_ov, file.path(OUT, "Check_PS_Distribution_overlap.csv"))
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(tab_ov)) |> theme_booktabs() |> autofit() |>
    set_caption(caption = "Propensity score distribution within the overlap region only.")
  print(read_docx() |> body_add_flextable(ft),
        target = file.path(OUT, "Check_PS_Distribution_overlap.docx"))
}, error = function(e) cat("[Word output failed, CSV saved]", conditionMessage(e), "\n"))

cat("\nKey numbers (deciding which population the main analysis should report):\n")
cat(sprintf("  Among %d Propofol procedures, %d (%.1f%%) have PS<0.034 (a range Ciprofol could never reach)\n",
            dw[Group == "Propofol_only", .N],
            dw[Group == "Propofol_only" & !in_overlap, .N],
            100 * dw[Group == "Propofol_only" & !in_overlap, .N] /
              dw[Group == "Propofol_only", .N]))
cat(sprintf("  Among %d Ciprofol procedures, only %d (%.1f%%) have PS<0.034\n",
            dw[Group == "Ciprofol_only", .N],
            dw[Group == "Ciprofol_only" & !in_overlap, .N],
            100 * dw[Group == "Ciprofol_only" & !in_overlap, .N] /
              dw[Group == "Ciprofol_only", .N]))

## ============================================================================
## Step 11: save the weight table
## ============================================================================
saveRDS(dw[, .(PATIENT_ID, VISIT_ID, OPER_ID, ps_hat, sw, sw_raw, in_overlap,
               Group, Hospital, year)],
        file.path(OUT, "IPTW_Weights.rds"))
