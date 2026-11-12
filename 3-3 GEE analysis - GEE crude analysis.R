## =============================================================================
## 3-3 GEE analysis - unweighted (crude) comparison
## =============================================================================
## ============================================================================
## Step 1: preparation
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # read and process tables
  library(geepack)      # run GEE models
})
# Set up the file output path
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## Three ID columns. Convert them all to character before merging
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 4 baseline variables adjusted for in the GEE model
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## Parse all timestamps as UTC
MY_TZ <- "UTC"

## ============================================================================
## Step 2: read the two tables and prepare the "unweighted analysis population" (weights all set to 1)
## ============================================================================
d <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
d <- d[!is.na(prolonged)]          # Keep only the surgeries with complete "prolonged hypotension" information

epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## Unweighted: the weight column sw is always 1 (so the subsequent model code is identical to the weighted version)
d[, sw := 1]

cat("Included:", nrow(d), "procedures |", uniqueN(d$PATIENT_ID), "patients\n\n")

## ============================================================================
## Step 3: model M1 - proportion of hypotension episodes lasting longer than 10 minutes, comparison between the two groups
##   Data: event-level (one row = one episode), weight = 1 means unweighted; report OR
## ============================================================================
## Attach the group, adjustment variables, and weights from the surgery table to every episode -> ep: one row = one episode
ep <- epi[d[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                 Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
          on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

## The geeglm id must be an integer, so convert the surgery key (patient + visit + surgery) into an integer id
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
result_crude <- data.table(
  analysis  = "Unweighted (crude)",
  model     = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
  unit      = "Episode (one hypotensive episode)",
  n_unit    = nrow(d1),                       # how many episodes
  n_patient = uniqueN(d1$PATIENT_ID),         # how many patients are involved
  estimate  = exp(est1), lo = exp(est1 - 1.96 * se1), hi = exp(est1 + 1.96 * se1),
  p         = p1,
  effect_type = "OR",
  effect    = sprintf("%.2f (%.2f-%.2f)", exp(est1), exp(est1 - 1.96 * se1), exp(est1 + 1.96 * se1)))

cat("M1 proportion of episodes >10 min  OR =", result_crude$effect,
    "  P =", sprintf("%.3g", p1), "\n")
cat("   （", nrow(d1), "episodes, involving", uniqueN(d1$PATIENT_ID), "patients)\n")


## ============================================================================
## Step 4: model M2 - severity score (AUC) of each hypotension episode, comparison between the two groups
##   AUC = the total "depth × duration" area of this episode, measuring how severe this episode is overall;
##   Take the log before comparing and convert the result back to "percent change". Report percent change
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
result_crude <- rbind(result_crude, data.table(
  analysis  = "Unweighted (crude)",
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

cat("M2 per-episode AUC change              ", result_crude[model %like% "M2", effect],
    "  P =", sprintf("%.3g", p2), "\n")

## ============================================================================
## Step 5: model M3 - duration of each hypotension episode in minutes, comparison between the two groups
##   Same approach as M2, except that the outcome is changed to episode duration; report results as percent change
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

result_crude <- rbind(result_crude, data.table(
  analysis  = "Unweighted (crude)",
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

cat("M3 per-episode duration change         ", result_crude[model %like% "M3", effect],
    "  P =", sprintf("%.3g", p3), "\n")


## ============================================================================
## Step 6: model M4 - number of hypotension episodes per procedure, comparison between the two groups
##   Data: surgery-level (one row = one surgery), outcome n_epi (number of episodes); add offset(log anesthesia duration),
##   The model actually compares "how many per hour", i.e. the "rate" of episodes. Report RR
## ============================================================================
adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  vv <- d[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]

d[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
d[, lg_dur  := log(anes_duration_mins)]   # Log-transform anesthesia duration and use it as an offset

## data.table with=FALSE only accepts a character vector, so here the column names are assembled into a character vector
d4 <- d[complete.cases(d[, c("n_epi", "Group", adj_ok, "lg_dur", "sw"), with = FALSE])]

f4 <- as.formula(paste("n_epi ~ Group +", paste(adj_ok, collapse = " + "), "+ offset(lg_dur)"))
m4 <- geeglm(f4, data = d4, id = sid_int, weights = d4$sw,
             family = poisson(), corstr = "exchangeable")

s4   <- summary(m4)$coefficients
est4 <- s4["GroupCiprofol_only", "Estimate"]
se4  <- s4["GroupCiprofol_only", "Std.err"]
p4   <- 2 * pnorm(-abs(est4 / se4))

result_crude <- rbind(result_crude, data.table(
  analysis  = "Unweighted (crude)",
  model     = "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)",
  unit      = "Procedure (one anesthesia)",
  n_unit    = nrow(d4), n_patient = uniqueN(d4$PATIENT_ID),
  estimate  = exp(est4), lo = exp(est4 - 1.96 * se4), hi = exp(est4 + 1.96 * se4),
  p         = p4,
  effect_type = "RR",
  effect    = sprintf("%.2f (%.2f-%.2f)", exp(est4), exp(est4 - 1.96 * se4), exp(est4 + 1.96 * se4))))

cat("M4 episode rate (per hour)         RR =", result_crude[model %like% "M4", effect],
    "  P =", sprintf("%.3g", p4), "\n")


## ============================================================================
## Step 7: save the results of the four unweighted models (csv + Word each) and print the descriptive statistics
## ============================================================================
fwrite(result_crude, file.path(OUT, "Crude_Table_Ciprofol.csv"))

## Save a Word version
save_docx_table <- function(df, path, caption = NULL) {
  ok <- tryCatch({
    suppressPackageStartupMessages({ library(flextable); library(officer) })
    ft <- flextable(as.data.frame(df)) |> theme_booktabs() |> autofit()
    if (!is.null(caption)) ft <- set_caption(ft, caption = caption)
    doc <- read_docx() |> body_add_flextable(ft)
    print(doc, target = path)
    TRUE
  }, error = function(e) {
    cat("  [Word output failed, CSV saved] ", conditionMessage(e), "\n")
    FALSE
  })
  invisible(ok)
}
save_docx_table(result_crude, file.path(OUT, "Crude_Table_Ciprofol.docx"),
                caption = "Unweighted (crude) GEE models for reference only. M1-M3: one hypotensive episode per row; M4: one surgical procedure per row; all clustered by surgical procedure.")

cat("\n==================== Unweighted (crude) results ====================\n")
print(result_crude[, .(model, unit, n_unit, n_patient, type = effect_type,
                        effect, p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])

## Group-wise description (surgery-level)
cat("\n---- Unweighted group description (procedure level)----\n")
print(d[, .(手术例次 = .N,
            患者数 = uniqueN(PATIENT_ID),
            长期低血压率 = sprintf("%.1f%%", 100 * mean(prolonged)),
            累积时长中位 = sprintf("%.1f min", median(MAP65_TIME)),
            累积AUC中位 = sprintf("%.1f", median(MAP65_AUC)),
            事件数每例次 = sprintf("%.2f", mean(n_epi))), by = Group])
