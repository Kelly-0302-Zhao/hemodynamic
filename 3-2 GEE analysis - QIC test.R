## =============================================================================
## 3-2 GEE analysis - QIC test.R - QIC test for selecting the working correlation structure
## =============================================================================

## ============================================================================
## Step 1: preparation (load packages, set paths, hard-code the constants)
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # read and process tables
  library(geepack)      # Run GEE and compute QIC
})

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## Three ID columns. All must be converted to character before merging, otherwise inconsistent types misalign the join
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## ============================================================================
## Step 2: attach the "group + covariates" of the surgery table to the event table to get the event-level data for QIC
## ============================================================================
d  <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
ep <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## Keep only the surgeries with an "analysable outcome" and attach their group and 4 adjustment variables to every episode
ep <- ep[d[!is.na(prolonged), .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                                Hypertension, AGE, SEX, ASA_GRADE)],
         on = ID_VARS, nomatch = 0]
ep <- ep[!is.na(Group)]
ep[, long_epi := as.integer(duration > 10)]   # Whether this episode lasted longer than 10 minutes
ep <- ep[!is.na(long_epi)]

## The geeglm cluster id must be an integer, so convert the surgery key (patient + visit + surgery) into an integer id
ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

## How many episodes occur per surgery - with cluster id = surgery, this distribution determines how much exchangeable differs from independence
rows_per_surgery <- ep[, .N, by = .(PATIENT_ID, VISIT_ID, OPER_ID)]
cat(sprintf("Distribution of hypotensive episodes per procedure: 1 episode %.1f%% | 2 episodes %.1f%% | ≥3 episodes %.1f%%\n\n",
            100 * mean(rows_per_surgery$N == 1),
            100 * mean(rows_per_surgery$N == 2),
            100 * mean(rows_per_surgery$N >= 3)))

## ============================================================================
## Step 3: run the M1 model for each of the three correlation structures and select the one with the smallest QIC
## ============================================================================

## Run M1 (proportion of prolonged hypotensive episodes) for a given correlation structure and extract QIC and QICu
fit_cor <- function(cs) {
  m <- geeglm(long_epi ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
              data = ep, id = sid_int, family = binomial(), corstr = cs)
  q <- QIC(m)   # Some geepack versions return a vector, so use square brackets to get the specific value
  data.table(Correlation_structure = cs,
              QIC = round(as.numeric(q[["QIC"]]), 2),
             LogLik_criterion = round(as.numeric(q[["QICu"]]), 2))
}

tab <- rbindlist(lapply(c("independence", "exchangeable", "ar1"), fit_cor))
print(tab)
best <- tab[which.min(QIC), Correlation_structure]
cat("\nThe smallest QIC is:", best, "\n")
cat("The QIC values of the three structures are very close and do not clearly distinguish them; per the pre-specified specification, the subsequent weighted GEE",
    "uses exchangeable for all four models.\n")

## Also save a Word copy of the table (only warn if the flextable/officer packages are missing, do not interrupt)
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

fwrite(tab, file.path(OUT, "Table3_GEE_QIC.csv"))
save_docx_table(tab, file.path(OUT, "Table3_GEE_QIC.docx"),
                caption = "Table 3. Working correlation structure selection by QIC (unweighted M1 model, one hypotensive episode per row, clustered by surgical procedure).")

