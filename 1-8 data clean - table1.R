## =============================================================================
## 1-8 data clean - Table 1 (baseline characteristics of the two groups + SMD + P values)
## =============================================================================
## ---- (i) Are the anesthesia start and end times usable ----
.t_s <- data1$ANES_START_TIME
.t_e <- data1$ANES_END_TIME
.bad_time <- is.na(.t_s) | is.na(.t_e)

## ---- (ii) Are there measurement points within the intraoperative window ----
## Split monitor1 into blocks of 5 million rows, match each block with the "start and end time of each procedure", and accumulate the number of measurement points within the window.
.bp <- data.table(
  PATIENT_ID = as.character(data1$PATIENT_ID),
  VISIT_ID   = as.character(data1$VISIT_ID),
  OPER_ID    = as.character(data1$OPER_ID),
  s = .t_s, e = .t_e,
  .row = seq_len(nrow(data1)),      # Used to map the results back in the original order
  n_win = 0L)
setkey(.bp, PATIENT_ID, VISIT_ID, OPER_ID)

.chunk <- 5e6L
.n_all <- nrow(monitor1)
for (.st in seq.int(1L, .n_all, by = .chunk)) {
  .en  <- min(.st + .chunk - 1L, .n_all)
  .blk <- data.table(
    PATIENT_ID = as.character(monitor1$PATIENT_ID[.st:.en]),
    VISIT_ID   = as.character(monitor1$VISIT_ID[.st:.en]),
    OPER_ID    = as.character(monitor1$OPER_ID[.st:.en]),
    t          = monitor1$TIME_POINT[.st:.en])
  setkey(.blk, PATIENT_ID, VISIT_ID, OPER_ID)
  ## How many measurement points of each procedure fall within the window in this block (by = .EACHI -> one row per procedure)
  .hit <- .blk[.bp, on = c("PATIENT_ID", "VISIT_ID", "OPER_ID"), nomatch = NULL,
               .(n = sum(t >= i.s & t <= i.e, na.rm = TRUE)), by = .EACHI]
  if (nrow(.hit)) {
    .bp[.hit, on = c("PATIENT_ID", "VISIT_ID", "OPER_ID"), n_win := n_win + i.n]
  }
}
setorder(.bp, .row)                 # Sort back into the original order of data1
.has_win <- .bp$n_win > 0L          # At least 1 measurement point within the window

## ---- Summary: who was excluded and why ----
.time_ok <- !.bad_time
.keep    <- .time_ok & .has_win

cat("\n---- Baseline table population screening----\n")
cat(sprintf("  Cleaned full cohort                        : %d procedures\n", nrow(data1)))
cat(sprintf("  Excluded: missing anesthesia times         : %d procedures\n", sum(!.time_ok)))
cat(sprintf("  Excluded: no measurement points in window  : %d procedures (among procedures with usable anesthesia times)\n",
            sum(.time_ok & !.has_win)))
cat(sprintf("  Total excluded                             : %d procedures\n", sum(!.keep)))
cat(sprintf("  Entered baseline table (denominator)       : %d procedures, %d patients\n",
            sum(.keep), length(unique(data1$PATIENT_ID[.keep]))))

## Write the exclusion list to disk for direct use in the manuscript flow diagram
.excl <- data.table(
  PATIENT_ID = as.character(data1$PATIENT_ID[!.keep]),
  VISIT_ID   = as.character(data1$VISIT_ID[!.keep]),
  OPER_ID    = as.character(data1$OPER_ID[!.keep]),
  Hospital   = as.character(data1$Hospital[!.keep]),
  Group      = as.character(data1$Group[!.keep]),
  exclude_reason = ifelse(!.time_ok[!.keep],
                          "anes_time_missing", "no_intraop_monitor"))
fwrite(.excl, "Check_Baseline_Excluded.csv", encoding = "UTF-8")

rm(.t_s, .t_e, .bad_time, .bp, .blk, .hit, .st, .en, .n_all, .chunk)
invisible(gc())

#Build the summary statistics table and output table1
data <- data1[.keep, ]
group_stats <- data %>% group_by(Group) %>% summarise(n = n(), .groups = "drop")
cipro_n <- group_stats$n[group_stats$Group == "Ciprofol_only"]
propo_n <- group_stats$n[group_stats$Group == "Propofol_only"]
#Variable definitions and label settings
continuous_vars <- c("AGE", "BMI","anes_duration_mins")
categorical_vars <- c("SEX", "ASA_GRADE", "operation_type",  "EMERGENCY_INDICATOR")
# Risk factor mapping (customized sublevel names)
risk_mapping <- c(
  "Diabetes" = "Diabetes mellitus (%)",
  "Hypertension" = "Chronic arterial hypertension (%)",
  "Coronary_artery_disease" = "Coronary artery disease (%)",
  "Heart_failure" = "Heart failure (%)")
# Labels for routine variables
var_labels <- c(
  "AGE" = "Age (SD), yr",
  "BMI" = "Body mass index, kg/m²",
  "SEX" = "Male sex (%)",
  "ASA_GRADE" = "ASA grade (%)",
  "operation_type" = "Operation type (%)",
  "anes_duration_mins" = "Anesthesia duration, min",
  "EMERGENCY_INDICATOR" = "Emergency (%)")
# Compute SMD (standardized mean difference)
## This function: custom SMD (standardized mean difference) - pooled SD for continuous variables, p-based for categorical variables, max/mean for multi-category variables
compute_smd <- function(var, data, group_var = "Group",
                        g1 = "Ciprofol_only", g2 = "Propofol_only",
                        multicat = c("max", "mean")) {
  multicat <- match.arg(multicat)
  x <- data[[var]]
  g <- trimws(as.character(data[[group_var]]))
  # Keep only the two groups to avoid out-of-range errors when a table row does not exist
  keep <- g %in% c(g1, g2)
  x <- x[keep]
  g <- g[keep]
  if (length(g) == 0 || any(!c(g1, g2) %in% unique(g))) return(NA_real_)
  if (is.numeric(x)) {
    # ===== SMD for continuous variables =====
    x1 <- x[g == g1]
    x2 <- x[g == g2]
    mean1 <- mean(x1, na.rm = TRUE); mean2 <- mean(x2, na.rm = TRUE)
    sd1 <- sd(x1, na.rm = TRUE);     sd2 <- sd(x2, na.rm = TRUE)
    pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
    return(if (is.na(pooled_sd) || pooled_sd == 0) NA_real_
           else abs(mean1 - mean2) / pooled_sd)
  } else {
    # ===== Standard SMD for categorical variables =====
    xf <- factor(x)             
    levs <- levels(xf)
    if (length(levs) < 2) return(NA_real_)
    x1 <- xf[g == g1]
    x2 <- xf[g == g2]
    # Binary: p-based SMD (standard definition)
    if (length(levs) == 2) {
      p1 <- mean(x1 == levs[2], na.rm = TRUE)
      p2 <- mean(x2 == levs[2], na.rm = TRUE)
      denom <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
      return(if (is.na(denom) || denom == 0) NA_real_
             else abs(p1 - p2) / denom)}
    # Multi-category: compute the SMD with each level treated as a dummy and return max (or mean)
    smds <- vapply(levs, function(L) {
      p1 <- mean(x1 == L, na.rm = TRUE)
      p2 <- mean(x2 == L, na.rm = TRUE)
      denom <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
      if (is.na(denom) || denom == 0) NA_real_ else abs(p1 - p2) / denom
    }, numeric(1))
    smds <- smds[is.finite(smds)]
    if (length(smds) == 0) return(NA_real_)
    if (multicat == "max") return(max(smds))
    return(mean(smds))}}
all_vars <- c(continuous_vars, categorical_vars, names(risk_mapping))
smd_df <- data.frame(
  Variable = all_vars,
  SMD = round(sapply(all_vars, compute_smd, data = data), 2))
# Build the Baseline Table data frame
baseline_table <- data.frame(
  Variable=character(), Total=character(), Ciprofol_only=character(), 
  Propofol_only=character(), SMD=character(), P_value=character(), stringsAsFactors=F)

## ============================================================================
## Step 1 (A): Continuous variables - mean±SD + Kruskal-Wallis P value, filled into the baseline table row by row
## ============================================================================
for(var in continuous_vars){
  vec <- data[[var]]
  # Use format(..., nsmall=2) to force 2 decimal places and keep the alignment neat
  total_m <- format(round(mean(vec, na.rm=T), 2), nsmall=2)
  total_s <- format(round(sd(vec, na.rm=T), 2), nsmall=2)
  total_stat <- paste0(total_m, " (", total_s, ")")
  stats <- data %>% group_by(Group) %>% 
    summarise(m=mean(!!sym(var), na.rm=T), s=sd(!!sym(var), na.rm=T))
  # Format the ciprofol group
  c_m <- format(round(stats$m[stats$Group=="Ciprofol_only"], 2), nsmall=2)
  c_s <- format(round(stats$s[stats$Group=="Ciprofol_only"], 2), nsmall=2)
  cipro <- paste0(c_m, " (", c_s, ")")
  # Format the propofol group
  p_m <- format(round(stats$m[stats$Group=="Propofol_only"], 2), nsmall=2)
  p_s <- format(round(stats$s[stats$Group=="Propofol_only"], 2), nsmall=2)
  propo <- paste0(p_m, " (", p_s, ")")
  # The P value calculation logic stays unchanged
  p_val_raw <- kruskal.test(vec ~ Group, data=data)$p.value
  p_val <- if(p_val_raw < 0.001) "<0.001" else format(round(p_val_raw, 3), nsmall=3)
  # Fill in following the original data frame structure
  baseline_table <- rbind(baseline_table, data.frame(
    Variable=var_labels[var], 
    Total=total_stat, 
    Ciprofol_only=cipro, 
    Propofol_only=propo, 
    SMD=as.character(smd_df$SMD[smd_df$Variable==var]), 
    P_value=p_val ))}

## ============================================================================
## Step 2 (B): Comorbidity block - bold the "Baseline risk factors" category header, then fill in the percentage + Fisher P value for each disease
## ============================================================================
baseline_table <- rbind(baseline_table, data.frame(
  Variable="Baseline risk factors", Total="", Ciprofol_only="", Propofol_only="", SMD="", P_value=""))
for(var in names(risk_mapping)){
  tbl_df <- data %>% count(!!sym(var), Group) %>% pivot_wider(names_from=Group, values_from=n, values_fill=0)
  row_1 <- tbl_df %>% filter(!!sym(var) == 1)
  # Compute the summary statistic
  smd_val <- smd_df$SMD[smd_df$Variable == var]
  tbl_matrix <- as.matrix(tbl_df[,c("Ciprofol_only","Propofol_only")])
  p_val_raw <- if(all(dim(tbl_matrix)==2)) fisher.test(tbl_matrix)$p.value else NA
  p_val <- if(is.na(p_val_raw)) "N/A" else if(p_val_raw < 0.001) "<0.001" else format(round(p_val_raw,3), nsmall=3)
  if(nrow(row_1) > 0){
    baseline_table <- rbind(baseline_table, data.frame(
      Variable = paste0("    ", risk_mapping[var]), 
      Total = as.character(row_1$Ciprofol_only + row_1$Propofol_only),
      Ciprofol_only = paste0(row_1$Ciprofol_only, " (", round(row_1$Ciprofol_only/cipro_n*100, 2), "%)"),
      Propofol_only = paste0(row_1$Propofol_only, " (", round(row_1$Propofol_only/propo_n*100, 2), "%)"),
      SMD = as.character(smd_val), P_value = p_val)) }}

## ============================================================================
## Step 3 (C): other categorical variables (ASA / surgery type / emergency) - category heading + percentage for each level + chi-square P value
## ============================================================================
for(var in categorical_vars){
  tbl_df <- data %>% count(!!sym(var), Group) %>% pivot_wider(names_from=Group, values_from=n, values_fill=0)
  tbl_df <- tbl_df %>% mutate(Total = Ciprofol_only + Propofol_only)
  # Compute the P value for the main row
  p_val_raw <- chisq.test(as.matrix(tbl_df[,c("Ciprofol_only","Propofol_only")]))$p.value
  p_val <- if(p_val_raw < 0.001) "<0.001" else format(round(p_val_raw, 3), nsmall = 3)
  baseline_table <- rbind(baseline_table, data.frame(
    Variable=var_labels[var], Total="", Ciprofol_only="", Propofol_only="", 
    SMD=as.character(smd_df$SMD[smd_df$Variable==var]), P_value=p_val ))
  for(i in 1:nrow(tbl_df)){
    baseline_table <- rbind(baseline_table, data.frame(
      Variable = paste0("    ", tbl_df[[var]][i]),
      Total = as.character(tbl_df$Total[i]),
      Ciprofol_only = paste0(tbl_df$Ciprofol_only[i], " (", round(tbl_df$Ciprofol_only[i]/cipro_n*100, 2), "%)"),
      Propofol_only = paste0(tbl_df$Propofol_only[i], " (", round(tbl_df$Propofol_only[i]/propo_n*100, 2), "%)"),
      SMD="", P_value="" )) }}
# Format the output and export
# Identify the rows to bold: rows not starting with four spaces
bold_rows <- which(!grepl("^    ", baseline_table$Variable))
## ============================================================================
## Format (bold category headings, italic sub-items) + export Word and CSV
## ============================================================================
ft <- flextable(baseline_table) %>%
  theme_booktabs() %>%
  bold(i = bold_rows, j = "Variable") %>% 
  italic(i = which(grepl("^    ", baseline_table$Variable)), j = "Variable") %>%
  set_header_labels(
    Ciprofol_only = paste0("Ciprofol only\n(n=", cipro_n, ")"),
    Propofol_only = paste0("Propofol only\n(n=", propo_n, ")"),
    SMD = "SMD", P_value = "P value") %>%
  align(align = "center", part = "header") %>%
  autofit()
# Export Word
doc <- read_docx() %>%
  body_add_par("Table 1. Baseline characteristics of patients", style = "heading 1") %>%
  body_add_flextable(ft)
print(doc, target = "Baseline_Table.docx")
# Output CSV
write.csv(baseline_table,"baseline_table1.csv", row.names=FALSE)
