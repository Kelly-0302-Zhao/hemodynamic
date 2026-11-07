## =============================================================================
## 1-8 data clean — Table 1（两组基线特征表 + SMD + P 值）
## =============================================================================
## ---- (i) 麻醉起止时间能不能用 ----
.t_s <- data1$ANES_START_TIME
.t_e <- data1$ANES_END_TIME
.bad_time <- is.na(.t_s) | is.na(.t_e)

## ---- (ii) 术中窗口内有没有测点 ----
## 把 monitor1 按 500 万行切开，逐块和"每台手术的起止时间"匹配，累加窗口内测点数。
.bp <- data.table(
  PATIENT_ID = as.character(data1$PATIENT_ID),
  VISIT_ID   = as.character(data1$VISIT_ID),
  OPER_ID    = as.character(data1$OPER_ID),
  s = .t_s, e = .t_e,
  .row = seq_len(nrow(data1)),      # 用来把结果按原顺序对回去
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
  ## 每台手术在这一块里有几个测点落在窗口内（by = .EACHI → 每台手术一行）
  .hit <- .blk[.bp, on = c("PATIENT_ID", "VISIT_ID", "OPER_ID"), nomatch = NULL,
               .(n = sum(t >= i.s & t <= i.e, na.rm = TRUE)), by = .EACHI]
  if (nrow(.hit)) {
    .bp[.hit, on = c("PATIENT_ID", "VISIT_ID", "OPER_ID"), n_win := n_win + i.n]
  }
}
setorder(.bp, .row)                 # 按 data1 的原顺序排回来
.has_win <- .bp$n_win > 0L          # 窗口内至少 1 个测点

## ---- 汇总：谁被筛掉、为什么 ----
.time_ok <- !.bad_time
.keep    <- .time_ok & .has_win

cat("\n---- 基线表人群筛选----\n")
cat(sprintf("  清洗后全队列              ：%d 台手术\n", nrow(data1)))
cat(sprintf("  剔除：麻醉起止时间缺失    ：%d 台\n", sum(!.time_ok)))
cat(sprintf("  剔除：窗口内无监测测点    ：%d 台（在上面时间可用的手术里）\n",
            sum(.time_ok & !.has_win)))
cat(sprintf("  合计剔除                  ：%d 台\n", sum(!.keep)))
cat(sprintf("  进入基线表（本表分母）    ：%d 台手术，%d 位患者\n",
            sum(.keep), length(unique(data1$PATIENT_ID[.keep]))))

## 剔除清单落盘，方便直接进论文流程图
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

#构建数据统计表，输出table1
data <- data1[.keep, ]
group_stats <- data %>% group_by(Group) %>% summarise(n = n(), .groups = "drop")
cipro_n <- group_stats$n[group_stats$Group == "Ciprofol_only"]
propo_n <- group_stats$n[group_stats$Group == "Propofol_only"]
#变量定义与标签设置
continuous_vars <- c("AGE", "BMI","anes_duration_mins")
categorical_vars <- c("SEX", "ASA_GRADE", "operation_type",  "EMERGENCY_INDICATOR")
# 风险因素映射关系（自定义子水平名称）
risk_mapping <- c(
  "Diabetes" = "Diabetes mellitus (%)",
  "Hypertension" = "Chronic arterial hypertension (%)",
  "Coronary_artery_disease" = "Coronary artery disease (%)",
  "Heart_failure" = "Heart failure (%)")
# 常规变量标签
var_labels <- c(
  "AGE" = "Age (SD), yr",
  "BMI" = "Body mass index, kg/m²",
  "SEX" = "Male sex (%)",
  "ASA_GRADE" = "ASA grade (%)",
  "operation_type" = "Operation type (%)",
  "anes_duration_mins" = "Anesthesia duration, min",
  "EMERGENCY_INDICATOR" = "Emergency (%)")
# 计算 SMD (标准化均数差)
## 这个函数：自定义算 SMD（标准化均数差）——连续变量用 pooled SD，分类变量用 p-based，多分类取 max/mean
compute_smd <- function(var, data, group_var = "Group",
                        g1 = "Ciprofol_only", g2 = "Propofol_only",
                        multicat = c("max", "mean")) {
  multicat <- match.arg(multicat)
  x <- data[[var]]
  g <- trimws(as.character(data[[group_var]]))
  # 只保留两组，避免 table 行不存在导致越界
  keep <- g %in% c(g1, g2)
  x <- x[keep]
  g <- g[keep]
  if (length(g) == 0 || any(!c(g1, g2) %in% unique(g))) return(NA_real_)
  if (is.numeric(x)) {
    # ===== 连续变量SMD =====
    x1 <- x[g == g1]
    x2 <- x[g == g2]
    mean1 <- mean(x1, na.rm = TRUE); mean2 <- mean(x2, na.rm = TRUE)
    sd1 <- sd(x1, na.rm = TRUE);     sd2 <- sd(x2, na.rm = TRUE)
    pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
    return(if (is.na(pooled_sd) || pooled_sd == 0) NA_real_
           else abs(mean1 - mean2) / pooled_sd)
  } else {
    # ===== 分类变量标准SMD =====
    xf <- factor(x)             
    levs <- levels(xf)
    if (length(levs) < 2) return(NA_real_)
    x1 <- xf[g == g1]
    x2 <- xf[g == g2]
    # 二分类：p-based SMD（标准定义）
    if (length(levs) == 2) {
      p1 <- mean(x1 == levs[2], na.rm = TRUE)
      p2 <- mean(x2 == levs[2], na.rm = TRUE)
      denom <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
      return(if (is.na(denom) || denom == 0) NA_real_
             else abs(p1 - p2) / denom)}
    # 多分类：每个level当dummy算SMD，返回max（或mean）
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
# 构建 Baseline Table 数据框
baseline_table <- data.frame(
  Variable=character(), Total=character(), Ciprofol_only=character(), 
  Propofol_only=character(), SMD=character(), P_value=character(), stringsAsFactors=F)

## ============================================================================
## 第 1 步（A）：连续变量 —— 均值±SD + Kruskal-Wallis P 值，逐行填入基线表
## ============================================================================
for(var in continuous_vars){
  vec <- data[[var]]
  # 使用 format(..., nsmall=2) 强制保留2位小数，确保对齐美观
  total_m <- format(round(mean(vec, na.rm=T), 2), nsmall=2)
  total_s <- format(round(sd(vec, na.rm=T), 2), nsmall=2)
  total_stat <- paste0(total_m, " (", total_s, ")")
  stats <- data %>% group_by(Group) %>% 
    summarise(m=mean(!!sym(var), na.rm=T), s=sd(!!sym(var), na.rm=T))
  # Ciprofol 组格式化
  c_m <- format(round(stats$m[stats$Group=="Ciprofol_only"], 2), nsmall=2)
  c_s <- format(round(stats$s[stats$Group=="Ciprofol_only"], 2), nsmall=2)
  cipro <- paste0(c_m, " (", c_s, ")")
  # Propofol 组格式化
  p_m <- format(round(stats$m[stats$Group=="Propofol_only"], 2), nsmall=2)
  p_s <- format(round(stats$s[stats$Group=="Propofol_only"], 2), nsmall=2)
  propo <- paste0(p_m, " (", p_s, ")")
  # P值计算逻辑保持不变
  p_val_raw <- kruskal.test(vec ~ Group, data=data)$p.value
  p_val <- if(p_val_raw < 0.001) "<0.001" else format(round(p_val_raw, 3), nsmall=3)
  # 按照原有数据框结构填充
  baseline_table <- rbind(baseline_table, data.frame(
    Variable=var_labels[var], 
    Total=total_stat, 
    Ciprofol_only=cipro, 
    Propofol_only=propo, 
    SMD=as.character(smd_df$SMD[smd_df$Variable==var]), 
    P_value=p_val ))}

## ============================================================================
## 第 2 步（B）：合并症块 —— 加粗"Baseline risk factors"大类标题，再逐病填百分比 + Fisher P 值
## ============================================================================
baseline_table <- rbind(baseline_table, data.frame(
  Variable="Baseline risk factors", Total="", Ciprofol_only="", Propofol_only="", SMD="", P_value=""))
for(var in names(risk_mapping)){
  tbl_df <- data %>% count(!!sym(var), Group) %>% pivot_wider(names_from=Group, values_from=n, values_fill=0)
  row_1 <- tbl_df %>% filter(!!sym(var) == 1)
  # 计算统计量
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
## 第 3 步（C）：其他分类变量（ASA / 手术类型 / 急诊）——大类标题 + 各水平百分比 + 卡方 P 值
## ============================================================================
for(var in categorical_vars){
  tbl_df <- data %>% count(!!sym(var), Group) %>% pivot_wider(names_from=Group, values_from=n, values_fill=0)
  tbl_df <- tbl_df %>% mutate(Total = Ciprofol_only + Propofol_only)
  # 计算主行 P 值
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
# 格式化输出与导出
# 确定需要加粗的行：不以四个空格开头的行
bold_rows <- which(!grepl("^    ", baseline_table$Variable))
## ============================================================================
## 格式化（加粗大标题、斜体子项）+ 导出 Word 与 CSV
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
# 导出 Word
doc <- read_docx() %>%
  body_add_par("Table 1. Baseline characteristics of patients", style = "heading 1") %>%
  body_add_flextable(ft)
print(doc, target = "Baseline_Table.docx")
# 输出 CSV
write.csv(baseline_table,"baseline_table1.csv", row.names=FALSE)
