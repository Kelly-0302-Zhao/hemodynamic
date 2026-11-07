## =============================================================================
## 3-6 GEE analysis - — 加权后的四个主模型（论文主结果）
## =============================================================================
## ============================================================================
## 第 1 步：准备工作
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # 读表、处理表
  library(geepack)      # 跑 GEE 模型
})

LEVEL_EN <- c("女" = "Female", "男" = "Male",
              "Ⅰ" = "I", "Ⅱ" = "II", "Ⅲ" = "III")
lvl_en <- function(x) {
  out <- as.character(x)
  hit <- LEVEL_EN[out]
  out[!is.na(hit)] <- hit[!is.na(hit)]
  out
}

## 根目录。所有输入输出都挂在这个目录下面
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
DATA <- file.path(D, "JCA_run")      # 清洗好的中间表放在这里（1-x 系列脚本的输出）
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## 三个 ID 列。
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## GEE 模型里调整的 4 个基线变量
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## 时间戳统一按 UTC 解析
MY_TZ <- "UTC"

## ============================================================================
## 第 2 步：读入三张表，拼出"加权分析人群"，报告纳入了多少例次、多少病人
## ============================================================================

## 手术表：一行 = 一次手术
d <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))

## 权重表：3-5 算好的倾向评分（ps_hat）和稳定化加权权重（sw）
w <- setDT(readRDS(file.path(OUT, "IPTW_Weights.rds")))

## 事件表：一行 = 一次低血压事件
epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## 三张表的 ID 列统一转成文字型（类型不一致时合并会错位，这是踩过的坑）
for (v in ID_VARS) {
  set(d,   j = v, value = as.character(d[[v]]))
  set(w,   j = v, value = as.character(w[[v]]))
  set(epi, j = v, value = as.character(epi[[v]]))
}

## 把权重挂到手术表上；再删掉"是否长期低血压"信息缺失的手术
dw <- d[w[, .(PATIENT_ID, VISIT_ID, OPER_ID, ps_hat, sw, in_overlap)],
        on = ID_VARS, nomatch = 0]
dw <- dw[!is.na(prolonged)]

cat("加权分析人群：手术", nrow(dw), "例次 | 患者", uniqueN(dw$PATIENT_ID), "人\n")
cat(sprintf("  （手术表全表 %d 例次；差额是倾向评分协变量缺失的手术）\n\n",
            d[!is.na(prolonged), .N]))


## ============================================================================
## 第 3 步：模型 M1 —— 每次低血压事件，持续超过 10 分钟的比例，两组比较
##   数据：事件表（一行 = 一次低血压事件），每个事件继承它所属手术的权重
##   模型：logistic 型 GEE，报告 OR（环泊酚 vs 丙泊酚）
## ============================================================================

## 把手术表的分组、权重、调整变量挂到每个事件上 → ep：一行 = 一次事件
ep <- epi[dw[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                 Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
          on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

## GEE 要求"聚类编号"必须是整数，这里把手术编号（患者+就诊+手术）转成整数编号
ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

## 因变量：这次事件持续超过 10 分钟记 1，否则记 0
ep[, long_epi := as.integer(duration > 10)]

## 删掉公式里用到但取值缺失的事件
d1 <- ep[complete.cases(ep[, .(long_epi, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]

## 跑模型
m1 <- geeglm(long_epi ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d1, id = sid_int, weights = d1$w,
             family = binomial(), corstr = "exchangeable")

## 从模型里取出"环泊酚 vs 丙泊酚"那个系数（对照 = 丙泊酚）
s1 <- summary(m1)$coefficients
est1 <- s1["GroupCiprofol_only", "Estimate"]   # 对数尺度系数
se1  <- s1["GroupCiprofol_only", "Std.err"]    # 标准误
p1   <- 2 * pnorm(-abs(est1 / se1))            # P 值

## OR 和置信区间 = 系数取指数
result_main <- data.table(
  analysis  = "IPTW-weighted (full cohort)",
  model     = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
  unit      = "Episode (one hypotensive episode)",
  n_unit    = nrow(d1),                       # 多少次事件
  n_patient = uniqueN(d1$PATIENT_ID),         # 涉及多少病人
  estimate  = exp(est1), lo = exp(est1 - 1.96 * se1), hi = exp(est1 + 1.96 * se1),
  p         = p1,
  effect_type = "OR",
  effect    = sprintf("%.2f (%.2f-%.2f)", exp(est1), exp(est1 - 1.96 * se1), exp(est1 + 1.96 * se1)))

cat("M1 每次事件 >10 min 的比例  OR =", result_main$effect,
    "  P =", sprintf("%.3g", p1), "\n")
cat("   （", nrow(d1), "次事件，涉及", uniqueN(d1$PATIENT_ID), "位病人）\n")


## ============================================================================
## 第 4 步：模型 M2 —— 每次低血压事件的严重程度总分（AUC），两组比较
##   AUC = 这次事件的"深度 × 时长"总面积，衡量这一次事件整体有多重
##   因为 AUC 偏态严重，先取对数再比较；结果换算回"百分比变化"
## ============================================================================

ep[, log_eauc := log(event_auc + 0.1)]   # +0.1 是为了避开 log(0)
d2 <- ep[complete.cases(ep[, .(log_eauc, Group, Hypertension, AGE, SEX, ASA_GRADE, w)])]

m2 <- geeglm(log_eauc ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
             data = d2, id = sid_int, weights = d2$w,
             family = gaussian(), corstr = "exchangeable")

s2   <- summary(m2)$coefficients
est2 <- s2["GroupCiprofol_only", "Estimate"]
se2  <- s2["GroupCiprofol_only", "Std.err"]
p2   <- 2 * pnorm(-abs(est2 / se2))

## 对数系数换算成百分比变化：100 * (exp(系数) - 1)
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

cat("M2 每次事件 AUC 变化      ", result_main[model %like% "M2", effect],
    "  P =", sprintf("%.3g", p2), "\n")


## ============================================================================
## 第 5 步：模型 M3 —— 每次低血压事件持续几分钟，两组比较
##   做法和 M2 一样，只是把因变量换成事件时长
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

cat("M3 每次事件时长变化       ", result_main[model %like% "M3", effect],
    "  P =", sprintf("%.3g", p3), "\n")


## ============================================================================
## 第 6 步：模型 M4 —— 每台手术发生几次低血压事件，两组比较
##   数据：手术表（一行 = 一次手术），因变量 = n_epi（这台手术的事件次数）
##   手术有长有短，时间长自然容易攒更多次，所以加 offset(log 麻醉时长)，
##   模型实际比较的是"每小时发生几次"——即事件次数的"率"
## ============================================================================
adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  vv <- dw[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]

dw[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
dw[, lg_dur  := log(anes_duration_mins)]   # 麻醉时长取对数，做 offset

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

cat("M4 事件次数率（每小时）   RR =", result_main[model %like% "M4", effect],
    "  P =", sprintf("%.3g", p4), "\n")

## ============================================================================
## 第 7 步：把四个模型的结果存盘（csv + Word 各一份）
##   预期结果：out\Analysis_Results_main.csv / .docx
## ============================================================================

fwrite(result_main, file.path(OUT, "Analysis_Results_main.csv"))

## 存一份 Word 版
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
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

cat("\n==================== 主结果（IPTW 加权后）====================\n")
print(result_main[, .(model, unit, n_unit, n_patient, type = effect_type,
                      effect, p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])

## ============================================================================
## 第 8 步：对照分析 —— 只在"两组特征有重叠"的范围内，重新加权、重做四个模型
##   安一院从没用过环泊酚，它的丙泊酚病人在"倾向评分"上和环泊酚病人差得太远，落在两组特征不重叠的区域
## ============================================================================

if (dw[in_overlap == FALSE, .N] > 0) {

  cat("\n==================== 对照：仅 PS 重叠区（重新 IPTW）====================\n")
  cat(sprintf("剔除不重叠 %d 例次（占 %.1f%%），其中安一院 %d 例次\n",
              dw[in_overlap == FALSE, .N], 100 * mean(!dw$in_overlap),
              dw[in_overlap == FALSE & Hospital == "AnYi", .N]))
  dw2 <- dw[in_overlap == TRUE]
  print(dw2[, .N, by = .(Hospital, Group)][order(Hospital, Group)])

  ## ---- 8a. 在重叠区内重新拟合倾向评分、重新算权重----
  dw2[, treat := as.integer(Group == "Ciprofol_only")]

  ## 倾向评分用的 10 个基线变量
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

  ## 稳定化权重：处理组 = 发生率/PS，对照组 = (1-发生率)/(1-PS)
  prev <- mean(dw2$treat)
  dw2[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

  ## 权重截断：在每组内部按 1% 和 99% 分位切掉极端权重
  dw2[, ql := as.numeric(quantile(sw_raw, 0.01)), by = Group]
  dw2[, qh := as.numeric(quantile(sw_raw, 0.99)), by = Group]
  dw2[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
  dw2[, c("ql", "qh") := NULL]

  ## 加权前后协变量平衡（SMD = 组间均值差 / 合并标准差，<0.1 算平衡）
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
  cat(sprintf("加权前最大 SMD %.3f → 加权后最大 SMD %.3f\n",
              max(balance_overlap$smd_before, na.rm = TRUE),
              max(balance_overlap$smd_after, na.rm = TRUE)))

  ## ---- 8b. 重叠区内重跑四个模型
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
  ## 在重叠区人群里重新检查一遍：4 个调整变量是否仍然人人有变化
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
  }, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))
}

## ============================================================================
## 第 9 步：两组各自的加权发生率，以及发生率差（绝对效应量）
##   每个结局算三样 —— 环泊酚加权发生率、丙泊酚加权发生率、两者的差
##   差的 95% 置信区间用 bootstrap 算：随机有放回地抽"病人"（不是手术）2000 次
## ============================================================================

cat("\n==================== 绝对效应量（加权风险 + bootstrap）====================\n")
ocs <- c("prolonged", "prolonged_cum", "multi_epi", "short_only", "any_hypo")

## 事件表也挂上权重，用来算"每次事件 >10 min"的发生比例（和 M1 同口径）
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

  ## 先把每台手术的总权重和加权结局求和（一次手术就是分析单位）
  per <- data.table(sid, isc, w, y)[, .(W = sum(w), Y = sum(w * y)), by = .(sid, isc)]
  CW <- numeric(np); CY <- numeric(np); PW <- numeric(np); PY <- numeric(np)
  CW[per$sid[per$isc == 1]] <- per$W[per$isc == 1]; CY[per$sid[per$isc == 1]] <- per$Y[per$isc == 1]
  PW[per$sid[per$isc == 0]] <- per$W[per$isc == 0]; PY[per$sid[per$isc == 0]] <- per$Y[per$isc == 0]

  ## 真实样本里的发生率差
  est <- sum(CY) / sum(CW) - sum(PY) / sum(PW)

  ## bootstrap：重抽病人 2000 次，每次重算发生率差
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

## 手术级结局（一行 = 一次手术）
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

## 事件级结局：每次事件 >10 min 的比例（分母 = 已发生的低血压事件，和手术级的发生率不是一回事，不能放一起比大小）
cw2 <- epw_abs$Group == "Ciprofol_only"; pw2 <- epw_abs$Group == "Propofol_only"
a_e <- boot_rd(epw_abs, "long_epi", wvar = "sw")
abs_epi <- data.table(
  Outcome = "long_epi", Unit = "Episode (one hypotensive episode)",
  Ciprofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(epw_abs$long_epi[cw2], epw_abs$sw[cw2])),
  Propofol_weighted_risk = sprintf("%.1f%%", 100 * weighted.mean(epw_abs$long_epi[pw2], epw_abs$sw[pw2])),
  Risk_difference      = sprintf("%+.1f pp", 100 * a_e["est"]),
  Risk_difference_95CI = sprintf("(%.1f to %.1f)", 100 * a_e["lo"], 100 * a_e["hi"]),
  Ciprofol_n = sum(cw2), Propofol_n = sum(pw2))
abs_tab <- rbind(abs_epi, abs_tab)   # 事件级那一行放最上面

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
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))


## ============================================================================
## 第 10 步：同第 9 步，但只在重叠区内做一遍
## ============================================================================

if (exists("dw2") && !is.null(dw2$sw)) {
  cat("\n==================== 绝对效应量（仅 PS 重叠区）====================\n")
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
  }, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))
}


## ============================================================================
## 第 11 步：E-value —— 需要多强的"未知混杂因素"才能把我们的结论推翻
##   越大说明结论越稳。算法是纯算术，两行就够
## ============================================================================

cat("\n==================== E-value（对 M1 主结局）====================\n")
## E-value 的计算公式：v + sqrt(v * (v - 1))；v<1 时先取倒数
evalue_of <- function(v) { if (v < 1) v <- 1 / v; v + sqrt(v * (v - 1)) }

m1_row <- result_main[model %like% "M1"]
or_pt <- m1_row$estimate                        # OR 点估计
or_ci <- if (or_pt < 1) m1_row$hi else m1_row$lo  # 置信区间里最接近 1 的那一端

cat(sprintf("  M1 点估计: OR = %.2f → RR ≈ %.3f → E-value = %.2f\n",
            or_pt, sqrt(or_pt), evalue_of(sqrt(or_pt))))
cat(sprintf("  M1 CI 端  : OR = %.2f → RR ≈ %.3f → E-value = %.2f\n",
            or_ci, sqrt(or_ci), evalue_of(sqrt(or_ci))))
cat(sprintf("  （若按 OR≈RR 的宽松口径，两个 E-value 分别为 %.2f / %.2f）\n",
            evalue_of(or_pt), evalue_of(or_ci)))


## ============================================================================
## 第 12 步：把主结果表和绝对风险表合成一份 Word
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
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

