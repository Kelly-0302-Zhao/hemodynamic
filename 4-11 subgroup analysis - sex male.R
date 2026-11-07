## =============================================================================
## 4-11 subgroup analysis - sex male.R
## -----------------------------------------------------------------------------
## 4-11 亚组分析：仅男性亚组
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
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## 三个 ID 列。三张表的 ID 列类型不一致，用之前都要先转成文字型
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 文件标签：决定所有输出文件名
lab <- "male"

PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR",
             "operation_type", "Hypertension", "Diabetes",
             "Coronary_artery_disease", "Heart_failure")

ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

cat("[准备] 已就绪 | 输出目录：", OUT, "\n")


## ============================================================================
## 第 2 步：读入手术级分析表，按"仅男性"取子集，报告人群规模
## ============================================================================

## 手术表：一行 = 一次手术
tab <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))

## 三张表的 ID 列统一转成文字型
for (v in ID_VARS) set(tab, j = v, value = as.character(tab[[v]]))


x <- tab[as.character(SEX) == "男"]   

cat("纳入敏感性队列（仅男性）：",
    sprintf("手术 %d 例次 | 患者 %d 人", nrow(x), uniqueN(x$PATIENT_ID)), "\n")
cat(sprintf("  分组：Ciprofol %d 例次 / Propofol %d 例次\n",
            sum(x$Group == "Ciprofol_only"), sum(x$Group == "Propofol_only")))


## ============================================================================
## 第 3 步：未加权（crude）粗对照 
## ============================================================================
crude <- x[, .(n_surgery = .N,
               n_patient = uniqueN(PATIENT_ID),
               prolonged_pct    = 100 * mean(prolonged, na.rm = TRUE),
               MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
               n_epi_mean = mean(n_epi, na.rm = TRUE)),
             by = Group]
print(crude)


## ============================================================================
## 第 4 步：在这个亚组内重新拟合倾向评分、重新算权重（IPTW）
## ============================================================================

## 给每个手术打上"是否环泊酚"的处理标志
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
  cat("   [PS] 本亚组内恒定、已自动剔除：", paste(dropped, collapse = ", "), "\n")

cc <- complete.cases(x[, ..use])
xp <- x[cc]
if (nrow(xp) < 100 || min(table(xp$treat)) < 10)
  stop("本亚组例数不足或某一组太少，无法拟合倾向评分")

ps_fit <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
              data = xp, family = binomial())
xp[, ps_hat := predict(ps_fit, type = "response")]

## 稳定化权重：处理组 = 发生率/PS，对照组 = (1-发生率)/(1-PS)
prev <- mean(xp$treat)
xp[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

xp[, `:=`(ql = as.numeric(quantile(sw_raw, 0.01)),
          qh = as.numeric(quantile(sw_raw, 0.99))), by = Group]
xp[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
xp[, c("ql", "qh") := NULL]

## 共同支持域（重叠区）标记：两组 PS 不重叠的手术标出来
lo_overlap <- max(min(xp$ps_hat[xp$treat == 1]), min(xp$ps_hat[xp$treat == 0]))
hi_overlap <- min(max(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 0]))
xp[, in_overlap := ps_hat >= lo_overlap & ps_hat <= hi_overlap]

## 加权前后平衡性（SMD = 组间均值差 / 合并标准差，<0.1 算平衡）
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

## 有效样本量（ESS）与权重分布，按组算
ess_tab <- xp[, .(Group = as.character(first(Group)),
                   n_procedures = .N,
                   n_patients   = uniqueN(PATIENT_ID),
                   weight_median = median(sw),
                   weight_max = max(sw),
                   ESS      = sum(sw)^2 / sum(sw^2)),
              by = Group][, Group := NULL]

cat(sprintf("   [PS] 参与加权手术 %d 例次（患者 %d 人）；PS 范围环泊酚 %.3f-%.3f，Propofol %.3f-%.3f\n",
            nrow(xp), uniqueN(xp$PATIENT_ID),
            min(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 1]),
            min(xp$ps_hat[xp$treat == 0]), max(xp$ps_hat[xp$treat == 0])))
cat(sprintf("   [平衡] 加权前最大 SMD %.3f（>0.1 的变量 %d 个）→ 加权后最大 SMD %.3f（%d 个）\n",
            max(balance$smd_before, na.rm = TRUE), sum(balance$smd_before > .1, na.rm = TRUE),
            max(balance$smd_after,  na.rm = TRUE), sum(balance$smd_after  > .1, na.rm = TRUE)))
e_c <- ess_tab[Group == "Ciprofol_only", ESS]
e_p <- ess_tab[Group == "Propofol_only", ESS]
cat(sprintf("   [权重] 组内截断分位 1%%/99%%；Ciprofol ESS %.0f / %d 例次，Propofol ESS %.0f / %d 例次\n",
            e_c, sum(xp$treat == 1), e_p, sum(xp$treat == 0)))
cat(sprintf("   [重叠] 两组 PS 重叠区间 = [%.3f, %.3f]；落在区间外的例次 %d（占 %.1f%%）\n",
            lo_overlap, hi_overlap, sum(!xp$in_overlap), 100 * mean(!xp$in_overlap)))
if (sum(!xp$in_overlap) > 0 && "Hospital" %in% names(xp)) {
  cat("          区间外例次的中心分布：\n")
  print(xp[in_overlap == FALSE, .N, by = .(Hospital, Group)][order(-N)])
}

## 加权分析人群：xp（含 ps_hat、sw）
xw <- xp


## ============================================================================
## 第 5 步：读入事件级表，只保留本队列里的事件，并挂上分组/权重/调整变量
## ============================================================================

epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))
for (v in ID_VARS) set(epi, j = v, value = as.character(epi[[v]]))

## 只保留"属于本队列手术"的事件
keep_ids <- xw[, .(PATIENT_ID, VISIT_ID, OPER_ID)]
ep <- epi[keep_ids, on = ID_VARS, nomatch = NULL]

## 把每个事件所属手术的分组、权重、调整变量挂上来 → ep：一行 = 一次事件
ep <- ep[xw[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
         on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

cat("事件表纳入：", nrow(ep), "次事件，涉及", uniqueN(ep$PATIENT_ID), "位病人\n")


## ============================================================================
## 第 6 步：模型 M1 —— 每次低血压事件，持续超过 10 分钟的比例，两组比较
## ============================================================================

adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  if (!v %in% names(xw)) return(FALSE)
  vv <- xw[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]
ADJ <- if (length(adj_ok)) paste(adj_ok, collapse = " + ") else "1"

## 事件表也要有整数型聚类编号（geeglm 的 id 必须整数）
ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

## 因变量：这次事件持续超过 10 分钟记 1，否则记 0
ep[, long_epi := as.integer(duration > 10)]

## 删掉公式里用到但取值缺失的事件
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
  cat("M1 每次事件 >10 min 的比例  OR =", result_sens[model %like% "M1", effect],
      "  P =", sprintf("%.3g", p1), "\n")
} else {
  cat("M1 事件数不足，跳过\n")
}


## ============================================================================
## 第 7 步：模型 M2 —— 每次低血压事件的严重程度总分（AUC），两组比较
## ============================================================================
ep[, log_eauc := log(event_auc + 0.1)]   # +0.1 是为了避开 log(0)
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
  cat("M2 每次事件 AUC 变化      ", result_sens[model %like% "M2", effect],
      "  P =", sprintf("%.3g", p2), "\n")
} else {
  cat("M2 事件数不足，跳过\n")
}


## ============================================================================
## 第 8 步：模型 M3 —— 每次低血压事件持续几分钟，两组比较
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
  cat("M3 每次事件时长变化       ", result_sens[model %like% "M3", effect],
      "  P =", sprintf("%.3g", p3), "\n")
} else {
  cat("M3 事件数不足，跳过\n")
}


## ============================================================================
## 第 9 步：模型 M4 —— 每台手术发生几次低血压事件，两组比较
## ============================================================================
xw[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
xw[, lg_dur  := log(anes_duration_mins)]   # 麻醉时长取对数，做 offset

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
cat("M4 事件次数率（每小时）   RR =", result_sens[model %like% "M4", effect],
    "  P =", sprintf("%.3g", p4), "\n")

cat("\n==================== 加权后四模型结果（仅男性）====================\n")
print(result_sens[, .(model, unit, n_unit, n_patient, type = effect_type, effect,
                      p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])


## ============================================================================
## 第 10 步：加权后的绝对指标 —— 两组各自的加权发生率、事件次数均值
## ============================================================================
abs_tab <- xw[, .(n_surgery = .N,
                  n_patient = uniqueN(PATIENT_ID),
                  prolonged_pct     = 100 * weighted.mean(prolonged, sw, na.rm = TRUE),
                  MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
                  n_epi_mean = weighted.mean(n_epi, sw, na.rm = TRUE)),
                by = Group]
print(abs_tab)


## ============================================================================
## 第 11 步：把结果存盘
## ============================================================================
fwrite(crude,       file.path(OUT, paste0("Sensitivity_Analysis_Crude_",    lab, ".csv")))
fwrite(abs_tab,     file.path(OUT, paste0("Sensitivity_Analysis_Absolute_", lab, ".csv")))
fwrite(balance,     file.path(OUT, paste0("Table1_IPTW_",                   lab, ".csv")))
fwrite(ess_tab,     file.path(OUT, paste0("Sensitivity_Analysis_ESS_",      lab, ".csv")))
fwrite(result_sens, file.path(OUT, paste0("Sensitivity_Analysis_Results_",  lab, ".csv")))

## docx：结果表
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
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

## docx：平衡表
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft2 <- flextable(as.data.frame(balance)) |> theme_booktabs() |> autofit() |>
        set_caption(caption = paste0("Table 1. Covariate balance before and after IPTW - ", lab, "."))
  print(read_docx() |> body_add_flextable(ft2),
        target = file.path(OUT, paste0("Table1_IPTW_", lab, ".docx")))
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

