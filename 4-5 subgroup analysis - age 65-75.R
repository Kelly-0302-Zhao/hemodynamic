## =============================================================================
## 4-5 subgroup analysis - age 65-75.R
## -----------------------------------------------------------------------------
## 4-5 亚组分析：只看年龄 65-74 岁的手术
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


## 根目录、输出目录
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")                 # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

## 三个 ID 列。三张表的 ID 列类型不一致，用之前都要先转成文字型
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 倾向评分用的 10 个基线变量（和全队列主分析一致；2026-09-30 起不含 PREOPERATIVE_MAP）
PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR",
             "operation_type", "Hypertension", "Diabetes", "Coronary_artery_disease",
             "Heart_failure")

ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## 本亚组标签
LABEL <- "AGE 65-75"

## 读入两张表
tab     <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))  # 手术表
epi_all <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))       # 事件表

## ID 列统一转成文字型
for (v in ID_VARS) {
  set(tab,     j = v, value = as.character(tab[[v]]))
  set(epi_all, j = v, value = as.character(epi_all[[v]]))
}

cat("已读入：手术表", nrow(tab), "例次 | 事件表", nrow(epi_all), "次事件\n")


## ============================================================================
## 第 2 步：取出这个亚组的手术，做未加权（crude）对照，并检查例数够不够
## ============================================================================

## 只挑年龄 65-74 岁的手术例次
x <- tab[AGE >= 65 & AGE < 75]
cat("纳入亚组：手术", nrow(x), "例次 | 患者", uniqueN(x$PATIENT_ID), "人\n")

if (nrow(x) < 200 || uniqueN(x$Group) < 2) {
  cat("  例数不足或只有一组，跳过。\n")
} else {

  cat(sprintf("  分组：Ciprofol %d 例次 / Propofol %d 例次\n",
              sum(x$Group == "Ciprofol_only"), sum(x$Group == "Propofol_only")))

  ## 未加权对照：只看两组在几个汇总指标上的原始差别（不加权，仅作对照）
  crude <- x[, .(n_surgery    = .N,
                 n_patient    = uniqueN(PATIENT_ID),
                 prolonged_pct = 100 * mean(prolonged, na.rm = TRUE),
                 MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
                 n_epi_mean   = mean(n_epi, na.rm = TRUE)),
              by = Group]
  cat("\n  ---- 第 2 步：未加权（crude）对照 ----\n")
  print(crude)


  ## ==========================================================================
  ## 第 3 步：在这个亚组内【重新】拟合倾向评分、算稳定化权重
  ## ==========================================================================
  xw <- copy(x)
  xw[, treat := as.integer(Group == "Ciprofol_only")]

  use <- character(0)
  for (v in PS_VARS) {
    if (!v %in% names(xw)) next
    vv <- xw[[v]]
    if (is.numeric(vv)) {
      if (is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0) use <- c(use, v)
    } else {
      if (uniqueN(as.character(vv[!is.na(vv)])) > 1) use <- c(use, v)
    }
  }
  dropped <- setdiff(PS_VARS, use)
  if (length(dropped))
    cat("   [PS] 本亚组内恒定、已自动剔除：", paste(dropped, collapse = ", "), "\n")

  ## 只在 11 个变量都填全的病人里拟合
  cc  <- complete.cases(xw[, ..use])
  xp  <- xw[cc]
  if (nrow(xp) < 100 || min(table(xp$treat)) < 10)
    stop("本亚组例数不足或某一组太少，无法拟合倾向评分")

  ## 逻辑回归算倾向评分，再算稳定化权重（处理组 = 发生率/PS，对照组 = (1-发生率)/(1-PS)）
  ps <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
            data = xp, family = binomial())
  xp[, ps_hat := predict(ps, type = "response")]
  prev <- mean(xp$treat)
  xp[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

  xp[, `:=`(ql = as.numeric(quantile(sw_raw, 0.01)),
            qh = as.numeric(quantile(sw_raw, 0.99))), by = Group]
  xp[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
  xp[, c("ql", "qh") := NULL]

  lo_overlap <- max(min(xp$ps_hat[xp$treat == 1]), min(xp$ps_hat[xp$treat == 0]))
  hi_overlap <- min(max(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 0]))
  xp[, in_overlap := ps_hat >= lo_overlap & ps_hat <= hi_overlap]

  ## 加权前后协变量平衡（SMD = 组间均值差 / 合并标准差，<0.1 算平衡）
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
  bal <- rbindlist(lapply(use, function(v) {
    vv <- xp[[v]]
    if (is.numeric(vv)) {
      data.table(variable   = v,
                 smd_before = wsmd(vv, rep(1, nrow(xp))),
                 smd_after  = wsmd(vv, xp$sw))
    } else {
      rbindlist(lapply(setdiff(unique(as.character(vv)), NA), function(l) {
        ind <- as.integer(as.character(vv) == l)
        data.table(variable   = paste0(v, ": ", lvl_en(l)),
                   smd_before = wsmd(ind, rep(1, nrow(xp))),
                   smd_after  = wsmd(ind, xp$sw))
      }))
    }
  }))

  ## 有效样本量（Kish）：权重失衡时判断"权重虚高、实际信息量有多少"
  ess <- function(w) sum(w)^2 / sum(w^2)
  ess_tab <- xp[, .(Group       = as.character(first(Group)),
                    n_procedures = .N,
                    n_patients   = uniqueN(PATIENT_ID),
                    weight_median = median(sw),
                    weight_max = max(sw),
                    ESS      = ess(sw)), by = Group][, Group := NULL]

  cat(sprintf("   [PS] 参与加权手术 %d 例次（患者 %d 人）；PS 范围环泊酚 %.3f-%.3f，Propofol %.3f-%.3f\n",
              nrow(xp), uniqueN(xp$PATIENT_ID),
              min(xp$ps_hat[xp$treat == 1]), max(xp$ps_hat[xp$treat == 1]),
              min(xp$ps_hat[xp$treat == 0]), max(xp$ps_hat[xp$treat == 0])))
  cat(sprintf("   [平衡] 加权前最大 SMD %.3f → 加权后最大 SMD %.3f\n",
              max(bal$smd_before, na.rm = TRUE), max(bal$smd_after, na.rm = TRUE)))
  e_c <- ess_tab[Group == "Ciprofol_only", ESS]
  e_p <- ess_tab[Group == "Propofol_only", ESS]
  cat(sprintf("   [权重] Ciprofol ESS %.0f / %d 例次，Propofol ESS %.0f / %d 例次\n",
              e_c, sum(xp$treat == 1), e_p, sum(xp$treat == 0)))
  cat(sprintf("   [重叠] 两组 PS 重叠区间 = [%.3f, %.3f]；区间外 %d 例次（%.1f%%）\n",
              lo_overlap, hi_overlap, sum(!xp$in_overlap), 100 * mean(!xp$in_overlap)))


  ## ==========================================================================
  ## 第 4 步：准备事件级数据 + 模型 M1
  ## ==========================================================================

  ## 把手术表的 分组、权重、调整变量 挂到每个事件上 → ep：一行 = 一次事件
  ep <- epi_all[xp[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                       Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
                on = ID_VARS, nomatch = NULL]
  ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

  ## GEE 要求"聚类编号"必须是整数，这里把手术编号（患者+就诊+手术）转成整数编号
  ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

  ## 调整变量：本亚组内人人一样的要剔除
  adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
    vv <- xp[[v]]
    if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
    else uniqueN(as.character(vv[!is.na(vv)])) > 1
  }, logical(1))]
  ADJ <- if (length(adj_ok)) paste(adj_ok, collapse = " + ") else "1"

  ## 因变量：这次事件持续超过 10 分钟记 1，否则记 0
  ep[, long_epi := as.integer(duration > 10)]
  keep1 <- c("long_epi", "Group", adj_ok, "sid_int")   # data.table 的 with=FALSE 只接受字符向量
  d1 <- ep[complete.cases(ep[, keep1, with = FALSE])]

  m1 <- geeglm(as.formula(paste("long_epi ~ Group +", ADJ)),
               data = d1, id = sid_int, weights = d1$w,
               family = binomial(), corstr = "exchangeable")

  ## 取出"环泊酚 vs 丙泊酚"那个系数（对照 = 丙泊酚）
  s1   <- summary(m1)$coefficients
  est1 <- s1["GroupCiprofol_only", "Estimate"]    # 对数尺度系数
  se1  <- s1["GroupCiprofol_only", "Std.err"]     # 标准误
  p1   <- 2 * pnorm(-abs(est1 / se1))            # P 值

  ## OR 和置信区间 = 系数取指数
  r1 <- data.table(
    analysis   = LABEL,
    model      = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d1),
    n_patient  = uniqueN(d1$PATIENT_ID),
    estimate   = exp(est1),
    lo         = exp(est1 - 1.96 * se1),
    hi         = exp(est1 + 1.96 * se1),
    p          = p1,
    effect_type = "OR",
    effect     = sprintf("%.2f (%.2f-%.2f)", exp(est1), exp(est1 - 1.96 * se1), exp(est1 + 1.96 * se1)))
  result <- r1

  cat("M1 每次事件 >10 min 的比例  OR =", r1$effect, "  P =", sprintf("%.3g", p1), "\n")
  cat("   （", nrow(d1), "次事件，涉及", uniqueN(d1$PATIENT_ID), "位病人）\n")


  ## ==========================================================================
  ## 第 5 步：模型 M2 —— 每次低血压事件的严重程度总分（AUC），两组比较
  ## ==========================================================================

  ep[, log_eauc := log(event_auc + 0.1)]   # +0.1 是为了避开 log(0)
  keep2 <- c("log_eauc", "Group", adj_ok, "sid_int")
  d2 <- ep[complete.cases(ep[, keep2, with = FALSE])]

  m2 <- geeglm(as.formula(paste("log_eauc ~ Group +", ADJ)),
               data = d2, id = sid_int, weights = d2$w,
               family = gaussian(), corstr = "exchangeable")

  s2   <- summary(m2)$coefficients
  est2 <- s2["GroupCiprofol_only", "Estimate"]
  se2  <- s2["GroupCiprofol_only", "Std.err"]
  p2   <- 2 * pnorm(-abs(est2 / se2))

  ## 对数系数换算成百分比变化：100 * (exp(系数) - 1)
  r2 <- data.table(
    analysis   = LABEL,
    model      = "M2 Per-episode hypotension burden (AUC)",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d2),
    n_patient  = uniqueN(d2$PATIENT_ID),
    estimate   = 100 * (exp(est2) - 1),
    lo         = 100 * (exp(est2 - 1.96 * se2) - 1),
    hi         = 100 * (exp(est2 + 1.96 * se2) - 1),
    p          = p2,
    effect_type = "Pct",
    effect     = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                         100 * (exp(est2) - 1),
                         100 * (exp(est2 - 1.96 * se2) - 1),
                         100 * (exp(est2 + 1.96 * se2) - 1)))
  result <- rbind(result, r2)

  cat("M2 每次事件 AUC 变化      ", r2$effect, "  P =", sprintf("%.3g", p2), "\n")


  ## ==========================================================================
  ## 第 6 步：模型 M3 —— 每次低血压事件持续几分钟，两组比较
  ## ==========================================================================

  ep[, log_duration := log(duration + 0.1)]
  keep3 <- c("log_duration", "Group", adj_ok, "sid_int")
  d3 <- ep[complete.cases(ep[, keep3, with = FALSE])]

  m3 <- geeglm(as.formula(paste("log_duration ~ Group +", ADJ)),
               data = d3, id = sid_int, weights = d3$w,
               family = gaussian(), corstr = "exchangeable")

  s3   <- summary(m3)$coefficients
  est3 <- s3["GroupCiprofol_only", "Estimate"]
  se3  <- s3["GroupCiprofol_only", "Std.err"]
  p3   <- 2 * pnorm(-abs(est3 / se3))

  r3 <- data.table(
    analysis   = LABEL,
    model      = "M3 Duration of each hypotensive episode",
    unit       = "Episode (one hypotensive episode)",
    n_unit     = nrow(d3),
    n_patient  = uniqueN(d3$PATIENT_ID),
    estimate   = 100 * (exp(est3) - 1),
    lo         = 100 * (exp(est3 - 1.96 * se3) - 1),
    hi         = 100 * (exp(est3 + 1.96 * se3) - 1),
    p          = p3,
    effect_type = "Pct",
    effect     = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)",
                         100 * (exp(est3) - 1),
                         100 * (exp(est3 - 1.96 * se3) - 1),
                         100 * (exp(est3 + 1.96 * se3) - 1)))
  result <- rbind(result, r3)

  cat("M3 每次事件时长变化       ", r3$effect, "  P =", sprintf("%.3g", p3), "\n")


  ## ==========================================================================
  ## 第 7 步：模型 M4 —— 每台手术发生几次低血压事件，两组比较
  ## ==========================================================================

  xp[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
  xp[, lg_dur  := log(anes_duration_mins)]   # 麻醉时长取对数，做 offset

  keep4 <- c("n_epi", "Group", adj_ok, "lg_dur", "sid_int")
  d4 <- xp[complete.cases(xp[, keep4, with = FALSE])]

  m4 <- geeglm(as.formula(paste("n_epi ~ Group +", ADJ, "+ offset(lg_dur)")),
               data = d4, id = sid_int, weights = d4$sw,
               family = poisson(), corstr = "exchangeable")

  s4   <- summary(m4)$coefficients
  est4 <- s4["GroupCiprofol_only", "Estimate"]
  se4  <- s4["GroupCiprofol_only", "Std.err"]
  p4   <- 2 * pnorm(-abs(est4 / se4))

  r4 <- data.table(
    analysis   = LABEL,
    model      = "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)",
    unit       = "Procedure (one anesthesia)",
    n_unit     = nrow(d4),
    n_patient  = uniqueN(d4$PATIENT_ID),
    estimate   = exp(est4),
    lo         = exp(est4 - 1.96 * se4),
    hi         = exp(est4 + 1.96 * se4),
    p          = p4,
    effect_type = "RR",
    effect     = sprintf("%.2f (%.2f-%.2f)", exp(est4), exp(est4 - 1.96 * se4), exp(est4 + 1.96 * se4)))
  result <- rbind(result, r4)

  cat("M4 事件次数率（每小时）   RR =", r4$effect, "  P =", sprintf("%.3g", p4), "\n")
  cat("\n  ---- 亚组内 IPTW 加权后的结果 ----\n")
  print(result[, .(model, unit, effect, effect_type, p)])


  ## ==========================================================================
  ## 第 8 步：加权后的绝对指标（两组各自的加权发生率等汇总）
  ## ==========================================================================

  abs_tab <- xp[, .(n_surgery    = .N,
                    n_patient    = uniqueN(PATIENT_ID),
                    prolonged_pct = 100 * weighted.mean(prolonged, sw, na.rm = TRUE),
                    MAP65_AUC_median = median(MAP65_AUC, na.rm = TRUE),
                    n_epi_mean   = weighted.mean(n_epi, sw, na.rm = TRUE)),
                 by = Group]
  cat("\n  ---- 加权后绝对指标 ----\n")
  print(abs_tab)


  ## ==========================================================================
  ## 第 9 步：把结果存盘
  ## ==========================================================================

  fwrite(crude,   file.path(OUT, paste0("Sensitivity_Analysis_Crude_",     LABEL, ".csv")))
  fwrite(abs_tab, file.path(OUT, paste0("Sensitivity_Analysis_Absolute_", LABEL, ".csv")))
  fwrite(bal,     file.path(OUT, paste0("Table1_IPTW_",                    LABEL, ".csv")))
  fwrite(ess_tab, file.path(OUT, paste0("Sensitivity_Analysis_ESS_",       LABEL, ".csv")))
  fwrite(result,  file.path(OUT, paste0("Sensitivity_Analysis_Results_",  LABEL, ".csv")))

  tryCatch({
    suppressPackageStartupMessages({ library(flextable); library(officer) })
    ft <- flextable(as.data.frame(result[, .(
             Model = model, Unit = unit, Sample_size = n_unit, n_patients = n_patient,
             Estimate = effect, Type = effect_type,
             P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])) |>
         theme_booktabs() |> autofit() |>
         set_caption(caption = sprintf(paste(
           "Sensitivity analysis - %s. Weighted GEE models after re-fitting IPTW within this subgroup.",
           "M1-M3: one hypotensive episode per row; M4: one surgical procedure per row;",
           "all clustered by surgical procedure."), LABEL))
    print(read_docx() |> body_add_flextable(ft),
          target = file.path(OUT, paste0("Sensitivity_Analysis_Results_", LABEL, ".docx")))

    ft2 <- flextable(as.data.frame(bal)) |> theme_booktabs() |> autofit() |>
           set_caption(caption = sprintf("Table 1. Covariate balance before and after IPTW - %s.", LABEL))
    print(read_docx() |> body_add_flextable(ft2),
          target = file.path(OUT, paste0("Table1_IPTW_", LABEL, ".docx")))
  }, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

}

