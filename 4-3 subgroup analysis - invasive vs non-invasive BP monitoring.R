## =============================================================================
## 4-3 亚组分析：有创 vs 无创血压监测（按结局血压实际来源归层）
## -----------------------------------------------------------------------------
## 分析单位：手术（procedure）。一台手术 = 一个观测 = 一个 cluster（cluster_id）。
## 分层：手术级唯一归组——每台手术只能进一层，两层互斥且并集 = 全队列。
##       同一患者的不同手术可以分属两层
##       不是重复计数，不做任何患者级剔除。
## =============================================================================

suppressPackageStartupMessages({ library(data.table); library(geepack) })

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
PS_VARS  <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR", "operation_type",
              "Hypertension", "Diabetes", "Coronary_artery_disease", "Heart_failure")
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## 主稿基准
REF <- list(M1 = c(0.51, 0.44, 0.58), M2 = c(-27.4, -30.9, -23.6),
            M3 = c(-27.2, -29.7, -24.5), M4 = c(1.22, 1.14, 1.30))

## ============================================================================
## 1. 读入 + 归层
## ============================================================================
tab <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))
lay <- fread(file.path(OUT, "Modality_Final_Layer_Assignment_v4.csv"), encoding = "UTF-8")

for (v in ID_VARS) {
  set(tab, j = v, value = as.character(tab[[v]]))
  set(epi, j = v, value = as.character(epi[[v]]))
  set(lay, j = v, value = as.character(lay[[v]]))
}

## 统一簇标识：cluster_id = "PATIENT-VISIT-OPER"（与归层文件同格式）
tab[, cluster_id := paste(PATIENT_ID, VISIT_ID, OPER_ID, sep = "-")]
if (!"cluster_id" %in% names(lay)) lay[, cluster_id := paste(PATIENT_ID, VISIT_ID, OPER_ID, sep = "-")]
lay[, cluster_id := as.character(cluster_id)]

stopifnot("手术级表 cluster_id 不唯一" = uniqueN(tab$cluster_id) == nrow(tab),
          "归层文件 cluster_id 不唯一" = uniqueN(lay$cluster_id) == nrow(lay))

tab <- merge(tab, unique(lay[, .SD, .SDcols = c("cluster_id", "layer_label", "layer_final",
                                                "sbp_src", "adjud_evid", "adjud_verdict",
                                                "n_cvp_rows")]),
             by = "cluster_id", all.x = TRUE)

stopifnot("归层合并后行数与手术数不符（归层文件有重复）" = nrow(tab) == 85397,
          "每台手术的 cluster_id 不唯一" = uniqueN(tab$cluster_id) == nrow(tab))

## ---- 规则链重建----
lay_keys <- lay[, .SD, .SDcols = c(ID_VARS, "sbp_src", "adjud_evid", "adjud_verdict", "layer_final")]
rule_layer <- function(sbp_src, evid, verdict) {
  fifelse(sbp_src == "arterial", "invasive",
  fifelse(sbp_src == "cuff", "non_invasive",
  fifelse(verdict == "arterial" |
            evid %in% c("结局值匹配动脉压", "术中窗动脉压证据",
                        "标签口径有创(MAP通道,袖带通道全空)+结局血压缺失", "CVP旁证"),
          "invasive",
  fifelse(verdict == "cuff" |
            evid %in% c("结局值匹配袖带压", "术中窗袖带压证据"),
          "non_invasive", "UNRESOLVED"))))
}
chk <- lay_keys[, .(layer_rule = rule_layer(sbp_src, adjud_evid, adjud_verdict),
                    layer_final)]
stopifnot("规则链与固化归层不一致" = all(chk$layer_rule == chk$layer_final),
          "归层存在缺失" = !any(is.na(tab$layer_final)),
          "归层存在非法值" = all(tab$layer_final %in% c("invasive", "non_invasive")),
          "两层台数之和 ≠ 全队列（存在重叠或遗漏）" =
            sum(tab$layer_final == "invasive") + sum(tab$layer_final == "non_invasive") == nrow(tab))

tab[, invasive := as.integer(layer_final == "invasive")]
stopifnot("同一台手术被归入两层" =
            tab[, uniqueN(layer_final), by = cluster_id][, max(V1)] == 1L)

## 聚类 id = cluster_id 的整数编码（geeglm 的 id 需整数）
tab[, sid_int := as.integer(factor(cluster_id))]
stopifnot("cluster_id 与 sid_int 非一一对应" =
            uniqueN(tab$cluster_id) == uniqueN(tab$sid_int) && nrow(tab) == uniqueN(tab$sid_int))

cat(sprintf("归层完成：%d 台 = 有创 %d + 无创 %d（零遗漏、零重叠）\n",
            nrow(tab), sum(tab$invasive == 1), sum(tab$invasive == 0)))
cat(sprintf("分析单位 = 手术：%d 台 / %d 名患者；每台手术一个簇（cluster_id 唯一，无一脚踏两层）\n",
            nrow(tab), uniqueN(tab$PATIENT_ID)))
cat("\n=== 归层证据分布===\n")
ev <- tab[, .N, by = .(sbp_src, adjud_evid, layer_final)][order(sbp_src, -N)]
print(ev)
fwrite(ev, file.path(OUT, "Subgroup_BPSourceMatched_LayerAssignment.csv"))
cat(sprintf("\nCVP>0 的手术：%d 台（有创 %d / 无创 %d；无创者结局值匹配袖带压，按结局来源归层）\n",
            sum(tab$n_cvp_rows > 0, na.rm = TRUE),
            sum(tab$n_cvp_rows > 0 & tab$invasive == 1, na.rm = TRUE),
            sum(tab$n_cvp_rows > 0 & tab$invasive == 0, na.rm = TRUE)))
cat(sprintf("标签口径 vs 结局口径不一致（双监测、结局实际取袖带压）：%d 台\n",
            sum(tab$layer_label == "invasive" & tab$layer_final == "non_invasive")))

## ============================================================================
## 2. 描述表
## ============================================================================
tab[, rate_map := n_map / anes_duration_mins]
desc <- tab[, .(n_surg = .N, n_patient = uniqueN(PATIENT_ID),
                n_map_med   = as.numeric(median(n_map, na.rm = TRUE)),
                rate_map_med = as.numeric(median(rate_map, na.rm = TRUE)),
                dur_med     = as.numeric(median(anes_duration_mins, na.rm = TRUE)),
                prolonged_pct = round(100 * mean(prolonged, na.rm = TRUE), 2),
                n_epi_mean  = round(mean(n_epi, na.rm = TRUE), 3)),
            by = .(layer_final, Group)]
cat("\n=== 描述：归层 × 组别 ===\n"); print(desc)

epi_dens <- merge(epi[, .SD, .SDcols = c(ID_VARS, "n_point", "duration", "event_auc")],
                  tab[, .SD, .SDcols = c(ID_VARS, "layer_final", "Group")],
                  by = ID_VARS)
dens <- epi_dens[, .(n_event = .N,
                     n_point_med = as.numeric(median(n_point, na.rm = TRUE)),
                     dur_med = as.numeric(median(duration, na.rm = TRUE))),
                 by = .(layer_final, Group)]
cat("\n=== 描述：事件级测压点数 ===\n"); print(dens)

fwrite(rbind(cbind(level = "surgery", desc), cbind(level = "event", dens), fill = TRUE),
       file.path(OUT, "Subgroup_BPSourceMatched_Descriptive.csv"))

## ============================================================================
## 3. 人群分析引擎：一个 -> 四模型全套（IPTW 组内重拟合 + 加权 GEE）
## ============================================================================
run_population <- function(surg, epi_all, lab) {

  x <- copy(surg)
  x[, treat := as.integer(Group == "Ciprofol_only")]

  ## 剔除本人群内恒定的 PS 变量
  use <- character(0)
  for (v in PS_VARS) {
    if (!v %in% names(x)) next
    vv <- x[[v]]
    if (is.numeric(vv)) { if (is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0) use <- c(use, v) }
    else if (uniqueN(as.character(vv[!is.na(vv)])) > 1) use <- c(use, v)
  }

  xp <- x[complete.cases(x[, .SD, .SDcols = use])]
  if (nrow(xp) < 100 || min(table(xp$treat)) < 10)
    stop(sprintf("[%s] 例数不足，无法拟合 PS", lab))
  xp <- droplevels(xp)

  ps_fit <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
                data = xp, family = binomial())
  xp[, ps_hat := predict(ps_fit, type = "response")]

  ## 稳定化权重 + 组内 1%/99% 截断
  prev <- mean(xp$treat)
  xp[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]
  xp[, `:=`(ql = as.numeric(quantile(sw_raw, 0.01)),
            qh = as.numeric(quantile(sw_raw, 0.99))), by = Group]
  xp[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
  xp[, c("ql", "qh") := NULL]

  ## 平衡性 SMD
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
  lvl_en <- function(x) { s <- c("女" = "Female", "男" = "Male", "Ⅰ" = "I", "Ⅱ" = "II", "Ⅲ" = "III")
                          r <- s[x]; ifelse(is.na(r), x, r) }
  bal_rows <- list()
  for (v in use) {
    vv <- xp[[v]]
    if (is.numeric(vv)) {
      bal_rows[[length(bal_rows) + 1]] <- data.table(
        population = lab, variable = v,
        smd_before = wsmd(vv, rep(1, nrow(xp))), smd_after = wsmd(vv, xp$sw))
    } else {
      for (lv in setdiff(unique(as.character(vv)), NA)) {
        ind <- as.integer(as.character(vv) == lv)
        bal_rows[[length(bal_rows) + 1]] <- data.table(
          population = lab, variable = paste0(v, ": ", lvl_en(lv)),
          smd_before = wsmd(ind, rep(1, nrow(xp))), smd_after = wsmd(ind, xp$sw))
      }
    }
  }
  balance <- rbindlist(bal_rows)

  ess_tab <- xp[, .(population = lab, Group = as.character(first(Group)),
                    n_procedures = .N, n_patients = uniqueN(PATIENT_ID),
                    weight_median = median(sw), weight_max = max(sw),
                    ESS = sum(sw)^2 / sum(sw^2)), by = Group][, Group := NULL]

  adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
    if (!v %in% names(xp)) return(FALSE)
    vv <- xp[[v]]
    if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
    else uniqueN(as.character(vv[!is.na(vv)])) > 1
  }, logical(1))]
  ADJ <- if (length(adj_ok)) paste(adj_ok, collapse = " + ") else "1"

  xw <- xp
  cat(sprintf("   [%s] %d 台；SMD %.3f -> %.3f；ESS 环 %.0f / 丙 %.0f\n",
              lab, nrow(xw), max(balance$smd_before, na.rm = TRUE),
              max(balance$smd_after,  na.rm = TRUE),
              ess_tab[Group == "Ciprofol_only", ESS], ess_tab[Group == "Propofol_only", ESS]))

  ## 事件表圈人 + 挂权重（sid_int 由手术级表继承：一台手术 = 一个簇）
  xw_small <- xw[, .SD, .SDcols = c(ID_VARS, "sid_int", "Group", adj_ok, "sw")]
  setnames(xw_small, "sw", "w")
  ep <- epi_all[xw[, .SD, .SDcols = ID_VARS], on = ID_VARS, nomatch = NULL]
  ep <- ep[xw_small, on = ID_VARS, nomatch = NULL]
  ep <- ep[!is.na(Group) & is.finite(w) & w > 0]
  stopifnot("事件级簇数 ≠ 手术数" =
              uniqueN(ep$sid_int) == uniqueN(ep[, paste(PATIENT_ID, VISIT_ID, OPER_ID)]))
  ep[, long_epi := as.integer(duration > 10)]

  result <- data.table(
    population = character(), model = character(), unit = character(),
    n_unit = integer(), n_patient = integer(),
    estimate = numeric(), lo = numeric(), hi = numeric(),
    p = numeric(), effect_type = character(), effect = character())

  ## M1：每次事件 >10min 比例（OR）
  d1 <- ep[complete.cases(ep[, c("long_epi", "Group", adj_ok, "sid_int"), with = FALSE])]
  if (nrow(d1) > 50 && uniqueN(d1$Group) > 1) {
    m1 <- geeglm(as.formula(paste("long_epi ~ Group +", ADJ)), data = d1,
                 id = sid_int, weights = d1$w, family = binomial(), corstr = "exchangeable")
    cm <- as.matrix(summary(m1)$coefficients)
    e <- cm["GroupCiprofol_only", "Estimate"]; s <- cm["GroupCiprofol_only", "Std.err"]
    or <- function(b) exp(b)
    result <- rbind(result, data.table(
      population = lab, model = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
      unit = "Episode", n_unit = nrow(d1), n_patient = uniqueN(d1$PATIENT_ID),
      estimate = or(e), lo = or(e - 1.96 * s), hi = or(e + 1.96 * s),
      p = 2 * pnorm(-abs(e / s)), effect_type = "OR",
      effect = sprintf("%.2f (%.2f-%.2f)", or(e), or(e - 1.96 * s), or(e + 1.96 * s))))
  }

  ## M2 / M3：每次事件 AUC、时长（% 差）
  for (mm in list(list("log_eauc", "event_auc", "M2 Per-episode hypotension burden (AUC)"),
                  list("log_dur", "duration",   "M3 Duration of each hypotensive episode"))) {
    vn <- mm[[1]]; raw <- mm[[2]]
    ep[, (vn) := log(get(raw) + 0.1)]
    dm <- ep[complete.cases(ep[, c(vn, "Group", adj_ok, "sid_int"), with = FALSE])]
    if (nrow(dm) > 50 && uniqueN(dm$Group) > 1) {
      gm <- geeglm(as.formula(paste(vn, "~ Group +", ADJ)), data = dm,
                   id = sid_int, weights = dm$w, family = gaussian(), corstr = "exchangeable")
      cm <- as.matrix(summary(gm)$coefficients)
      e <- cm["GroupCiprofol_only", "Estimate"]; s <- cm["GroupCiprofol_only", "Std.err"]
      pc <- function(b) 100 * (exp(b) - 1)
      result <- rbind(result, data.table(
        population = lab, model = mm[[3]], unit = "Episode", n_unit = nrow(dm),
        n_patient = uniqueN(dm$PATIENT_ID),
        estimate = pc(e), lo = pc(e - 1.96 * s), hi = pc(e + 1.96 * s),
        p = 2 * pnorm(-abs(e / s)), effect_type = "Pct",
        effect = sprintf("%+.1f%% (%+.1f%% to %+.1f%%)", pc(e), pc(e - 1.96 * s), pc(e + 1.96 * s))))
    }
  }

  ## M4：每台事件次数率（RR，offset log 时长；簇 = 手术）
  xw[, lg_dur := log(anes_duration_mins)]
  d4 <- xw[complete.cases(xw[, c("n_epi", "Group", adj_ok, "lg_dur", "sid_int"), with = FALSE])]
  if (nrow(d4) > 50 && uniqueN(d4$Group) > 1) {
    m4 <- geeglm(as.formula(paste("n_epi ~ Group +", ADJ, "+ offset(lg_dur)")), data = d4,
                 id = sid_int, weights = d4$sw, family = poisson(), corstr = "exchangeable")
    cm <- as.matrix(summary(m4)$coefficients)
    e <- cm["GroupCiprofol_only", "Estimate"]; s <- cm["GroupCiprofol_only", "Std.err"]
    result <- rbind(result, data.table(
      population = lab, model = "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)",
      unit = "Procedure", n_unit = nrow(d4), n_patient = uniqueN(d4$PATIENT_ID),
      estimate = exp(e), lo = exp(e - 1.96 * s), hi = exp(e + 1.96 * s),
      p = 2 * pnorm(-abs(e / s)), effect_type = "RR",
      effect = sprintf("%.2f (%.2f-%.2f)", exp(e), exp(e - 1.96 * s), exp(e + 1.96 * s))))
  }

  crude <- x[, .(population = lab, n_surgery = .N, n_patient = uniqueN(PATIENT_ID),
                 prolonged_pct = 100 * mean(prolonged, na.rm = TRUE),
                 MAP65_AUC_median = as.numeric(median(MAP65_AUC, na.rm = TRUE)),
                 n_epi_mean = mean(n_epi, na.rm = TRUE)), by = Group]
  abs_tab <- xw[, .(population = lab, n_surgery = .N, n_patient = uniqueN(PATIENT_ID),
                    prolonged_pct = 100 * weighted.mean(prolonged, sw, na.rm = TRUE),
                    MAP65_AUC_median = as.numeric(median(MAP65_AUC, na.rm = TRUE)),
                    n_epi_mean = weighted.mean(n_epi, sw, na.rm = TRUE)), by = Group]

  cat(sprintf("   [%s] %s\n", lab,
              paste(result[, effect], collapse = " | ")))
  list(crude = crude, balance = balance, ess = ess_tab, result = result,
       abs_tab = abs_tab, xw = xw, adj_ok = adj_ok)
}

## ============================================================================
## 4. 基线复现（整体队列，必须复现主稿 M1-M4）
## ============================================================================
cat("\n==================== 基线复现：整体队列 ====================\n")
run_overall <- run_population(tab, epi, "overall")

got <- run_overall$result
ref_v <- c(REF$M1, REF$M2, REF$M3, REF$M4)
chk <- data.table(
  model  = rep(c("M1", "M2", "M3", "M4"), each = 3),
  metric = rep(c("est", "lo", "hi"), 4),
  got = c(got[model == "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode", estimate],
          got[model == "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode", lo],
          got[model == "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode", hi],
          got[model == "M2 Per-episode hypotension burden (AUC)", estimate],
          got[model == "M2 Per-episode hypotension burden (AUC)", lo],
          got[model == "M2 Per-episode hypotension burden (AUC)", hi],
          got[model == "M3 Duration of each hypotensive episode", estimate],
          got[model == "M3 Duration of each hypotensive episode", lo],
          got[model == "M3 Duration of each hypotensive episode", hi],
          got[model == "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)", estimate],
          got[model == "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)", lo],
          got[model == "M4 Number of hypotensive episodes per procedure (offset: anesthesia duration)", hi]),
  ref = ref_v)
chk[, ok := abs(got - ref) <= 0.051]
print(chk)
if (!all(chk$ok)) stop("基线未能复现主稿数字，停止。")
cat("基线复现通过。\n")

## ============================================================================
## 5. 有创 / 无创亚组（各自组内重拟合 IPTW）
## ============================================================================
cat("\n==================== 有创亚组 ====================\n")
run_inv <- run_population(tab[invasive == 1], epi, "invasive")
cat("\n==================== 无创亚组 ====================\n")
run_non <- run_population(tab[invasive == 0], epi, "non_invasive")

result_all <- rbind(run_overall$result, run_inv$result, run_non$result, fill = TRUE)
crude_all  <- rbind(run_overall$crude,  run_inv$crude,  run_non$crude,  fill = TRUE)
abs_all    <- rbind(run_overall$abs_tab, run_inv$abs_tab, run_non$abs_tab, fill = TRUE)
bal_all    <- rbind(run_overall$balance, run_inv$balance, run_non$balance, fill = TRUE)
ess_all    <- rbind(run_overall$ess,    run_inv$ess,    run_non$ess,    fill = TRUE)

## ============================================================================
## 6. 交互检验（整体队列权重，加 Group x invasive）
## ============================================================================
cat("\n==================== 交互检验 ====================\n")
xw <- run_overall$xw[!is.na(invasive)]
## invasive 已在 tab 上按结局口径赋值，随 run_overall 继承
ADJ <- if (length(run_overall$adj_ok)) paste(run_overall$adj_ok, collapse = " + ") else "1"
xw_small <- xw[, .SD, .SDcols = c(ID_VARS, "sid_int", "Group", run_overall$adj_ok, "invasive", "sw")]
setnames(xw_small, "sw", "w")
ep <- epi[xw[, .SD, .SDcols = ID_VARS], on = ID_VARS, nomatch = NULL]
ep <- ep[xw_small, on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]
ep[, long_epi := as.integer(duration > 10)]
xw[, lg_dur := log(anes_duration_mins)]

inter_rows <- list()
d1 <- ep[complete.cases(ep[, c("long_epi", "Group", "invasive", run_overall$adj_ok, "sid_int"), with = FALSE])]
if (nrow(d1) > 50) {
  mi <- geeglm(as.formula(paste("long_epi ~ Group * invasive +", ADJ)), data = d1,
               id = sid_int, weights = d1$w, family = binomial(), corstr = "exchangeable")
  cm <- as.matrix(summary(mi)$coefficients)
  k <- "GroupCiprofol_only:invasive"
  inter_rows[[1]] <- data.table(
    model = "M1 prolonged episode (>10 min), episode-level", interaction_term = k,
    coef_log = cm[k, "Estimate"], se = cm[k, "Std.err"], p = 2 * pnorm(-abs(cm[k, "Estimate"] / cm[k, "Std.err"])))
}
d4 <- xw[complete.cases(xw[, c("n_epi", "Group", "invasive", run_overall$adj_ok, "lg_dur"), with = FALSE])]
if (nrow(d4) > 50) {
  mi4 <- geeglm(as.formula(paste("n_epi ~ Group * invasive +", ADJ, "+ offset(lg_dur)")),
                data = d4, id = sid_int, weights = d4$sw, family = poisson(), corstr = "exchangeable")
  cm <- as.matrix(summary(mi4)$coefficients)
  k <- "GroupCiprofol_only:invasive"
  inter_rows[[2]] <- data.table(
    model = "M4 episodes per procedure (offset log duration), procedure-level", interaction_term = k,
    coef_log = cm[k, "Estimate"], se = cm[k, "Std.err"], p = 2 * pnorm(-abs(cm[k, "Estimate"] / cm[k, "Std.err"])))
}
interaction <- rbindlist(inter_rows)
interaction[, ratio := exp(coef_log)]
print(interaction)

## ============================================================================
## 7. 患者级跨层：仅描述，不做任何剔除（分析单位 = 手术）
## ============================================================================
pl <- tab[, .(has_inv = any(invasive == 1), has_non = any(invasive == 0)), by = PATIENT_ID]
cross_pat <- pl[has_inv & has_non]$PATIENT_ID
cat(sprintf("患者级跨层（仅描述，不剔除）：%d 人 / %d 台 = %.2f%%\n",
            length(cross_pat), sum(tab$PATIENT_ID %in% cross_pat),
            100 * mean(tab$PATIENT_ID %in% cross_pat)))

## ============================================================================
## 8. 保存输出
## ============================================================================
fwrite(result_all,    file.path(OUT, "Subgroup_BPSourceMatched_Results.csv"))
fwrite(crude_all,     file.path(OUT, "Subgroup_BPSourceMatched_Crude.csv"))
fwrite(abs_all,       file.path(OUT, "Subgroup_BPSourceMatched_Absolute.csv"))
fwrite(bal_all,       file.path(OUT, "Subgroup_BPSourceMatched_Balance.csv"))
fwrite(ess_all,       file.path(OUT, "Subgroup_BPSourceMatched_ESS.csv"))
fwrite(interaction,   file.path(OUT, "Subgroup_BPSourceMatched_Interaction.csv"))

tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(result_all[, .(
    Population = population, Model = model, Unit = unit,
    Sample_size = n_unit, n_patients = n_patient,
    Estimate = effect, Type = effect_type,
    P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])) |>
    theme_booktabs() |> autofit() |>
    set_caption(caption = paste(
      "Subgroup analysis by intraoperative monitoring modality (invasive vs non-invasive).",
      "Layers were assigned by the actual source of the outcome blood-pressure values",
      "(arterial vs non-invasive cuff channel); all 85,397 procedures were classified.",
      "IPTW re-fitted within each subgroup; weighted GEE clustered by procedure.",
      "M1-M3: episode-level; M4: procedure-level with offset of log anesthesia duration."))
  print(read_docx() |> body_add_flextable(ft),
        target = file.path(OUT, "Subgroup_BPSourceMatched_Results.docx"))
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft2 <- flextable(as.data.frame(interaction[, .(
    Model = model, Interaction = interaction_term,
    Coef = round(coef_log, 4), SE = round(se, 4),
    P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    Ratio_invasive_vs_noninvasive = round(ratio, 3))])) |>
    theme_booktabs() |> autofit() |>
    set_caption(caption = paste(
      "Interaction test: does the ciprofol-vs-propofol effect differ by monitoring modality?",
      "Fitted on the overall cohort with primary-analysis IPTW weights.",
      "Modality is strongly associated with centre, era and case mix; the interaction term",
      "therefore carries those confounders as well. Subgroup-specific estimates are the primary readout."))
  print(read_docx() |> body_add_flextable(ft2),
        target = file.path(OUT, "Subgroup_BPSourceMatched_Interaction.docx"))
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

cat("\n==================== 汇总 ====================\n")
print(result_all[, .(population, model = substr(model, 1, 28), n_unit, effect,
                     p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])
cat("\n交互检验：\n"); print(interaction)
