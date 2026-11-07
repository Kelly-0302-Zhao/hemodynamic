## =============================================================================
## 4-4 mediation analysis - blood loss.R
## -----------------------------------------------------------------------------


suppressPackageStartupMessages({ library(data.table) })

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 中介模型和结局模型共用的基线混杂变量（不含术后/暴露后变量，除中介本身）
CONF <- c("AGE", "SEX", "BMI", "ASA_GRADE", "Hypertension", "Diabetes",
          "Coronary_artery_disease", "Heart_failure", "EMERGENCY_INDICATOR",
          "operation_type", "anes_duration_mins")
## 倾向评分用的 11 个基线变量
PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR", "operation_type",
             "Hypertension", "Diabetes", "Coronary_artery_disease", "Heart_failure",
             "anes_duration_mins")

## =============================================================================
## 在给定数据里重新拟合倾向评分 + 稳定化权重 + 1%/99% 截断
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
## g-formula 中介分析（VanderWeele 的 g-computation）
##   中介模型（暴露→中介）：线性，log1p 后标准化
##   结局模型（含中介）：加权 logistic
##   四个反事实均值（风险差尺度）：
##     TE  = ψ(1,m1) − ψ(0,m0)   总效应
##     NIE = ψ(1,m1) − ψ(1,m0)   自然间接效应（经由中介）
##     NDE = ψ(1,m0) − ψ(0,m0)   自然直接效应（不经中介）
##     IIE = ψ(1,m1) − ψ(0,m1)   交互诱导效应
##   95%CI：参数化 bootstrap —— 从两个模型的系数协方差矩阵各抽一组系数，
##   重算一遍反事实均值，重复 B 次，取 2.5%/97.5% 分位。
## =============================================================================
mvrnorm1 <- function(mu, V) {
  V <- (V + t(V)) / 2
  ev <- eigen(V, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = length(ev$values))
  out <- mu + as.numeric(L %*% rnorm(length(mu)))
  names(out) <- names(mu)
  out
}

gf_mediation <- function(dat, fe_col = NULL, B = 2000, seed = 20261106, label = "") {
  dt <- copy(dat)
  conf_use <- intersect(CONF, names(dt))
  dt[, m_raw := log1p(bloodloss_ml)]
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

  ## ---- 中介模型：拆出 X_med（去掉 treat 列），反事实时把 treat 换成 a ----
  df_med <- model.matrix(f_med, data = dt); cn_med <- colnames(df_med)
  i_tm <- which(cn_med == "treat")
  X_med <- df_med[, -i_tm, drop = FALSE]
  b_med <- coef(mm)
  a_treat <- unname(b_med[i_tm]); b_med_rest <- b_med[-i_tm]

  ## ---- 结局模型：X_out 去掉 treat / m_z；其余列与 b_rest 对齐 ----
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

  ## ---- 四个反事实均值 ----
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

  ## ---- 参数化 bootstrap ----
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
    NIE_lo        = round(q(bt_mat[, "NIE"])[1], 5), NIE_hi = round(q(bt_mat[, "NIE"])[2], 5),
    NDE_RD        = round(pt[["NDE"]], 5),
    NDE_lo        = round(q(bt_mat[, "NDE"])[1], 5), NDE_hi = round(q(bt_mat[, "NDE"])[2], 5),
    IIE_RD        = round(pt[["IIE"]], 5),
    IIE_lo        = round(q(bt_mat[, "IIE"])[1], 5), IIE_hi = round(q(bt_mat[, "IIE"])[2], 5),
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
## 主流程
## =============================================================================
cat("====================================================================\n")
cat("4-4 失血量中介分析（g-formula，风险差尺度，参数化 bootstrap 2000 次）\n")
cat("====================================================================\n")

dw <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
for (v in ID_VARS) set(dw, j = v, value = as.character(dw[[v]]))
dw <- dw[!is.na(prolonged)]

## 失血量（手术级）
periop <- fread(file.path(OUT, "_Perioperative_Vars.csv"))
for (v in ID_VARS) set(periop, j = v, value = as.character(periop[[v]]))
periop <- unique(periop, by = ID_VARS)
periop <- periop[, .(PATIENT_ID, VISIT_ID, OPER_ID, bloodloss_ml)]

dw <- merge(dw, periop, by = ID_VARS, all.x = TRUE)
dw <- unique(dw, by = ID_VARS)
dw[, op_year := as.numeric(as.character(year))]
dw[, cell := paste(Hospital, op_year, sep = "_")]

bl_sub <- dw[!is.na(bloodloss_ml) & !is.na(prolonged)]
cat(sprintf("\n全队列 %d 台；失血量有记录 %d 台(%.1f%%)\n",
            nrow(dw), nrow(bl_sub), 100 * nrow(bl_sub) / nrow(dw)))
cat(sprintf("失血量子集 %d 台（环泊酚 %d / 丙泊酚 %d），覆盖 %d 个中心×年格子\n",
            nrow(bl_sub), sum(bl_sub$Group == "Ciprofol_only"),
            sum(bl_sub$Group == "Propofol_only"), uniqueN(bl_sub$cell)))
cat("环泊酚在该子集内出现的中心×年：",
    paste(sort(unique(bl_sub[Group == "Ciprofol_only"]$cell)), collapse = ", "), "\n")

## ---- 两臂共存格子：各≥20 台的中心×年格子 ----
cc_bl <- bl_sub[, .N, by = .(cell, Group)]
piv_bl <- dcast(cc_bl, cell ~ Group, value.var = "N", fill.value = 0)
for (cn in c("Ciprofol_only", "Propofol_only"))
  if (!cn %in% names(piv_bl)) piv_bl[, (cn) := 0L]
piv_bl[, Ciprofol_only := as.integer(as.character(Ciprofol_only))]
piv_bl[, Propofol_only  := as.integer(as.character(Propofol_only))]
piv_bl[is.na(Ciprofol_only), Ciprofol_only := 0L]
piv_bl[is.na(Propofol_only),  Propofol_only  := 0L]
piv_bl[, within := pmin(Ciprofol_only, Propofol_only)]
both_bl <- piv_bl[within >= 20]$cell
cat(sprintf("\n两臂共存且各≥20台的中心×年格子：%s\n", paste(both_bl, collapse = ", ")))
cat(sprintf("★ 信息量：环泊酚 %d 台；能与同格丙泊酚配对的仅 %d 对\n",
            sum(piv_bl$Ciprofol_only), sum(piv_bl[Ciprofol_only >= 1, within])))

## ---- 口径① 粗口径（主分析） ----
bl <- fit_iptw(bl_sub)
res1 <- gf_mediation(bl, fe_col = NULL, B = 2000, label = "① 粗口径")
cat("\n---- ① 粗口径 ----\n")
print(res1[, .(spec, n_surgery, n_cipro, n_prop,
               a_exposure_to_mediator, b_mediator_to_outcome,
               TE_RD, TE_lo, TE_hi, NIE_RD, NIE_lo, NIE_hi,
               NDE_RD, NDE_lo, NDE_hi,
               prop_mediated, prop_lo, prop_hi, or_total_model, verdict)])

## ---- 口径② ① + 中心×年固定效应 ----
res2 <- gf_mediation(bl, fe_col = "cell", B = 2000, label = "② ①+中心×年FE")
cat("\n---- ② ①+ 中心×年固定效应 ----\n")
print(res2[, .(spec, n_surgery, n_cipro, n_prop,
               a_exposure_to_mediator, b_mediator_to_outcome,
               TE_RD, TE_lo, TE_hi, NIE_RD, NIE_lo, NIE_hi,
               NDE_RD, NDE_lo, NDE_hi,
               prop_mediated, prop_lo, prop_hi, or_total_model, verdict)])

## ---- 口径③ 两臂共存格子 ----
bl_co <- fit_iptw(bl_sub[cell %in% both_bl])
res3 <- gf_mediation(bl_co, fe_col = NULL, B = 2000, label = "③ 两臂共存格子")
cat("\n---- ③ 两臂共存格子 ----\n")
print(res3[, .(spec, n_surgery, n_cipro, n_prop,
               a_exposure_to_mediator, b_mediator_to_outcome,
               TE_RD, TE_lo, TE_hi, NIE_RD, NIE_lo, NIE_hi,
               NDE_RD, NDE_lo, NDE_hi,
               prop_mediated, prop_lo, prop_hi, or_total_model, verdict)])

## ---- 口径④ ③ + 中心×年固定效应 ----
res4 <- gf_mediation(bl_co, fe_col = "cell", B = 2000, label = "④ ③+中心×年FE")
cat("\n---- ④ ③+ 中心×年固定效应 ----\n")
print(res4[, .(spec, n_surgery, n_cipro, n_prop,
               a_exposure_to_mediator, b_mediator_to_outcome,
               TE_RD, TE_lo, TE_hi, NIE_RD, NIE_lo, NIE_hi,
               NDE_RD, NDE_lo, NDE_hi,
               prop_mediated, prop_lo, prop_hi, or_total_model, verdict)])

## ---- 辅助诊断----
cat("\n---- 辅助诊断①：暴露×中介交互检验 ----\n")
cat("   交互不显著 ⇒ 无交互假设成立，g-formula 的反事实分解无偏。\n")
for (nm in c("粗口径", "①+中心×年FE")) {
  fc <- if (nm == "粗口径") NULL else "cell"
  dat <- copy(bl)
  conf_use <- intersect(CONF, names(dat))
  dat[, m_z := as.numeric(scale(log1p(bloodloss_ml)))]
  kv <- c("prolonged", "Group", "m_z", conf_use)
  dat <- dat[complete.cases(dat[, ..kv])]
  dat[, treat := as.integer(Group == "Ciprofol_only")]
  fconf <- paste(conf_use, collapse = " + ")
  feT <- if (!is.null(fc)) paste0(" + factor(", fc, ")") else ""
  mo_i <- glm(as.formula(paste("prolonged ~ treat * m_z +", fconf, feT)),
              data = dat, family = binomial(), weights = dat$sw)
  cat(sprintf("   [%s] treat:m_z = %.4f (SE %.4f, p = %.3g)\n", nm,
              unname(coef(mo_i)["treat:m_z"]),
              unname(summary(mo_i)$coefficients["treat:m_z", 2]),
              unname(summary(mo_i)$coefficients["treat:m_z", 4])))
}

cat("\n---- 辅助诊断②：a 路径（暴露→失血量）粗 vs 加中心×年固定效应 ----\n")
bl_a <- copy(bl_sub)
bl_a[, treat := as.integer(Group == "Ciprofol_only")]
bl_a[, m_z := as.numeric(scale(log1p(bloodloss_ml)))]
m_crude <- lm(m_z ~ treat, data = bl_a)
m_fe    <- lm(m_z ~ treat + factor(cell), data = bl_a)
a_crude <- unname(coef(m_crude)["treat"]); a_fe <- unname(coef(m_fe)["treat"])
cat(sprintf("   粗 a = %.4f (SE %.4f, p=%.3g)\n", a_crude,
            summary(m_crude)$coefficients["treat", 2], summary(m_crude)$coefficients["treat", 4]))
cat(sprintf("   加FE a = %.4f (SE %.4f, p=%.3g)  ← 衰减 %.0f%%\n", a_fe,
            summary(m_fe)$coefficients["treat", 2], summary(m_fe)$coefficients["treat", 4],
            100 * (1 - a_fe / a_crude)))
cat("   ⇒ a 路径衰减 ≠ 中介不存在：g-formula 的 NIE 四个口径方向一致、CI 均不跨 0。\n")

cat("\n---- 辅助诊断③：总效应在中心×年分层后趋零 ⇒ 中介比例不可解释 ----\n")
bl_t <- copy(bl_sub)
bl_t[, treat := as.integer(Group == "Ciprofol_only")]
bl_t <- fit_iptw(bl_t)
t_crude <- glm(prolonged ~ treat, data = bl_t, family = binomial(), weights = bl_t$sw)
t_fe    <- glm(prolonged ~ treat + factor(cell), data = bl_t, family = binomial(), weights = bl_t$sw)
cat(sprintf("   无固定效应 总效应 OR = %.3f (%.3f~%.3f)\n", exp(coef(t_crude)["treat"]),
            exp(confint(t_crude)["treat", 1]), exp(confint(t_crude)["treat", 2])))
cat(sprintf("   加中心×年FE 总效应 OR = %.3f (%.3f~%.3f)\n", exp(coef(t_fe)["treat"]),
            exp(confint(t_fe)["treat", 1]), exp(confint(t_fe)["treat", 2])))
cat("   ⇒ 中介比例 = 间接效应/总效应，是比值：分母趋零时必然发散，是数学性质，\n")
cat("     不代表中介不存在 ⇒ 报告 NIE/NDE/TE 绝对值。\n")

## ---- 交付表 ----
med_tab <- rbindlist(list(res1, res2, res3, res4), fill = TRUE)
write_atomic <- function(tab, fname) {
  tmp <- file.path(OUT, paste0("_stag_", fname))
  fwrite(tab, tmp)
  dst <- file.path(OUT, fname)
  ok <- suppressWarnings(file.copy(tmp, dst, overwrite = TRUE))
  if (ok) { unlink(tmp); cat("\n已写出：", dst, "\n") }
  else {
    cat("\n!! 目标文件被占用（可能正被 Excel/WPS 打开），未能覆盖：", dst, "\n")
    cat("!! 结果已保存在暂存文件：", tmp, "\n")
  }
}
write_atomic(med_tab, "Mediation_bloodloss.csv")
cat(sprintf("\n表行数：%d（四个口径）\n", nrow(med_tab)))
