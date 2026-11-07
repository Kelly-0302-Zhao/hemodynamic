## =============================================================================
## 3-5 GEE analysis -  算倾向评分（PS）和 IPTW 加权权重
## =============================================================================

## ============================================================================
## 第 1 步：准备工作（加载包、设定路径、写死几个全文件共用的常量）
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # 读表、处理表
  library(geepack)      # 本脚本只造权重；geeglm 在 3-6 才会用到，这里先加载备用
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
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")   # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

## 三个 ID 列（患者 / 就诊 / 手术）。
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 时间戳统一按 UTC 解析（原始数据带 Z 后缀）
MY_TZ <- "UTC"

## 倾向评分用的 10 个基线变量。
## 主分析不调中心(Hospital)与年代(year)：中心/年代的结构性差异由敏感性分析验证。
PS_VARS <- c("AGE", "SEX", "BMI", "ASA_GRADE", "EMERGENCY_INDICATOR",
             "operation_type", "Hypertension", "Diabetes", "Coronary_artery_disease",
             "Heart_failure")

## GEE 模型里统一的调整变量
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## 权重截断分位：在每组内部把极端权重按 1% / 99% 分位切掉（见第 4 步坑提示）
TRUNC <- c(.01, .99)

## ============================================================================
## 第 2 步：读入手术级分析表，报告全队列纳入了多少例次、多少病人
## ============================================================================

## 读 3-1 产出的手术级分析表（一行 = 一次手术）
p_tab <- file.path(OUT, "Final_Analysis_Long_Table.rds")
if (!file.exists(p_tab)) stop("找不到 Final_Analysis_Long_Table.rds，请先跑 3-1。")
d <- setDT(readRDS(p_tab))

## 删掉"是否长期低血压"信息缺失的手术
n_before <- nrow(d)
d <- d[!is.na(prolonged)]

cat("====================================================================\n")
cat("3-5 倾向评分与 IPTW 权重\n")
cat(sprintf("   读入手术表 %d 行 → 剔除 prolonged 缺失 %d 行 → 进入分析 %d 行\n",
            n_before, n_before - nrow(d), nrow(d)))
cat(sprintf("   全队列分析人群：手术 %d 例次 | 患者 %d 人\n",
            nrow(d), uniqueN(d$PATIENT_ID)))

## ============================================================================
## 第 3 步：给每个病人打"是否处理组"标记，并选出能进倾向评分的变量
## ============================================================================

## 处理组 = 1（环泊酚），对照组 = 0（丙泊酚）
d[, treat := as.integer(Group == "Ciprofol_only")]

## 逐个变量检查：数值型看标准差是否 > 0；分类型看是否多于一个取值
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

cat(sprintf("   选入倾向评分的变量（%d 个）：%s\n", length(use), paste(use, collapse = ", ")))
if (length(dropped))
  cat(sprintf("   本人群里恒定、已自动剔除：%s\n", paste(dropped, collapse = ", ")))
cat(sprintf("   处理组（Ciprofol）%d 例次 / 对照组（Propofol）%d 例次\n",
            sum(d$treat == 1), sum(d$treat == 0)))

## ============================================================================
## 第 4 步：拟合倾向评分、算稳定化权重、组内截断、标记重叠区
## ============================================================================

## (1) 完整观测子集
cc <- complete.cases(d[, ..use])
dw <- d[cc]
if (nrow(dw) < 100 || min(table(dw$treat)) < 10)
  stop("例数不足或某一组太少，无法拟合倾向评分")

## (2) 逻辑回归拟合 PS
ps_fit <- glm(as.formula(paste("treat ~", paste(use, collapse = " + "))),
              data = dw, family = binomial())
dw[, ps_hat := predict(ps_fit, type = "response")]

## (3) 稳定化权重：当前人群里用环泊酚的真实比例
prev <- mean(dw$treat)
dw[, sw_raw := ifelse(treat == 1, prev / ps_hat, (1 - prev) / (1 - ps_hat))]

## (4) 截断【分组内】做（坑：按总体分位会把处理组压平）
dw[, `:=`(ql = as.numeric(quantile(sw_raw, TRUNC[1])),
          qh = as.numeric(quantile(sw_raw, TRUNC[2]))), by = Group]
dw[, sw := as.numeric(pmin(pmax(sw_raw, ql), qh))]
dw[, c("ql", "qh") := NULL]

## (5) 共同支持域（重叠区）标记
## 安一院从未使用环泊酚，其丙泊酚病例的 PS 会低到 0 附近，IPTW 只能把它们的权重做到约等于 1，无法降权剔除——这是 IPTW 的固有局限，这里仅作标记。
lo_overlap <- max(min(dw$ps_hat[dw$treat == 1]), min(dw$ps_hat[dw$treat == 0]))
hi_overlap <- min(max(dw$ps_hat[dw$treat == 1]), max(dw$ps_hat[dw$treat == 0]))
dw[, in_overlap := ps_hat >= lo_overlap & ps_hat <= hi_overlap]

cat(sprintf("   拟合完成：纳入加权 %d 例次（患者 %d 人）\n",
            nrow(dw), uniqueN(dw$PATIENT_ID)))
cat(sprintf("   PS 范围  Propofol %.3f-%.3f | Ciprofol %.3f-%.3f\n",
            min(dw$ps_hat[dw$treat == 0]), max(dw$ps_hat[dw$treat == 0]),
            min(dw$ps_hat[dw$treat == 1]), max(dw$ps_hat[dw$treat == 1])))
cat(sprintf("   稳定化权重（组内截断 %.0f%%/%.0f%%）范围：%.3f - %.1f\n",
            100 * TRUNC[1], 100 * TRUNC[2], min(dw$sw), max(dw$sw)))


## ============================================================================
## 第 5 步：算加权前后平衡表（每个协变量的标准化均值差 SMD）
## ============================================================================

## 单个变量的加权 SMD。注意 v 必须是数值型；分类型会先变成 0/1 指示变量再算
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

## 看选入的变量，数值型直接算；分类型对每个水平算一次
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

cat("\n---- 协变量平衡表（SMD：加权前 → 加权后）----\n")
print(bal)
cat(sprintf("   加权前最大 SMD %.3f（>0.1 的变量 %d 个）→ 加权后最大 SMD %.3f（%d 个）\n",
            max(bal$smd_before, na.rm = TRUE), sum(bal$smd_before > .1, na.rm = TRUE),
            max(bal$smd_after,  na.rm = TRUE), sum(bal$smd_after  > .1, na.rm = TRUE)))

fwrite(bal, file.path(OUT, "Table1_IPTW_Balance_SMD.csv"))


## ============================================================================
## 第 6 步：算有效样本量（ESS）表，并打印全部诊断信息
## ============================================================================

## Kish 有效样本量
ess <- function(w) sum(w)^2 / sum(w^2)

## ESS 表：按组算 手术例次 / 患者数 / 权重中位 / 权重最大 / ESS
ess_tab <- dw[, .(Group = as.character(first(Group)),
                   n_procedures = .N,
                   n_patients = uniqueN(PATIENT_ID),
                   weight_median = median(sw),
                   weight_max = max(sw),
                   ESS = ess(sw)), by = Group][, Group := NULL]

## ess_tab 行序取决于组出现顺序，取数必须按组名取，不能按位置取
e_c <- ess_tab[Group == "Ciprofol_only", ESS]
e_p <- ess_tab[Group == "Propofol_only", ESS]

cat(sprintf("   [PS]   参与加权手术 %d 例次（患者 %d 人）；PS 范围环泊酚 %.3f-%.3f，Propofol %.3f-%.3f\n",
            nrow(dw), uniqueN(dw$PATIENT_ID),
            min(dw$ps_hat[dw$treat == 1]), max(dw$ps_hat[dw$treat == 1]),
            min(dw$ps_hat[dw$treat == 0]), max(dw$ps_hat[dw$treat == 0])))
cat(sprintf("   [平衡] 加权前最大 SMD %.3f（>0.1 的变量 %d 个）→ 加权后最大 SMD %.3f（%d 个）\n",
            max(bal$smd_before, na.rm = TRUE), sum(bal$smd_before > .1, na.rm = TRUE),
            max(bal$smd_after,  na.rm = TRUE), sum(bal$smd_after  > .1, na.rm = TRUE)))
cat(sprintf("   [权重] 组内截断分位 %.0f%%/%.0f%%；Ciprofol ESS %.0f / %d 例次，Propofol ESS %.0f / %d 例次\n",
            100 * TRUNC[1], 100 * TRUNC[2],
            e_c, sum(dw$treat == 1), e_p, sum(dw$treat == 0)))
cat(sprintf("   [重叠] 两组 PS 重叠区间 = [%.3f, %.3f]；落在区间外的例次 %d（占 %.1f%%）\n",
            lo_overlap, hi_overlap, sum(!dw$in_overlap),
            100 * mean(!dw$in_overlap)))
if (sum(!dw$in_overlap) > 0 && "Hospital" %in% names(dw)) {
  cat("          区间外例次的中心分布：\n")
  print(dw[in_overlap == FALSE, .N, by = .(Hospital, Group)][order(-N)])
}

cat("\n---- 权重描述（分中心）----\n")
print(dw[, .(n_procedures = .N, n_patients = uniqueN(PATIENT_ID),
             weight_median = round(median(sw), 2), weight_max = round(max(sw), 1),
             ESS = round(ess(sw))), by = .(Hospital, Group)][order(Hospital, Group)])

cat("\n---- 各中心对加权的贡献（ESS/例次）----\n")
print(dw[, .(n_procedures = .N, ESS = round(ess(sw)),
             ESS占比 = sprintf("%.1f%%", 100 * ess(sw) / .N)), by = Hospital][order(-n_procedures)])

cat("\n提示：安一院只有丙泊酚组，加权后 ESS 很有限，主结果主要由中大与南医两院决定；这点正文要讲清。\n")

cat("\n---- 有效样本量（ESS）表 ----\n")
print(ess_tab)
fwrite(ess_tab, file.path(OUT, "Table1_IPTW_ESS.csv"))


## ============================================================================
## 第 7 步：生成 Table 1（加权前 / 加权后，两组基线特征对比）
## ============================================================================

## 连续变量和分类变量清单（与论文 Table 1 一致）
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

## 加权后的 均值(标准差) 字符串（本步内部用）
wmean_sd <- function(v, w) {
  m <- weighted.mean(v, w, na.rm = TRUE)
  s <- sqrt(sum(w * (v - m)^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
  sprintf("%.1f (%.1f)", m, s)
}

## 取 SMD 
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
      ## 加权后只给百分比：IPTW 之后"加权例次"是个小数，摆出来会被误读成真有那么多例
      Ciprofol_weighted   = sprintf("%.1f%%", 100 * weighted.mean(ind_a, a$sw)),
      Propofol_weighted   = sprintf("%.1f%%", 100 * weighted.mean(ind_b, b$sw)),
      ## 先按 "变量: 水平" 查（真分类变量），查不到再按变量名查（0/1 的合并症）
      SMD_before = round(pick_smd(c(key, v), "smd_before"), 3),
      SMD_after = round(pick_smd(c(key, v), "smd_after"), 3))
  }
}
tab1 <- rbindlist(rows, fill = TRUE)

cat(sprintf("\n表头：Ciprofol %d 例次 / %d 人；Propofol %d 例次 / %d 人\n",
            dw[Group == "Ciprofol_only", .N], uniqueN(dw[Group == "Ciprofol_only"]$PATIENT_ID),
            dw[Group == "Propofol_only", .N], uniqueN(dw[Group == "Propofol_only"]$PATIENT_ID)))
cat("\n---- Table 1：加权前 / 加权后基线特征 ---- \n")
print(tab1)
fwrite(tab1, file.path(OUT, "Table1_IPTW.csv"))


## ============================================================================
## 第 8 步：把 Table 1 另存一份 Word 版
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
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))


## ============================================================================
## 第 9 步：画 PS 重叠诊断图（两组 PS 分布叠加 + 权重分布箱线图）
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
## 第 10 步：PS 分布分箱（量化重叠区）
## ============================================================================

## 分箱断点
br <- c(0, 0.0338, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001)
lab <- c("<0.034", "0.034-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4",
         "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
dw[, ps_bin := cut(ps_hat, breaks = br, right = FALSE, labels = lab)]

## 把 PS 按箱数出两组各多少例次、各占本组百分比
bin_tab <- function(x) {
  t <- dcast(x[, .N, by = .(ps_bin, Group)], ps_bin ~ Group, value.var = "N", fill = 0)
  setnames(t, c("Propofol_only", "Ciprofol_only"), c("Propofol_n", "Ciprofol_n"))
  t[, `:=`(Propofol_pct = round(100 * Propofol_n / sum(Propofol_n), 1),
           Ciprofol_pct = round(100 * Ciprofol_n / sum(Ciprofol_n), 1))]
  t[]
}

cat("\n---- 全队列 PS 分箱（例次 + 组内占比 %）----\n")
tab_all <- bin_tab(dw)
print(tab_all)
fwrite(tab_all, file.path(OUT, "Check_PS_Distribution.csv"))
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(tab_all)) |> theme_booktabs() |> autofit() |>
    set_caption(caption = "Propensity score distribution by treatment group (full cohort).")
  print(read_docx() |> body_add_flextable(ft),
        target = file.path(OUT, "Check_PS_Distribution.docx"))
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

cat("\n---- 仅重叠区（in_overlap == TRUE）内的 PS 分箱 ----\n")
tab_ov <- bin_tab(dw[in_overlap == TRUE])
print(tab_ov)
fwrite(tab_ov, file.path(OUT, "Check_PS_Distribution_overlap.csv"))
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft <- flextable(as.data.frame(tab_ov)) |> theme_booktabs() |> autofit() |>
    set_caption(caption = "Propensity score distribution within the overlap region only.")
  print(read_docx() |> body_add_flextable(ft),
        target = file.path(OUT, "Check_PS_Distribution_overlap.docx"))
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

cat("\n核心数字（决定主分析该报哪个人群）：\n")
cat(sprintf("  Propofol %d 例次中，PS<0.034（环泊酚根本用不到的区间）有 %d 例次（%.1f%%）\n",
            dw[Group == "Propofol_only", .N],
            dw[Group == "Propofol_only" & !in_overlap, .N],
            100 * dw[Group == "Propofol_only" & !in_overlap, .N] /
              dw[Group == "Propofol_only", .N]))
cat(sprintf("  Ciprofol %d 例次中，PS<0.034 的只有 %d 例次（%.1f%%）\n",
            dw[Group == "Ciprofol_only", .N],
            dw[Group == "Ciprofol_only" & !in_overlap, .N],
            100 * dw[Group == "Ciprofol_only" & !in_overlap, .N] /
              dw[Group == "Ciprofol_only", .N]))

## ============================================================================
## 第 11 步：保存权重表
## ============================================================================
saveRDS(dw[, .(PATIENT_ID, VISIT_ID, OPER_ID, ps_hat, sw, sw_raw, in_overlap,
               Group, Hospital, year)],
        file.path(OUT, "IPTW_Weights.rds"))
