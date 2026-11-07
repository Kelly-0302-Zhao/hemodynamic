## =============================================================================
## 3-3 GEE analysis - 未加权（crude）对照
## =============================================================================
## ============================================================================
## 第 1 步：准备工作
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # 读表、处理表
  library(geepack)      # 跑 GEE 模型
})
# 建立一个文件产出路径
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## 三个 ID 列。合并前都要先转成文字型
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## GEE 模型里调整的 4 个基线变量
ADJ_VARS <- c("Hypertension", "AGE", "SEX", "ASA_GRADE")

## 时间戳统一按 UTC 解析
MY_TZ <- "UTC"

## ============================================================================
## 第 2 步：读入两张表，准备"未加权分析人群"（权重统一设成 1）
## ============================================================================
d <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
d <- d[!is.na(prolonged)]          # 只保留"是否长期低血压"信息齐全的手术

epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## 未加权：权重列 sw 恒等于 1（这样后面的模型代码和加权版完全一致）
d[, sw := 1]

cat("纳入：手术", nrow(d), "例次 | 患者", uniqueN(d$PATIENT_ID), "人\n\n")

## ============================================================================
## 第 3 步：模型 M1 —— 每次低血压事件，持续超过 10 分钟的比例，两组比较
##   数据：事件级（一行 = 一次事件），权重=1 即未加权；报告 OR
## ============================================================================
## 把手术表的分组、调整变量、权重挂到每个事件上 → ep：一行 = 一次事件
ep <- epi[d[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                 Hypertension, AGE, SEX, ASA_GRADE, w = sw)],
          on = ID_VARS, nomatch = NULL]
ep <- ep[!is.na(Group) & is.finite(w) & w > 0]

## geeglm 的 id 必须整数，把手术编号（患者+就诊+手术）转成整数编号
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
result_crude <- data.table(
  analysis  = "Unweighted (crude)",
  model     = "M1 Proportion of prolonged hypotensive episodes (>10 min) per episode",
  unit      = "Episode (one hypotensive episode)",
  n_unit    = nrow(d1),                       # 多少次事件
  n_patient = uniqueN(d1$PATIENT_ID),         # 涉及多少病人
  estimate  = exp(est1), lo = exp(est1 - 1.96 * se1), hi = exp(est1 + 1.96 * se1),
  p         = p1,
  effect_type = "OR",
  effect    = sprintf("%.2f (%.2f-%.2f)", exp(est1), exp(est1 - 1.96 * se1), exp(est1 + 1.96 * se1)))

cat("M1 每次事件 >10 min 的比例  OR =", result_crude$effect,
    "  P =", sprintf("%.3g", p1), "\n")
cat("   （", nrow(d1), "次事件，涉及", uniqueN(d1$PATIENT_ID), "位病人）\n")


## ============================================================================
## 第 4 步：模型 M2 —— 每次低血压事件的严重程度总分（AUC），两组比较
##   AUC = 这次事件的"深度 × 时长"总面积，衡量这一次事件整体有多重；
##   先取对数再比较，结果换算回"百分比变化"。报告百分比变化
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

cat("M2 每次事件 AUC 变化      ", result_crude[model %like% "M2", effect],
    "  P =", sprintf("%.3g", p2), "\n")

## ============================================================================
## 第 5 步：模型 M3 —— 每次低血压事件持续几分钟，两组比较
##   做法和 M2 一样，只是把因变量换成事件时长；结果用百分比变化
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

cat("M3 每次事件时长变化       ", result_crude[model %like% "M3", effect],
    "  P =", sprintf("%.3g", p3), "\n")


## ============================================================================
## 第 6 步：模型 M4 —— 每台手术发生几次低血压事件，两组比较
##   数据：手术级（一行 = 一次手术），因变量 n_epi（事件次数）；加 offset(log 麻醉时长)，
##   模型实际比较"每小时发生几次"——即事件次数的"率"。报告 RR
## ============================================================================
adj_ok <- ADJ_VARS[vapply(ADJ_VARS, function(v) {
  vv <- d[[v]]
  if (is.numeric(vv)) is.finite(sd(vv, na.rm = TRUE)) && sd(vv, na.rm = TRUE) > 0
  else uniqueN(as.character(vv[!is.na(vv)])) > 1
}, logical(1))]

d[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
d[, lg_dur  := log(anes_duration_mins)]   # 麻醉时长取对数，做 offset

## data.table 的 with=FALSE 只接受字符向量，所以这里把列名拼成字符向量
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

cat("M4 事件次数率（每小时）   RR =", result_crude[model %like% "M4", effect],
    "  P =", sprintf("%.3g", p4), "\n")


## ============================================================================
## 第 7 步：把四个未加权模型的结果存盘（csv + Word 各一份），并打印描述性统计
## ============================================================================
fwrite(result_crude, file.path(OUT, "Crude_Table_Ciprofol.csv"))

## 存一份 Word 版
save_docx_table <- function(df, path, caption = NULL) {
  ok <- tryCatch({
    suppressPackageStartupMessages({ library(flextable); library(officer) })
    ft <- flextable(as.data.frame(df)) |> theme_booktabs() |> autofit()
    if (!is.null(caption)) ft <- set_caption(ft, caption = caption)
    doc <- read_docx() |> body_add_flextable(ft)
    print(doc, target = path)
    TRUE
  }, error = function(e) {
    cat("  [Word 输出失败，CSV 已保存] ", conditionMessage(e), "\n")
    FALSE
  })
  invisible(ok)
}
save_docx_table(result_crude, file.path(OUT, "Crude_Table_Ciprofol.docx"),
                caption = "Unweighted (crude) GEE models for reference only. M1-M3: one hypotensive episode per row; M4: one surgical procedure per row; all clustered by surgical procedure.")

cat("\n==================== 未加权（crude）结果 ====================\n")
print(result_crude[, .(model, unit, n_unit, n_patient, type = effect_type,
                        effect, p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))])

## 分组描述（手术级）
cat("\n---- 未加权分组描述（手术级）----\n")
print(d[, .(手术例次 = .N,
            患者数 = uniqueN(PATIENT_ID),
            长期低血压率 = sprintf("%.1f%%", 100 * mean(prolonged)),
            累积时长中位 = sprintf("%.1f min", median(MAP65_TIME)),
            累积AUC中位 = sprintf("%.1f", median(MAP65_AUC)),
            事件数每例次 = sprintf("%.2f", mean(n_epi))), by = Group])
