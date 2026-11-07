## =============================================================================
## 2-3 description - — 低血压情况描述性统计（Table 2）
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)   # 读表、处理表
  library(geepack)      # 跑 GEE 模型
})

## 根目录。所有输入输出都挂在这个目录下面
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
DATA <- file.path(D, "JCA_run")   # 清洗好的中间表（master / monitor / event / lab）的家
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")     # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

## 三个 ID 列。三张表 ID 列类型不一致，合并前都要先转成文字型（否则合并会错位）
ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
MY_TZ    <- "UTC"      # 原始时间戳带 Z 后缀，全程按 UTC 解析，保证内部一致
GAP_LIMIT <- 10        # 分钟：相邻测点超过 10 分钟视为"监测中断"，不跨中断累计低血压时长

## ：把 "2024-01-01T08:30:00Z" 这种 ISO 时间字符串转成 POSIXct
parse_iso <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  as.POSIXct(gsub("T", " ", sub("Z$", "", as.character(x))),
             tz = MY_TZ, format = "%Y-%m-%d %H:%M:%S")
}
## 把三列 ID 统一转成文字型
to_chr_ids <- function(DT) {
  DT[, (ID_VARS) := lapply(.SD, as.character), .SDcols = ID_VARS]
  DT
}

cat("[公共设置] 已就绪 | 输出目录：", OUT, "\n")


## ============================================================================
## 第 1 步：读入手术表 / 事件表，并补上分组、中心、年份
## ============================================================================

d   <- setDT(readRDS(file.path(OUT, "Final_Analysis_Surgery.rds")))
epi <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## 读手术基础表（一行 = 一次手术），取分组 / 中心 / 年份用于拆分两组做描述
master <- fread(file.path(DATA, "master_clean_85397obs.csv"), encoding = "UTF-8")
master <- to_chr_ids(master)
master[, anes_start := parse_iso(ANES_START_TIME)]
master[, anes_end   := parse_iso(ANES_END_TIME)]
master[, year       := as.integer(format(as.Date(anes_start), "%Y"))]
master[, Group      := relevel(factor(as.character(Group)), ref = "Propofol_only")]
for (v in c("ASA_GRADE", "operation_type", "Hospital", "SEX")) {
  if (v %in% names(master)) master[[v]] <- factor(master[[v]])
}
info <- master[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, Hospital, year)]

## 把分组 / 中心 / 年份挂到手术表上，再删掉"是否长期低血压"信息缺失的手术
d   <- info[d, on = ID_VARS]
d   <- d[!is.na(prolonged)]
d[, any_hypo := as.integer(n_epi >= 1L)]

cat("手术表：", nrow(d), "次手术（已剔除 prolonged 缺失）\n")
cat("事件表：", nrow(epi), "次低血压 episode\n\n")


## ============================================================================
## 第 2 步：算手术级描述表
## ============================================================================

by_surg <- d[, .(
  `Surgical procedures, n`                = .N,
  `Patients, n`                    = uniqueN(PATIENT_ID),
  `Any IOH (MAP <65 mm Hg), n (%)`        = sprintf("%d (%.1f%%)", sum(any_hypo),
                                          100 * mean(any_hypo)),
  `Prolonged IOH (>10 min), n (%)`     = sprintf("%d (%.1f%%)", sum(prolonged),
                                          100 * mean(prolonged)),
  `Cumulative IOH (>10 min), n (%)`           = sprintf("%d (%.1f%%)", sum(prolonged_cum),
                                          100 * mean(prolonged_cum)),
  `>=2 episodes, n (%)`             = sprintf("%d (%.1f%%)", sum(multi_epi),
                                          100 * mean(multi_epi)),
  `Short episodes only, n (%)`           = sprintf("%d (%.1f%%)", sum(short_only),
                                          100 * mean(short_only)),
  `Episodes per procedure, mean`          = sprintf("%.2f", mean(n_epi)),
  `Cumulative IOH duration, min, median (IQR)` = sprintf("%.1f (%.1f, %.1f)",
                                          median(MAP65_TIME),
                                          quantile(MAP65_TIME, .25),
                                          quantile(MAP65_TIME, .75)),
  `Cumulative IOH burden (AUC), median (IQR)`     = sprintf("%.1f (%.1f, %.1f)",
                                          median(MAP65_AUC),
                                          quantile(MAP65_AUC, .25),
                                          quantile(MAP65_AUC, .75)),
  `Intraoperative MAP readings, median`       = sprintf("%.0f", median(n_map)),
  `Anesthesia duration, min, median`       = sprintf("%.0f", median(anes_duration_mins)),
  `Monitoring density, readings/min, median`    = sprintf("%.2f", median(rate_map))
), by = Group]

cat("---- 表2A：手术级指标（分母 = 手术例次）----\n")
print(by_surg)
cat(sprintf("\n单位提示：手术 %d 例次 | 患者 %d 人\n\n", nrow(d), uniqueN(d$PATIENT_ID)))


## ============================================================================
## 第 3 步：算事件级描述表（表 2B，分母 = 一次低血压 episode）
## ============================================================================

de <- merge(epi, d[, .(PATIENT_ID, VISIT_ID, OPER_ID, Group)], by = ID_VARS)

by_epi <- de[, .(
  `Hypotensive episodes, n`            = .N,
  `Procedures involved, n`            = uniqueN(paste(PATIENT_ID, VISIT_ID, OPER_ID)),
  `Duration per episode, min, median (IQR)`  = sprintf("%.1f (%.1f, %.1f)", median(duration),
                                          quantile(duration, .25), quantile(duration, .75)),
  `Episodes >10 min, n (%)`        = sprintf("%d (%.1f%%)", sum(duration > 10),
                                          100 * mean(duration > 10)),
  `Burden per episode (AUC), median (IQR)`      = sprintf("%.1f (%.1f, %.1f)", median(event_auc),
                                          quantile(event_auc, .25), quantile(event_auc, .75)),
  `Maximum depth, mm Hg, median`      = sprintf("%.1f", median(max_map_depth))
), by = Group]

cat("---- 表2B：事件级指标（分母 = 一次低血压 episode）----\n")
print(by_epi)


## ============================================================================
## 第 4 步：把两张表存盘（csv + Word 各一份）
## ============================================================================

fwrite(by_surg, file.path(OUT, "Table2_SurgeryLevel.csv"))
fwrite(by_epi,  file.path(OUT, "Table2_EpisodeLevel.csv"))

## 存 Word 版
tryCatch({
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  ft1 <- flextable(as.data.frame(by_surg)) |> theme_booktabs() |> autofit() |>
         set_caption(caption = "Table 2A. Intraoperative hypotension by treatment group. Unit of analysis: surgical procedure (one anaesthesia episode).")
  print(read_docx() |> body_add_flextable(ft1),
        target = file.path(OUT, "Table2_SurgeryLevel.docx"))
  ft2 <- flextable(as.data.frame(by_epi)) |> theme_booktabs() |> autofit() |>
         set_caption(caption = "Table 2B. Episode-level characteristics of intraoperative hypotension. Unit of analysis: one hypotensive episode.")
  print(read_docx() |> body_add_flextable(ft2),
        target = file.path(OUT, "Table2_EpisodeLevel.docx"))
}, error = function(e) cat("[Word 输出失败，CSV 已保存]", conditionMessage(e), "\n"))

