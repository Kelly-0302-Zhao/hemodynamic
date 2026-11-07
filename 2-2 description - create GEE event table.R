## =============================================================================
## 2-2 description - create GEE event table.R —— 重建术中低血压事件表（事件级 + 手术级）
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

## 三个 ID 列。三张表 ID 列类型不一致，合并前都要先转成文字型
ID_VARS  <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")
MY_TZ    <- "UTC"      # 原始时间戳带 Z 后缀，全程按 UTC 解析，保证内部一致
GAP_LIMIT <- 10        # 分钟：相邻测点超过 10 分钟视为"监测中断"，不跨中断累计低血压时长


## 把 "2024-01-01T08:30:00Z" 这种 ISO 时间字符串转成 POSIXct
parse_iso <- function(x) {
  ## 已是 POSIXct 就直接返回，不再过 as.character()：
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
## 第 1 步：读入两张基础表（手术表 + 术中监测表）
## ============================================================================

## 手术表：一行 = 一次手术，已含分组、协变量、麻醉起止时间、麻醉时长
master <- fread(file.path(DATA, "master_clean_85397obs.csv"), encoding = "UTF-8") 
master <- to_chr_ids(master)
master[, anes_start := parse_iso(ANES_START_TIME)]
master[, anes_end   := parse_iso(ANES_END_TIME)]
master[, year       := as.integer(format(as.Date(anes_start), "%Y"))]
master[, Group      := relevel(factor(as.character(Group)), ref = "Propofol_only")]
for (v in c("ASA_GRADE", "operation_type", "Hospital", "SEX")) {
  if (v %in% names(master)) master[[v]] <- factor(master[[v]])
}
cat("手术表：", nrow(master), "次手术，", uniqueN(master$PATIENT_ID), "位患者\n")

## 术中监测表：约 7246 万条原始记录（同一分钟多条，全部保留不去重）。只读需要的列（三列 ID + 时间 + MAP）
bp <- fread(file.path(DATA, "monitor_clean_72460818obs.csv"), encoding = "UTF-8",
            select = c(ID_VARS, "TIME_POINT", "MAP"))
bp <- to_chr_ids(bp)
setnames(bp, "TIME_POINT", "t_raw")
bp[, t := parse_iso(t_raw)]
bp[, t_raw := NULL]
cat("监测表：", nrow(bp), "个测点\n\n")


## ============================================================================
## 第 2 步：先做时间基准自检，再施加术中时间窗过滤
## ============================================================================

## 把每台手术的麻醉起止时间挂到每个测点上
bp <- master[, .(PATIENT_ID, VISIT_ID, OPER_ID, anes_start, anes_end)][bp, on = ID_VARS]
n_all <- nrow(bp)
bp[, in_window    := t >= anes_start & t <= anes_end]
bp[, pos_min      := as.numeric(difftime(t, anes_start, units = "mins"))]
bp[, before_start := t < anes_start]
bp[, after_end    := t > anes_end]

cat("---- 时间基准自检（测点相对麻醉开始时间的位置）----\n")
cat(sprintf("master 里 ANES_START_TIME 缺失：%d 次手术；ANES_END_TIME 缺失：%d 次\n",
            master[is.na(anes_start), .N], master[is.na(anes_end), .N]))
cat(sprintf("测点位置：麻醉开始前 %.2f%% | 术中窗口内 %.2f%% | 麻醉结束后 %.2f%%\n",
            100 * mean(bp$before_start, na.rm = TRUE),
            100 * mean(bp$in_window,    na.rm = TRUE),
            100 * mean(bp$after_end,    na.rm = TRUE)))
cat("测点相对麻醉开始时间的分位数（分钟，负值 = 麻醉前）：\n")
print(round(quantile(bp$pos_min, c(.001, .01, .05, .25, .50, .75, .95, .99, .999),
                     na.rm = TRUE), 1))

## 逐台手术看术中窗口内到底有没有测点
ov <- bp[!is.na(anes_start) & !is.na(anes_end),
         .(n_point = .N, n_in = sum(in_window, na.rm = TRUE)), by = ID_VARS]
pct_cov <- 100 * ov[n_in > 0, .N] / nrow(ov)
cat(sprintf("有麻醉起止时间的手术中，术中窗口内至少 1 个测点的比例：%.2f%%\n", pct_cov))

if (!is.finite(pct_cov) || pct_cov < 50) {
  stop("术中窗口内几乎没有测点 —— monitor 与 master 的时间基准很可能不一致（时区差/单位差）。\n",
       "       请先核对两表的时间列定义，不要直接套用术中窗口过滤。")
}

bp[, bin_lo := pmin(pmax(floor(pos_min / 10) * 10, -180), 480)]
chk <- bp[!is.na(bin_lo), .(n_point = .N), by = bin_lo][order(bin_lo)]
chk[, pct := round(100 * n_point / sum(n_point), 3)]
fwrite(chk, file.path(OUT, "Check_Window_Position.csv"))

cat(sprintf("\n时间窗过滤：%d → %d 个测点（剔除 %.1f%%）\n",
            n_all, sum(bp$in_window, na.rm = TRUE),
            100 * (1 - sum(bp$in_window, na.rm = TRUE) / n_all)))
bp <- bp[in_window == TRUE]
bp[, `:=`(in_window = NULL, pos_min = NULL, before_start = NULL,
          after_end = NULL, bin_lo = NULL)]


## ============================================================================
## 第 3 步：同一分钟有多条记录 -> 折叠成"一分钟一个测点"（取该分钟内的最低 MAP）
## ============================================================================
bp[, MAP := {
  v <- MAP[!is.na(MAP)]
  if (length(v)) min(v) else NA_real_
}, by = c(ID_VARS, "t")]
bp[, n_raw := .N, by = c(ID_VARS, "t")]

n_before_agg <- nrow(bp)
bp <- unique(bp, by = c(ID_VARS, "t"))
cat(sprintf("同一分钟多条记录折叠：%d → %d 个测点（平均每点 %.2f 条原始记录）；MAP 缺失 %d 个\n\n",
            n_before_agg, nrow(bp), n_before_agg / nrow(bp), sum(is.na(bp$MAP))))


## ============================================================================
## 第 4 步：算采样间隔，给每个测点定"覆盖时长"
## ============================================================================

setorder(bp, PATIENT_ID, VISIT_ID, OPER_ID, t)
bp[, idx   := seq_len(.N), by = ID_VARS]
bp[, next_t := shift(t, type = "lead"), by = ID_VARS]
bp[, delta  := as.numeric(difftime(next_t, t, units = "mins"))]

## 该次手术的"典型采样间隔"（只统计正常范围内 0<间隔<=10 分钟的那些点）
bp[, med_delta := as.numeric(median(delta[delta > 0 & delta <= GAP_LIMIT], na.rm = TRUE)),
   by = ID_VARS]
bp[!is.finite(med_delta) | med_delta <= 0, med_delta := 5]

## 单点覆盖时长：正常间隔用真实间隔；缺失或超过中断阈值，只用典型间隔
bp[, cover := ifelse(is.na(delta) | delta <= 0 | delta > GAP_LIMIT, med_delta, delta)]
bp[cover > GAP_LIMIT, cover := GAP_LIMIT]      # 兜底封顶，任何单点最多代表 10 分钟

cat(sprintf("典型采样间隔中位数：%.1f 分钟\n", median(bp$med_delta, na.rm = TRUE)))
cat(sprintf("测点间隔 >%d 分钟（判定为监测中断、会切断事件）的比例：%.3f%%；MAP 缺失测点 %d 个\n",
            GAP_LIMIT, 100 * mean(bp$delta > GAP_LIMIT, na.rm = TRUE), sum(is.na(bp$MAP))))

## ============================================================================
## 第 5 步：划分"连续片段"，切出每个低血压 episode（事件级表）
##   事件定义：MAP<65 的连续测点算 1 次事件（单点也算 1 次）；只有 MAP 回升到 ≥65、或相邻测点间隔 >10 分钟（监测中断）才另起一次。
## ============================================================================

bp[, is_low := as.integer(MAP < 65)]
bp[, depth  := pmax(0, 65 - MAP)]
bp[, step_auc := depth * cover]

## 开新片段的条件：本次手术第一个点 / 低血压状态发生变化 / 真实的监测中断（间隔 > 10 分钟）
bp[, brk := as.integer(idx == 1L | is_low != shift(is_low) | (!is.na(delta) & delta > GAP_LIMIT))]
bp[, seg := cumsum(brk), by = ID_VARS]

## 事件级：一行 = 一次低血压 episode
epi_episodes <- bp[is_low == 1L, .(
  start_time = min(t),
  end_time   = max(t),
  duration   = sum(cover, na.rm = TRUE),        # 逐点覆盖时长求和（不跨中断）
  event_auc  = sum(step_auc, na.rm = TRUE),
  n_point    = .N,
  max_map_depth = max(depth, na.rm = TRUE)
), by = c(ID_VARS, "seg")]
epi_episodes[, seg := NULL]

cat("\n事件级表（一行 = 一次低血压 episode）：", nrow(epi_episodes), "个 episode\n")
cat(sprintf("  episode 时长：中位 %.1f 分钟，>10 分钟的占 %.1f%%\n",
            median(epi_episodes$duration), 100 * mean(epi_episodes$duration > 10)))
cat(sprintf("  其中只采到 1 个测点的单点事件 %d 个（%.1f%%）—— 它们同样各算 1 次事件\n",
            sum(epi_episodes$n_point == 1L), 100 * mean(epi_episodes$n_point == 1L)))


## ============================================================================
## 第 6 步：汇总手术级表
## ============================================================================

surg <- bp[, .(
  n_epi        = uniqueN(seg[is_low %in% 1L]),                 # ★ 事件次数（M4 的因变量）
  MAP65_TIME   = sum(cover * is_low, na.rm = TRUE),            # 累积低于 65 的分钟数
  MAP65_AUC    = sum(step_auc, na.rm = TRUE),                  # 累积低血压面积（负荷）
  n_map        = .N,                                           # 术中有效测点数
  first_low_t  = suppressWarnings(min(t[is_low == 1L]))
), by = ID_VARS]

## 补上"完全没有低血压"的手术（episode 数 = 0）
surg <- merge(master[, .(PATIENT_ID, VISIT_ID, OPER_ID, anes_start, anes_end,
                         anes_duration_mins)], surg, by = ID_VARS, all.x = TRUE)

surg[, no_monitor := is.na(n_map)]

surg[is.na(n_epi) & !is.na(anes_start), `:=`(n_epi = 0L, MAP65_TIME = 0,
                                             MAP65_AUC = 0, n_map = 0L)]
surg[, rate_map := as.numeric(NA)]
surg[!no_monitor & !is.na(anes_duration_mins),
     rate_map := n_map / anes_duration_mins]

## 结局变量（全部定义在"一次手术"上）
surg[, n_long     := NA_integer_]
surg[, prolonged  := NA_integer_]     # 该次手术是否存在 >10 分钟的低血压
surg[, prolonged_cum := NA_integer_]  # 该次手术累积低于 65 的时长是否 >10 分钟
surg[, multi_epi  := NA_integer_]     # 该次手术是否发生 ≥2 次低血压
surg[, short_only := NA_integer_]     # 该次手术只有短事件、没有长事件

agg <- epi_episodes[, .(n_long = sum(duration > 10, na.rm = TRUE)), by = ID_VARS]
surg <- merge(surg, agg, by = ID_VARS, all.x = TRUE, suffixes = c("", ".epi"))
surg[, n_long := as.integer(fcoalesce(n_long.epi, 0L))][, n_long.epi := NULL]

surg[!is.na(n_epi), `:=`(
  prolonged     = as.integer(n_long >= 1L),
  prolonged_cum = as.integer(MAP65_TIME > 10),
  multi_epi     = as.integer(n_epi >= 2L),
  short_only    = as.integer(n_epi >= 1L & n_long == 0L))]

## 麻醉时间窗缺失、或术中窗口内无任何测点的手术：结局无法定义，置为 NA （后续所有模型会自动剔除这类手术，并如实报告例数）
bad <- surg[is.na(anes_start) | is.na(anes_end) | no_monitor]
if (nrow(bad) > 0) {
  cat(sprintf("\n[注意] %d 次手术无法定义术中结局（缺麻醉起止时间 %d 例 + 术中无测点 %d 例）：\n",
              nrow(bad), surg[is.na(anes_start) | is.na(anes_end), .N], surg[no_monitor == TRUE, .N]))
  print(master[bad[, .(PATIENT_ID, VISIT_ID, OPER_ID)], on = ID_VARS,
               nomatch = NULL][, .N, by = .(Hospital, Group)][order(-N)])
}
surg[is.na(anes_start) | is.na(anes_end) | no_monitor, `:=`(
  n_epi = NA_integer_, n_long = NA_integer_, prolonged = NA_integer_,
  prolonged_cum = NA_integer_, multi_epi = NA_integer_, short_only = NA_integer_)]


## ============================================================================
## 第 7 步：和旧流程（data_ioh.csv）做一次手术级对照，确认数字合理
## ============================================================================

old_path <- file.path(D, "data_ioh.csv")
if (file.exists(old_path)) {
  old <- fread(old_path, encoding = "UTF-8")
  old <- to_chr_ids(old)
  cmp <- merge(surg[, .(PATIENT_ID, VISIT_ID, OPER_ID, new_time = MAP65_TIME,
                        new_auc = MAP65_AUC, new_n = n_epi)],
               old[, .(PATIENT_ID, VISIT_ID, OPER_ID, old_time = MAP65_TIME,
                       old_auc = MAP65_AUC, old_ioh = MAP65_IOH)],
               by = ID_VARS, all.x = TRUE)
  cat("\n---- 新旧口径对照（手术级，中位数）----\n")
  cat(sprintf("  累积低于 65 时长：旧 %.1f 分钟 → 新 %.1f 分钟\n",
              median(cmp$old_time, na.rm = TRUE), median(cmp$new_time, na.rm = TRUE)))
  cat(sprintf("  累积低血压 AUC  ：旧 %.1f → 新 %.1f\n",
              median(cmp$old_auc, na.rm = TRUE), median(cmp$new_auc, na.rm = TRUE)))
  cat(sprintf("  旧口径有 IOH 的手术：%d（%.1f%%）\n",
              sum(cmp$old_ioh == 1, na.rm = TRUE), 100 * mean(cmp$old_ioh == 1, na.rm = TRUE)))
  cat(sprintf("  新口径有 IOH 的手术：%d（%.1f%%）\n",
              sum(cmp$new_n > 0, na.rm = TRUE), 100 * mean(cmp$new_n > 0, na.rm = TRUE)))
  fwrite(cmp, file.path(OUT, "Check_EventTable_Old_vs_New.csv"))
} else {
  cat("\n[注意] 未找到旧流程对照表：", old_path, "\n")
  cat("       最近一份旧流程 data_ioh.csv 在 20260623/ 下；因其队列快照与本轮不同，\n")
  cat("       这里不做强行对照，跳过新旧口径对照（不影响任何主结果）。\n")
  cat("       如需这项 QC，把对应快照的 data_ioh.csv 复制到 ", D, " 再重跑即可。\n")
}

cat("\n=========== 手术级结局分布（一行 = 一次手术）===========\n")
print(surg[!is.na(prolonged), .(n_surgery = .N,
                                IOH率 = sprintf("%.1f%%", 100 * mean(n_epi > 0)),
                                `长期IOH率>10min` = sprintf("%.1f%%", 100 * mean(prolonged == 1)),
                                `≥2次事件率` = sprintf("%.1f%%", 100 * mean(multi_epi == 1)),
                                `事件次数/例` = sprintf("%.2f", mean(n_epi)),
                                `累积时长中位` = sprintf("%.1f", median(MAP65_TIME)))])
cat("\n患者数一并给出：手术", nrow(surg[!is.na(prolonged)]), "例次，患者",
    uniqueN(surg[!is.na(prolonged)]$PATIENT_ID), "人\n")


## ============================================================================
## 第 8 步：把两张主表存盘
## ============================================================================

saveRDS(epi_episodes, file.path(OUT, "Final_Analysis_Events.rds"))
saveRDS(surg,         file.path(OUT, "Final_Analysis_Surgery.rds"))

