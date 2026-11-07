## =============================================================================
## forest.R —— 主结果森林图（未加权对照 / 全队列 IPTW / 重叠区对照，四个结局）
## =============================================================================

## ============================================================================
## 第 1 步：准备工作（加载包、设定路径、确认输出目录）
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # 读 csv、做表
})

## 根目录。所有输入输出都挂在这个目录下面
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")   # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

cat("输出目录：", OUT, "\n")


## ============================================================================
## 第 2 步：读入三份结果表（未加权 / 全队列 IPTW / 重叠区对照）
## ============================================================================
crude <- fread(file.path(OUT, "Crude_Table_Ciprofol.csv"))
main  <- fread(file.path(OUT, "Analysis_Results_main.csv"))
ov_path <- file.path(OUT, "Analysis_Results_overlap.csv")
ov <- if (file.exists(ov_path)) fread(ov_path) else NULL

cat(sprintf("读入：crude %d 行 | 全队列主结果 %d 行 | 重叠区对照 %s\n",
            nrow(crude), nrow(main),
            if (is.null(ov)) "无（不画这一条）" else paste0(nrow(ov), " 行")))


## ============================================================================
## 第 3 步：把"百分比变化"行换算成"比值"尺度，拼成统一绘图数据
## ============================================================================

## Pct 行换算成比值；OR/RR 行比值就是 estimate 本身
crude <- crude[, `:=`(
  ratio = ifelse(effect_type == "Pct", 1 + estimate / 100, estimate),
  rlo   = ifelse(effect_type == "Pct", 1 + lo / 100, lo),
  rhi   = ifelse(effect_type == "Pct", 1 + hi / 100, hi))]
crude <- crude[!is.finite(ratio) | !is.finite(rlo) | !is.finite(rhi),
               `:=`(ratio = NA_real_, rlo = NA_real_, rhi = NA_real_)]

main <- main[, `:=`(
  ratio = ifelse(effect_type == "Pct", 1 + estimate / 100, estimate),
  rlo   = ifelse(effect_type == "Pct", 1 + lo / 100, lo),
  rhi   = ifelse(effect_type == "Pct", 1 + hi / 100, hi))]
main <- main[!is.finite(ratio) | !is.finite(rlo) | !is.finite(rhi),
             `:=`(ratio = NA_real_, rlo = NA_real_, rhi = NA_real_)]

if (!is.null(ov)) {
  ov <- ov[, `:=`(
    ratio = ifelse(effect_type == "Pct", 1 + estimate / 100, estimate),
    rlo   = ifelse(effect_type == "Pct", 1 + lo / 100, lo),
    rhi   = ifelse(effect_type == "Pct", 1 + hi / 100, hi))]
  ov <- ov[!is.finite(ratio) | !is.finite(rlo) | !is.finite(rhi),
           `:=`(ratio = NA_real_, rlo = NA_real_, rhi = NA_real_)]
}

## 四个模型按 M1→M2→M3→M4 顺序画
MODELS <- unique(main$model)
cat(sprintf("要画的四个模型：%s\n", paste(MODELS, collapse = "  |  ")))
if (length(MODELS) < 4L) stop("主结果表里只读到 ", length(MODELS), " 个模型，期望 4 个（M1~M4）")

MT <- c(
  "M1 长期低血压事件比例（每次事件 >10 min）"         = "M1 Prolonged hypotensive episodes (>10 min)",
  "M2 单次事件低血压负荷（事件 AUC）"                 = "M2 Burden of a single episode (AUC)",
  "M3 每次低血压持续时间"                             = "M3 Duration per hypotensive episode",
  "M4 低血压事件次数（事件/例次，含麻醉时长 offset）" = "M4 Hypotensive episodes per procedure")
mt <- function(m) if (m %in% names(MT)) MT[[m]] else m


## ============================================================================
## 第 4 步：打开图片画布，画 M1 面板（每次事件 >10 分钟的长期低血压比例，OR）
## ============================================================================

png(file.path(OUT, "ForestPlot_main.png"), width = 1700, height = 1100, res = 150)
par(mfrow = c(2, 2), mar = c(4, 15, 3, 2), oma = c(3.5, 0, 0, 0))

mm <- MODELS[1]
rows <- list(list(lab = "Unweighted (crude)", r = crude[model == mm]),
             list(lab = "IPTW, full cohort", r = main[model == mm]))
if (!is.null(ov)) rows[[length(rows) + 1]] <- list(lab = "IPTW, PS-overlap only", r = ov[model == mm])
rows <- rows[vapply(rows, function(z) nrow(z$r) > 0, logical(1))]
est <- sapply(rows, function(z) z$r$ratio)
lo  <- sapply(rows, function(z) z$r$rlo)
hi  <- sapply(rows, function(z) z$r$rhi)
ok  <- is.finite(est) & is.finite(lo) & is.finite(hi)
if (!any(ok)) {
  plot.new(); title(main = paste(mt(mm), "(no valid estimate)"))
  cat(sprintf("  %s：没有可画的点\n", mm))
} else {
  xs <- range(c(lo[ok], hi[ok], 1))
  n  <- length(rows)
  plot(NA, xlim = log(xs), ylim = c(0.4, n + 0.6), yaxt = "n", xaxt = "n",
       xlab = "OR / RR (log scale)", ylab = "", main = mt(mm))
  abline(v = 0, lty = 2, col = "grey50")
  ticks <- c(0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
  ticks <- ticks[ticks >= xs[1] * 0.9 & ticks <= xs[2] * 1.1]
  axis(1, at = log(ticks), labels = as.character(ticks))
  axis(2, at = n:1, labels = sapply(rows, function(z) z$lab), las = 1, cex.axis = 0.85)
  for (i in seq_len(n)) {
    if (!ok[i]) next
    y <- n - i + 1
    col <- c("grey35", "#B2182B", "#2166AC")[i]
    segments(log(lo[i]), y, log(hi[i]), y, lwd = 2, col = col)
    points(log(est[i]), y, pch = 18, cex = 1.6, col = col)
  }
  cat(sprintf("  M1 面板：画了 %d 条（crude / 全队列%s）\n",
              n, if (is.null(ov)) "" else " / 重叠区"))
}


## ============================================================================
## 第 5 步：画 M2 面板
## ============================================================================
mm <- MODELS[2]
rows <- list(list(lab = "Unweighted (crude)", r = crude[model == mm]),
             list(lab = "IPTW, full cohort", r = main[model == mm]))
if (!is.null(ov)) rows[[length(rows) + 1]] <- list(lab = "IPTW, PS-overlap only", r = ov[model == mm])
rows <- rows[vapply(rows, function(z) nrow(z$r) > 0, logical(1))]
est <- sapply(rows, function(z) z$r$ratio)
lo  <- sapply(rows, function(z) z$r$rlo)
hi  <- sapply(rows, function(z) z$r$rhi)
ok  <- is.finite(est) & is.finite(lo) & is.finite(hi)
if (!any(ok)) {
  plot.new(); title(main = paste(mt(mm), "(no valid estimate)"))
  cat(sprintf("  %s：没有可画的点\n", mm))
} else {
  xs <- range(c(lo[ok], hi[ok], 1))
  n  <- length(rows)
  plot(NA, xlim = log(xs), ylim = c(0.4, n + 0.6), yaxt = "n", xaxt = "n",
       xlab = "Ratio of geometric mean (log scale)", ylab = "", main = mt(mm))
  abline(v = 0, lty = 2, col = "grey50")
  ticks <- c(0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
  ticks <- ticks[ticks >= xs[1] * 0.9 & ticks <= xs[2] * 1.1]
  axis(1, at = log(ticks), labels = as.character(ticks))
  axis(2, at = n:1, labels = sapply(rows, function(z) z$lab), las = 1, cex.axis = 0.85)
  for (i in seq_len(n)) {
    if (!ok[i]) next
    y <- n - i + 1
    col <- c("grey35", "#B2182B", "#2166AC")[i]
    segments(log(lo[i]), y, log(hi[i]), y, lwd = 2, col = col)
    points(log(est[i]), y, pch = 18, cex = 1.6, col = col)
  }
  cat(sprintf("  M2 面板：画了 %d 条\n", n))
}


## ============================================================================
## 第 6 步：画 M3 面板
## ============================================================================
mm <- MODELS[3]
rows <- list(list(lab = "Unweighted (crude)", r = crude[model == mm]),
             list(lab = "IPTW, full cohort", r = main[model == mm]))
if (!is.null(ov)) rows[[length(rows) + 1]] <- list(lab = "IPTW, PS-overlap only", r = ov[model == mm])
rows <- rows[vapply(rows, function(z) nrow(z$r) > 0, logical(1))]
est <- sapply(rows, function(z) z$r$ratio)
lo  <- sapply(rows, function(z) z$r$rlo)
hi  <- sapply(rows, function(z) z$r$rhi)
ok  <- is.finite(est) & is.finite(lo) & is.finite(hi)
if (!any(ok)) {
  plot.new(); title(main = paste(mt(mm), "(no valid estimate)"))
  cat(sprintf("  %s：没有可画的点\n", mm))
} else {
  xs <- range(c(lo[ok], hi[ok], 1))
  n  <- length(rows)
  plot(NA, xlim = log(xs), ylim = c(0.4, n + 0.6), yaxt = "n", xaxt = "n",
       xlab = "Ratio of geometric mean (log scale)", ylab = "", main = mt(mm))
  abline(v = 0, lty = 2, col = "grey50")
  ticks <- c(0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
  ticks <- ticks[ticks >= xs[1] * 0.9 & ticks <= xs[2] * 1.1]
  axis(1, at = log(ticks), labels = as.character(ticks))
  axis(2, at = n:1, labels = sapply(rows, function(z) z$lab), las = 1, cex.axis = 0.85)
  for (i in seq_len(n)) {
    if (!ok[i]) next
    y <- n - i + 1
    col <- c("grey35", "#B2182B", "#2166AC")[i]
    segments(log(lo[i]), y, log(hi[i]), y, lwd = 2, col = col)
    points(log(est[i]), y, pch = 18, cex = 1.6, col = col)
  }
  cat(sprintf("  M3 面板：画了 %d 条\n", n))
}


## ============================================================================
## 第 7 步：画 M4 面板
## ============================================================================
mm <- MODELS[4]
rows <- list(list(lab = "Unweighted (crude)", r = crude[model == mm]),
             list(lab = "IPTW, full cohort", r = main[model == mm]))
if (!is.null(ov)) rows[[length(rows) + 1]] <- list(lab = "IPTW, PS-overlap only", r = ov[model == mm])
rows <- rows[vapply(rows, function(z) nrow(z$r) > 0, logical(1))]
est <- sapply(rows, function(z) z$r$ratio)
lo  <- sapply(rows, function(z) z$r$rlo)
hi  <- sapply(rows, function(z) z$r$rhi)
ok  <- is.finite(est) & is.finite(lo) & is.finite(hi)
if (!any(ok)) {
  plot.new(); title(main = paste(mt(mm), "(no valid estimate)"))
  cat(sprintf("  %s：没有可画的点\n", mm))
} else {
  xs <- range(c(lo[ok], hi[ok], 1))
  n  <- length(rows)
  plot(NA, xlim = log(xs), ylim = c(0.4, n + 0.6), yaxt = "n", xaxt = "n",
       xlab = "OR / RR (log scale)", ylab = "", main = mt(mm))
  abline(v = 0, lty = 2, col = "grey50")
  ticks <- c(0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
  ticks <- ticks[ticks >= xs[1] * 0.9 & ticks <= xs[2] * 1.1]
  axis(1, at = log(ticks), labels = as.character(ticks))
  axis(2, at = n:1, labels = sapply(rows, function(z) z$lab), las = 1, cex.axis = 0.85)
  for (i in seq_len(n)) {
    if (!ok[i]) next
    y <- n - i + 1
    col <- c("grey35", "#B2182B", "#2166AC")[i]
    segments(log(lo[i]), y, log(hi[i]), y, lwd = 2, col = col)
    points(log(est[i]), y, pch = 18, cex = 1.6, col = col)
  }
  cat(sprintf("  M4 面板：画了 %d 条\n", n))
}

## ============================================================================
## 第 8 步：加上整张图的说明文字
## ============================================================================
mtext("Unit of analysis: one hypotensive episode for M1-M3 (M1 = episode-level incidence, >10 min), one surgical procedure for M4. GEE with an exchangeable working correlation structure, clustered by surgical procedure.",
      side = 1, line = 0.4, outer = TRUE, cex = 0.62, col = "grey30")
dev.off()
