## =============================================================================
## forest.R - forest plot of the main results (unweighted comparison / whole-cohort IPTW / overlap-region comparison, four outcomes)
## =============================================================================

## ============================================================================
## Step 1: preparation (load packages, set paths, confirm the output directory)
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # Read the csv and build the tables
})

## Root directory. All inputs and outputs are placed under this directory
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
CODE <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11")
OUT  <- file.path(CODE, "out")   # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(CODE)

cat("Output directory: ", OUT, "\n")


## ============================================================================
## Step 2: read the three result tables (unweighted / whole-cohort IPTW / overlap-region comparison)
## ============================================================================
crude <- fread(file.path(OUT, "Crude_Table_Ciprofol.csv"))
main  <- fread(file.path(OUT, "Analysis_Results_main.csv"))
ov_path <- file.path(OUT, "Analysis_Results_overlap.csv")
ov <- if (file.exists(ov_path)) fread(ov_path) else NULL

cat(sprintf("Read in: crude %d rows | full-cohort main results %d rows | overlap-region control %s\n",
            nrow(crude), nrow(main),
            if (is.null(ov)) "none (this panel is not drawn)" else paste0(nrow(ov), " rows")))


## ============================================================================
## Step 3: convert the "percent change" rows to the "ratio" scale and assemble a unified plotting dataset
## ============================================================================

## Convert Pct rows to ratios; for OR/RR rows the ratio is the estimate itself
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

## Draw the four models in the order M1->M2->M3->M4
MODELS <- unique(main$model)
cat(sprintf("The four models to be drawn: %s\n", paste(MODELS, collapse = "  |  ")))
if (length(MODELS) < 4L) stop("Only ", length(MODELS), " models read from the main results table; 4 expected (M1~M4)")

MT <- c(
  "M1 长期低血压事件比例（每次事件 >10 min）"         = "M1 Prolonged hypotensive episodes (>10 min)",
  "M2 单次事件低血压负荷（事件 AUC）"                 = "M2 Burden of a single episode (AUC)",
  "M3 每次低血压持续时间"                             = "M3 Duration per hypotensive episode",
  "M4 低血压事件次数（事件/例次，含麻醉时长 offset）" = "M4 Hypotensive episodes per procedure")
mt <- function(m) if (m %in% names(MT)) MT[[m]] else m


## ============================================================================
## Step 4: open the graphics device and draw the M1 panel (proportion of episodes with prolonged hypotension lasting >10 minutes, OR)
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
  cat(sprintf("  %s: no plottable point\n", mm))
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
  cat(sprintf("  M1 panel: %d points drawn (crude / full cohort%s)\n",
              n, if (is.null(ov)) "" else " / overlap region"))
}


## ============================================================================
## Step 5: draw the M2 panel
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
  cat(sprintf("  %s: no plottable point\n", mm))
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
  cat(sprintf("  M2 panel: %d points drawn\n", n))
}


## ============================================================================
## Step 6: draw the M3 panel
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
  cat(sprintf("  %s: no plottable point\n", mm))
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
  cat(sprintf("  M3 panel: %d points drawn\n", n))
}


## ============================================================================
## Step 7: draw the M4 panel
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
  cat(sprintf("  %s: no plottable point\n", mm))
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
  cat(sprintf("  M4 panel: %d points drawn\n", n))
}

## ============================================================================
## Step 8: add the caption text for the whole figure
## ============================================================================
mtext("Unit of analysis: one hypotensive episode for M1-M3 (M1 = episode-level incidence, >10 min), one surgical procedure for M4. GEE with an exchangeable working correlation structure, clustered by surgical procedure.",
      side = 1, line = 0.4, outer = TRUE, cex = 0.62, col = "grey30")
dev.off()
