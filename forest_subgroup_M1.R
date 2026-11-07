suppressMessages({
  library(forestplot)
  library(grid)
  library(showtext)
})

font_add("Arial",
         regular    = "C:/Windows/Fonts/arial.ttf",
         bold       = "C:/Windows/Fonts/arialbd.ttf",
         italic     = "C:/Windows/Fonts/ariali.ttf",
         bolditalic = "C:/Windows/Fonts/arialbi.ttf")
showtext_opts(dpisfamily = "Arial")


OUT_DIR  <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827/投稿材料/JCA/投稿包_20261009/10-11/out"
DEST_DIR <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827/投稿材料/JCA/投稿包_20261009/乱七八糟/Figures"

get_est <- function(tag, model_id) {
  f <- file.path(OUT_DIR, paste0("Sensitivity_Analysis_Results_", tag, ".csv"))
  if (!file.exists(f)) stop("missing: ", f)
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  row <- d[grepl(paste0("^", model_id, " "), d$model), , drop = FALSE]
  if (nrow(row) != 1) stop("expect 1 ", model_id, " row in ", f)
  list(est = as.numeric(row$estimate),
       lo  = as.numeric(row$lo),
       hi  = as.numeric(row$hi),
       p   = as.numeric(row$p),
       n   = as.numeric(row$n_unit))
}

LAYERS <- list(
  list(label = "  65-75",  tag = "AGE 65-75"),
  list(label = "  75-85",  tag = "AGE 75-85"),
  list(label = "  \u2265 85", tag = "AGE 85+"),
  list(label = "  Female", tag = "female"),
  list(label = "  Male",   tag = "male")
)

fmt_p <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)

res <- lapply(LAYERS, function(x) c(get_est(x$tag, "M1"), list(label = x$label)))

cat("=== M1 subgroup values ===\n")
for (r in res) {
  cat(sprintf("  %-9s %8.2f (%8.2f to %8.2f)  P=%-8s n=%d\n",
              trimws(r$label), r$est, r$lo, r$hi, fmt_p(r$p), as.integer(r$n)))
}

# 统一用 "to" 分隔 CI 区间：负值区间若用连字符会写成 "−33.1−−24.9" 造成歧义，
# 全部指标（OR/RR/%）一律用 to，正文表格亦同，全篇一致。
fmt_eff <- function(r) sprintf("%.2f (%.2f to %.2f)", r$est, r$lo, r$hi)

labels <- c("Subgroup", "", "Age (years)",
            res[[1]]$label, res[[2]]$label, res[[3]]$label,
            "", "Sex",
            res[[4]]$label, res[[5]]$label)

is_summary <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)

effect_text <- c("Effect (95% CI)", "", "",
                 fmt_eff(res[[1]]), fmt_eff(res[[2]]), fmt_eff(res[[3]]),
                 "", "",
                 fmt_eff(res[[4]]), fmt_eff(res[[5]]))

p_text <- c("P-value", "", "",
            fmt_p(res[[1]]$p), fmt_p(res[[2]]$p), fmt_p(res[[3]]$p),
            "", "",
            fmt_p(res[[4]]$p), fmt_p(res[[5]]$p))

m <- c(NA, NA, NA, res[[1]]$est, res[[2]]$est, res[[3]]$est, NA, NA,
       res[[4]]$est, res[[5]]$est)
l <- c(NA, NA, NA, res[[1]]$lo,  res[[2]]$lo,  res[[3]]$lo,  NA, NA,
       res[[4]]$lo,  res[[5]]$lo)
u <- c(NA, NA, NA, res[[1]]$hi,  res[[2]]$hi,  res[[3]]$hi,  NA, NA,
       res[[4]]$hi,  res[[5]]$hi)

draw <- function() {
  forestplot(
    labeltext = cbind(labels, effect_text, p_text),
    mean = m, lower = l, upper = u,
    is.summary = is_summary,
    title = "Prolonged Hypotension (>10 min): Odds Ratio",
    graph.pos = 2,
    zero = 1,
    xticks = c(0.25, 0.50, 0.75, 1.00),
    xlab = "Ratio",
    graphwidth = unit(125, "mm"),
    colgap = unit(7, "mm"),
    line.margin = unit(2.0, "mm"),
    boxsize = 0.42,
    fn.ci_norm = fpDrawCircleCI,
    ci.vertices = TRUE,
    ci.vertices.height = 0.14,
    lwd.ci = 2.6,
    col = fpColors(box = "#1F4E79", line = "#1F4E79",
                   summary = "#1F4E79", zero = "#555555"),
    txt_gp = fpTxtGp(
      label = gpar(fontsize = 17, fontfamily = "Arial"),
      ticks = gpar(fontsize = 17, fontfamily = "Arial"),
      xlab  = gpar(fontsize = 19, fontface = "bold", fontfamily = "Arial"),
      title = gpar(fontsize = 22, fontface = "bold", fontfamily = "Arial")
    ),
    hrzl_lines = list(
      "1" = gpar(lwd = 1.5, col = "black"),
      "2" = gpar(lwd = 0.9, col = "black"),
      "7" = gpar(lwd = 2.8, col = "black")
    ),
    align = c("l", "l", "l")
  )
}

BASE <- "Supplemental Figure S1_subgroup_M1"

tiff(file.path(DEST_DIR, paste0(BASE, ".tiff")),
     width = 5600, height = 2500, units = "px", res = 400, compression = "lzw")
grid.newpage(); draw(); dev.off()

png(file.path(DEST_DIR, paste0(BASE, ".png")),
    width = 5600, height = 2500, units = "px", res = 400)
grid.newpage(); draw(); dev.off()
