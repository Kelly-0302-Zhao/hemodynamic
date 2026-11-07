library(forestplot)
library(grid)

# ============================================================
# M2 (总低血压负荷 AUC) 分层亚组森林图 —— 自动读 CSV 取数版
# ============================================================

OUT_DIR <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827/投稿材料/JCA/投稿包_20261009/10-11/out"

# ---------- 自动取数 ----------
get_m2 <- function(tag) {
  f <- file.path(OUT_DIR, paste0("Sensitivity_Analysis_Results_", tag, ".csv"))
  if (!file.exists(f)) stop("missing: ", f)
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  row <- d[grepl("^M2 ", d$model), , drop = FALSE]
  if (nrow(row) != 1) stop("expect 1 M2 row in ", f, " got ", nrow(row))
  list(
    est = as.numeric(row$estimate),
    lo  = as.numeric(row$lo),
    hi  = as.numeric(row$hi),
    p   = as.numeric(row$p),
    n   = as.numeric(row$n_unit)
  )
}

# 顺序与原脚本一致：65-75 / 75-85 / >=85 / Female / Male
LAYERS <- list(
  list(label = "  65-75", tag = "AGE 65-75"),
  list(label = "  75-85", tag = "AGE 75-85"),
  list(label = "  \u226585", tag = "AGE 85+"),
  list(label = "  Female", tag = "female"),
  list(label = "  Male",   tag = "male")
)

res <- lapply(LAYERS, function(x) c(get_m2(x$tag), list(label = x$label)))

# 百分数 >=1 取 1 位；负号用 ASCII '-'（forestplot 内为纯文本）
fmt_pct <- function(v) sprintf("%.1f%%", v)
fmt_eff <- function(r) {
  sprintf("%s (%s to %s)", fmt_pct(r$est), fmt_pct(r$lo), fmt_pct(r$hi))
}
fmt_p <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)

cat("=== M2 subgroup values read from CSV ===\n")
for (r in res) {
  cat(sprintf("  %-10s %-34s P=%-8s n_epi=%d\n",
              trimws(r$label), fmt_eff(r), fmt_p(r$p), as.integer(r$n)))
}

# ---------- 组装向量（10 行：表头2 + Age标题1 + Age 3层 + Sex标题1 + Sex 2层 + 空1）----------
labels <- c("Subgroup", "",
            "Age (years)",
            res[[1]]$label, res[[2]]$label, res[[3]]$label,
            "",
            "Sex",
            res[[4]]$label, res[[5]]$label)

is_summary <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)

effect_text <- c("Effect (95% CI)", "", "",
                 fmt_eff(res[[1]]), fmt_eff(res[[2]]), fmt_eff(res[[3]]),
                 "",
                 "",
                 fmt_eff(res[[4]]), fmt_eff(res[[5]]))

p_text <- c("P value", "", "",
            fmt_p(res[[1]]$p), fmt_p(res[[2]]$p), fmt_p(res[[3]]$p),
            "",
            "",
            fmt_p(res[[4]]$p), fmt_p(res[[5]]$p))

text <- cbind(labels, effect_text, p_text)

m <- c(NA, NA, NA,
       res[[1]]$est, res[[2]]$est, res[[3]]$est,
       NA, NA,
       res[[4]]$est, res[[5]]$est)

l <- c(NA, NA, NA,
       res[[1]]$lo, res[[2]]$lo, res[[3]]$lo,
       NA, NA,
       res[[4]]$lo, res[[5]]$lo)

u <- c(NA, NA, NA,
       res[[1]]$hi, res[[2]]$hi, res[[3]]$hi,
       NA, NA,
       res[[4]]$hi, res[[5]]$hi)

# ---------- 绘图----------
fp <- function() {
  forestplot(
    labeltext = text,
    mean = m,
    lower = l,
    upper = u,
    is.summary = is_summary,

    title = "Total Hypotension Burden (AUC)",

    graph.pos = 2,
    zero = 0,

    xticks = c(-40, -20, 0, 20),
    xlab = "Percentage change (%)",

    graphwidth = unit(125, "mm"),

    colgap = unit(7, "mm"),
    line.margin = unit(2.0, "mm"),

    boxsize = 0.42,
    fn.ci_norm = fpDrawCircleCI,
    ci.vertices = TRUE,
    ci.vertices.height = 0.14,
    lwd.ci = 2.6,

    col = fpColors(
      box = "#1F4E79",
      line = "#1F4E79",
      summary = "#1F4E79",
      zero = "#555555"
    ),

    txt_gp = fpTxtGp(
      label = gpar(fontsize = 17, fontfamily = "serif"),
      ticks = gpar(fontsize = 17, fontfamily = "serif"),
      xlab = gpar(fontsize = 19, fontface = "bold", fontfamily = "serif"),
      title = gpar(fontsize = 22, fontface = "bold", fontfamily = "serif")
    ),

    hrzl_lines = list(
      "1" = gpar(lwd = 1.5, col = "black"),
      "2" = gpar(lwd = 0.9, col = "black"),
      "7" = gpar(lwd = 2.8, col = "black")
    ),

    align = c("l", "l", "l")
  )
}

DEST <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827/投稿材料/JCA/投稿包_20261009/乱七八糟/Figures"

tiff(
  filename = file.path(DEST, "Supplemental Figure S1_Subgroup_burden.tiff"),
  width = 5600, height = 2500, units = "px", res = 400, compression = "lzw"
)
grid.newpage(); fp(); dev.off()

pdf(
  file = file.path(DEST, "Supplemental Figure S1_Subgroup_burden.pdf"),
  width = 18, height = 8
)
grid.newpage(); fp(); dev.off()

png(
  filename = file.path(DEST, "Supplemental Figure S1_Subgroup_burden.png"),
  width = 5600, height = 2500, units = "px", res = 400
)
grid.newpage(); fp(); dev.off()
