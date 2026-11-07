## =============================================================================
## 3-4 GEE analysis -— 低血压发生率的时间趋势
## =============================================================================
## ============================================================================
## 第 1 步：准备工作
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # 读表、处理表
  library(geepack)      # 跑 GEE
})
# 设定路径
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## 三个 ID 列。合并前都要先转成文字型，类型不一致会合并错位
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 时间戳统一按 UTC 解析
MY_TZ <- "UTC"


## ============================================================================
## 第 2 步：读入手术级分析表，报告环泊酚最早使用日期
##   预期结果：屏幕上打印环泊酚最早上线的日期
## ============================================================================
d <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
d <- d[!is.na(prolonged)]   # 只保留结局齐全的手术

cat("环泊酚最早使用日期：", as.character(min(d[Group == "Ciprofol_only", as.Date(anes_start)])), "\n\n")


## ============================================================================
## 第 3 步：逐年 × 分组发生率表（分母 = 当年手术例次），存盘
## ============================================================================
trend <- d[, .(
  n_procedures   = .N,
  n_patients     = uniqueN(PATIENT_ID),
  Any_IOH_pct= round(100 * mean(any_hypo), 1),
  Prolonged_IOH_pct= round(100 * mean(prolonged), 1),
  Cum_duration_median= round(median(MAP65_TIME), 1),
  Cum_AUC_median= round(median(MAP65_AUC), 1),
  Episodes_per_procedure= round(mean(n_epi), 2)
), by = .(year, Group)][order(year, Group)]

cat("---- 逐年 × 分组（分母 = 当年手术例次）----\n")
print(trend)
fwrite(trend, file.path(OUT, "Table4_IncidenceTrend.csv"))

## 只保留两药同时存在的年份，避免"单臂年份"造成假趋势
both_years <- trend[, .N, by = year][N == 2, year]
cat("\n两药同时存在的年份：", paste(both_years, collapse = ", "), "\n")
fwrite(trend[year %in% both_years], file.path(OUT, "Table4_IncidenceTrend_BothArms.csv"))


## ============================================================================
## 第 4 步：年份 × 分组交互检验（GEE，事件级，按手术聚类）
## ============================================================================
## 把事件表挂上所属手术的年份、分组、协变量、中心
epi4 <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))
dd <- epi4[d[year %in% both_years, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, year,
                                     Hypertension, AGE, SEX, ASA_GRADE, Hospital)],
           on = ID_VARS, nomatch = 0]
dd <- dd[!is.na(Group)]
dd[, long_epi := as.integer(duration > 10)]   # 这次事件是否持续超过 10 分钟
dd <- dd[!is.na(long_epi)]

## geeglm 的 id 必须整数，把手术编号（患者+就诊+手术）转成整数编号；
## year_c = 年份减去最早年份（让交互项系数更好解释）
dd[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
dd[, year_c  := as.numeric(year) - min(as.numeric(year))]

m <- geeglm(long_epi ~ Group * year_c + Hypertension + AGE + SEX + ASA_GRADE + Hospital,
            data = dd, id = sid_int, family = binomial(), corstr = "exchangeable")
s <- summary(m)$coefficients
print(round(s, 4))

ix <- grep(":", rownames(s), value = TRUE)
if (length(ix)) {
  cat(sprintf("\n交互项 %s：OR=%.2f (%.2f-%.2f)，p=%s\n", ix[1],
              exp(s[ix[1], "Estimate"]),
              exp(s[ix[1], "Estimate"] - 1.96 * s[ix[1], "Std.err"]),
              exp(s[ix[1], "Estimate"] + 1.96 * s[ix[1], "Std.err"]),
              ifelse(s[ix[1], ncol(s)] < 0.001, "<0.001",
                     sprintf("%.3f", s[ix[1], ncol(s)]))))
  cat("解释：交互项不显著，说明「两组随时间的相对差异」没有统计学证据，\n")
  cat("      即观察到的组间差异不能归因于「环泊酚越用越好」。\n")
}

## ============================================================================
## 第 5 步：首次 IOH 的"时间分布条形图"
## ============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

bin_dur      <- 5    # 每 5 分钟一个时间窗
max_obs_time <- 60   # 只看诱导后 60 分钟

## first_low_t 在"从未发生 IOH"的手术里是 Inf/NA；剔除结局不可定义的那 1 台，与主分析人群一致
s5 <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
s5 <- s5[!is.na(n_epi)]
s5[, rel_min := as.numeric(difftime(first_low_t, anes_start, units = "mins"))]

group_totals <- s5[, .(total_n = .N), by = Group]
cat("\n---- 第 5 步：首次 IOH 时间分布（分母 = 各组全部手术例次）----\n")
print(group_totals)

## 先把"箱 × 组"的完整网格列出来
bins <- seq(0, max_obs_time - bin_dur, by = bin_dur)
grid <- data.table::CJ(time_bin = bins, Group = c("Ciprofol_only", "Propofol_only"))
ev <- s5[!is.na(rel_min) & rel_min >= 0 & rel_min < max_obs_time,
         .(events = .N), by = .(Group, time_bin = floor(rel_min / bin_dur) * bin_dur)]

tr5 <- merge(grid, ev, by = c("time_bin", "Group"), all.x = TRUE)
tr5[is.na(events), events := 0L]
tr5 <- merge(tr5, group_totals, by = "Group")
tr5[, incidence_pct := 100 * events / total_n]
tr5[, time_label := paste0(time_bin, "-", time_bin + bin_dur)]
data.table::setorder(tr5, time_bin, Group)

## 宽表：每箱两组的比例与差值，直接看"有没有早期优势"
wide <- data.table::dcast(tr5, time_bin + time_label ~ Group, value.var = "incidence_pct")
wide[, diff_ciprofol_minus_propofol := Ciprofol_only - Propofol_only]
cat("\n每 5 分钟箱的首次 IOH 发生比例（%）与组间差（环泊酚 − 丙泊酚，负值 = 环泊酚更低）：\n")
wide_show <- copy(wide)
num_cols  <- setdiff(names(wide_show), c("time_bin", "time_label"))
wide_show[, (num_cols) := lapply(.SD, function(x) round(x, 3)), .SDcols = num_cols]
print(wide_show)
fwrite(tr5,  file.path(OUT, "Table_TimeToFirst_IOH.csv"))
fwrite(wide, file.path(OUT, "Table_TimeToFirst_IOH_wide.csv"))

## 作图：并列条形图（配色沿用原脚本：环泊酚 #E18476 / 丙泊酚 #00B2B2）
tr5[, Group := factor(Group, levels = c("Ciprofol_only", "Propofol_only"))]
p_bar <- ggplot(tr5, aes(x = factor(time_bin), y = incidence_pct, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = paste0(bins, "-", bins + bin_dur)) +
  scale_fill_manual(values = c("Ciprofol_only" = "#E18476", "Propofol_only" = "#00B2B2"),
                    labels = c("Ciprofol", "Propofol")) +
  labs(title = "Distribution of Time to First IOH Episode",
       subtitle = "Denominator: total population in each group",
       x = "Minutes since induction",
       y = "Incidence of first IOH event (%)") +
  theme_bw() +
  theme(text = element_text(family = "Times New Roman"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT, "Figure_TimeToFirst_IOH.pdf"), p_bar,
       width = 9, height = 5.5, device = cairo_pdf)
ggsave(file.path(OUT, "Figure_TimeToFirst_IOH.png"), p_bar,
       width = 9, height = 5.5, dpi = 600, device = "png", type = "cairo")
tiff(file.path(OUT, "Figure_TimeToFirst_IOH.tiff"),
     width = 9 * 600, height = 5.5 * 600, res = 600, compression = "lzw")
print(p_bar); dev.off()
