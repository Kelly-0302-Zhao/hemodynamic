## =============================================================================
## 3-4 GEE analysis - temporal trend of hypotension incidence
## =============================================================================
## ============================================================================
## Step 1: preparation
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # read and process tables
  library(geepack)      # Run GEE
})
# Set up paths
D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # All result files are saved here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## Three ID columns. All must be converted to character before merging, otherwise inconsistent types misalign the join
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## Parse all timestamps as UTC
MY_TZ <- "UTC"


## ============================================================================
## Step 2: read the procedure-level analysis table and report the earliest date of ciprofol use
##   Expected output: the date when ciprofol first became available is printed to the screen
## ============================================================================
d <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
d <- d[!is.na(prolonged)]   # Keep only procedures with complete outcome data

cat("Earliest Ciprofol use date:", as.character(min(d[Group == "Ciprofol_only", as.Date(anes_start)])), "\n\n")


## ============================================================================
## Step 3: annual x group incidence table (denominator = procedures in that year), save to disk
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

cat("---- Year × group (denominator = procedures in that year)----\n")
print(trend)
fwrite(trend, file.path(OUT, "Table4_IncidenceTrend.csv"))

## Keep only years in which both drugs were used, to avoid spurious trends from "single-arm years"
both_years <- trend[, .N, by = year][N == 2, year]
cat("\nYears with both drugs present:", paste(both_years, collapse = ", "), "\n")
fwrite(trend[year %in% both_years], file.path(OUT, "Table4_IncidenceTrend_BothArms.csv"))


## ============================================================================
## Step 4: year x group interaction test (GEE, episode level, clustered by procedure)
## ============================================================================
## Attach the year, group, covariates and center of the parent procedure to the episode table
epi4 <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))
dd <- epi4[d[year %in% both_years, .(PATIENT_ID, VISIT_ID, OPER_ID, Group, year,
                                     Hypertension, AGE, SEX, ASA_GRADE, Hospital)],
           on = ID_VARS, nomatch = 0]
dd <- dd[!is.na(Group)]
dd[, long_epi := as.integer(duration > 10)]   # Whether this episode lasted longer than 10 minutes
dd <- dd[!is.na(long_epi)]

## geeglm requires an integer id; convert the procedure identifier (patient + visit + surgery) into an integer index;
## year_c = year minus the earliest year (makes the interaction coefficient easier to interpret)
dd[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
dd[, year_c  := as.numeric(year) - min(as.numeric(year))]

m <- geeglm(long_epi ~ Group * year_c + Hypertension + AGE + SEX + ASA_GRADE + Hospital,
            data = dd, id = sid_int, family = binomial(), corstr = "exchangeable")
s <- summary(m)$coefficients
print(round(s, 4))

ix <- grep(":", rownames(s), value = TRUE)
if (length(ix)) {
  cat(sprintf("\nInteraction term %s: OR=%.2f (%.2f-%.2f), p=%s\n", ix[1],
              exp(s[ix[1], "Estimate"]),
              exp(s[ix[1], "Estimate"] - 1.96 * s[ix[1], "Std.err"]),
              exp(s[ix[1], "Estimate"] + 1.96 * s[ix[1], "Std.err"]),
              ifelse(s[ix[1], ncol(s)] < 0.001, "<0.001",
                     sprintf("%.3f", s[ix[1], ncol(s)]))))
  cat("Interpretation: the interaction term is not significant, so there is no statistical evidence of a \"relative difference between the two groups over time\",\n")
  cat("      i.e. the observed between-group difference cannot be attributed to \"the more Ciprofol is used, the better\".\n")
}

## ============================================================================
## Step 5: bar chart of the timing distribution of first IOH
## ============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

bin_dur      <- 5    # One time bin every 5 minutes
max_obs_time <- 60   # Look only at the 60 minutes after induction

## first_low_t is Inf/NA in procedures with no IOH episode; drop the 1 procedure with an undefined outcome, consistent with the main analysis population
s5 <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
s5 <- s5[!is.na(n_epi)]
s5[, rel_min := as.numeric(difftime(first_low_t, anes_start, units = "mins"))]

group_totals <- s5[, .(total_n = .N), by = Group]
cat("\n---- Step 5: distribution of time to first IOH (denominator = all procedures in each group)----\n")
print(group_totals)

## First build the complete grid of "bin x group"
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

## Wide table: proportions and difference between the two groups in each bin, to see directly whether there is an early advantage
wide <- data.table::dcast(tr5, time_bin + time_label ~ Group, value.var = "incidence_pct")
wide[, diff_ciprofol_minus_propofol := Ciprofol_only - Propofol_only]
cat("\nProportion of first IOH occurrences per 5-minute bin (%) and between-group difference (Ciprofol − Propofol; negative = lower with Ciprofol):\n")
wide_show <- copy(wide)
num_cols  <- setdiff(names(wide_show), c("time_bin", "time_label"))
wide_show[, (num_cols) := lapply(.SD, function(x) round(x, 3)), .SDcols = num_cols]
print(wide_show)
fwrite(tr5,  file.path(OUT, "Table_TimeToFirst_IOH.csv"))
fwrite(wide, file.path(OUT, "Table_TimeToFirst_IOH_wide.csv"))

## Plot: side-by-side bar chart (colors follow the original script: ciprofol #E18476 / propofol #00B2B2)
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
