## =============================================================================
## 1-7 data clean - group of susceptible.R —— 造"易感人群"亚组（给药前 MAP<65）
## =============================================================================
#构建给药前即发生低血压的人群，为易感人群，用于后续亚组分析
monitor_susceptible <- monitor1 %>%
  semi_join(bp.exclude1, by = c("PATIENT_ID", "VISIT_ID", "OPER_ID"))

## ============================================================================
## 第 1 步：把易感监测点关联上"手术时间轴 + 药物组别"，并打上 cluster_id
## ============================================================================
monitor_susceptible_final <- monitor_susceptible %>%
  left_join(data1 %>% select(PATIENT_ID, VISIT_ID, OPER_ID, 
                             Group, Hospital, 
                             START_DATE_TIME, END_DATE_TIME), 
            by = c("PATIENT_ID", "VISIT_ID", "OPER_ID"))
# 统计一下这部分人的数据量 ---
message("易感人群总监测点数量: ", nrow(monitor_susceptible_final)) # 25338320 obs
monitor_susceptible_final <- monitor_susceptible_final %>%
  mutate(cluster_id   = paste(PATIENT_ID, VISIT_ID, OPER_ID, sep = "_") )
# 保存为独立亚组分析文件 ---
# 输出 CSV
fwrite(monitor_susceptible_final, "monitor_below65_before_start.csv")
#获取易感人群基线数据表
## ============================================================================
## 第 2 步：取易感人群的唯一手术清单，按每台手术首条记录造基线表，并写出两个 CSV
## ============================================================================
susceptible_patients_list <- monitor_susceptible_final %>%
  distinct(cluster_id)
message("易感人群中的唯一手术数: ", nrow(susceptible_patients_list)) # 27594 obs
# 2. 从 data1 提取这些人的基线数据并去重，对于同一个病人多次手术的情况，我们以手术为单位
baseline_data_susceptible <- data1 %>%
  semi_join(susceptible_patients_list, by = "cluster_id") %>%
  arrange(cluster_id, START_DATE_TIME) %>%
  group_by(cluster_id) %>%
  slice(1) %>%
  ungroup()
# 3. 统计该亚组的人数分布
message("易感人群基线表 (N=", nrow(baseline_data_susceptible), ") 构建完成。") # n=27449
print(table(baseline_data_susceptible$Group)) # Ciprofol_only 1051 : Propofol_only 26398 
# 保存，方便后续亚组分析
fwrite(baseline_data_susceptible, "monitor_below65_before_start_baseline.csv")

## ============================================================================
## 第 3 步：选择偏倚检查 —— 术前就低血压的"易感人群"在两组各占多少百分比
## ============================================================================
## 3.1 给每个手术打上"是否易感"标记和年份
sb <- data1 %>%
  mutate(
    susceptible = as.integer(cluster_id %in% baseline_data_susceptible$cluster_id),
    year        = as.integer(substr(as.character(START_DATE_TIME), 1, 4))
  )

## 3.2 两函数：一个算每组的 n 与百分比，一个算组间 RR（环泊酚 vs 丙泊酚）
sb_count <- function(dt, level, hosp = NA_character_, yr = NA_integer_) {
  dt %>%
    filter(Group %in% c("Ciprofol_only", "Propofol_only")) %>%
    group_by(Group) %>%
    summarise(n_surgery       = n(),
              n_susceptible   = sum(susceptible),
              pct_susceptible = round(100 * mean(susceptible), 2),
              .groups = "drop") %>%
    mutate(level = level, Hospital = hosp, year = yr) %>%
    select(level, Hospital, year, Group, n_surgery, n_susceptible, pct_susceptible)
}
sb_rr <- function(dt, level, hosp = NA_character_, yr = NA_integer_) {
  dt <- dt %>% filter(Group %in% c("Ciprofol_only", "Propofol_only"))
  n1 <- sum(dt$Group == "Ciprofol_only"); a <- sum(dt$Group == "Ciprofol_only" & dt$susceptible == 1L)
  n0 <- sum(dt$Group == "Propofol_only");  b <- sum(dt$Group == "Propofol_only"  & dt$susceptible == 1L)
  if (n1 == 0 || n0 == 0) return(NULL)                     # 该层只有一种药，算不出 RR
  p1 <- a / n1; p0 <- b / n0
  if (a == 0 || b == 0) {                                  # 分子为 0 时用 +0.5 校正，不给 CI
    rr <- ((a + 0.5) / (n1 + 1)) / ((b + 0.5) / (n0 + 1)); se <- NA_real_
  } else {
    rr <- p1 / p0; se <- sqrt(1/a - 1/n1 + 1/b - 1/n0)
  }
  data.frame(level = level, Hospital = hosp, year = yr,
             n_ciprofol   = n1, pct_ciprofol = round(100 * p1, 2),
             n_propofol   = n0, pct_propofol = round(100 * p0, 2),
             RR = round(rr, 3),
             CI_low  = round(exp(log(rr) - 1.96 * se), 3),
             CI_high = round(exp(log(rr) + 1.96 * se), 3))
}

## 3.3 组装三张表：全部三院 / 按医院分层 / 中大按年份
sb_year_list <- sort(unique(sb$year[sb$Hospital == "ZhongDa" & !is.na(sb$year)]))
sb_share <- bind_rows(
  sb_count(sb, "全部三院"),
  sb_count(sb %>% filter(Hospital == "ZhongDa"), "仅 ZhongDa（两药都有）", "ZhongDa"),
  sb_count(sb %>% filter(Hospital == "NanYi"),   "仅 NanYi（环泊酚仅 2024 年）", "NanYi"),
  sb_count(sb %>% filter(Hospital == "AnYi"),    "仅 AnYi（零环泊酚）", "AnYi"),
  bind_rows(lapply(sb_year_list, function(y)
    sb_count(sb %>% filter(Hospital == "ZhongDa", year == y), "ZhongDa 按年份", "ZhongDa", y)))
)
sb_rr_tab <- bind_rows(
  sb_rr(sb, "全部三院"),
  sb_rr(sb %>% filter(Hospital == "ZhongDa"), "仅 ZhongDa（两药都有）", "ZhongDa"),
  sb_rr(sb %>% filter(Hospital == "ZhongDa", year %in% 2022:2024), "ZhongDa 2022-2024（同期）", "ZhongDa"),
  bind_rows(lapply(sb_year_list, function(y)
    sb_rr(sb %>% filter(Hospital == "ZhongDa", year == y), "ZhongDa 按年份", "ZhongDa", y)))
)

## 3.4 打印
cat("\n【选择偏倚检查】易感人群（术前 MAP<65）在两组各占该组手术的百分比\n")
print(as.data.frame(sb_share), row.names = FALSE)
cat("\n【组间风险比】RR = 环泊酚组易感占比 / 丙泊酚组易感占比（RR>1 表示环泊酚组里重的病人更多）\n")
print(as.data.frame(sb_rr_tab), row.names = FALSE) #全部三院或是仅zhongda中显著，且RR均小于1，提示丙泊酚里反而重的病人多
fwrite(sb_share,  "Susceptible_Share_By_Group.csv")
fwrite(sb_rr_tab, "Susceptible_Share_RR.csv")
