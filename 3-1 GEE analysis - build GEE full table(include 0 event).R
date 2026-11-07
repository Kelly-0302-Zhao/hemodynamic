## =============================================================================
## 3-1 GEE analysis - 构建手术级分析表
## =============================================================================
## ============================================================================
## 第 1 步：准备工作（加载包、设定路径、写死常量）
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # 读表、处理表
  library(geepack)      # 后面跑 GEE 模型用
})

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
DATA <- file.path(D, "JCA_run")      # 清洗好的中间表
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## 三个 ID 列。三张表的 ID 列类型可能不一致，用之前都要先转成文字型
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## 时间戳统一按 UTC 解析（原始数据带 Z 后缀）
MY_TZ <- "UTC"

## ============================================================================
## 第 2 步：读入手术基础表，解析麻醉起止时间，整理成"一行 = 一次手术"
## ============================================================================
parse_iso <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  as.POSIXct(gsub("T", " ", sub("Z$", "", as.character(x))),
             tz = MY_TZ, format = "%Y-%m-%d %H:%M:%S")
}
to_chr_ids <- function(DT) {
  DT[, (ID_VARS) := lapply(.SD, as.character), .SDcols = ID_VARS]
  DT
}

master <- fread(file.path(DATA, "master_clean_85397obs.csv"), encoding = "UTF-8")
master <- to_chr_ids(master)
master[, anes_start := parse_iso(ANES_START_TIME)]
master[, anes_end   := parse_iso(ANES_END_TIME)]
master[, year       := as.integer(format(as.Date(anes_start), "%Y"))]
## 分组因子：以丙泊酚为参照组（环泊酚 vs 丙泊酚）
master[, Group := relevel(factor(as.character(Group)), ref = "Propofol_only")]
## 以下几个分类变量转成因子，方便后续模型处理
for (v in c("ASA_GRADE", "operation_type", "Hospital", "SEX")) {
  if (v %in% names(master)) master[[v]] <- factor(master[[v]])
}

## ============================================================================
## 第 3 步：读入每台手术的低血压结局表，按 ID 拼到基础表上
## ============================================================================
surg <- setDT(readRDS(file.path(OUT, "Final_Analysis_Surgery.rds")))

## 基础表挑出需要的列
keep <- c(ID_VARS, "Group", "Hospital", "year", "anes_start", "anes_end",
          "anes_duration_mins", "AGE", "SEX", "BMI", "ASA_GRADE",
          "EMERGENCY_INDICATOR", "operation_type",
          "Hypertension", "Diabetes", "Coronary_artery_disease",
          "Heart_failure", "IN_ICU_INDICATOR")

d <- master[, ..keep][surg[, !c("anes_start", "anes_end", "anes_duration_mins"),
                           with = FALSE], on = ID_VARS]



## ============================================================================
## 第 4 步：给总表补上分析要用的派生列
## ============================================================================
d[, lg_dur := log(anes_duration_mins)]
d[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]
d[, year_f  := factor(year)]
d[, any_hypo := as.integer(n_epi >= 1L)]

cat("已补上派生列：lg_dur / sid_int / year_f / any_hypo\n")


## ============================================================================
## 第 5 步：单位核对（打印"例次 vs 人数"，并看分组 × 中心 / 年份）
## ============================================================================
cat("=============== 单位核对（很重要）===============\n")
cat("分析表的行数（= 手术例次）：", nrow(d), "\n")
cat("其中患者人数              ：", uniqueN(d$PATIENT_ID), "\n")
cat("同一患者做过多次手术的人  ：",
    d[, .N, by = PATIENT_ID][N > 1, .N], "人（共",
    d[, .N, by = PATIENT_ID][N > 1, sum(N)], "例次）\n")
cat("可分析结局的手术例次      ：", d[!is.na(prolonged), .N], "\n")
cat("可分析结局的患者人数      ：", uniqueN(d[!is.na(prolonged)]$PATIENT_ID), "\n\n")

cat("---- 分组 × 中心（例次）----\n")
print(dcast(d[!is.na(prolonged)], Hospital ~ Group, value.var = "prolonged", fun.aggregate = length))
cat("\n---- 分组 × 年份（例次）----\n")
print(dcast(d[!is.na(prolonged)], year ~ Group, value.var = "prolonged", fun.aggregate = length))


## ============================================================================
## 第 6 步：保存手术级分析总表
## ============================================================================
saveRDS(d, file.path(OUT, "Final_Analysis_Long_Table.rds"))
fwrite(d[, .(PATIENT_ID, VISIT_ID, OPER_ID, Hospital, year, Group, n_epi, n_long,
             prolonged, prolonged_cum, multi_epi, short_only, any_hypo,
             MAP65_TIME, MAP65_AUC, n_map, rate_map, anes_duration_mins,
             AGE, SEX, BMI, ASA_GRADE, Hypertension, Diabetes,
             Coronary_artery_disease, Heart_failure)],
       file.path(OUT, "Final_Analysis_Long_Table.csv"))
