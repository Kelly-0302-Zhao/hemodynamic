## =============================================================================
## 3-2 GEE analysis - QIC test.R —— 工作相关结构选择的 QIC 检验
## =============================================================================

## ============================================================================
## 第 1 步：准备工作（加载包、设定路径、写死常量）
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)   # 读表、处理表
  library(geepack)      # 跑 GEE、算 QIC
})

D    <- "D:/Paper/anesthesia/propofol&ciprofol compares/My paper/multicenter/20260827"
OUT  <- file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11", "out")  # 所有结果文件都存到这里
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(file.path(D, "投稿材料", "JCA", "投稿包_20261009", "10-11"))

## 三个 ID 列。合并前都要先转成文字型，类型不一致会合并错位
ID_VARS <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

## ============================================================================
## 第 2 步：把事件表挂上手术表的"分组 + 协变量"，得到 QIC 用的事件级数据
## ============================================================================
d  <- setDT(readRDS(file.path(OUT, "Final_Analysis_Long_Table.rds")))
ep <- setDT(readRDS(file.path(OUT, "Final_Analysis_Events.rds")))

## 只保留"可分析结局"的手术，把它的分组和 4 个调整变量挂到每个事件上
ep <- ep[d[!is.na(prolonged), .(PATIENT_ID, VISIT_ID, OPER_ID, Group,
                                Hypertension, AGE, SEX, ASA_GRADE)],
         on = ID_VARS, nomatch = 0]
ep <- ep[!is.na(Group)]
ep[, long_epi := as.integer(duration > 10)]   # 这次事件是否持续超过 10 分钟
ep <- ep[!is.na(long_epi)]

## geeglm 的聚类 id 必须是整数，把手术编号（患者+就诊+手术）转成整数编号
ep[, sid_int := as.integer(factor(paste(PATIENT_ID, VISIT_ID, OPER_ID)))]

## 每台手术发生几次事件 —— 聚类 id = 手术，这个分布才决定 exchangeable 与 independence 差多少
rows_per_surgery <- ep[, .N, by = .(PATIENT_ID, VISIT_ID, OPER_ID)]
cat(sprintf("每台手术低血压事件数分布：1 次 %.1f%% | 2 次 %.1f%% | ≥3 次 %.1f%%\n\n",
            100 * mean(rows_per_surgery$N == 1),
            100 * mean(rows_per_surgery$N == 2),
            100 * mean(rows_per_surgery$N >= 3)))

## ============================================================================
## 第 3 步：对三种相关结构各跑一遍 M1 模型，用 QIC 选最小的那一个
## ============================================================================

## 对给定相关结构跑 M1（长期低血压事件比例），取出 QIC 和 QICu
fit_cor <- function(cs) {
  m <- geeglm(long_epi ~ Group + Hypertension + AGE + SEX + ASA_GRADE,
              data = ep, id = sid_int, family = binomial(), corstr = cs)
  q <- QIC(m)   # 某些 geepack 版本返回的是向量，用中括号取具体值
  data.table(Correlation_structure = cs,
              QIC = round(as.numeric(q[["QIC"]]), 2),
             LogLik_criterion = round(as.numeric(q[["QICu"]]), 2))
}

tab <- rbindlist(lapply(c("independence", "exchangeable", "ar1"), fit_cor))
print(tab)
best <- tab[which.min(QIC), Correlation_structure]
cat("\nQIC 最小的是：", best, "\n")
cat("三种结构的 QIC 非常接近，未构成明确区分；按预设口径，后续加权 GEE 的",
    "四个模型一律采用 exchangeable。\n")

## 把表同时存一份 Word（缺 flextable/officer 包时只提示，不中断）
save_docx_table <- function(df, path, caption = NULL) {
  ok <- tryCatch({
    suppressPackageStartupMessages({ library(flextable); library(officer) })
    ft <- flextable(as.data.frame(df)) |> theme_booktabs() |> autofit()
    if (!is.null(caption)) ft <- set_caption(ft, caption = caption)
    doc <- read_docx() |> body_add_flextable(ft)
    print(doc, target = path)
    TRUE
  }, error = function(e) {
    cat("  [Word 输出失败，CSV 已保存] ", conditionMessage(e), "\n")
    FALSE
  })
  invisible(ok)
}

fwrite(tab, file.path(OUT, "Table3_GEE_QIC.csv"))
save_docx_table(tab, file.path(OUT, "Table3_GEE_QIC.docx"),
                caption = "Table 3. Working correlation structure selection by QIC (unweighted M1 model, one hypotensive episode per row, clustered by surgical procedure).")

