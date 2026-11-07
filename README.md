# hemodynamic
R code about hemodynamics

# 代码使用说明

本目录是论文的数据分析代码，共 **34 个 R 脚本**，按编号顺序执行即可从清洗后的数据跑出全部表格、图形和敏感性/亚组/中介分析。

---

[简体中文](./README_zh.md) | [English](./README.md)

**脚本命名规则**：`编号 - 类型 - 具体内容.R`

| 编号          | 段         | 含义                                       |
| ----------- | --------- | ---------------------------------------- |
| `1-x`       | 数据清洗      | 原始表                             |
| `2-x`       | 描述统计      | 事件表重建、Table 2                            |
| `3-x`       | 主分析       | IPTW + 四个 GEE 主模型（论文主结果）                 |
| `4-x`       | 敏感性/亚组/中介 | 剂量中介、时长分层、监测方式分层、失血量中介，以及年龄/时间/中心/性别七个亚组 |
| `forest*.R` | 出图        | 森林图                                      |

`4-x` 段共 11 个脚本，按分析域排列：

| 编号              | 分析域    | 脚本                                                                            | 中介/分层口径                                          |
| --------------- | ------ | ----------------------------------------------------------------------------- | ------------------------------------------------ |
| `4-1`           | 剂量中介   | `4-1 mediation analysis - dose.R`                                             | 四个口径：① 粗口径 ② ①+中心×年FE ③ 两臂共存格子 ④ ③+中心×年FE        |
| `4-2`           | 时长分层   | `4-2 subgroup analysis - operation duration.R`                                | 三分位 + 临床切点（<120 / 120–240 / ≥240 min）+ 交互检验      |
| `4-3`           | 监测方式分层 | `4-3 subgroup analysis - invasive vs non-invasive BP monitoring.R`            | 按**结局血压实际来源**归层（有创 37,913 台 / 无创 47,484 台）+ 交互检验 |
| `4-4`           | 失血量中介  | `4-4 mediation analysis - blood loss.R`                                       | 同 `4-1` 的四个口径                                    |
| `4-5` ~ `4-7`   | 年龄亚组   | `4-5 subgroup analysis - age 65-75.R` 等 3 个                                   | 65–74 / 75–84 / ≥85 岁                            |
| `4-8`           | 时间亚组   | `4-8 subgroup analysis - year 2022-2024.R`                                    | 2022 年及以后                                        |
| `4-9`           | 中心亚组   | `4-9 subgroup analysis - hospital without AnYi.R`                             | 剔除安一院                                            |
| `4-10` ~ `4-11` | 性别亚组   | `4-10 subgroup analysis - sex female.R`、`4-11 subgroup analysis - sex male.R` | 女性 / 男性                                          |

---

## 2. 运行环境

| 项   | 要求                                     |
| --- | -------------------------------------- |
| R   | **4.6.1**（建议 ≥ 4.5）                    |
| 编码  | 命令行跑必须加 `LC_ALL=zh_CN.UTF-8`（否则中文路径报错） |
| 依赖包 | 见下（27 个，按用途分组）                         |

| 用途      | 包                                                                                                                             |
| ------- | ----------------------------------------------------------------------------------------------------------------------------- |
| 数据处理    | `data.table` `dplyr` `tidyr` `purrr` `tidyverse` `readr` `stringr` `stringi` `stringdist` `lubridate` `scales` `broom` `grid` |
| 建模与加权   | `geepack` `survey` `MatchIt` `MuMIn`                                                                                          |
| 表格与文档   | `flextable` `officer` `kableExtra` `openxlsx`                                                                                 |
| 描述统计与作图 | `tableone` `cobalt` `Hmisc` `ggplot2` `forestplot` `showtext`                                                                 |

> `WeightIt` 不需要装（倾向评分用的是 `MatchIt` + 自写的稳定化权重）。

每个脚本开头都写死了这三个路径常量（一般在第 7–25 行），换电脑时改这里即可：

```r
D    <- ".../multicenter/20260827"                       # 根目录
DATA <- file.path(D, "JCA_run")                          # 中间大表
OUT  <- file.path(D, "投稿材料/JCA/投稿包_20261009/10-11/out")   # 结果
```

---

## 3. 怎么跑

### 3.1 跑单个脚本
缺什么补什么，脚本之间只通过 `out/` 里的文件传递结果：

```bash
LC_ALL=zh_CN.UTF-8 Rscript "全套代码/3-6 GEE analysis - GEE after IPTW.R"
```

### 3.2 哪些脚本要单独跑

`run_all.R` 覆盖 26 个脚本（`1-1 ~ 3-6`、`4-5 ~ 4-11`、`forest.R`、`forest_subgroup.R`）。
以下 8 个**不在**主控里，需单独跑：

| 脚本                                                 | 原因                                         |
| -------------------------------------------------- | ------------------------------------------ |
| `4-1 mediation analysis - dose.R`                  | 2000 次参数化 bootstrap，耗时较长；只依赖 `out/` 里的既有结果 |
| `4-4 mediation analysis - blood loss.R`            | 同上                                         |
| `4-2 subgroup analysis - operation duration.R`     | 时长分层只读 `out/` 里的事件表与手术表，不进主流程              |
| `4-3 ... invasive vs non-invasive BP monitoring.R` | 依赖 `out/` 里的归层中间表（由本管线之外的核验脚本产出）           |
| `forest_subgroup_M1.R` ~ `M4.R`                    | 亚组森林图，读 `4-5 ~ 4-11` 的输出，单独出图              |

### 3.3 RStudio

打开脚本 → 点右上角 `Source`。清洗段（1-x）建议用这种方式，便于一步一步检查中间结果。

---

## 4. 脚本清单

### 4.1 数据清洗（`1-1` ~ `1-8`）

⚠️ 这一段**必须在同一个 R 会话里顺序执行**（后一个脚本依赖前一个留在内存里的 `data1`）

| 脚本                                                     | 干什么                          | 产出                                                                                  |
| ------------------------------------------------------ | ---------------------------- | ----------------------------------------------------------------------------------- |
| `1-1 data clean - operation type.R`                    | 三家基础表统一格式、合并、只留单药、打手术类型      | 内存 `data1`                                                                          |
| `1-2 ... ASA&GENERAL&DURATION.R`                       | ASA 统一分级、只保留全麻、补麻醉起止时间       | 内存 `data1`                                                                          |
| `1-3 ... fill operation time by event.R`               | 与 event 表对齐 ID、补麻醉结束时间、算麻醉时长 | 内存 `data1`                                                                          |
| `1-4 ... baseline risk factors by diagnosis.R`         | 从术前诊断生成四大合并症 0/1             | 内存 `data1`                                                                          |
| `1-5 ... baseline MAP by preoperative&monitor&data5.R` | 清洗术中监测大表、生成术前低血压排除名单         | 内存 `monitor1`                                                                       |
| `1-6 data clean - AGE&Emergency.R`                     | 最终筛选、清洗年龄性别、**写出 clean 大表**  | `JCA_run/master_clean_*.csv`、`monitor_clean_*.csv`、`event_clean_*.csv`、`术前诊断_*.csv` |
| `1-7 data clean - group of susceptible.R`              | 造易感人群亚组                      | `JCA_run/monitor_below65_before_start*.csv`、`Susceptible_Share_*.csv`               |
| `1-8 data clean - table1.R`                            | Table 1 基线表 + 排除清单           | `JCA_run/baseline_table1.csv`、`Baseline_Table.docx`、`Check_Baseline_Excluded.csv`   |

### 4.2 描述统计（`2-1` ~ `2-3`）

| 脚本                                           | 干什么                             | 产出                                                           |
| -------------------------------------------- | ------------------------------- | ------------------------------------------------------------ |
| `2-1 description - lab data clean.R`         | 三家检验表合并 + 数据质量摘要                | `JCA_run/lab_clean_*.csv`                                    |
| `2-2 description - create GEE event table.R` | **重建低血压事件表**（事件级 + 手术级），含时间基准自检 | `out/Final_Analysis_Events.rds`、`Final_Analysis_Surgery.rds` |
| `2-3 description - table2.R`                 | Table 2 描述性统计                   | `out/Table2_EpisodeLevel.*`、`Table2_SurgeryLevel.*`          |

> `2-2` 是全流程的地基，后面所有模型都读这两个 `.rds`。

### 4.3 主分析（`3-1` ~ `3-6`）

| 脚本                                                | 干什么                                       | 产出                                                                                            |
| ------------------------------------------------- | ----------------------------------------- | --------------------------------------------------------------------------------------------- |
| `3-1 ... build GEE full table(include 0 event).R` | 生成手术级分析表（含 0 事件的手术）                       | `out/Final_Analysis_Long_Table.rds/.csv`                                                      |
| `3-2 ... QIC test.R`                              | 工作相关结构选择（QIC）                             | `out/Table3_GEE_QIC.csv/.docx`                                                                |
| `3-3 ... GEE crude analysis.R`                    | 未加权（crude）对照                              | `out/Crude_Table_Ciprofol.csv/.docx`                                                          |
| `3-4 ... incidence trend time.R`                  | 发生率时间趋势、首次低血压时间                           | `out/Table4_IncidenceTrend*.csv`、`Figure_TimeToFirst_IOH.*`                                   |
| `3-5 ... IPTW.R`                                  | 倾向评分 + 稳定化权重 + 加权后 Table 1                | `out/IPTW_Weights.rds`、`Table1_IPTW*`、`Figure_PS_Overlap.png`                                 |
| `3-6 ... GEE after IPTW.R`                        | ★ **主结果**：四个加权 GEE + 绝对效应 + 重叠区 + E-value | `out/Analysis_Results_main.*`、`Analysis_Results_absolute_main.*`、`Analysis_Results_overlap.*` |

### 4.4 敏感性 / 亚组 / 中介：`4-1` ~ `4-4`（需单独跑）

| 脚本                                                                 | 干什么                                   | 产出                                                       |
| ------------------------------------------------------------------ | ------------------------------------- | -------------------------------------------------------- |
| `4-1 mediation analysis - dose.R`                                  | 诱导剂量中介分析（丙泊酚当量 `ind_equiv`，g-formula） | `out/Mediation_dose.csv`                                 |
| `4-2 subgroup analysis - operation duration.R`                     | 按麻醉时长分层 + 交互检验                        | `out/Subgroup_Duration_*.csv`（3 个）                       |
| `4-3 subgroup analysis - invasive vs non-invasive BP monitoring.R` | 有创/无创血压监测分层 + 交互                      | `out/Subgroup_BPSourceMatched_*.csv`（8 个 CSV + 2 个 docx） |
| `4-4 mediation analysis - blood loss.R`                            | 失血量中介分析（g-formula）                    | `out/Mediation_bloodloss.csv`                            |

**依赖的外部文件**（不在本目录，已放在 `out/`，缺了会报错）：

| 文件                                       | 给谁用            |
| ---------------------------------------- | -------------- |
| `_dose_metrics_v12.csv`                  | 4-1（剂量指标）      |
| `_Perioperative_Vars.csv`                | 4-4（失血量等围术期变量） |
| `Modality_Final_Layer_Assignment_v4.csv` | 4-3（有创/无创归层结果） |


### 4.5 七个亚组分析：`4-5` ~ `4-11`

七个脚本结构完全相同，只改筛选条件和输出标签：

| 脚本                                                | 筛选        | 输出标签           |
| ------------------------------------------------- | --------- | -------------- |
| `4-5 subgroup analysis - age 65-75.R`             | 65–74 岁   | `AGE 65-75`    |
| `4-6 subgroup analysis - age 75-85.R`             | 75–84 岁   | `AGE 75-85`    |
| `4-7 subgroup analysis - age 85+.R`               | ≥85 岁     | `AGE 85+`      |
| `4-8 subgroup analysis - year 2022-2024.R`        | 2022 年及以后 | `2022-2024`    |
| `4-9 subgroup analysis - hospital without AnYi.R` | 剔除安一院     | `without AnYi` |
| `4-10 subgroup analysis - sex female.R`           | 女性        | `female`       |
| `4-11 subgroup analysis - sex male.R`             | 男性        | `male`         |

每个输出 5 个 CSV：`Sensitivity_Analysis_Results_<标签>.csv`、`_Crude_`、`_Absolute_`、`_ESS_`、`Table1_IPTW_<标签>.csv`。

> **输出标签不要改。** `forest_subgroup.R` 与 `forest_subgroup_M1~M4.R` 是按标签字面量去读这些 CSV 的。

### 4.6 出图

| 脚本                              | 产出                                                                              |
| ------------------------------- | ------------------------------------------------------------------------------- |
| `forest.R`                      | 主结果森林图 → `out/ForestPlot_main.png`                                              |
| `forest_subgroup.R`             | M2 亚组森林图 → `Figures/Supplemental Figure S1_Subgroup_burden.*`                   |
| `forest_subgroup_M1.R` ~ `M4.R` | 每个指标一张亚组森林图 → `Figures/Supplemental Figure S1~S4_subgroup_M{1,2,3,4}.png/.tiff` |



## 5. 结果去哪找

全部在 `out/`，按前缀认：

| 前缀                                    | 是什么                                             |
| ------------------------------------- | ----------------------------------------------- |
| `Table1_IPTW*`                        | 加权后基线平衡表                                        |
| `Table2_*` / `Table3_*` / `Table4_*`  | 论文 Table 2 / 3 / 4                              |
| `Analysis_Results_*`                  | ★ 主结果（`main`=全队列，`overlap`=重叠区；`absolute`=绝对效应） |
| `Crude_Table_*`                       | 未加权对照                                           |
| `Sensitivity_Analysis_*`              | 七个亚组分析（4-5 ~ 4-11）                              |
| `Subgroup_*`                          | 时长分层（4-2）、监测方式分层（4-3）                           |
| `Mediation_*`                         | 中介分析：剂量（4-1）、失血量（4-4）                           |
| `Diag_* / Check_* / Explore_* / QA_*` | 过程自检文件，**不是论文结果，但请保留**（核对数据质量时要用）               |

多数表格同时输出 `.csv`（便于再分析）和 `.docx`（便于直接放进稿件）。

---

## 6. 常见问题

**报错 `cannot open the connection` / 找不到文件**
前置脚本没跑。按编号顺序补齐即可，例如跑 4-3 前必须先有 `Final_Analysis_*.rds` 和归层 CSV。

**报错 `file name conversion problem`**
没加 `LC_ALL=zh_CN.UTF-8`。命令行运行必须带这个前缀。

**报错 `not enough contiguous virtual memory`**
内存里还留着清洗段的中间表。重启 R 会话，或先跑 `rm(list = ls()); gc()` 再继续。

**只想重跑某一个表/图**
直接跑对应脚本；它们只读 `out/` 里已有的文件，不会重算上游。

**改了某个公共函数（如 `fit_iptw()`）**
这段代码在多个脚本里各有一份（每个脚本自包含）。**必须同步修改所有含该函数的脚本**，否则不同分析的口径会悄悄不一致。


# Code User Guide (English)

This directory contains the data-analysis code for the manuscript: **34 R scripts**. Executed in numerical order, they rebuild every table, figure, and sensitivity / subgroup / mediation analysis from the cleaned data.
---

[English](./README.md) | [简体中文](./README_zh.md)



**Script naming convention**: `<number> - <type> - <specific content>.R`

| Prefix      | Stage                              | Meaning                                                                                                                                                               |
| ----------- | ---------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `1-x`       | Data cleaning                      | raw tables → clean analysis tables                                                                                                                                    |
| `2-x`       | Descriptive statistics             | rebuild the event table, Table 2                                                                                                                                      |
| `3-x`       | Main analysis                      | IPTW + the four primary GEE models (main results of the manuscript)                                                                                                   |
| `4-x`       | Sensitivity / subgroup / mediation | dose mediation, duration stratification, monitoring-modality stratification, blood-loss mediation, plus seven subgroup analyses by age / calendar year / centre / sex |
| `forest*.R` | Figures                            | forest plots                                                                                                                                                          |

The `4-x` stage holds 11 scripts, ordered by analysis domain:

| Number          | Domain                             | Script                                                                         | Mediation / stratification spec                                                                                                        |
| --------------- | ---------------------------------- | ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------- |
| `4-1`           | Dose mediation                     | `4-1 mediation analysis - dose.R`                                              | four specs: ① crude ② ① + centre×year FE ③ two-arm-coexisting cells ④ ③ + centre×year FE                                               |
| `4-2`           | Duration stratification            | `4-2 subgroup analysis - operation duration.R`                                 | tertiles + clinical cut-points (<120 / 120–240 / ≥240 min) + interaction test                                                          |
| `4-3`           | Monitoring-modality stratification | `4-3 subgroup analysis - invasive vs non-invasive BP monitoring.R`             | stratified by the **actual source of the outcome blood pressure** (invasive 37,913 / non-invasive 47,484 surgeries) + interaction test |
| `4-4`           | Blood-loss mediation               | `4-4 mediation analysis - blood loss.R`                                        | same four specs as `4-1`                                                                                                               |
| `4-5` – `4-7`   | Age subgroups                      | `4-5 subgroup analysis - age 65-75.R` and two siblings                         | 65–74 / 75–84 / ≥85 years                                                                                                              |
| `4-8`           | Calendar-year subgroup             | `4-8 subgroup analysis - year 2022-2024.R`                                     | 2022 and later                                                                                                                         |
| `4-9`           | Centre subgroup                    | `4-9 subgroup analysis - hospital without AnYi.R`                              | excluding AnYi                                                                                                                         |
| `4-10` – `4-11` | Sex subgroups                      | `4-10 subgroup analysis - sex female.R`, `4-11 subgroup analysis - sex male.R` | female / male                                                                                                                          |

---

## 2. Runtime environment

| Item     | Requirement                                                                                |
| -------- | ------------------------------------------------------------------------------------------ |
| R        | **4.6.1** (≥ 4.5 recommended)                                                              |
| Encoding | command-line runs must set `LC_ALL=zh_CN.UTF-8` (otherwise non-ASCII paths raise an error) |
| Packages | see below (27, grouped by purpose)                                                         |

| Purpose                          | Packages                                                                                                                      |
| -------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Data wrangling                   | `data.table` `dplyr` `tidyr` `purrr` `tidyverse` `readr` `stringr` `stringi` `stringdist` `lubridate` `scales` `broom` `grid` |
| Modelling and weighting          | `geepack` `survey` `MatchIt` `MuMIn`                                                                                          |
| Tables and documents             | `flextable` `officer` `kableExtra` `openxlsx`                                                                                 |
| Descriptive statistics and plots | `tableone` `cobalt` `Hmisc` `ggplot2` `forestplot` `showtext`                                                                 |

> `WeightIt` is not required (the propensity score uses `MatchIt` plus a hand-written stabilisation step).

Each script hard-codes three path constants (typically lines 7–25); change them there when moving to another machine:

```r
D    <- ".../multicenter/20260827"                       # root directory
DATA <- file.path(D, "JCA_run")                          # intermediate tables
OUT  <- file.path(D, "投稿材料/JCA/投稿包_20261009/10-11/out")   # results
```

## 3. How to run

### 3.1 Running a single script 

Run only what is missing; scripts communicate exclusively through files in `out/`:

```bash
LC_ALL=zh_CN.UTF-8 Rscript "全套代码/3-6 GEE analysis - GEE after IPTW.R"
```

### 3.2 Scripts that must be run separately

`run_all.R` covers 26 scripts (`1-1` – `3-6`, `4-5` – `4-11`, `forest.R`, `forest_subgroup.R`).
The following eight are **not** in the runner and must be executed individually:

| Script                                             | Why                                                                                                  |
| -------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `4-1 mediation analysis - dose.R`                  | 2000 parametric bootstrap replicates — slow, and it only reads files already present in `out/`       |
| `4-4 mediation analysis - blood loss.R`            | as above                                                                                             |
| `4-2 subgroup analysis - operation duration.R`     | duration stratification reads only the event and surgery tables in `out/`; not part of the main flow |
| `4-3 ... invasive vs non-invasive BP monitoring.R` | depends on the layer-assignment tables in `out/`, which are produced outside this pipeline           |
| `forest_subgroup_M1.R` – `M4.R`                    | subgroup forest plots; they read the output of `4-5` – `4-11`                                        |

### 3.3 RStudio

Open the script and click `Source` in the top-right corner. For the cleaning stage (`1-x`) this is the recommended route, so that intermediate objects can be inspected step by step.

---

## 4. Script inventory

### 4.1 Data cleaning (`1-1` to `1-8`)

⚠️ This stage **must be executed sequentially within a single R session** (each script depends on the `data1` object left in memory by the previous one). 

| Script                                                 | What it does                                                                                                          | Output                                                                                 |
| ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| `1-1 data clean - operation type.R`                    | harmonise the format of the three hospital tables, merge them, keep single-drug cases, derive the surgical speciality | in-memory `data1`                                                                      |
| `1-2 ... ASA&GENERAL&DURATION.R`                       | unify ASA grading, restrict to general anaesthesia, add anaesthesia start/end times                                   | in-memory `data1`                                                                      |
| `1-3 ... fill operation time by event.R`               | align identifiers with the event table, backfill anaesthesia end time, compute anaesthesia duration                   | in-memory `data1`                                                                      |
| `1-4 ... baseline risk factors by diagnosis.R`         | derive the four binary comorbidity indicators from the preoperative diagnoses                                         | in-memory `data1`                                                                      |
| `1-5 ... baseline MAP by preoperative&monitor&data5.R` | clean the large intraoperative monitoring table, build the exclusion list for preoperative hypotension                | in-memory `monitor1`                                                                   |
| `1-6 data clean - AGE&Emergency.R`                     | final cohort filtering, clean age and sex, **write the clean analysis tables**                                        | `JCA_run/master_clean_*.csv`, `monitor_clean_*.csv`, `event_clean_*.csv`, `术前诊断_*.csv` |
| `1-7 data clean - group of susceptible.R`              | build the susceptible subgroup                                                                                        | `JCA_run/monitor_below65_before_start*.csv`, `Susceptible_Share_*.csv`                 |
| `1-8 data clean - table1.R`                            | baseline Table 1 plus the exclusion list                                                                              | `JCA_run/baseline_table1.csv`, `Baseline_Table.docx`, `Check_Baseline_Excluded.csv`    |

### 4.2 Descriptive statistics (`2-1` to `2-3`)

| Script                                       | What it does                                                                                                               | Output                                                        |
| -------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| `2-1 description - lab data clean.R`         | merge the laboratory tables of the three centres, summarise data quality                                                   | `JCA_run/lab_clean_*.csv`                                     |
| `2-2 description - create GEE event table.R` | **rebuild the intraoperative hypotension event table** (episode-level and surgery-level), with a time-reference self-check | `out/Final_Analysis_Events.rds`, `Final_Analysis_Surgery.rds` |
| `2-3 description - table2.R`                 | Table 2 descriptive statistics                                                                                             | `out/Table2_EpisodeLevel.*`, `Table2_SurgeryLevel.*`          |

> `2-2` is the foundation of the whole pipeline: every downstream model reads these two `.rds` files.

### 4.3 Main analysis (`3-1` to `3-6`)

| Script                                            | What it does                                                                                | Output                                                                                          |
| ------------------------------------------------- | ------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `3-1 ... build GEE full table(include 0 event).R` | build the surgery-level analysis table (including operations with zero events)              | `out/Final_Analysis_Long_Table.rds/.csv`                                                        |
| `3-2 ... QIC test.R`                              | select the working correlation structure (QIC)                                              | `out/Table3_GEE_QIC.csv/.docx`                                                                  |
| `3-3 ... GEE crude analysis.R`                    | unweighted (crude) reference analysis                                                       | `out/Crude_Table_Ciprofol.csv/.docx`                                                            |
| `3-4 ... incidence trend time.R`                  | incidence trend over time, time to first hypotension episode                                | `out/Table4_IncidenceTrend*.csv`, `Figure_TimeToFirst_IOH.*`                                    |
| `3-5 ... IPTW.R`                                  | propensity score, stabilised weights, weighted Table 1                                      | `out/IPTW_Weights.rds`, `Table1_IPTW*`, `Figure_PS_Overlap.png`                                 |
| `3-6 ... GEE after IPTW.R`                        | ★ **main results**: the four weighted GEE models, absolute effects, overlap region, E-value | `out/Analysis_Results_main.*`, `Analysis_Results_absolute_main.*`, `Analysis_Results_overlap.*` |

### 4.4 Sensitivity / subgroup / mediation: `4-1` to `4-4` (run separately)

| Script                                                             | What it does                                                                           | Output                                                |
| ------------------------------------------------------------------ | -------------------------------------------------------------------------------------- | ----------------------------------------------------- |
| `4-1 mediation analysis - dose.R`                                  | mediation analysis for the induction dose (propofol-equivalent `ind_equiv`, g-formula) | `out/Mediation_dose.csv`                              |
| `4-2 subgroup analysis - operation duration.R`                     | stratify by anaesthesia duration, test for interaction                                 | `out/Subgroup_Duration_*.csv` (3 files)               |
| `4-3 subgroup analysis - invasive vs non-invasive BP monitoring.R` | stratify by invasive vs non-invasive blood-pressure monitoring, interaction            | `out/Subgroup_BPSourceMatched_*.csv` (8 CSV + 2 docx) |
| `4-4 mediation analysis - blood loss.R`                            | mediation analysis for blood loss (g-formula)                                          | `out/Mediation_bloodloss.csv`                         |

**External files required** (not in this directory; kept in `out/`; the scripts fail if they are missing):

| File                                     | Used by                                            |
| ---------------------------------------- | -------------------------------------------------- |
| `_dose_metrics_v12.csv`                  | 4-1 (dose metrics)                                 |
| `_Perioperative_Vars.csv`                | 4-4 (perioperative variables including blood loss) |
| `Modality_Final_Layer_Assignment_v4.csv` | 4-3 (monitoring-modality layer assignment)         |


### 4.5 Seven subgroup analyses: `4-5` to `4-11` 

The seven scripts are structurally identical; only the filtering condition and the output label change:

| Script                                            | Filter         | Output label   |
| ------------------------------------------------- | -------------- | -------------- |
| `4-5 subgroup analysis - age 65-75.R`             | 65–74 years    | `AGE 65-75`    |
| `4-6 subgroup analysis - age 75-85.R`             | 75–84 years    | `AGE 75-85`    |
| `4-7 subgroup analysis - age 85+.R`               | ≥85 years      | `AGE 85+`      |
| `4-8 subgroup analysis - year 2022-2024.R`        | 2022 and later | `2022-2024`    |
| `4-9 subgroup analysis - hospital without AnYi.R` | excluding AnYi | `without AnYi` |
| `4-10 subgroup analysis - sex female.R`           | female         | `female`       |
| `4-11 subgroup analysis - sex male.R`             | male           | `male`         |

Each writes five CSV files: `Sensitivity_Analysis_Results_<label>.csv`, `_Crude_`, `_Absolute_`, `_ESS_`, and `Table1_IPTW_<label>.csv`.

> **Do not change the output labels.** `forest_subgroup.R` and `forest_subgroup_M1~M4.R` look these CSV files up by their literal label.

### 4.6 Figures

| Script                           | Output                                                                                                  |
| -------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `forest.R`                       | forest plot of the main results → `out/ForestPlot_main.png`                                             |
| `forest_subgroup.R`              | M2 subgroup forest plot → `Figures/Supplemental Figure S1_Subgroup_burden.*`                            |
| `forest_subgroup_M1.R` to `M4.R` | one subgroup forest plot per metric → `Figures/Supplemental Figure S1~S4_subgroup_M{1,2,3,4}.png/.tiff` |


---

## 5. Where to find the results

Everything is in `out/`; identify files by prefix:

| Prefix                                | What it is                                                                                                  |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `Table1_IPTW*`                        | post-weighting baseline balance table                                                                       |
| `Table2_*` / `Table3_*` / `Table4_*`  | Tables 2 / 3 / 4 of the manuscript                                                                          |
| `Analysis_Results_*`                  | ★ main results (`main` = full cohort, `overlap` = overlap region, `absolute` = absolute effects)            |
| `Crude_Table_*`                       | unweighted reference analysis                                                                               |
| `Sensitivity_Analysis_*`              | the seven subgroup analyses (4-5 – 4-11)                                                                    |
| `Subgroup_*`                          | duration stratification (4-2), monitoring-modality stratification (4-3)                                     |
| `Mediation_*`                         | mediation analyses: dose (4-1), blood loss (4-4)                                                            |
| `Diag_* / Check_* / Explore_* / QA_*` | process self-check files — **not manuscript results, but please keep them** (needed to verify data quality) |

Most tables are written in both `.csv` (for re-analysis) and `.docx` (ready to paste into the manuscript).

---

## 6. Common problems

**Error `cannot open the connection` / file not found**
A prerequisite script has not been run. Complete the sequence in numerical order; for example, `Final_Analysis_*.rds` and the layer-assignment CSV must exist before `4-3` can run.

**Error `file name conversion problem`**
`LC_ALL=zh_CN.UTF-8` was not set. This prefix is mandatory for command-line runs.

**Error `not enough contiguous virtual memory`**
Intermediate tables from the cleaning stage are still held in memory. Restart the R session, or run `rm(list = ls()); gc()` before continuing.

**Re-running only one table or figure**
Run the corresponding script directly; it only reads files already present in `out/` and never recomputes upstream stages.

**After editing a shared helper function (e.g. `fit_iptw()`)**
This function exists as a separate copy in each script (every script is self-contained). **All scripts containing that function must be edited in parallel**, otherwise the analytic conventions of different analyses will silently diverge.


