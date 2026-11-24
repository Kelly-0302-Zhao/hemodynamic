# hemodynamic
R code about hemodynamics
# Code User Guide 

This directory contains the data-analysis code for the manuscript: **36 R scripts**, each self-contained and runnable individually. Run them one by one in numerical order to rebuild every table, figure, and sensitivity / subgroup / mediation analysis from the cleaned data.

## 1. What is in this directory

**Script naming convention**: `<number> - <type> - <specific content>.R`

| Prefix      | Stage                              | Meaning                                                                                                                                                                                                                              |
| ----------- | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `1-x`       | Data cleaning                      | raw tables → clean analysis tables                                                                                                                                                                                                   |
| `2-x`       | Descriptive statistics             | rebuild the event table, Table 2                                                                                                                                                                                                     |
| `3-x`       | Main analysis                      | IPTW + the four primary GEE models (main results of the manuscript)                                                                                                                                                                  |
| `4-x`       | Sensitivity / subgroup / mediation | dose mediation, duration stratification, monitoring-modality stratification, concurrent-cohort sensitivity, absolute-outcomes reporting, harmonized downsampling, plus seven subgroup analyses by age / calendar year / centre / sex |
| `forest*.R` | Figures                            | forest plots                                                                                                                                                                                                                         |

The `4-x` stage holds 13 scripts, ordered by analysis domain:

| Number         | Domain                                  | Script                                                                        | Mediation / stratification spec                                                                                                                                                                                  |
| -------------- | --------------------------------------- | ----------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `4-1`          | Dose mediation                          | `4-1 mediation analysis - dose.R`                                             | four specs: ① crude ② ① + centre×year FE ③ two-arm-coexisting cells ④ ③ + centre×year FE                                                                                                                         |
| `4-2`          | Duration stratification                 | `4-2 subgroup analysis - operation duration.R`                                | tertiles + clinical cut-points (<120 / 120–240 / ≥240 min) + interaction test                                                                                                                                    |
| `4-3`          | Monitoring-modality stratification      | `4-3 subgroup analysis - invasive vs non-invasive BP monitoring.R`            | stratified by the **actual source of the outcome blood pressure** (invasive 37,913 / non-invasive 47,484 surgeries) + interaction test                                                                           |
| `4-4` – `4-6`  | Age subgroups                           | `4-4 subgroup analysis - age 65-75.R` and two siblings                        | 65–74 / 75–84 / ≥85 years                                                                                                                                                                                        |
| `4-7`          | Calendar-year subgroup                  | `4-7 subgroup analysis - year 2022-2024.R`                                    | 2022 and later                                                                                                                                                                                                   |
| `4-8`          | Centre subgroup                         | `4-8 subgroup analysis - hospital without AnYi.R`                             | excluding AnYi                                                                                                                                                                                                   |
| `4-9` – `4-10` | Sex subgroups                           | `4-9 subgroup analysis - sex female.R`, `4-10 subgroup analysis - sex male.R` | female / male                                                                                                                                                                                                    |
| `4-11`         | Concurrent centre-specific cohort       | `4-11 sensitivity analysis - concurrent centre-specific cohort.R`             | subset to centre×year cells where ciprofol was actually available (ZhongDa 2022–2024 + NanYi 2024); refit PS + IPTW + four GEE within the subset to address positivity                                           |
| `4-12`         | Absolute outcomes + log-transform check | `4-12 sensitivity analysis - M2M3 absolute outcomes and log transform.R`      | weighted geometric / arithmetic mean AUC, median / mean duration with 2000-replicate bootstrap CI; compares `log(AUC+0.1)` vs `log(AUC)` (AUC min = 0.20, no zeros)                                              |
| `4-13`         | Harmonized downsampling                 | `4-13 sensitivity analysis - harmonized downsampling.R`                       | downsample invasive-layer MAP readings to 5-min intervals (simulating non-invasive sampling frequency), rebuild the event table, refit GEE — tests whether the invasive-layer effect is driven by detection bias |

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

---

## 3. How to run — one script at a time

Every script is self-contained: it reads the files it needs from `out/` (or `JCA_run/` for the cleaning stage) and writes its results back to `out/`. There is no master runner — execute each script individually, in numerical order, and only run the ones you need.

### 3.1 From the command line

```bash
LC_ALL=zh_CN.UTF-8 Rscript "全套代码/3-6 GEE analysis - GEE after IPTW.R"
```

The `LC_ALL=zh_CN.UTF-8` prefix is mandatory for command-line runs (the project paths contain non-ASCII characters).

### 3.2 From RStudio

Open the script and click `Source` in the top-right corner. This is the recommended route for the cleaning stage (`1-x`), so intermediate objects can be inspected step by step.

### 3.3 Execution order and prerequisites

Scripts communicate only through files on disk. Run them in numerical order so that each script's inputs are already present in `out/`. The dependency chain is:

```
1-x  →  JCA_run/*.csv           (cleaned tables)
2-2  →  out/Final_Analysis_*.rds (event + surgery tables — the foundation)
3-1  →  out/Final_Analysis_Long_Table.rds
3-5  →  out/IPTW_Weights.rds
3-6  →  out/Analysis_Results_main.*          ★ main results
4-x  →  out/<domain-specific outputs>        (each reads the .rds above)
forest*.R → out/Figures/                     (reads 4-x outputs)
```

> **Cleaning stage (`1-1` to `1-8`)** must be executed sequentially within a single R session — each script depends on the `data1` object left in memory by the previous one. Only `1-6` writes to disk. If you already have the `JCA_run/*.csv` files, you can skip this stage entirely and start from `2-2`.

### 3.4 Scripts that take longer

These scripts are correct but slow (bootstrap or table rebuild); budget extra time for them:

| Script                                                                   | Why it takes longer                                                                          |
| ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------- |
| `4-1 mediation analysis - dose.R`                                        | 2000 parametric bootstrap replicates                                                         |
| `4-12 sensitivity analysis - M2M3 absolute outcomes and log transform.R` | 2000-replicate bootstrap for absolute-outcome CIs                                            |
| `4-13 sensitivity analysis - harmonized downsampling.R`                  | rebuilds the event table from the raw monitoring CSV (~72 million rows) before refitting GEE |

---

## 4. Script inventory

### 4.1 Data cleaning (`1-1` to `1-8`)

⚠️ This stage **must be executed sequentially within a single R session** (each script depends on the `data1` object left in memory by the previous one). Only `1-6` writes to disk.

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

### 4.4 Sensitivity / mediation / supplementary analyses: `4-1` to `4-3`, `4-11` to `4-13`

| Script                                                                   | What it does                                                                                                                                  | Output                                                                                                                |
| ------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `4-1 mediation analysis - dose.R`                                        | mediation analysis for the induction dose (propofol-equivalent `ind_equiv`, g-formula)                                                        | `out/Mediation_dose.csv`                                                                                              |
| `4-2 subgroup analysis - operation duration.R`                           | stratify by anaesthesia duration, test for interaction                                                                                        | `out/Subgroup_Duration_*.csv` (3 files)                                                                               |
| `4-3 subgroup analysis - invasive vs non-invasive BP monitoring.R`       | stratify by invasive vs non-invasive blood-pressure monitoring, interaction                                                                   | `out/Subgroup_BPSourceMatched_*.csv` (8 CSV + 2 docx)                                                                 |
| `4-11 sensitivity analysis - concurrent centre-specific cohort.R`        | restrict to centre×year cells where ciprofol was actually used (ZhongDa 2022–2024 + NanYi 2024), refit PS + IPTW + four GEE within the subset | `out/Sensitivity_Analysis_Results_Concurrent_CentreSpecific.csv` (+ `_Crude_`, `_Absolute_`, `_ESS_`, `Table1_IPTW_`) |
| `4-12 sensitivity analysis - M2M3 absolute outcomes and log transform.R` | report weighted geometric / arithmetic mean AUC, median / mean duration with bootstrap CI; compare `log(AUC+0.1)` vs `log(AUC)` for M2/M3     | `out/M2M3_Absolute_Outcomes.csv`, `out/M2M3_LogTransform_Comparison.csv`                                              |
| `4-13 sensitivity analysis - harmonized downsampling.R`                  | downsample invasive-layer MAP to 5-min intervals, rebuild event table, refit GEE; compare with original invasive and non-invasive results     | `out/Harmonized_Downsampling_Invasive_5min.csv`, `out/Harmonized_Downsampling_Comparison.csv`                         |

**External files required** (not in this directory; kept in `out/` or `JCA_run/`; the scripts fail if they are missing):

| File                                       | Location   | Used by                                                            |
| ------------------------------------------ | ---------- | ------------------------------------------------------------------ |
| `_dose_metrics_v12.csv`                    | `out/`     | 4-1 (dose metrics)                                                 |
| `Modality_Final_Layer_Assignment_v4.csv`   | `out/`     | 4-3, 4-13 (monitoring-modality layer assignment)                   |
| `Modality_Layer_For_Subgroup_v4.csv`       | `out/`     | 4-13 (layer assignment for subgroup linkage)                       |
| `Subgroup_Monitoring_Modality_Results.csv` | `out/`     | 4-13 (reference results for comparison)                            |
| `monitor_clean_72460818obs.csv`            | `JCA_run/` | 4-13 (raw intraoperative monitoring table — output of `1-5`/`1-6`) |
| `master_clean_85397obs.csv`                | `JCA_run/` | 4-13 (raw master table — output of `1-6`)                          |

> **The dose domain contains mediation only.** The separate dose-as-covariate attenuation analysis is not part of this delivered script set; its script and outputs are archived under `.workbuddy/backup/已归档_剂量协变量调整_20261107/`. Its role in this project is corroborative: after adding dose to M1 the OR moves only from 0.564 to 0.567 (change-in-estimate 0.9%, below the 10% threshold), i.e. the outcome association is not driven by dose differences.

### 4.5 Seven subgroup analyses: `4-4` to `4-10`

The seven scripts are structurally identical; only the filtering condition and the output label change:

| Script                                            | Filter         | Output label   |
| ------------------------------------------------- | -------------- | -------------- |
| `4-4 subgroup analysis - age 65-75.R`             | 65–74 years    | `AGE 65-75`    |
| `4-5 subgroup analysis - age 75-85.R`             | 75–84 years    | `AGE 75-85`    |
| `4-6 subgroup analysis - age 85+.R`               | ≥85 years      | `AGE 85+`      |
| `4-7 subgroup analysis - year 2022-2024.R`        | 2022 and later | `2022-2024`    |
| `4-8 subgroup analysis - hospital without AnYi.R` | excluding AnYi | `without AnYi` |
| `4-9 subgroup analysis - sex female.R`            | female         | `female`       |
| `4-10 subgroup analysis - sex male.R`             | male           | `male`         |

Each writes five CSV files: `Sensitivity_Analysis_Results_<label>.csv`, `_Crude_`, `_Absolute_`, `_ESS_`, and `Table1_IPTW_<label>.csv`.

> **Do not change the output labels.** `forest_subgroup.R` and `forest_subgroup_M1~M4.R` look these CSV files up by their literal label.

### 4.6 Figures

| Script                           | Output                                                                                                  |
| -------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `forest.R`                       | forest plot of the main results → `out/ForestPlot_main.png`                                             |
| `forest_subgroup.R`              | M2 subgroup forest plot → `Figures/Supplemental Figure S1_Subgroup_burden.*`                            |
| `forest_subgroup_M1.R` to `M4.R` | one subgroup forest plot per metric → `Figures/Supplemental Figure S1~S4_subgroup_M{1,2,3,4}.png/.tiff` |

> The forest plots are written to `投稿包_20261009/乱七八糟/Figures`, **not** to `out/`.

---

## 5. Where to find the results

Everything is in `out/`; identify files by prefix:

| Prefix                                | What it is                                                                                                  |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `Table1_IPTW*`                        | post-weighting baseline balance table                                                                       |
| `Table2_*` / `Table3_*` / `Table4_*`  | Tables 2 / 3 / 4 of the manuscript                                                                          |
| `Analysis_Results_*`                  | ★ main results (`main` = full cohort, `overlap` = overlap region, `absolute` = absolute effects)            |
| `Crude_Table_*`                       | unweighted reference analysis                                                                               |
| `Sensitivity_Analysis_*`              | the seven subgroup analyses (4-4 – 4-10) plus the concurrent centre-specific cohort (4-11)                  |
| `Subgroup_*`                          | duration stratification (4-2), monitoring-modality stratification (4-3)                                     |
| `Mediation_*`                         | mediation analysis: dose (4-1)                                                                              |
| `M2M3_*`                              | absolute outcomes and log-transform comparison (4-12)                                                       |
| `Harmonized_Downsampling_*`           | invasive-layer 5-min downsampling sensitivity (4-13)                                                        |
| `Diag_* / Check_* / Explore_* / QA_*` | process self-check files — **not manuscript results, but please keep them** (needed to verify data quality) |
