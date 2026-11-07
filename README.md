# hemodynamic
R code about hemodynamics
（Chinese）
脚本命名规则**：`编号 - 类型 - 具体内容.R`

| 编号                         | 段         | 含义                       |
| -------------------------- | --------- | ------------------------ |
| `1-x`                      | 数据清洗      | 原始表 → 干净的分析表             |
| `2-x`                      | 描述统计      | 描述性统计分析，建立事件级及手术级表格          |
| `3-x`                      | 主分析       | IPTW + 四个 GEE 主模型 |
| `4-x`                      | 敏感性/亚组/中介 | 剂量中介分析，时长分层分析，监测方式分层分析，失血量中介分析，年龄亚组分析，时间亚组分析，中心亚组分析，性别亚组分析        |
| `forest*.R`                | 出图        | 森林图                      |

###1 数据清洗`1-x` 

⚠️ 这一段**必须在同一个 R 会话里顺序执行**（后一个脚本依赖前一个留在内存里的 `data1`）

| 脚本                                                     | 干什么                          | 产出                                                                                  |
| ------------------------------------------------------ | ---------------------------- | ----------------------------------------------------------------------------------- |
| `1-1 data clean - operation type.R`                    | 三家基础表统一格式、合并、只留单药、打手术类型      | 内存 `data1`                                                                          |
| `1-2 ... ASA&GENERAL&DURATION.R`                       | ASA 统一分级、只保留全麻、补麻醉起止时间       | 内存 `data1`                                                                          |
| `1-3 ... fill operation time by event.R`               | 与 event 表对齐 ID、补麻醉结束时间、算麻醉时长 | 内存 `data1`                                                                          |
| `1-4 ... baseline risk factors by diagnosis.R`         | 从术前诊断生成四大合并症 0/1             | 内存 `data1`                                                                          |
| `1-5 ... baseline MAP by preoperative&monitor&data5.R` | 清洗术中监测大表、生成术前低血压排除名单         | 内存 `monitor1`                                                                       |
| `1-6 data clean - AGE&Emergency.R`                     | 最终筛选、清洗年龄性别、**写出 clean 大表**  | `master_clean_*.csv`、`monitor_clean_*.csv`、`event_clean_*.csv`、`术前诊断_*.csv` |
| `1-7 data clean - group of susceptible.R`              | 造易感人群亚组                      | `monitor_below65_before_start*.csv`、`Susceptible_Share_*.csv`                       |
| `1-8 data clean - table1.R`                            | Table 1 基线表 + 排除清单           | `baseline_table1.csv`、`Baseline_Table.docx`、`Check_Baseline_Excluded.csv`       |

###2 描述统计`2-X` 

| 脚本                                           | 干什么                             | 产出                                                           |
| -------------------------------------------- | ------------------------------- | ------------------------------------------------------------ |
| `2-1 description - lab data clean.R`         | 三家检验表合并 + 数据质量摘要                | `lab_clean_*.csv`                                    |
| `2-2 description - create GEE event table.R` | **重建低血压事件表**（事件级 + 手术级），含时间基准自检 | `Final_Analysis_Events.rds`、`Final_Analysis_Surgery.rds` |
| `2-3 description - table2.R`                 | Table 2 描述性统计                   | `Table2_EpisodeLevel.*`、`Table2_SurgeryLevel.*`              |

###3 主分析`3-X`

| 脚本                                                | 干什么                                       | 产出                                                                                        |
| ------------------------------------------------- | ----------------------------------------- | ----------------------------------------------------------------------------------------- |
| `3-1 ... build GEE full table(include 0 event).R` | 生成手术级分析表（含 0 事件的手术）                       | `Final_Analysis_Long_Table.rds/.csv`                                                      |
| `3-2 ... QIC test.R`                              | 工作相关结构选择（QIC）                             | `Table3_GEE_QIC.csv/.docx`                                                                |
| `3-3 ... GEE crude analysis.R`                    | 未加权（crude）对照                              | `Crude_Table_Ciprofol.csv/.docx`                                                          |
| `3-4 ... incidence trend time.R`                  | 发生率时间趋势、首次低血压时间                           | `Table4_IncidenceTrend*.csv`、`Figure_TimeToFirst_IOH.*`                                   |
| `3-5 ... IPTW.R`                                  | 倾向评分 + 稳定化权重 + 加权后 Table 1                | `IPTW_Weights.rds`、`Table1_IPTW*`、`Figure_PS_Overlap.png`                                 |
| `3-6 ... GEE after IPTW.R`                        | ★ **主结果**：四个加权 GEE + 绝对效应 + 重叠区 + E-value | `Analysis_Results_main.*`、`Analysis_Results_absolute_main.*`、`Analysis_Results_overlap.*` |

