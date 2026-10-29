# BRIWECS-ML Phase 2B-Closure
## Association–Prediction Decoupling & Novelty Assessment
### Codex Scientific Execution Specification v2.4

**项目状态：** Phase 2B-Pilot CONDITIONAL GO → Phase 2B-Closure AUTHORIZED

**本阶段性质：** 科学假设验证、预测误差分解、病害关联稳健性分析及论文创新性判定

**核心病害：** Wheat stripe rust

**数据版本：** BRIWECS Official Figshare v1

**资源原则：** 本地CPU优先，严格限制新模型训练

**终点：** Manuscript GO / CONDITIONAL GO / NO-GO

**禁止：** 新湿实验、数据版本混用、无必要的模型扩展、提前进入完整Phase 2C

---

# 0. 研究目标与科学边界

本阶段是BRIWECS-ML项目的一次正式科学关闭节点。

前期已经完成：

- Phase 1：官方数据审计与队列可行性分析。
- Phase 2A：跨环境产量预测与严格分组验证。
- Phase 2B-Pilot：条锈病和Septoria的病害信息增量预测分析。

当前核心结果：

- 条锈病完整队列：8,539条、13个地点×年份。
- 自然感染子集：5,899条、10个地点×年份。
- 完整队列中，条锈病与产量在12/13个环境中呈负相关。
- 自然感染子集中，9/10个环境呈负相关。
- CatBoost F1加入条锈病后，汇总ΔMAE为+1.581 dt/ha。
- 环境等权ΔMAE为+0.434 dt/ha，95%描述性区间跨零。
- 去掉QLB 2015后，环境等权ΔMAE转为−0.143 dt/ha。
- 自然感染子集CatBoost F1的汇总ΔMAE仅+0.073 dt/ha。
- Septoria不支持作为稳定跨环境预测增量的主要研究对象。

以上结果不代表本论文已经具有充分创新性。

本阶段必须回答：

**Q1. 条锈病与产量的负向关联，在充分考虑品种身份、环境和管理后是否依然存在？**

**Q2. 加入病害信息后的预测改善或恶化，主要体现为环境平均产量估计，还是环境内部样地与品种的预测差异？**

**Q3. 当测试样本严格固定为相同自然感染环境时，训练集中包含人工接种相关数据是否改变病害特征的增量预测价值？**

**Q4. 以上现象相对于已有小麦病害–产量及机器学习研究，是否真正形成具有独立学术价值的新结论？**

上述四项完成后，必须签发明确的论文路线裁定。

不要为了追求正向结果继续增加新的分析任务。

---

# 1. 工程继承与执行范围

项目根目录：

D:\机器学习\BRIWECS_ML_Phase1

优先读取：

Phase2A/reports/BRIWECS_Phase2A_Modeling_Report.md

Phase2A/reports/BRIWECS_Phase2A_Scientific_Assessment.md

Phase2A/results/predictions/out_of_fold_predictions.parquet

Phase2B_Pilot/reports/BRIWECS_Phase2B_Pilot_Modeling_Report.md

Phase2B_Pilot/reports/BRIWECS_Phase2B_Pilot_Scientific_Assessment.md

Phase2B_Pilot/reports/BRIWECS_Phase2B_Pilot_GO_NO_GO.md

Phase2B_Pilot/results/metrics/disease_incremental_performance.csv

Phase2B_Pilot/results/metrics/environment_incremental_metrics.csv

Phase2B_Pilot/results/metrics/within_environment_association.csv

Phase2B_Pilot/results/metrics/inoculation_sensitivity.csv

Phase2B_Pilot/results/predictions/disease_pilot_oof_predictions.parquet

同时读取前期队列冻结文件、模型配置、数据字典、实际训练脚本和SHA256清单。

实际名称不同则依据README定位。

创建：

Phase2B_Closure/

下设：

- data/
- config/
- scripts/
- results/tables/
- results/figures/
- reports/
- logs/

不能覆盖任何既有交付文件。

工作完成后生成独立的复现入口和交付清单。

---

# 2. Work Package 1 — Prediction Error Decomposition

## 2.1 研究目的

解释此前发现的跨环境预测不稳定现象。

重点检查：

条锈病信息是否主要改变了测试环境的平均产量预测偏差，而没有稳定改善环境内部不同样地的预测差异？

本部分原则上完全复用既有OOF预测，无须重新训练模型。

## 2.2 分析范围

主要分析：

Stripe rust，CatBoost，M0/M1。

其中：

- M0：品种＋管理。
- M1：品种＋管理＋条锈病严重度。

辅助分析：

Stripe rust，Ridge，M0/M1。

如果现有文件完整、无需额外计算，可增加M2/M3作为天气增强的对照。

只采用真实LOSYO外层测试预测。

## 2.3 必须完成的误差分解

对每个地点×年份环境，基于相同测试样本计算：

1. 实际平均产量。
2. 预测平均产量。
3. 平均预测偏差。
4. 环境内部真实产量标准差。
5. 原始MAE与RMSE。
6. 环境内部中心化误差。
7. 环境内部预测排序Spearman相关。
8. 同环境下M0与M1的配对误差差异。

对平方误差使用如下精确恒等式：

MSE_environment = (Mean prediction error)² + Variance of prediction error

其中误差方差使用与MSE一致的总体方差定义。

进一步计算：

- 病害特征对环境平均偏差平方项的改变。
- 病害特征对中心化残差方差项的改变。
- 两部分分别对总MSE差异的贡献。

注意：MAE没有同样的简单加性分解。MAE应分别报告，但不要将MSE分解结果直接等同于MAE分解。

## 2.4 环境内部品种表现

在具有足够品种和产量变异的环境内，检查：

- M0和M1的样地级排序相关。
- 按品种汇总后的排序相关（仅在重复结构支持时）。
- 病害特征是否改善高产和低产品种的相对辨别能力。

环境中心化只能作为事后诊断。

**严禁将测试集的真实平均产量用于模型预测校准，再把调整后的预测误差报告为实际外部泛化性能。**

## 2.5 结果解释

区分三种情况：

A. 病害信息改善环境平均水平估计。

B. 病害信息改善环境内部样地或品种差异的预测。

C. 病害信息在两方面均无稳定增益。

如果无法明确分解此前的性能异质性，应标记PARTIAL_EXPLANATION。

不得为了形成科学故事而强制选择某一种解释。

## 2.6 输出

- `01_environment_error_decomposition.csv`
- `02_within_environment_ranking.csv`
- `03_disease_delta_error_components.csv`
- `01_Error_Decomposition_Findings.md`

---

# 3. Work Package 2 — Adjusted Disease–Yield Association

## 3.1 研究目的

判断前期12/13环境的条锈病–产量负相关是否只是表面现象。

必须区分：

- 未调整的病害–产量关联。
- 环境控制后的关联。
- 环境与管理控制后的关联。
- 进一步控制品种身份后的关联。

分析对象为审计后的条锈病严重度和Seedyield。

不能使用原作者清洗过程错误生成的零值。

重新确认原始小数评分恢复记录及最终零值比例，不要把原始零值比例误当成纠正后的最终比例。

## 3.2 逐层模型

建议采用四层分析。

**A0 — Crude association**

Yield ~ Stripe rust severity

**A1 — Environment-adjusted**

Yield ~ Stripe rust severity + Site-year fixed effects

**A2 — Environment and management-adjusted**

Yield ~ Stripe rust severity + Site-year × Management fixed effects

**A3 — Genotype-adjusted**

Yield ~ Stripe rust severity + Site-year × Management fixed effects + Cultivar fixed effects

使用实际观测数据可支持的固定效应结构。

如果出现完美共线性、秩亏或病害变异不足，应记录不可识别的参数，不能人为生成系数。

管理处理和品种固定效应仅用于调整可观察的混杂结构，不代表完全消除因果混杂。

## 3.3 主要评价

报告：

- 条锈病严重度回归系数。
- 每增加10个病害评分百分点时的产量关联变化（dt/ha）。
- 95%置信区间。
- 调整前后方向及幅度。
- 模型使用样本数与环境数。
- 模型可识别性和残差诊断。

严重度评分必须确认其单位允许按10个百分点解释；如跨地点调查评分语义不同，应按有证据的可比子集报告，不能无条件合并解释。

## 3.4 推断不确定性

由于只有约13个地点×年份环境，环境聚类稳健标准误的大样本近似可能不可靠。

优先结合：

- 环境级效应方向。
- 环境级Bootstrap或合适的小样本聚类推断方法。
- Leave-One-Site-Year-Out关联敏感性。
- 固定效应估计与环境内部关联的一致性。

如果采用wild cluster bootstrap，应记录实际方法、重采样数和适用条件。

不得使用普通样地级Bootstrap制造过窄的置信区间。

## 3.5 自然感染敏感性

在5,899条自然感染记录内重复主要关联分析。

比较完整队列与自然感染子集的系数方向、变化幅度和环境覆盖。

两者样本构成不同，不能将系数差异解释为人工接种的因果作用。

## 3.6 研究结论边界

允许：

“Stripe rust severity was negatively associated with grain yield after adjustment for observed environmental, management, and cultivar factors.”

前提是计算结果确实支持。

不允许：

“Stripe rust caused a specific yield loss.”

也不能用本研究的回归系数直接推导喷药经济收益。

## 3.7 输出

- `04_adjusted_association_models.csv`
- `05_within_environment_effects.csv`
- `06_cultivar_management_sensitivity.csv`
- `02_Adjusted_Association_Findings.md`

---

# 4. Work Package 3 — Matched Natural-Domain Transfer Test

## 4.1 研究目的

前期发现：

完整条锈病队列有一定汇总预测增量，自然感染子集增量极小。

但是两个队列的环境构成不同，不能直接比较。

因此，本阶段要严格固定测试域，检查不同训练数据构成是否影响病害信息的预测增量。

## 4.2 测试集定义

测试域限定为：

有可靠证据属于自然感染条件的条锈病样地。

保持原始LOSYO环境分组。

对每个留出的Site-year：

- 只在该环境内预先确认的自然感染样地上评价。
- 同一个测试环境中的M0和M1必须使用完全相同的测试样本。
- 两种训练策略必须使用相同的自然感染测试样本。

如果个别测试环境不能合法评价，记录UNSCORABLE。

## 4.3 两种训练策略

**Strategy A — Natural-only training**

只使用除测试环境之外、已确认自然感染的训练样本。

**Strategy B — Full eligible training**

使用除测试环境之外的全部合格条锈病样本，包括自然感染和人工接种制度相关样本。

两种策略均不得包含任何来自留出Site-year的训练记录。

特别注意：

Strategy B比Strategy A通常有更多训练样本，而且训练分布不同。

因此两种策略的差异只能解释为训练数据构成策略的敏感性，不能归因于人工接种本身。

## 4.4 配对实验

优先开展：

- CatBoost M0/M1
- Ridge M0/M1作为敏感性验证

保持Phase 2B-Pilot的配置及特征定义。

不重新大规模调参。

使用完全相同的自然感染测试样本，对每种策略计算：

ΔMAE = MAE_without_disease − MAE_with_disease

同时报告：

- 环境等权ΔMAE。
- Pooled ΔMAE。
- 每个自然感染测试环境的ΔMAE。
- 两种训练策略之间的配对增量差异。
- 训练样本数量与环境覆盖。
- 训练集中人工接种相关样本比例。
- 测试环境删除敏感性。

如已有OOF预测恰好满足同域比较条件，可直接复用；否则仅重训练必要的折。

## 4.5 计算预算

优先限于CatBoost和Ridge的F1特征组合。

若按10个测试环境、2种训练策略、2个病害特征方案、2种算法计算，预期约80次拟合。

实际数量以可评价环境及已有模型复用情况为准。

不开展新的Bayesian Optimization。

不加入神经网络、XGBoost扩展或新数据源。

若有限计算结束后结论仍不确定，保留不确定结论，不再追加大规模重训练。

## 4.6 输出

- `07_matched_test_domain_manifest.csv`
- `08_matched_training_strategy_metrics.csv`
- `09_environment_level_strategy_comparison.csv`
- `03_Natural_Domain_Transfer_Findings.md`

---

# 5. Work Package 4 — Scientific Novelty Assessment

## 5.1 研究目的

判断这项工作能否在已有农业科学和机器学习文献基础上提出真正独立的贡献。

不能仅仅因为采用了新模型或新数据处理方式，就认定具有创新性。

需要同时核查三类研究：

A. 小麦条锈病严重度、抗性与产量关系。

B. 多环境小麦产量机器学习预测及跨域泛化。

C. 作物病害特征的增量预测信息价值及其环境稳定性。

## 5.2 必须核查的文献

**P1. Zetzsche et al. (2020)**

*Breeding progress for pathogen resistance is a second major driver for yield increase in German winter wheat at contrasting N levels.*

Scientific Reports.

DOI: 10.1038/s41598-020-77200-0

重点检查Figure 2、相关统计模型、试验数据与BRIWECS的重叠关系。

**P2. Wang et al. (2025)**

*Multi-environment field trials for wheat yield, stability and breeding progress in Germany.*

Scientific Data.

DOI: 10.1038/s41597-024-04332-7

重点核查原始数据论文已有的病害–产量及G×E×M分析结论。

**P3. Wang et al. (2026)**

*Plasticity of source-sink dynamics contributes to wheat yield stability.*

Nature Communications.

DOI: 10.1038/s41467-026-72330-x

重点检查其复用BRIWECS的研究范围、病害相关性分析和已有生物学结论。

**P4. Maatougui et al. (2026)**

*Wheat yield prediction using machine learning, proximal sensing, nitrogen and climate data for Oklahoma.*

Precision Agriculture.

DOI: 10.1007/s11119-026-10462-2

重点对比跨环境验证、预测特征、农业管理应用场景和传统统计基线。

**P5. Bebie & Kyparissis (2026)**

*Durum Wheat Yield Prediction: A Machine Learning Framework Integrating Sentinel-2 Imagery and Meteorological Data.*

Remote Sensing.

DOI: 10.3390/rs18162694

重点比较跨季节泛化性能下降的既有证据，避免将泛化失败现象本身作为首次发现。

另外针对2023—2026年相关研究进行有边界的补充检索。

优先选择最直接相关的约6—12篇原始研究，不需要进行无止境的综述式检索。

检索可使用Google Scholar、Crossref、Web of Science（如可访问）、PubMed及出版商官网。

对于有直接创新性冲突的文献，应尽可能核查全文Methods和Results，而不仅是摘要。

## 5.3 Novelty Matrix

建立：

`10_claim_evidence_novelty_matrix.csv`

字段至少包含：

- Candidate claim
- Our supporting evidence
- Prior study
- DOI
- Dataset overlap
- Method overlap
- Scientific finding overlap
- What remains different
- Evidence strength
- Novelty status
- Recommended claim wording

重点核查三个候选创新点：

**Claim A**

条锈病与产量的负向关联具有跨环境重复性。

预期：既有相关研究较多，独立创新性可能较弱。

**Claim B**

病害严重度与产量存在可重复关联，但病害特征的跨环境预测增量并不稳定。

核查是否有直接研究过这一“关联–预测脱钩”问题。

**Claim C**

通过环境平均偏差与环境内部样地差异分解，能够识别病害特征预测增量受限的具体层面。

核查相关分析是否提供了超出已有研究的实质性贡献。

不得预设Claim B或Claim C一定成立。

## 5.4 论文贡献等级

对每个Claim给出：

- STRONG
- MODERATE
- WEAK
- NOT NOVEL
- UNSUPPORTED

并解释依据。

如果只有“首次在BRIWECS上用CatBoost测试某种特征”这类方法应用差异，不应视为强创新。

---

# 6. Work Package 5 — Integration & Scientific Decision

整合前三项统计结果与文献核查，回答：

1. 条锈病与产量的负向关联在控制品种与管理后是否仍存在？
2. 关联结果在自然感染条件下是否保持相同方向？
3. 加入病害信息的收益主要来自何种误差组成？
4. 环境平均偏差能否解释跨环境泛化下降的部分现象？
5. 相同自然感染测试域下的结果是否改变此前的判断？
6. 人工接种相关训练样本的作用能否被可靠限定？
7. 哪些结果与2020、2025和2026年的已有论文实质性重叠？
8. 是否形成真正有区别的研究问题？
9. 论文现有证据适合何种学科定位？
10. 是否应继续投入科研计算与论文生产？

**不能以“结果复杂”“需要进一步提高模型性能”作为自动延续项目的理由。**

---

# 7. 最终裁定规则

## Manuscript GO

应同时满足：

- 数据与预测评估不存在阻断性错误。
- 在控制主要可观察混杂结构后，条锈病–产量关联仍具有明确且可解释的证据。
- 误差分解或共同测试域分析提供了实质性的新信息，而不仅是重复此前的性能表。
- 主要论文结论不依赖未经承认的因果假设。
- 文献对照显示至少一项可清晰表述、受到结果支持的独立贡献。

满足上述要求，批准精简版Phase 2C，完成有限的最终解释分析与论文整合。

## Manuscript CONDITIONAL GO

如果科学现象真实，但：

- 关联稳健性有限；
- 误差分解不能充分说明预测失效；
- 自然感染测试域仍存在明显不确定性；
- 或创新性仅达到方法学边界研究水平；

则判断是否值得转向较窄的SCI文章。

必须明确需要删减的研究主张和后续工作的实际范围。

不批准再次进行广泛模型探索。

## Manuscript NO-GO

如果：

- 控制品种和环境后主要关联不再稳定；
- 病害预测增量缺乏可靠性；
- 误差分解不能提供新的可解释证据；
- 且主要结论与现有发表研究高度重复；

则停止当前论文方向的新增计算。

保留完整研究成果、代码和统计报告，建议将资源转向其他农业机器学习问题。

不得强行撰写一篇仅由模型比较和SHAP图构成的文章。

---

# 8. 最终交付文件

所有新文件保存在：

`Phase2B_Closure/`

## 计算结果

1. `01_environment_error_decomposition.csv`
2. `02_within_environment_ranking.csv`
3. `03_disease_delta_error_components.csv`
4. `04_adjusted_association_models.csv`
5. `05_within_environment_effects.csv`
6. `06_cultivar_management_sensitivity.csv`
7. `07_matched_test_domain_manifest.csv`
8. `08_matched_training_strategy_metrics.csv`
9. `09_environment_level_strategy_comparison.csv`
10. `10_claim_evidence_novelty_matrix.csv`

## 科学报告

- `reports/01_Error_Decomposition_Findings.md`
- `reports/02_Adjusted_Association_Findings.md`
- `reports/03_Natural_Domain_Transfer_Findings.md`
- `reports/04_Literature_Novelty_Assessment.md`
- `reports/BRIWECS_Phase2B_Closure_Scientific_Integration.md`
- `reports/BRIWECS_Phase2B_Closure_GO_NO_GO.md`

## 图表

只生成四类必要的科学诊断图：

**C1.** 条锈病加入前后的环境平均误差与环境内部误差分解。

**C2.** 调整前后条锈病–产量关联估计及不确定性。

**C3.** 相同自然感染测试域下两种训练策略的病害预测增量。

**C4.** 关联强度与跨环境预测增量之间的环境级对比，展示二者是否一致。

图表采用简洁学术风格，输出PNG/PDF。

如果某项分析结果不支持有效作图，明确说明，不制造示意性结果图。

## 复现交付

- README
- 分析脚本
- 配置
- 实际重训练日志
- 模型拟合状态
- 版本与SHA256清单
- 输出文件清单

---

# 9. 任务停止规则

本阶段明确禁止：

1. 再次清洗整个BRIWECS数据库。
2. 扩大到全部六种病害。
3. 大规模机器学习超参数搜索。
4. 引入新的深度学习网络。
5. 因为某个环境结果不理想就删除样本。
6. 在测试结果基础上挑选最有利的病害调查时间点。
7. 进行超出原始实验设计支持的因果推断。
8. 使用测试环境产量进行模型校准。
9. 将不确定的人工接种状态归为自然感染。
10. 为形成论文创新而追加没有明确科学必要性的分析。
11. 未经批准直接进入Phase 2C或正式论文写作。

如有技术问题，优先进行一次合理修复。

非阻断性问题记录在限制清单中，不应导致无限期延迟科学裁定。

---

# 10. Codex最终中文汇报格式

完成实际计算后，依次报告：

**一、误差分解**

条锈病特征对环境平均预测和环境内部差异预测分别产生什么影响？

**二、病害关联**

调整环境、管理、品种之后，条锈病–产量关联的方向、系数及置信区间是什么？

**三、自然感染测试域**

两种训练策略在完全相同测试样本上的ΔMAE有何变化？

**四、创新性**

与已有文献相比，哪些贡献真正成立，哪些已经被前人完成？

**五、最合适的论文故事**

提供一个推荐英文题目、三项核心研究贡献、三项主要局限。

不得夸大预测精度、因果作用或实际农药管理能力。

**六、投稿潜力**

分别判断：

- 农业机器学习方向
- 植物病理学方向
- 农学与作物管理方向

指出最适合的定位及其证据基础。

**七、正式裁定**

Manuscript GO / CONDITIONAL GO / NO-GO。

必须给出清楚的理由，而不是笼统地说“需要继续完善”。

**八、交付**

列出实际文件路径，未完成内容标记PARTIAL或BLOCKED。

完成后停止执行。

**本轮的最终目的不是继续优化一个模型，而是决定这篇论文究竟值不值得写。**

现在开始Phase 2B-Closure的实际执行。