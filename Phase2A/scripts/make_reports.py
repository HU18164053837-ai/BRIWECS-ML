from common import *
def make_reports():
 metric=pd.read_csv(ROOT/'results/metrics/model_performance_summary.csv');env=pd.read_csv(ROOT/'results/metrics/environment_level_metrics.csv');paired=pd.read_csv(ROOT/'results/metrics/model_comparison_paired.csv');status=pd.read_csv(ROOT/'results/metrics/model_fit_status.csv');disease=pd.read_csv(ROOT/'results/metrics/disease_track_prescreen.csv');frozen=json.loads((ROOT/'config/cohort_freeze.json').read_text());c=pd.read_parquet(ROOT/'data/cohort/phaseI_yield_cohort.parquet');diag=pd.read_csv(ROOT/'results/metrics/error_stratification.csv')
 primary=metric[metric.validation_strategy.eq('LOSYO')];macro=primary[primary.aggregation.eq('ENVIRONMENT_MACRO')].sort_values('MAE');best=macro.iloc[0];bestname=f'{best.model_name} {best.feature_set}'
 bestpool=primary[primary.aggregation.eq('POOLED')&primary.model_name.eq(best.model_name)&primary.feature_set.eq(best.feature_set)].iloc[0]
 null=primary[primary.aggregation.eq('ENVIRONMENT_MACRO')&primary.model_name.eq('TrainingMean')].iloc[0]
 be=env[env.validation_strategy.eq('LOSYO')&env.model_name.eq(best.model_name)&env.feature_set.eq(best.feature_set)].sort_values('MAE',ascending=False)
 pairprim=paired[paired.validation_strategy.eq('LOSYO')];best_vs_null=paired[paired.validation_strategy.eq('LOSYO')&paired.model_a.eq(best.model_name)&paired.feature_a.eq(best.feature_set)&paired.model_b.eq('TrainingMean')]
 fail=status[status.status.eq('FAILED')];warn=status[status.status.eq('SUCCESS_WITH_WARNINGS')]
 mainmetric=metric[metric.validation_strategy.isin(['LOSYO','LOYO','LOLO','LOEO'])&metric.aggregation.eq('POOLED')][['validation_strategy','feature_set','model_name','n','MAE','RMSE','R2','MBE']]
 macmetric=metric[metric.validation_strategy.isin(['LOSYO','LOYO','LOLO','LOEO'])&metric.aggregation.eq('ENVIRONMENT_MACRO')][['validation_strategy','feature_set','model_name','environments','MAE','RMSE','R2','MBE','MAE_env_bootstrap_lo','MAE_env_bootstrap_hi']]
 sensitivity=metric[metric.validation_strategy.eq('LOSYO_NO_QLB')][['aggregation','feature_set','model_name','n','environments','MAE','RMSE','R2','MBE']]
 catweather=metric[metric.model_name.eq('CatBoost')&metric.feature_set.eq('F2')&metric.aggregation.eq('POOLED')].set_index('validation_strategy')
 wd=pd.read_csv(ROOT/'results/metrics/environment_weather_diagnostics.csv');weatherbest=wd[wd.model_name.eq(best.model_name)&wd.feature_set.eq(best.feature_set)]
 model_report=f'''# BRIWECS Phase 2A Modeling Report

状态：实际计算COMPLETED；个别模型失败以技术原因保留，见model_fit_status。总体裁定 **CONDITIONAL GO**。此报告不等于模型在所有未见环境均有稳定优势。

## 1. 冻结队列与预测任务

官方Figshare v1；{frozen['samples']:,}条独立可追溯样地记录，220品种，5地点，2015–2017三年，15地点×年份、51地点×年份×管理。单位dt/ha、100%干物质整样地产量。并非论文假设的45环境，也没有加入Phase II。全部所选记录有完整坐标。142条Phase1争议记录仍排除；本范围另有308条缺失/无效产量排除，3,040条RHH及6,202条范围外年份记录按互斥台账记录。343条IQR极端值只标记，未删除。队列及配置在训练前冻结，哈希见config/cohort_freeze.json。

任务是已知品种向新环境迁移，不是未见品种预测。所有测试折的品种覆盖见fold_integrity_check；原始地点、年份、坐标和区组ID均不是F1/F2输入。

## 2. 特征时间和模型

F0仅训练均值。F1为品种ID、release year（1,433行缺失，训练折内中位数填补）、计划氮/杀菌剂/水分处理。GGE 2016/2017水分设UNKNOWN_V1，保留原v1管理标签用于分组；绝不导入更新Github标签。未加入土壤等来源/测量时点不一致的变量。

F2只增加当年1月1日–3月31日的平均温度、温度SD、累计降水和累计辐射。所有15环境均有晚于截止日的真实抽穗记录（最早123 DOY，窗口截止90/91），提供环境季节级时间排序证据；未臆造播种/收获日，也未使用整季或收获后天气。F2标记RETROSPECTIVE_ENVIRONMENTAL_MODEL，属于固定春季截止的历史天气产量估计，不能称播前预测或早期病害预报。精确每样地收获日仍不完整。该窗口仅表示已知的冬季/早春环境，不能解释成整生长季干旱/热害。

所有中位数、尺度和编码在外层训练折内拟合；无全局目标编码、BLUE、特征筛选或测试集early stopping。250轮CatBoost和XGBoost、固定Ridge alpha=10、CPU、seed=42，无调参。训练配置、各折日志、主LOSYO模型权重均保存，未进行全数据拟合后伪装OOF。

MixedLM为较低维结构基线：release year和管理（F2加四个天气变量）作固定效应，site-year及其内management×replicate作随机截距；它不具有Ridge/CatBoost的品种ID固定效应信息，不能把其差异纯粹归结算法优劣。所有测试均只使用固定效应边际预测，未知环境随机效应为零；LOEO也使用同一规则，不估测试产量的随机效应。

## 3. 验证方案

LOSYO 15折为主，地点×年份的全部处理同时留出；LOYO 3折、LOLO 5折、LOEO 51折为补充。这里**LOEO专指地点×年份×管理留出**，与Phase1报告中名称用法不同，已明确映射。LOEO训练集仍见相同site-year其他处理，不能代替主要泛化结论。LOYO是跨年留出、可含比测试年更晚训练年份，不等于严格前向时间预测。

参数预先固定，inner_cv_manifest明确NOT APPLICABLE；没有虚构内层调参。所有训练均排除对应测试ID；各策略样本只测试一次。敏感性检查预定义排除QLB整地点，使用其余12个site-year的LOSYO，不能解释为独立消除接种的因果实验。

## 4. Pooled指标（完整列出）

{table(mainmetric)}
MBE=预测−观测，正值高估、负值低估。R²保留负值；低方差不可定义时保留UNDEFINED状态，不人为补值。若模型某折失败，n不足24,991，必须按不完整覆盖解释，不与全队列得分直接排名。

## 5. 环境等权指标

先在每site-year计算，再对环境等权平均；macro RMSE是环境RMSE的算术均值，不是pooled RMSE。macro R²是定义良好环境R²的均值，与pooled R²不同。

{table(macmetric)}
以主验证macro MAE排名，最低为**{bestname}**：macro MAE={best.MAE:.3f}，训练均值={null.MAE:.3f}；相应pooled MAE={bestpool.MAE:.3f}、RMSE={bestpool.RMSE:.3f}、R²={bestpool.R2:.3f}，但macro R²={best.R2:.3f}。多数环境的平均绝对误差较均值改善，不等于全部环境内预测有效。排名仅是这15环境上的描述性结果，不能用同一组测试环境排名后再声称选中模型已获独立确认。summary另提供TEST_GROUP_MACRO，分别对15/3/5/51个外层测试组等权；主环境macro始终以site-year为单位，避免把共享天气的管理组当独立环境。

## 6. 配对比较与不确定性

负delta_MAE表示A优于B。每项在完全相同sample_id+fold_id交集上比较；失败模型的缺失折不会伪造预测，原覆盖和匹配数量列出。

{table(pairprim[['comparison','matched_samples','environments','macro_delta_MAE','delta_MAE_CI95_lo','delta_MAE_CI95_hi','environments_improved','direction_consistency','exact_same_test_samples']])}

Bootstrap为2,000次site-year级重抽样，不把2万多样地当独立观测。15环境仍共享5地点/3年份、训练集相互重叠，因此区间为描述性的条件不确定性，不能当15个完全独立未来环境的严格置信保证。不以单一p值裁定；没有逐行bootstrap或正式SHAP。

## 7. 失败折与警告

共{len(status)}个预定模型×特征×外层任务，FAILED={len(fail)}、SUCCESS_WITH_WARNINGS={len(warn)}。失败没有从状态表删除，不通过增加搜索或看测试效果来重试。

{table(fail[['validation_strategy','fold_id','feature_set','model_name','error']]) if len(fail) else '无失败折。'}

混合模型不收敛、边界方差或Hessian警告保留在逐折JSON/错误日志；收敛且预测有限的警告折可用于预测描述，但不解释其参数p值或方差分量显著性。

## 8. 误差诊断

主验证最佳macro模型误差最大的环境如下，其余全部环境亦保留：

{table(be[['site_year','n','MAE','RMSE','R2','MBE']].head(5))}

{table(diag[diag.model_name.eq(best.model_name)&diag.feature_set.eq(best.feature_set)][['stratum_type','stratum','n','MAE','RMSE','R2','MBE']])}

高/低产区间按冻结队列四分位数事后描述，仅用于诊断，不反馈筛选或调参。CatBoost F2对低产Q1平均高估14.39、高产Q4平均低估18.35 dt/ha，显示明显向中间收缩。每环境实际与预测分布、极端产量、冬季天气范围和外推表现分别见error_stratification、extreme_yield_diagnostics与environment_weather_diagnostics；不能从单一小环境误差推出病害或气候机制。最小环境亦未删除。

冬季窗口的气象外推诊断：每个测试环境分别与其余14环境的原始天气范围比较，不使用测试产量定义气象极端。

{table(weatherbest[['site_year','weather_outside_training_range','outside_variables','MAE','MBE']])}

GGE 2016冬季均温4.58°C为本队列最高，对应CatBoost F2 MAE27.89；但最冷KIE 2016（2.78°C）的MAE只有8.93。最大失败KAL 2017（MAE36.45）不在上述四变量训练范围之外。因此失败不只集中在这些冬季气候极端，不能归因于某个未经检验的气候机制，也不能用前三个月代替整季胁迫。

## 9. QLB敏感性

{table(sensitivity)}
排除QLB改变训练及测试环境组成，不与原全队列pooled R²直接做因果差值。重点看相同12环境内各模型相对各自训练均值/Ridge的配对增益，详model_comparison_paired中LOSYO_NO_QLB。

## 10. 图与原始预测

Figure 2A-1覆盖热图；2A-2全部主验证环境误差；2A-3所有合法LOSYO预测的密度散点；2A-4所有site-year偏差；2A-5同一测试样本F2−F1增量。PNG/PDF均保存，未挑选好折。主要产物out_of_fold_predictions.parquet逐条含样本、真值、预测、模型、特征组、验证策略、fold、site-year、地点、年份、管理。

Phase2B/2C和正式论文撰写均未启动。
'''
 (ROOT/'reports/BRIWECS_Phase2A_Modeling_Report.md').write_text(model_report,encoding='utf-8')
 assessment=f'''# BRIWECS Phase 2A Scientific Assessment

## Overall decision — CONDITIONAL GO

合法队列和严格验证成立，存在一定跨环境预测信息，但不足以把这套基线包装为普遍稳定、显著优于统计模型的机器学习成果。后续只值得进行有界的、同队列同验证折的病害增量分析；不应无限调参追逐测试分数，也不能承诺SCI接收。

| Readiness | Decision | Evidence |
| --- | --- | --- |
| Data readiness | GO | 24,991真实样地、220品种、15环境；v1哈希、排除台账、稳定ID和冻结配置 |
| Modeling readiness | GO WITH DOCUMENTED FAILURES | 五类模型实际尝试；固定配置、折内预处理；失败原因与预测覆盖完整保存 |
| Cross-environment generalization readiness | CONDITIONAL GO | 主LOSYO可评价，但仅5地点×3年，改进与环境/模型有关，F2并非一致提升 |
| Disease track readiness | CONDITIONAL GO | 条锈病/Septoria可作配对增量研究；接种、评分定义与缺失域限制仍在 |
| Overall Phase 2A | CONDITIONAL GO | 保留严格范围，等待独立科学评审，不自动进入2B |

## G1–G4验收

G1 PASS：来源可信、目标清晰、样本ID唯一、无人工复制。G2 PASS WITH SCOPE CONDITIONS：排除产量派生/成熟性状，训练折内编码，LOSYO不拆处理；F2固定冬季窗口的时间证据在site-year层面，不作播前/早期预报结论，管理标签不确定性有掩码。G3 PASS：规定基线和两种树模型均实际运行，任何失败保留可复现技术原因。G4 EVIDENCE AVAILABLE / STABILITY LIMITED：已获得完整环境留出证据，不能把“可以评价”误写成“具有稳定优势”。

## 核心量化证据

主验证macro MAE最低为{bestname}，{best.MAE:.3f} dt/ha；训练均值{null.MAE:.3f} dt/ha。最佳模型对应pooled R²={bestpool.R2:.3f}。应同时查阅各环境R²、偏差和配对区间，不能凭该总体R²断言大多数环境均能准确预测。

{table(best_vs_null[['comparison','macro_delta_MAE','delta_MAE_CI95_lo','delta_MAE_CI95_hi','environments_improved','environments']]) if len(best_vs_null) else '最佳模型与均值的配对差异见完整结果表。'}

全部树模型对Ridge的配对结果及F2增量均在建模报告列出。CatBoost F2相对均值在12/15环境改善，但delta MAE的描述性区间跨零；F1 CatBoost与Ridge的macro差异只有+0.003 dt/ha，不能宣称树模型有稳定算法优势。CatBoost的F2−F1 macro MAE差异约−0.499，区间也跨零；XGBoost加入天气反而增加误差。

同一CatBoost F2，LOEO pooled R²={catweather.loc['LOEO','R2']:.3f}，LOSYO仅{catweather.loc['LOSYO','R2']:.3f}，LOLO为{catweather.loc['LOLO','R2']:.3f}。补充LOEO仍可从其他处理见到相同site-year天气向量，这个差异直接说明验证边界会显著改变结论，不能把LOEO的较高分数当作真正新环境成功。

## 可以支持怎样的论文？

可以支持以“严格环境留出揭示品种与管理模型的外推边界”为核心的研究基础，或在经过独立评审后开展病害信息是否增加跨环境产量估计价值的限定实验。当前仅为基线，不足以单独证明新方法优越、新机制、管理因果效应或外部有效性。SCI可发表性还需新颖性、独立复核和后续有限增量证据，不由本阶段最高R²裁定。

## 主要局限与下一步边界

三个年份不能充分代表未来气候；LOYO非前向预测；少数地点的管理处理与地点绑定。F2只有冬季窗口而无整季预测，15个气候向量无法提供广泛气候外推保证。品种ID编码适用于已知品种，不能称未见材料预测。GGE水分标签争议、QLB病害调查汇总语义、整数截断、病害观测选择偏倚仍需单独控制。无独立外部测试。

建议只允许一次预先注册的病害增量研究；始终在同一病害配对队列内重新拟合无病害/有病害两套模型，锁定LOSYO。不能把全队列基线的分数与较容易的病害子集分数直接相减。除QLB预定义敏感性外，本阶段没有依据误差反向删环境。

## 完成状态

COMPLETED：队列冻结、F1/F2、四种分组验证、五类模型、OOF预测、环境等权与配对区间、诊断图、QLB敏感性和病害预筛选。

PARTIAL：失败/警告的混合模型折、环境层面有限不确定性、精确逐样地收获时间与全季时序定义。不存在将失败估计伪装成功。

BLOCKED：声称稳定普遍外推、正式外部验证、疾病早期预报或管理因果解释。

NOT APPLICABLE：大规模搜索、深度学习、SHAP、Phase2B实际训练、湿实验、正式论文撰写。
'''
 (ROOT/'reports/BRIWECS_Phase2A_Scientific_Assessment.md').write_text(assessment,encoding='utf-8')
 disease_report=f'''# BRIWECS Phase 2B Readiness — 仅预筛选

建议优先 **Stripe_rust + Septoria**，CONDITIONAL GO；Powdery_mildew暂为备选。没有开展F3或病害预测模型训练，没有自动进入Phase2B。

{table(disease[['disease','paired_samples','site_years','cultivars','managements','raw_zero','raw_positive','raw_zero_proportion','raw_positive_proportion','artificial_samples','natural_samples','largest_environment_share','raw_positive_truncated_to_published_zero']])}

条锈病覆盖最广，但2,640条人工接种占{2640/8539:.1%}，来自QLB三年NF处理；自然记录5,899条。对其增量应分层/排除QLB敏感性，不能当成完全自然发生预测。Septoria有10个有效site-year、5,722条配对记录，较适合作第二病种；白粉病5,070条、9环境，覆盖更窄且正值占比高。完整环境和管理覆盖见两个prescreen表。

本阶段再次按原始source_file/source_row提取真实小数，输出disease_prescreen_raw_value_provenance.csv；表内raw_zero/raw_positive为原始数值口径，不把作者整数截断后的0当成真实未发病。published_supported_zero/positive另列以复核Phase1，缺失转0的可疑值仍排除。原始小数恢复只是预筛选证据，不是新训练特征。

未来必须在每个病害的相同配对样本内，重拟合F1/F2与加病害模型，训练和测试样本完全一致；复用已冻结LOSYO group定义，缺乏病害测试样本的环境明确UNSCORABLE，不补齐。不能直接拿全24,991样地产量模型与病害子集结果比较，避免选择偏倚假增益。

Phase1的简单病害邻接日期提取在QLB会碰到Inoculation Date字段，不能当病害调查日期；原xlsx实际还含1st–4th measure date。Phase2A未据此生成任何病害时间特征。未来如要用时点模型，需按原始表头区分接种和调查日期，而不是复用该邻接字段。AO/即时严重度的跨站定义仍待统一。

本建议用于独立科学评审；不授权实际进入Phase2B。
'''
 (ROOT/'reports/BRIWECS_Phase2A_Phase2B_Readiness.md').write_text(disease_report,encoding='utf-8')
 decisions=dict(Data_readiness='GO',Modeling_readiness='GO_WITH_DOCUMENTED_FAILURES',Cross_environment_generalization_readiness='CONDITIONAL GO',Disease_track_readiness='CONDITIONAL GO',Overall_Phase2A='CONDITIONAL GO',best_primary_macro_model=bestname,best_primary_macro_MAE=float(best.MAE),best_primary_pooled_RMSE=float(bestpool.RMSE),best_primary_pooled_R2=float(bestpool.R2),failed_fits=len(fail),warning_fits=len(warn),prescreen_recommendation=['Stripe_rust','Septoria'],phase2b_started=False)
 write_json(decisions,'reports/stage_decision.json')
 print(json.dumps(decisions,indent=2),flush=True)
if __name__=='__main__':make_reports()
