from common import *
def reports():
 p=ROOT/'results/metrics';co=pd.read_csv(ROOT/'data/disease_cohort_summary.csv');perf=pd.read_csv(p/'disease_incremental_performance.csv');u=pd.read_csv(p/'paired_uncertainty_analysis.csv');env=pd.read_csv(p/'environment_incremental_metrics.csv');a=pd.read_csv(p/'between_within_association_summary.csv');ae=pd.read_csv(p/'within_environment_association.csv');status=pd.read_csv(p/'model_fit_status.csv');drop=pd.read_csv(p/'drop_one_environment_sensitivity.csv');contrib=pd.read_csv(p/'inoculation_regime_contribution.csv');qc=pd.read_csv(ROOT/'data/disease_label_qc.csv');fold=pd.read_csv(ROOT/'splits/disease_outer_fold_manifest.csv')
 primary=perf[perf.validation_strategy.eq('LOSYO')];full=primary[primary.scope.eq('FULL')];nat=primary[primary.scope.eq('NATURAL')];pu=u[u.validation_strategy.eq('LOSYO')];puf=pu[pu.scope.eq('FULL')]
 gaincols=['disease','scope','model_name','comparison','aggregation','n','environments','MAE_without','MAE_with','delta_MAE','relative_MAE_gain_percent','RMSE_without','RMSE_with','delta_RMSE','R2_without','R2_with','delta_R2','MBE_without','MBE_with','improved_environments','worsened_environments']
 acols=['disease','scope','model_name','comparison','macro_delta_MAE','macro_CI95_lo','macro_CI95_hi','pooled_delta_MAE','pooled_CI95_lo','pooled_CI95_hi','improved_environments','environments','two_largest_positive_gain_share']
 mainpivot=env[env.validation_strategy.eq('LOSYO')&env.scope.eq('FULL')].assign(label=lambda d:d.model_name+' '+d.comparison).pivot(index=['disease','site_year','n'],columns='label',values='delta_MAE').reset_index()
 dropsummary=drop.groupby(['disease','scope','model_name','comparison']).agg(macro_sign_flips=('macro_sign_flip','sum'),pooled_sign_flips=('pooled_sign_flip','sum'),minimum_remaining_macro_delta=('macro_delta_MAE','min'),maximum_remaining_macro_delta=('macro_delta_MAE','max')).reset_index()
 secondary=perf[perf.validation_strategy.ne('LOSYO')&perf.aggregation.eq('ENVIRONMENT_MACRO')][['disease','scope','validation_strategy','model_name','comparison','outer_test_groups','delta_MAE','delta_RMSE','delta_R2','improved_environments','environments']]
 failed=status[status.status.eq('FAILED')];nfit=int(status.status.str.startswith('SUCCESS').sum());uns=int(status.status.eq('UNSCORABLE').sum())
 report=f'''# BRIWECS Phase 2B-Pilot Modeling Report

状态：COMPLETED。仅Ridge与CatBoost，M0–M3，官方Figshare v1；未调参、未扩模、未进入Phase2C。所有增量均来自相同配对队列中重训的对照，不使用Phase2A全队列成绩代替M0/M2。

## 1. 冻结队列与病害标签

{table(co)}
全队列沿用Phase2A五地点、2015–2017年合法样地；未重新选择品种、极端产量或误差较大的环境。sample_id/原始联合键均唯一。条锈病8,539，Septoria5,722，经原始source_file/source_row独立复核后与预筛选数相同，不是照抄预设数字。条锈病真实零453（5.31%）、正值8,086；Septoria真实零264（4.61%）、正值5,458。原始小数恢复分别纠正199和619个被作者截断为零的正值。缺失转0且无原始支持的值仍排除。

病害单位来自官方Unit.xlsx，为叶部感染面积百分比。原始单列样地评分是本次终点；QLB条锈病存在AO历史汇总语义，其他地点不保证同一调查时点。没有选择最大值、最相关日期或补造重复调查。measurement_calendar_evidence.csv记录98项原管理日历证据，区分Inoculation Date与1st–4th measure date；这些处理级日历不等于每样地标签的完整时间戳。调查部位和终点限制见data/disease_data_dictionary.csv。

条锈病自然感染子集5,899、10环境，仅纳入明确NATURAL。UNKNOWN/MIXED不自动转自然，本批合格样本中这两类为零。Septoria全部合格值属于该病害特异性自然感染，所以其NATURAL与FULL完全相同，复用同一拟合而不重复训练；QLB部分样地的其他病原人工接种背景另标记，不能解读为整个试验没有人工接种。

## 2. 预测任务、特征与严格配对

M0=F1（品种ID、release year、氮/杀菌剂/水分管理）；M1=M0+原始单病害严重度；M2=F2（F1加冬季天气）；M3=M2+同一病害严重度。除病害列外，成对特征完全相同；训练/测试ID、目标与外层折相同。release year缺失填补、尺度和类别编码均只在训练折拟合。严重度不缺失，不生成新的病害插补值。

沿用Phase2A Ridge alpha=10、lsqr；CatBoost250轮、depth6、learning_rate0.05、l2=5、CPU、seed42，不使用外层测试eval_set、early stopping或新搜索。F2只使用目标年份1月1日–3月31日的平均温度、温度SD、累计降水、累计辐射，保持Phase2A矩阵不变。已有site-year抽穗证据晚于窗口，但病害调查日期未完全同步，因此所有M0–M3统一定位 **RETROSPECTIVE ASSOCIATION-ORIENTED ESTIMATION**，不称播种前、发病前或已证实的调查后实时预测。

## 3. 分组、实际拟合与完整性

复用Phase2A全部15个LOSYO环境身份；缺乏病害配对测试记录的环境保留UNSCORABLE，不补样本。LOYO按整个年份留出，LOLO按整个地点留出；同site-year全部处理同时离开主训练集。LOYO只有3个测试年，LOLO只有4–5地点，补充证据有限。

{table(fold.groupby(['disease','scope','validation_strategy','status']).size().rename('folds').reset_index())}
共有{len(status)}个计划任务状态：实际成功拟合{nfit}，UNSCORABLE={uns}，模型失败={len(failed)}。112个不可评价状态来自14个没有配对测试观测的预定组乘以8模型配置，不能称失败模型或因效果差删折。440个成功拟合构成全部OOF预测；主LOSYO模型及其预处理保存models。两CPU进程仅调度独立固定任务，不改变科学方案。

## 4. 主验证：绝对表现与增量

**ΔMAE=MAE_without−MAE_with，正值改善**；ΔRMSE同方向；ΔR²=R²_with−R²_without。MBE=预测−观测。不得与Phase2A部分表使用的相反差值符号混淆。

{table(full[gaincols])}

环境macro是各site-year指标的算术平均；pooled对样地加权，两者不同。所有负R²保留。KIE 2016的Septoria只有7条配对、产量变异有限：Ridge M2/M3的该环境R²分别约−332.99/−285.72，数学上有定义但极不稳定，造成宏ΔR²约+4.36；这绝不代表模型获得良好解释能力。该小环境不因结果异常删除，判断以MAE为主，同时公开分母风险。

## 5. 配对不确定性与环境方向

{table(puf[acols])}
区间采用2,000次配对site-year簇bootstrap，避免数千样地IID重抽样。只有10–13环境且共享地点、年份和训练集，属于有限环境下的描述性条件区间，不能视为大量独立环境的强保证。不做单一p值裁定。

全部主验证环境的ΔMAE（包括恶化与小环境）：

{table(mainpivot)}
条锈病CatBoost M1/M0在9/13环境改善，但HAN三年和KIE 2016均恶化；KIE 2016为−6.17，QLB 2015为+7.36 dt/ha。总体增益+1.581（约11.0%）与环境macro仅+0.434形成明显差别。加入天气后CatBoost条锈病增益缩至pooled+0.185、macro+0.078，且RMSE和pooled R²并不改善。

Septoria依赖算法和背景：CatBoost M1/M0略恶化，M3/M2改善约+0.394；Ridge M1/M0改善约+0.369，但M3/M2 pooled恶化−0.807。不能按最好一个组合宣布通用增益。

## 6. 自然感染敏感性与增量来源

{table(nat[gaincols])}

{table(pu[pu.scope.eq('NATURAL')][acols])}

条锈病自然子集CatBoost M1/M0 pooled仅+0.073（0.56%），macro+0.445；M3/M2 pooled+0.193、macro+0.181。Ridge M1/M0 macro+0.343的区间约[0.005,0.645]，但pooled区间跨零，且M3/M2自然子集转为恶化。因此有有限F1信号，但尚未证实独立于天气的稳定价值。

完整队列和自然子集的环境、训练样本与管理构成不同，不能将两者得分之差解释为接种因果效应。以下是在**同一完整队列模型的既有OOF预测内部**分解误差改善贡献，与自然子集重新训练结果分开：

{table(contrib)}
条锈病接种制度局限于QLB部分NF，与地点/管理混杂，标CONFOUNDED。人工子集仅描述三个QLB环境的分布和相关，不另做仅三环境的正式外推模型。

## 7. 删除单个环境后的稳定性（不重训）

{table(dropsummary)}
完整条锈病CatBoost M1/M0去掉QLB 2015后，macro增量由+0.434变为−0.143，说明环境平均结论可受单环境改变；pooled仍为正，不能混为同一结论。条锈病CatBoost M3/M2两个环境的删除可改变macro方向。最大两个正向贡献环境占比等见区间表，分母为正增益总和，不是净增益，以免正负抵消时产生误导比例。

## 8. 补充LOYO与LOLO

{table(secondary)}
这些验证并非新增大量独立环境。条锈病CatBoost M1/M0在LOLO的macro增量约−0.895，进一步限制新地点稳定性；不因主LOSYO有正增量而隐藏补充恶化。

## 9. 环境内/间关联及管理背景

{table(a)}
条锈病13环境中12个负相关，自然10环境中9个负相关；环境内中心化及进一步环境×管理内中心化的相关仍约−0.420/−0.445，说明关联不只由环境均值或粗管理均值差异构成。但品种抗性、潜在病害、测量时点和选择观测仍可混杂，不能推断病害因果减产。

Septoria pooled相关约+0.100、between约+0.285，而within-env约−0.021、within-env×management约−0.007，环境内信号接近零。其7/10环境负向不等于有一致强关系。每环境相关CI为500次品种簇bootstrap；考虑同品种多样地相关，但未完全排除品种亲缘或管理混杂。真实产量中心化仅用于此诊断，从未进入预测。

WF/NF与氮处理的严重度、产量、相关和OOF误差见management_association.csv与management_error_diagnostics.csv。不对覆盖不平衡的处理作农药推荐或因果效应估计。

## 10. 交付与阶段边界

B1–B5均提供PNG/PDF，全部由实际结果生成。逐样本预测、模型状态、逐折样本/目标分布、标签QC、原始时间/单位证据和全部配对检查保存。PARTIAL：病害终点的跨地点时间同步与AO语义统一；不是训练未执行。NOT ESTIMABLE/UNSCORABLE：只有预定测试组无配对记录。Phase2C、SHAP、正式论文与新湿实验均未启动。
'''
 (ROOT/'reports/BRIWECS_Phase2B_Pilot_Modeling_Report.md').write_text(report,encoding='utf-8')
 assessment=f'''# BRIWECS Phase 2B-Pilot Scientific Assessment

## 主要判断

**总体CONDITIONAL GO，研究范围须收缩。** 条锈病存在可重复的环境内负向关联和有限预测信息；“病害信息在未知环境普遍稳定地增加产量预测价值，且独立于管理/冬季天气”的强主张未获支持。Septoria不足以承担稳定增量预测的主病害路线。

## 核心证据

{table(puf[['disease','model_name','comparison','macro_delta_MAE','macro_CI95_lo','macro_CI95_hi','pooled_delta_MAE','improved_environments','environments']])}
全队列所有这些macro增量区间均跨零；这不是单凭不显著否定价值，而是结合方向异质性、pooled/macro分歧、新地点验证及自然子集衰减判断。没有更换参数、删坏环境或只选择最好组合。

## 关联与预测不能互相替代

条锈病12/13环境负相关，自然子集9/10；环境/管理内中心化相关仍为负，支持植物保护方向的关联问题。但自然感染下CatBoost F1 pooled增益仅0.073 dt/ha，加入天气后Ridge的病害增量反而为负。强相关不保证在不同时点/环境中的条件增量可迁移，更不证明因果。

Septoria汇总相关与环境内相关方向/强度不一致，within-env×management约−0.007；CatBoost F2下约0.394 dt/ha增益是有限、情景依赖信号，不能将其当作已解释的病害减产机制。其Ridge F2出现严重外推和小样本R²分母问题，所有预测保留。

## 人工接种和少数环境

完整条锈病CatBoost F1 pooled增益1.581，macro只有0.434；去掉QLB 2015即改变macro方向。自然子集重新拟合后仍有小增益，故不能说一切信息只来自人工接种；但数据不支持强自然场景普遍优势。接种与QLB/NF混杂，完整与自然模型差异不是人工接种效应。

## 信息时间与测量可比性

原始小数可靠恢复、缺失转0已排除；单列终点未事后挑选日期。QLB AO历史汇总与其他地点调查终点、调查日期和部位/程序仍未完全统一，可能影响跨环境迁移。F2只有Jan1–Mar31，不能代表整季病害气象，也不证明可在播前/发病前预报；全部结论仅限回顾性估计。

## 论文创新性评价

现有证据可支持“病害关联具有农业意义，但其预测增量受环境、接种背景与验证边界限制”的限定植物保护数据科学研究基础。尚不足以支撑以稳定高精度增益或一般化防治决策为主张的论文。条锈病比Septoria更有理由保留，但核心应是外推边界和条件价值，而非换算法制造更高分。

本Pilot没有开展系统性文献新颖性检索，不能断言首次发现或保证SCI发表。下一步若经独立科学评审决定继续，应仅针对条锈病已有信号预先制定有限稳健性目标，不自动扩展到全面Phase2C或重复调参。正式论文未撰写。

## 完成与限制

COMPLETED：两病独立冻结、原始小数QC、440次真实配对拟合、LOSYO/LOYO/LOLO、配对CI、环境内/间关联、自然子集、管理与删环境诊断、五类图及复现材料。

PARTIAL：完整逐样地同步调查时点、跨地点AO/评分程序统一、外部独立确认。不得把这些限制称作模型训练失败。

UNSCORABLE：原预定组缺病害配对测试观测（明确保留112个模型状态）。

NOT APPLICABLE：算法搜索、SHAP、农药因果推荐、Phase2C、正式论文或湿实验。
'''
 (ROOT/'reports/BRIWECS_Phase2B_Pilot_Scientific_Assessment.md').write_text(assessment,encoding='utf-8')
 decision='''# BRIWECS Phase 2B-Pilot GO / NO-GO

## Overall Decision — CONDITIONAL GO

本阶段完成，停止新增病害预测模型。不得把结果解读为已证实普遍稳定的病害增量，不直接全面进入Phase2C。

| Track / claim | Decision | Evidence |
| --- | --- | --- |
| Data and paired validation | GO | 官方v1、原始小数、唯一样地、匹配重训、无测试标签预处理、440次实际成功拟合 |
| Stripe rust retrospective association / limited F1 value | CONDITIONAL GO | 多环境负向关系重复；部分配对增益存在，但自然场景增益小、天气调整后依赖模型、少环境可改变结论 |
| Septoria as principal stable-incremental predictor | NO-GO | 环境/管理内相关接近零，算法/背景下增益方向改变；不宜继续为其扩模寻找正结果 |
| Broad stable increment independent of management/weather | NO-GO AS A CLAIM | 全队列macro区间跨零、LOLO和自然子集结果不稳定，接种/终点语义限制未消除 |
| Overall Pilot research direction | CONDITIONAL GO | 仅保留条锈病与外推边界的限定问题，等待独立评审，不保证论文发表 |

## Evidence Supporting Limited Continuation

条锈病完整队列CatBoost M1−M0的pooled ΔMAE为+1.581 dt/ha（约11.0%），但macro仅+0.434；自然子集F1模型有小的正向信号，环境内负向关联并非完全由环境均值差异造成。

## Critical Counter-Evidence

完整条锈病CatBoost F1去掉QLB 2015后macro转−0.143；自然子集CatBoost pooled只+0.073（0.56%）；加入天气后并无跨算法一致增量。Septoria在CatBoost F2下有约+0.394的情景信号，但其他设置不稳定。不能只引用最好组合。

## Recommended Scope

保留条锈病作为主要关联终点，明确自然/人工接种、终点AO与日期限制；暂停Septoria作为预测主线。若独立评审认为值得继续，只允许有界、预先冻结的稳健性检查，不能通过反复调整队列、模型或指标追求GO。不得推断病害减产因果、人工接种独立效应或推荐农药方案。

## Phase Boundary

Phase2B-Pilot结束。Phase2C未启动，未自动撰写正式论文。真实文件与SHA256交付清单见本目录reports。
'''
 (ROOT/'reports/BRIWECS_Phase2B_Pilot_GO_NO_GO.md').write_text(decision,encoding='utf-8')
 js(dict(overall='CONDITIONAL GO',stripe_rust='CONDITIONAL GO',septoria_primary_stable_increment='NO-GO',broad_stable_increment_claim='NO-GO',data_and_validation='GO',successful_fits=nfit,failed_fits=len(failed),unscorable_states=uns,phase2c_started=False,manuscript_started=False),'reports/stage_decision.json')
 print('Three scientific reports and stage decision written; no new models fitted.',flush=True)
if __name__=='__main__':reports()
