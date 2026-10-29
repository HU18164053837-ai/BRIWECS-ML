from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import platform, importlib.metadata as md

def read(n):return pd.read_csv(ROOT/'results/tables'/n)
def report(name,s):(ROOT/'reports'/name).write_text(s,encoding='utf-8')
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(ROOT/'results/figures'/f'{name}.{ext}',dpi=220,bbox_inches='tight')
 plt.close(fig)

def literature():
 sources=[
 ('P1','10.1038/s41598-020-77200-0','Zetzsche et al. (2020)','QLB 2015–17 overlapping trial; 178 cultivars, adjusted means; not independent replication','Figure 2; Materials and methods / statistical analysis','FULLTEXT_LOCAL','混合模型处理后品种均值的病害—产量关系已被分析；Figure 2明确包含条锈病。','原始地块级、真正留出site-year的病害增量及误差组成'),
 ('P2','10.1038/s41597-024-04332-7','Wang et al. (2025)','Same parent Figshare BRIWECS dataset; technical-validation subset 220 cultivars, 45 growing conditions','Methods; Technical Validation I–IV, especially TFI and stability regression','FULLTEXT_LOCAL','已研究G×E×M、总真菌感染TFI与育种进展/产量稳定性；不是独立新数据。','单病害同折增量与同自然测试域的训练构成对照'),
 ('P3','10.1038/s41467-026-72330-x','Wang et al. (2026)','BRIWECS1/2 reused with four additional physiological datasets; substantial parent-trial overlap','Results Fig1/2 and Supplementary Fig3; Methods Statistics; data availability','FULLTEXT_LOCAL','已综合61性状、抗病育种进展、源库与稳定性；不能再将抗病维持源能力当作新机制。','文中所核查正文未见本研究的条锈病配对LOSYO增量和误差分解；不是全补充材料排他性检索'),
 ('P4','10.1007/s11119-026-10462-2','Maatougui et al. (2026)','No BRIWECS overlap; Oklahoma 57 site-years','Materials and methods: GroupShuffleSplit and GroupKFold; Results','FULLTEXT_LOCAL','已经采用site-year分组测试、CatBoost及INSEY/OLS/LMM/GAM基线，限定施氮前信息。','本研究病害评分无同步日期，不能复制其施氮决策场景；关注病害增量而非算法优势'),
 ('P5','10.3390/rs18162694','Bebie & Kyparissis (2026)','No overlap; Greece durum wheat, 2018–25','Publisher indexed full-text sections 1.4,2.5.3,3.1/3.3,4.1','PUBLISHER_INDEXED_FULLTEXT; LOCAL_HTML/PDF_403','已有随机划分到LOYO性能下降、空间差异与年际气候基线分离；泛化失败和分层概念不是首次。','条锈病配对OOF的偏差平方/误差方差精确核算及固定自然测试域；不是新数学恒等式'),
 ('P6','10.1038/s41477-023-01516-8','Sabir et al. (2023)','Overlapping BRIWECS cultivar/trial lineage; climate and yield-component reanalysis','Publisher abstract, figure captions and data availability; subscription Methods unavailable','PARTIAL_FULLTEXT_ACCESS','220品种、81发育窗口的环境敏感性已研究，环境异质性与阶段效应已有先例。','本研究不选择病害时间窗口、不重复生理机制；受限全文不用于排他性新颖判定'),
 ('P7','10.1002/tpg2.20532','McBreen et al. (online 2024; issue 2025)','No overlap; Sungrains southeastern US, 11 environments','Fulltext XML sections 2.6,2.7 and Results; CV0-S/CV0-Y','FULLTEXT_EUROPE_PMC_XML','已区分新材料/新地点/新年份任务并以环境内部相关评价辅助表型信息。','条锈病单特征的增量、混杂调整与自然域匹配；不能把环境内排序评价当新方法')]
 registry=[];matrix=[]
 statuses={'A':'NOT NOVEL','B':'WEAK','C':'MODERATE'}
 support={'A':'Crude negative 12/13, natural 9/10; A3 full -1.001 [-2.820,0.186], natural -0.606 [-3.149,0.514] per10 nominal score units; robust adjusted claim unsupported',
 'B':'Association attenuates after cultivar FE; natural CB macro increment 0.445 [-0.321,1.169]; no stable incremental advantage. Strong independent adjusted association premise not established',
 'C':'FULL CB macro deltaMSE +27.923 = bias-square +45.227 + centered variance -17.304; plot/cultivar ranks improve only5/13 and4/13. Empirical scope-specific explanation, not new decomposition identity'}
 wording={'A':'Raw within-environment negative associations recur, but cultivar-adjusted evidence is uncertain; this is a replication, not a discovery.',
 'B':'Observed disease-yield associations do not establish a stable incremental cross-environment prediction benefit in these trials; an independently robust association is not claimed.',
 'C':'Paired OOF diagnostics localize heterogeneous disease-feature gains to mean bias and within-environment error components; ranking gains are inconsistent.'}
 for ident,doi,author,overlap,sections,access,finding,difference in sources:
  fn=ROOT/'data/literature'/f'{ident}_crossref.json'
  if not fn.exists():fn=ROOT/'data/literature'/f'{ident}_crossref_retry.json'
  message=json.loads(fn.read_text(encoding='utf-8'))['message'];title=message['title'][0]
  registry.append(dict(id=ident,study=author,title=title,doi=doi,url='https://doi.org/'+doi,metadata_verified='Crossref + primary publisher/PMC',dataset_overlap=overlap,inspected_sections=sections,access_level=access,finding=finding,remaining_difference=difference))
  for claim in ['A','B','C']:
   matrix.append({'Candidate claim':claim,'Our supporting evidence':support[claim],'Prior study':author,'DOI':doi,'Dataset overlap':overlap,'Method overlap':sections,'Scientific finding overlap':finding,'What remains different':difference,'Evidence strength':'Descriptive/limited; adjusted robustness not established' if claim!='C' else 'Exact algebra and matched OOF verified; limited external scope','Novelty status':statuses[claim],'Recommended claim wording':wording[claim],'Access limitation':access})
 csv(registry,'literature_evidence_registry.csv');csv(matrix,'10_claim_evidence_novelty_matrix.csv')
 report('04_Literature_Novelty_Assessment.md','''# Bounded novelty audit — 26 October 2026

结论：A = NOT NOVEL；B = WEAK（原始“独立稳健关联—预测脱钩”强版本为UNSUPPORTED）；C = MODERATE，仅限这套数据的实证误差定位。没有STRONG贡献，不声称首次。

检索边界为7篇原始研究：指定P1–P5，加2023年同数据谱系P6和2024在线/2025卷期的多环境辅助表型预测P7。使用nature-academic-search技能；专用academic-search MCP未挂载，改用公共Crossref、出版商、PubMed/Europe PMC及web索引。未使用Web of Science、未声称系统综述或独立他引认证。日期和DOI均以Crossref记录与出版商/PMC互核，P2 Crossref一次重试成功。查询与访问限制见logs/literature_search_log.md，原始响应及SHA保留。

'''+ '\n\n'.join(f"## {r['id']} — {r['study']}\n\n[{r['title']}]({r['url']})。\n\n{r['finding']}\n\n数据关系：{r['dataset_overlap']}。查阅范围：{r['inspected_sections']}；证据层级：{r['access_level']}。区别：{r['remaining_difference']}。" for r in registry)+'''

## 创新主张裁定

A不能重复宣称“条锈病负关联的新发现”。P1 Figure 2使用处理特异、跨年重复调整均值；其基因型随机效应与本研究A3品种固定效应的条件斜率并非同一估计目标，不能因本研究区间跨零就推翻原研究。

B仅支持一项受限警示：描述性关联不能代替增量预测检验。品种调整后的全队列及自然域区间均跨零，不能宣称已经证明独立关联稳健、同时预测信息失效；也不能由不显著推出精确零效应。所核查7篇中未找到完全同设计的条锈病配对LOSYO实验，但有限检索无法证明世界首次。

C的数学恒等式没有新颖性，尺度分离、组间泛化与环境内部相关评价也有先例。保留的中等贡献是具体审计证据：原始总体收益如何依赖均值偏差、环境权重和人工相关训练组成，且排序并未普遍改善。C是已有工具产生的限定实证新信息，不是新算法或生物机制。

完整21行claim×paper对照见10_claim_evidence_novelty_matrix.csv。P5原站HTML/PDF访问受限，Methods/Results来自出版商索引全文片段；P6订阅正文受限。这两项是明确的证据限制，因此禁止强排他性创新表述。没有因访问受限追加无止境检索。
''')

def main():
 d=read('03_disease_delta_error_components.csv');r=read('02_within_environment_ranking.csv');a=read('04_adjusted_association_models.csv');w=read('05_within_environment_effects.csv');s=read('06_cultivar_management_sensitivity.csv');m=read('08_matched_training_strategy_metrics.csv');e=read('matched_strategy_environment_metrics.csv');co=pd.read_parquet(ROOT/'data/stripe_rust_inherited.parquet')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 # Summaries explicitly distinguish averaging weights and model scope.
 ag=[]
 for keys,g in d.groupby(['scope','model_name','comparison']):
  for weighting in ['ENVIRONMENT_EQUAL','PLOT_WEIGHTED']:
   wt=np.ones(len(g)) if weighting=='ENVIRONMENT_EQUAL' else g.n
   ag.append(dict(zip(['scope','model_name','comparison'],keys),weighting=weighting,environments=len(g),**{k:np.average(g[k],weights=wt) for k in ['delta_MSE','delta_bias_squared','delta_centered_variance','delta_MAE']},bias_improved=int(g.delta_bias_squared.gt(0).sum()),centered_variance_improved=int(g.delta_centered_variance.gt(0).sum())))
 ag=csv(ag,'error_component_aggregation.csv')
 rank=r[r.context.eq('M0')].merge(r[r.context.eq('M1')],on=['scope','model_name','site_year','level'],suffixes=('_M0','_M1'));rank['delta_spearman']=rank.spearman_M1-rank.spearman_M0;csv(rank,'paired_ranking_changes.csv')
 fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
 for ax,scope in zip(axs,['FULL','NATURAL']):
  g=d[d.scope.eq(scope)&d.model_name.eq('CatBoost')&d.comparison.eq('M1_vs_M0')];x=np.arange(len(g));ax.bar(x-.18,g.delta_bias_squared,.36,label='Bias squared reduction',color='#39749e');ax.bar(x+.18,g.delta_centered_variance,.36,label='Centered variance reduction',color='#d88a44');ax.axhline(0,color='black',lw=.7);ax.set_xticks(x,g.site_year,rotation=75);ax.set_title(scope+' training / test cohort');ax.set_ylabel('M0 minus M1 MSE component (dt/ha)$^2$')
 axs[0].legend(fontsize=8);fig.suptitle('C1 | CatBoost LOSYO disease increment: exact error decomposition');save(fig,'C1_error_decomposition')
 fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for ax,scope in zip(axs,['FULL','NATURAL']):
  g=a[a.scope.eq(scope)];v=g.association_per_10_score_units.to_numpy();ax.errorbar(v,np.arange(4),xerr=np.vstack([v-g.per10_CI95_lo,g.per10_CI95_hi-v]),fmt='o',capsize=3,color='#39749e');ax.axvline(0,color='gray',ls='--');ax.set_yticks(range(4),['A0 crude','A1 + site-year','A2 + environment-management','A3 + cultivar']);ax.invert_yaxis();ax.set_title(scope+f' ({int(g.environments.iloc[0])} site-years)');ax.set_xlabel('dt/ha per 10 nominal score units\n499 environment bootstrap resamples\n95% descriptive interval')
 fig.suptitle('C2 | Association attenuates after cultivar adjustment; endpoints differ by site');save(fig,'C2_adjusted_association')
 fig,axs=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
 for ax,model in zip(axs,['CatBoost','Ridge']):
  g=e[e.model_name.eq(model)];pv=g.pivot(index='site_year',columns='strategy',values='delta_MAE');x=np.arange(len(pv));ax.plot(x,pv.A_NATURAL_ONLY,'o-',label='A natural-only training',color='#39749e');ax.plot(x,pv.B_FULL_ELIGIBLE,'s--',label='B full eligible training',color='#d88a44');ax.axhline(0,color='gray',lw=.7);ax.set_xticks(x,pv.index,rotation=75);ax.set_ylabel('Disease increment: MAE(M0) - MAE(M1), dt/ha');z=m[m.model_name.eq(model)&m.strategy.eq('B_MINUS_A_INCREMENT')].iloc[0];ax.set_title(f'{model}: macro B-A {z.macro_delta_MAE:.2f}\n95% [{z.macro_CI95_lo:.2f}, {z.macro_CI95_hi:.2f}]')
 axs[0].legend(fontsize=8);fig.suptitle('C3 | Identical 5,899 natural test plots; positive increment means improvement');save(fig,'C3_matched_natural_domain')
 assoc=[]
 for scope,c in [('FULL',co),('NATURAL',co[co.inoculation_status.eq('NATURAL')])]:
  for sy,g in c.groupby('site_year'):assoc.append(dict(scope=scope,site_year=sy,raw_spearman=spearmanr(g.disease_severity,g.observed_yield).statistic))
 assoc=pd.DataFrame(assoc);j=assoc.merge(d[d.model_name.eq('CatBoost')&d.comparison.eq('M1_vs_M0')],on=['scope','site_year']);csv(j,'association_prediction_environment_comparison.csv')
 fig,axs=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
 for ax,scope in zip(axs,['FULL','NATURAL']):
  g=j[j.scope.eq(scope)];ax.scatter(g.raw_spearman,g.delta_MAE,c=np.where(g.site_year.str.startswith('QLB'),'#d88a44','#39749e'),s=35)
  for z in g.itertuples():ax.annotate(z.site_year,(z.raw_spearman,z.delta_MAE),fontsize=6,xytext=(3,3),textcoords='offset points')
  ax.axhline(0,color='gray',lw=.7);ax.axvline(0,color='gray',lw=.7);ax.set_xlabel('Within-site-year raw disease-yield Spearman');ax.set_ylabel('CatBoost disease increment (dt/ha)');ax.set_title(scope+' fitted cohort; descriptive, not causal')
 fig.suptitle('C4 | Negative association does not ensure positive predictive increment');save(fig,'C4_association_prediction')
 # Scientific reports use computed tables, not hand-entered favorable subsets.
 primary=ag[ag.model_name.eq('CatBoost')&ag.comparison.eq('M1_vs_M0')]
 report('01_Error_Decomposition_Findings.md','''# Error decomposition findings

完成：真实LOSYO预测，同地块M0/M1；Ridge及M2/M3辅助结果也全部保留。正Δ表示加入病害后误差减少。每个环境严格满足MSE = MBE² + Var(error)，使用总体方差ddof=0；最大恒等式误差见prediction_verification.json（约2.27e-13）。单位：MAE/RMSE为dt/ha，平方项为(dt/ha)²。误差定义predicted−observed。

'''+table(primary)+'''
主结果不是“病害只改善环境均值”的普遍规律。FULL环境等权总ΔMSE=27.923，其中偏差平方改善45.227，环境内部误差方差反而恶化17.304；样本量加权时两部分分别改善51.152和11.236。NATURAL单独训练的环境等权两部分为+16.943与−5.019。平均方式可改变总体叙事，不能只给pooled指标。

FULL的13环境中，偏差平方改善8个、中心化方差改善6个。QLB2015 ΔMSE=239.650，其中204.718（85.4%）来自偏差平方，34.932来自中心化方差。QLB2017偏差项改善65.913、方差项恶化3.160。HAN2015/16/17的方差项分别恶化53.069/70.306/107.883，不能用环境均值偏移单独解释；KIE2015原始MAE改善2.637，中心化MAE反而恶化2.263，且仅35地块。

FULL样地Spearman改善5/13，品种均值Spearman改善4/13；平均变化分别−0.0455/−0.0642。NATURAL为4/10与3/10，平均−0.0127/−0.0118。因此没有稳定的田间或品种排序收益。按真实产量上下四分位定义的高低产品种MAE已输出于02表；这是诊断分组，不是可部署的品种选择规则。品种均值未按管理平衡、不是空间BLUE，重复支持数逐环境保留，不将相关地块当独立重复。

裁定：PARTIAL_EXPLANATION。精确分解定位了各环境的误差组成，却没有识别造成分布迁移的唯一原因。MAE没有此加法恒等式，不能把MSE贡献比例当MAE贡献比例。中心化使用测试真均值仅作事后诊断，绝未校准预测后报告为实际泛化精度。

辅助Ridge/M2/M3及所有环境数据见01–03表和error_component_aggregation.csv；未按结果删除任何环境。
''')
 report('02_Adjusted_Association_Findings.md','''# Adjusted disease-yield association

复用审计后的真实小数。全队列8539、自然域5899；最终真实零值均为453，分别5.305%及7.679%。Pilot已从原始source_file/source_row恢复199条被整数截断为零的正小数；本次未重新清洗、未将缺失补零。来源和哈希继承见input_integrity与inheritance_receipt。

A0为截距+严重度；A1吸收site-year；A2吸收site-year×management；A3再吸收品种固定效应。采用加权交替去均值FWL估计，每次环境重采样重新吸收所有固定效应；A3另用独立稠密哑变量投影核验到1e-7。不是预测模型拟合。499次site-year pairs bootstrap、种子42，百分位95%描述性区间；不是普通地块bootstrap或wild bootstrap。只有13/10环境，且同地点、品种、年份有依赖，区间不能当精确确认性推断。

'''+table(a[['scope','model','n','environments','association_per_10_score_units','per10_CI95_lo','per10_CI95_hi','residual_severity_sd','residual_df']])+'''
所有主模型可识别且吸收收敛；A3残差病害SD为6.017（FULL）/6.612（NATURAL），残差df8289/5655。不是因为病害完全无变异而得到零信息。相对A0，A3点估计绝对值衰减约90.5%/91.9%；CI跨零，不能表述为“控制品种后依然存在明确的独立负关联”，也不能证明没有作用。品种固定效应会吸收遗传抗性及其与产量的共同变化，其目标是条件关联，不是疾病总因果损失。

逐环境A3负向10/13，自然域7/10；比原始负相关12/13、9/10减弱。删除任一环境后A3点值仍负：FULL每10名义分值−1.565至−0.609，NATURAL−1.233至−0.031；方向稳于精度，但不替代不确定区间。05表列环境内A0/A2/A3，06表列所有删除环境与同地点敏感性，无人为生成不可识别系数。残差偏度、异方差比、设计秩与收敛也随表保存；中位数两侧少于2例时异方差比为NA，不表示没有异方差。

**单位边界**：全/自然合并系数仅以每10个名义报告分值缩放。QLB有平均病程AO证据，其他地点调查时间、汇总端点不完全明确；不能无条件解释为统一调查日的10个百分点，更不是因果减产值。地点内估计也不应跨地点比较绝对幅度或把低范围评分外推10分。

有文献支持的QLB可比AO子集单独给出A3：每10个AO叶面积百分点评分，关联−3.017 dt/ha，描述性区间[−8.124, −0.362]；仅3年份，用27种有序环境pairs样本穷举，区间高度离散、缺乏独立地点。其负向证据限于同一人工接种试验，不能覆盖自然感染条件。见QLB_comparable_AO_association.csv。QLB的AO是时间汇总端点，不是同步调查日严重度。
''')
 manifest=read('07_matched_test_domain_manifest.csv');drop=read('matched_domain_delete_environment.csv')
 report('03_Natural_Domain_Transfer_Findings.md','''# Matched natural-domain transfer

完成：两策略×两算法×M0/M1×10留出环境，80项既有模型/OOF精确复用；新机器学习拟合0。所有测试均为同5899个自然感染ID，训练集完全排除留出site-year。按sample_id做一对一配对，逐模型train_id_sha256与原始Pilot状态一致，测试ID哈希也一致；07表为实际模型复用状态及域清单，不是拟合计划。

A_NATURAL_ONLY每折4579–5864训练地块、9环境、人工相关比例0；B_FULL_ELIGIBLE为7219–8504地块、12环境、人工相关比例31.04%–36.57%。B同时增加训练量并改变地点/年份/管理与病害分布，结果不能归因人工接种本身。自然域不含UNKNOWN归类记录。

下表正值表示病害特征改善；B_MINUS_A_INCREMENT比较的是两策略的“病害增量差”，不是直接比较两模型绝对MAE。两策略绝对MAE/RMSE/R²逐环境见matched_strategy_environment_metrics.csv；2000次配对环境重采样区间仍属描述性。

'''+table(m)+'''
CatBoost在共同测试域，A的pooled/macro ΔMAE=0.073/0.445，B=0.550/−0.602 dt/ha。B−A的macro=−1.047，95%[−3.804,1.383]，pooled=+0.477，95%[−2.719,2.309]；5/10环境病害增量更好。权重方向相反，不能选pooled得出人工数据有益。Ridge也给出B−A macro−1.033，区间[−4.280,1.048]。

'''+table(drop.groupby('model_name')[['macro_A','macro_B','macro_B_minus_A','pooled_B_minus_A']].agg(['min','max']).reset_index().pipe(lambda q:q.set_axis(['model_name']+['_'.join(c).strip('_') for c in q.columns[1:]],axis=1)))+'''
全部删除单环境结果保留；不通过删掉不利环境生成主结果。裁定：训练组成敏感，未证明B稳定提升；同域比较支持收窄而非推翻Pilot的谨慎判断。复用而非重新拟合完全符合预定预算；模型配置和折内预处理继承Pilot，无调参。
''')
 literature()
 report('BRIWECS_Phase2B_Closure_Scientific_Integration.md','''# Phase 2B Closure scientific integration

## Q1 — 调整后关联

全队列和自然域A3点估计都负，但描述性95%区间跨零；只能说观察关联显著衰减、独立稳健性不足。QLB AO同端点负向结果局限于3年人工相关试验，不能替代自然域证据。品种抗性与品种产量共同变化是可能解释，不能由固定效应衰减证明全部因果混杂已经识别。

## Q2 — 误差层面

偏差平方部分解释了总体收益，特别是QLB2015；FULL环境等权中心化误差方差恶化，排序多数未改善。HAN的失败主要在环境内误差，故“仅均值偏差导致失效”不成立。结论PARTIAL_EXPLANATION；没有使用测试真均值校准实际预测。

## Q3 — 固定自然测试域

共同5899样本上，B−A CatBoost macro增量−1.047 [−3.804,1.383]，pooled+0.477，方向不一致；训练人工相关比例与训练量等共同改变，不能归因人工接种。没有稳定预测增益，也没有证明严格零增益。

## Q4 — 贡献与文献重叠

2020同QLB试验已给出病害—产量关系；2025同数据论文已有G×E×M、TFI与稳定性；2026源库论文复用BRIWECS并包含抗病育种进展。跨环境验证、性能下降和尺度分离在外部文献中均已存在。A不新颖，B强版本证据不足、修正版WEAK，C限域实证MODERATE。独立价值是可复核的病害特征信息价值边界审计，而不是新抗病机制、首次关联或高精度管理工具。

## 十项整合答复与定位

1/2：调整后的全/自然负点值相同方向，区间均跨零。3/4：均值偏差提供部分解释，田间差异也会恶化，非唯一机制。5：共同测试域使组成混杂更清楚，并未支持稳定优势。6：人工相关训练效应只能限定为训练策略敏感性。7：P1/P2/P3同谱系重叠已明确。8：具体的同折误差审计与共同域对照有区别，但没有强新机制。9：最适合农业机器学习的负结果/评估边界案例，植物病理机制与农药管理证据不足。10：不再投入模型搜索；仅有条件考虑简短边界论文，未授权写稿或2C。

推荐英文题目：**Limits of stripe-rust severity as an incremental predictor of wheat yield across environments: a BRIWECS reanalysis**。

三项可保留贡献：恢复审计小数后区分粗关联和品种调整关联；以配对OOF精确定位均值偏差与田间误差、核查排序；固定自然测试域量化训练组成敏感性。三项主要局限：仅5地点3年且环境依赖、区间宽；评分端点和调查时间异质、同数据文献重叠；非因果设计、无独立外部验证且增益不稳定。

农业ML：有限的应用评估/负结果短文潜力。植物病理：独立机制贡献弱，不适合作强病害损失研究。农学与管理：仅回顾性信息价值边界，无农药收益或施药阈值证据。这里是证据定位，不保证期刊接收或分区。
''')
 report('BRIWECS_Phase2B_Closure_GO_NO_GO.md','''# Formal scientific closure decision

**Manuscript CONDITIONAL GO — 仅限较窄的农业机器学习评估边界案例。**

原主线“控制品种后仍稳健的病害负关联，却不能跨环境预测”的强版本为NO-GO：A3全/自然区间跨零，不能作为独立稳健生物关联前提。稳定预测改进、人工接种训练优越性、新抗病机制等主张均不获支持。

不裁定全面Manuscript GO，因调整关联不够明确、自然域策略比较区间宽、创新没有STRONG等级。也不把所有工作判为没有信息：误差分解发现环境等权中心化方差恶化、pooled与macro反向，并在固定自然测试域识别组成策略不稳定，已超出此前性能表。因此选择CONDITIONAL GO而非自动追加计算或强行造正结果。

| 验收项 | 状态 | 后果 |
| --- | --- | --- |
| 数据/OOF继承、同域配对、误差恒等式 | PASS | 可复核分析成立 |
| 全/自然品种调整负关联 | LIMITED | 点值负、区间跨零；不支持强独立关联 |
| QLB可比AO负向关联 | SCOPE LIMITED | 仅3年人工相关试验，不能推广 |
| 病害跨环境预测增益 | UNSTABLE | 宏平均与pooled权重敏感 |
| 误差/同域分析的新信息 | MODERATE EMPIRICAL | 部分解释，非新算法/机制 |
| 文献新颖性 | LIMITED | A不新颖；B弱；C限域中等 |
| 进入Phase2C或正式写作 | NOT AUTHORIZED / NOT STARTED | 本阶段完成后停止 |

后续可接受范围只可能是：经用户另行明确授权，以现有冻结表图完成短文结构与措辞收敛，核对有限引用、附录和可复现性。不是本轮待办，不要求扩大模型/病害/数据、调参、SHAP、追加时间窗口或湿实验。如果目标仍是强植物病理机制论文或高水平新预测方法文章，则应停止该方向投入；现有证据达不到目标。

必须删除的主张：稳健独立负关联已证实；病害导致精确减产；模型普遍提升；人工接种数据有因果益处；首次发现泛化下降/误差分解；可指导施药经济收益。

本次新增机器学习拟合0；复用80项主要模型预测，关联统计不是新ML训练。四类图、十项核心表、六份报告均已完成；文献全文访问限制已记录，属非阻断性，不伪报全文完整。禁止自动进入2C，至此停止新增科研计算。
''')
 versions={'python':platform.python_version()}
 for package in ['numpy','pandas','scipy','scikit-learn','statsmodels','catboost','matplotlib','pyarrow']:versions[package]=md.version(package)
 js(versions,'config/runtime_versions.json');(ROOT/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in versions.items() if k!='python')+'\n',encoding='utf-8')
 (ROOT/'README.md').write_text('''# BRIWECS-ML Phase 2B Closure

状态：COMPLETE WITH DOCUMENTED LITERATURE ACCESS LIMITS。正式裁定为Manuscript CONDITIONAL GO，仅限边界案例；没有启动Phase2C或正式写作。

数据：Figshare官方v1，继承Phase2A与Phase2B_Pilot；Stripe rust FULL8539/13SY、NATURAL5899/10SY。只写本目录。新ML拟合0，复用80项既有模型/OOF；全部核心输入按既有SHA清单验证。

先读reports/BRIWECS_Phase2B_Closure_GO_NO_GO.md及Scientific_Integration.md；各WP报告给方法与数字。results/tables/01–10为核心交付，另附汇总、敏感性、文献来源与复核表。results/figures/C1–C4各有PNG/PDF。data/stripe_rust_inherited.parquet为值不变的队列快照，data/literature含公开来源响应；部分HTTP200页面仅反机器人或订阅预览，实际证据层级以literature_evidence_registry.csv为准。

## 本地离线复现

在此目录运行PowerShell：

```powershell
$env:PYTHONIOENCODING='utf-8'
$env:OPENBLAS_NUM_THREADS='2'
$env:MKL_NUM_THREADS='2'
& '..\\Phase2A\\.venv\\Scripts\\python.exe' scripts\\run_closure.py
```

需要相邻Phase2A、Phase2B_Pilot和本次已归档Crossref记录。requirements.txt记录实际版本；无需下载新模型。run_closure依次重算预测诊断/关联与图表报告，然后校验并生成SHA清单；不联网，不新拟合ML。关联499次环境pairs bootstrap非普通地块bootstrap；QLB3环境另用27种pairs穷举。data/literature的来源下载是可选独立脚本literature_sources.py，会联网，离线复现不调用。

config/request.md为本轮完整授权规范，analysis.json为冻结分析设置。logs/retraining_status.json与07表记录真实复用状态；logs记录实际执行。统计模型与新ML模型计数分开。输出manifest排除自身及动态缓存，不具有对第三方独立审核的替代意义。既有阶段文件全部保持不变。
''',encoding='utf-8')
 print('Four PNG/PDF figure pairs and six scientific reports generated. Manuscript CONDITIONAL GO, narrow scope only.')

if __name__=='__main__':main()
