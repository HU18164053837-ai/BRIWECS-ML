from pathlib import Path
import hashlib,json,datetime as dt
import pandas as pd
def md(df):
 df=df.copy().fillna('')
 def fmt(v):
  if isinstance(v,float):return f'{v:.4g}'
  return str(v).replace('|',' / ').replace('\n',' ')
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(fmt,row))+' |' for row in df.itertuples(index=False,name=None))+'\n'
def make_reports(root):
 a=root/'audit';r=root/'reports';s=json.loads((a/'summary.json').read_text());m=pd.read_csv(a/'source_file_manifest.csv');cov=pd.read_csv(a/'03_trait_coverage.csv');ov=pd.read_csv(a/'05_disease_yield_overlap.csv');comp=pd.read_csv(a/'07_management_comparison_feasibility.csv');versions=pd.read_csv(a/'version_population_counts.csv');prov=pd.read_csv(a/'disease_missingness_provenance.csv');ph=pd.read_csv(a/'phase_actual_samples.csv');v=pd.read_csv(a/'11_validation_feasibility.csv');dates=pd.read_csv(a/'disease_date_evidence.csv');git=json.loads((root/'sources/git_versions.json').read_text());metadata=json.loads((root/'sources/figshare_article.json').read_text())
 disease=cov[(cov.basis=='source_supported')&cov.trait.isin(['Stripe_rust','Leaf_rust','Powdery_mildew','Septoria','DTR','Fusarium'])].copy()
 disease=disease.set_index('trait').loc[['Stripe_rust','Leaf_rust','Powdery_mildew','Septoria','DTR','Fusarium']].reset_index()
 dc=disease[['trait','valid','positive','years','sites','environments','managements','minimum','median','maximum']]
 pairs=ov[(ov.basis=='source_supported')&ov.Location.eq('ALL')][['trait','paired_units','environments']]
 cp=comp[comp.common_cultivars.gt(0)].groupby(['treatment_a','treatment_b','status']).agg(comparable_location_years=('Year','size'),cultivar_environment_intersections=('common_cultivars','sum'),units_a=('units_a','sum'),units_b=('units_b','sum')).reset_index()
 filecounts=m.groupby('source_repository').agg(files=('relative_path','size'),bytes=('file_size_bytes','sum')).reset_index()
 sources=f'''# BRIWECS Phase 1：官方来源与版本核验

状态：COMPLETED。冻结主分析版本为Figshare **v1**（{metadata['published_date']}），DOI [10.6084/m9.figshare.27910269.v1](https://doi.org/10.6084/m9.figshare.27910269.v1)。当前GitHub仅用于差异核查，不与v1拼接。

## 获取与完整性

官方文件列表由[Figshare API](https://api.figshare.com/v2/articles/27910269)实际解析：一个Briwecs_data.zip，11,329,195 bytes；原包保留，官方MD5为fb6543627c67a368b4fe5f078757f8b5，本地核验一致。SHA256逐文件见`audit/source_file_manifest.csv`，包括原压缩包及解压文件。

{md(filecounts)}
以上为实际保存的非.git文件，共{len(m)}个；按SHA256去重为{m.sha256.nunique()}种内容。原始网络下载包1个，不能把解压文件数当作独立下载请求数。完整仓库克隆等待期间又成功浅克隆一次，保留为冗余副本并标记同哈希；不参与主分析、不增加样本量。Git内部对象不计入清单。

{md(pd.DataFrame(git)[['repository','branch','commit','status']])}

Git status为空，作者代码未修改。raw数据未写回；所有派生审计在audit。压缩包使用下载落地时间；逐文件清单注明该时间定义。公开候选BASF小表另存sources并单独哈希。

## 526,751与561,700差异

{md(versions)}
按Unit.xlsx的24项性状逐单元格计数，Figshare精确复现论文526,751；GitHub精确得到561,722，与README约561.7K一致。因此均非样地数。差额34,971个性状单元格、210行。具体环境差异见`version_environment_counts.csv`；同路径文件逐字节对比见`figshare_github_comparison.csv`。

已确认：GitHub增加RHH 2018两种WF处理共200行，GGE 2019五种处理各增加2行，共10行；GGE 2016/2017存在RF→IR重新标记。`version_delta_reconciliation.csv`逐地点×年份精确加总34,971差额：KAL 2018/2019/2020分别+13,169/+13,248/+7,956；RHH 2018 +1,778；GGE 2016/2017/2018/2019分别−880/−885/+517/+68。

代码证据进一步表明：旧版写出屏蔽KAL第二阶段的res；当前版虽仍计算res，却写出未屏蔽的export_dat，解释KAL新增大量性状单元格。当前版还为KAL将BRISONr_214–220映射到222–228，修补HAN 2018坐标并借助补丁表补GGE 2019身份。RHH 2018当前单独用逗号读取，旧脚本统一read_csv2。每个具体值变化的完整历史因果链未重跑R闭合，仍标 **UNRESOLVED_VERSION_DIFFERENCE**；不能把新旧文件当等价替代。仓库描述2015–2019与实际两个版本均有2020不一致；2020仅KAL。README KAL Phase II整理警告与这些脚本规则对应，版本选择因此是研究限制而非单纯总数不一致。

## 特别文件核查

{md(pd.read_csv(a/'expected_file_check.csv'))}
品种信息实际名为`metadata/BRIWECs_cultivar_info.csv`，不是论文写出的文件名。`disease_record.xlsx`缺失有脚本证据：extract_management.R的第三项条件、disease提取和导出均被注释；原始管理xlsx仍保留部分调查时间。因此不视为核心性状表缺失，但早期预报时间门未通过。

## 来源与获取限制

最初沙箱命令行网络被拒，联网授权后curl下载成功；web读取Nature一度重定向失败，改读同文PMC官方全文。第二仓库首次探测main返回422，真实分支master已克隆并冻结。不存在来源不明镜像替代。没有重建作者全部旧R依赖或重跑清洗。Scientific Research Design Specification v1.0的另一个完整文本未提供，本次依据用户Phase 1执行规格，未臆造额外设计要求。

来源：[原论文](https://doi.org/10.1038/s41597-024-04332-7)、[PMC全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11732995/)、[作者预处理仓库](https://github.com/tillrose/BRIWECS_Data_Publication)、[作者验证仓库](https://github.com/Illustratien/Scientific_Data_Analyis)。本报告中的规模由本地冻结数据直接计算。
'''
 (r/'BRIWECS_Phase1_Source_Audit.md').write_text(sources,encoding='utf-8')
 science=f'''# BRIWECS Phase 1：科学可行性审计

整体：**CONDITIONAL GO**。数据允许开展有限范围的跨环境产量预测与回顾性单病害–产量关联；当前不能支持全六病完整病例跨环境主线，也没有完成无泄漏的早期病害预报数据集或正式独立外部测试。审计已执行，未训练模型。

## 1. 下载到了哪些官方文件、使用哪个版本？

主分析采用Figshare v1；原始数据包、解压结果、两个作者Git仓库及commit均保存，核心表完整。获取数量与哈希见Source Audit及清单。不存在的合并病害时间表没有伪造。

## 2. 真正统计单位和合法样本规模

发布表34,683行、33列，是混合整样地与50 cm取样性状的宽表；非长表、非52万独立样地。共228品种、6地点、6年份、28地点×年份、98地点×年份×管理、9管理情景。论文所谓100环境与v1实际98也应保留为计数差异。

候选键为Year×Location×Treatment×Block×Row×Column×BRISONr。排除GGE 2019的22行重复身份，并进一步排除HAN/KAL的120行坐标冲突后，有**{s['reliable_units']:,}个可唯一追溯的记录单位**。其中**{s['coordinate_complete_units']:,}条同时有完整且不冲突坐标**，499条无完整坐标但品种×区组×环境键可唯一追溯。后者可用于记录级配对，但不能声称空间样地身份已完全证实；主分析优先坐标完整子集，499条作敏感性扩展。142条争议记录单列，未任意选一条或平均。作物试验的区组/环境内相关仍存在，这些单位并非统计独立同分布。

`02_join_key_integrity.csv`显示简单品种×环境或再加Block仍非唯一，不允许据此多对多并表。所有病害与产量来自同一发布行，通过唯一原始键检查；只做many-to-one原始来源关联并断言行数不扩张。KIE 2017作者曾平均一组重复，源文件保持不变。本次无法从发布表恢复原始逐日重复调查次数，记NOT_IDENTIFIABLE，绝不将其作为新样本。

Seedyield为整样地dt/ha、100%干物质；Seedyield_bio为成熟期50 cm取样g/m²，二者不得互换。后者不能另算一块独立田间样地。

## 3. 六病标签及与产量共同覆盖

以下采用“原始非缺失非负病害值支持+键无争议”的保守口径，数值仍是官方清洗后的整数分值；并非声称恢复全部原始精度。

{md(dc)}
有效产量{ s['strict_yield']:,}条（发布表原有{s['published_yield']:,}条）；至少一病+产量{s['any_pair']:,}条，覆盖22地点×年份。各病配对如下：

{md(pairs)}
六病全部完整+产量仅416条，全部HAN 2019；TFI完整值最大121，符合加总指标可大于100的含义。不能用这一单环境子集验证跨环境六病联合预测，也不把缺失病害补0计算TFI。没有发现六病只有零的全局标签，但DTR仅3地点、7环境、4管理，适合作探索或次要终点；Fusarium仅8环境且有接种混杂。条锈病覆盖最广；Septoria、叶锈病、白粉病可作次要单病任务，须分别限定目标定义和可用验证环境。

![病害覆盖](disease_coverage.svg)

## 4. 最关键的标签审计问题

作者`data_cleaning.R`明确在Y/L/M内非全缺失的病害列中将NA或负值替0；另以as.integer将原始小数向零截断。

{md(prov)}
前三列为实际清洗影响，不是模型插补。对于发生/不发生分类，原始0到1之间的正值被截为0也会改变标签；故不能直接用发布表0训练发生分类。Phase 2必须重新冻结基于原始病害值的目标规则，同时保留日期/评分定义；本阶段不替作者重做整套清洗。`03_trait_coverage.csv`分别提供published与source_supported口径，缺失率分母为34,683行，零比例分母为对应有效值。

原始文献QLB方法以重复调查的AO汇总锈病等，FHB为组合指标，而公开整合表是单列数值；混合即时严重度、AO与整数截断不能未经核查当相同时间点标签。这个问题限制前瞻性疾病任务，不阻断回顾性关联，但应分层或限制至定义兼容的数据。

## 5. 哪些管理处理可比？

{md(cp)}
表内cultivar_environment_intersections是每个地点×年份共同品种数的合计，不是全项目独立品种数；逐环境共同品种及各病、产量覆盖见07。COMPARABLE只表示同环境有处理重叠，不代表随机化设计已经逐地重建或支持因果效应。加总units是处理的全部可靠单位，不全是共同品种子集；不得作为一对一配对数。品种支持逐项见management_cultivar_support.csv。

优先同环境HN_NF_RF/LN_NF_RF氮水平比较和非QLB的WF/NF描述性比较。灌溉对比仅GGE三年，遮雨仅KIE三年，水分处理与地点完全绑定，不能推断一般化水分因果效应。GGE 2016/2017处理标签在GitHub变更为IR，因此使用v1 RF氮/杀菌剂支持域仅是版本内暂定结论，**Phase 2须先核查原始水分管理、必要时排除这两年相关比较**。同Block命名仅支持区组信息存在，不证明跨处理同名区组或重复样地能一对一匹配。

## 6. 人工接种与两阶段

按病害特异性标记，QLB 2015–2017 HN_NF_RF/LN_NF_RF的条锈病、叶锈病与Fusarium为ARTIFICIAL；其余病害/场景依据主论文自然侵染说明标NATURAL，标记证据层级为论文，不是逐样地接种日记。不存在将QLB全部病害全部年份标人工接种的操作。NF与WF杀菌剂比较在该地点同时改变接种条件，标签LIMITED；氮比较可在同接种制度内进行。人工接种日期仅有季节或生育期范围，不可伪造精确日期。

{md(ph)}
按2015–2017/2018–2020划分，实际ID并集为Phase I 228、Phase II 221，阶段共同ID为{s['cross_phase_cultivars']}，并非全局52。除KAL后第二阶段及跨阶段共同ID均为52。Phase I多数单个地点×年份为220，但跨环境集合并不完全相同，因此并集228不能解释成每环境228。非KAL第二阶段多为52品种，GGE 2019实际51；RHH 2017仅50。2020仅KAL 884行、221品种，发布版病害均缺失。当前作者脚本还对KAL第二阶段7个品种编号做替换，因此上述为观测ID口径，不证明生物材料身份全已复核；详phase_panel_overlap.csv。

Phase I→II可以作组成变化下迁移敏感性分析，不能说是相同试验设计下的纯时间泛化。主论文对QLB/RHH 2017及GGE/KAL/RHH 2018列有设计例外，08中逐环境标记。Phase I子集对主要终点保留足够跨环境记录，适合单独分析。

## 7. 天气和农事时间是否合法？

6份DWD插值逐日天气，每地点3,696日，2014-01-01至2024-02-13；温度均/最小/最大、降水、湿度、辐射、风速共7变量，日期完整无重复，各变量有效值详09。文件为分号分隔、小数逗号；本次正确解析而非把一行误作单列。气象完整性不等于预测合法性，也不等于测站实测无误。

播种/出苗在性状表为年内日序，需确定冬麦播种所属前一年；BBCH59/87有缺失。管理表有238条施肥、388条植保记录及日期，但作者Excel日期origin需用原xlsx复核。原管理工作簿找到{len(dates)}个病害模板单元格，其中205个邻接字段非空，只有{int(dates.date_parsed.notna().sum())}个是本阶段严格解析的单一ISO日期；其他含日期范围、多日期、BBCH、0、No/none。GGE 2018还有重复病害模板列，不能当重复实测调查。逐单元格证据保存坐标与原字符串，未用部分日期推断全部样地调查时间。

因此：pre-season只能用品种、可确认的基线土壤和计划管理；实际施肥/杀菌剂、未来整季天气不可用。in-season forecasting目前BLOCKED，缺少逐样地一致的目标调查日期/定义。retrospective association可用明确生长季窗口的季内天气，但须标回顾性。收获/灌溉信息存在于部分原始管理表，覆盖和字符串不统一，已保留原始证据，尚不适合宣称全体可用的精确截点。

产量预测必须排除Grain、Protein_yield、Biomass、Straw等由产量推导的指标；Seedyield_bio、成熟期质量/生物量也不作为早期预测输入。若使用病害分值预测产量，仅能表述季内/回顾性关联，不能冒称播前预测。10列出13项泄漏风险及控制。

## 8. 哪种严格验证最有意义？

11逐目标、逐折给出LOYO、LOLO、LOEO、阶段与管理迁移的训练/测试样本、品种覆盖、正/零类别、管理交集及未见管理。产量有6年份、6地点、28环境可留出；病害仅2015–2019有实测支持，2020病害折无测试样本，不能计入有效病害验证。每病有效环境见上表。

建议主评估LOEO（同地点×年份下所有管理与重复同时留出），辅以LOYO/LOLO。LOEO是新地点-年份组合检验，不等于全新地点；LOYO/LOLO分别报告才有清晰的外推边界。病害按各折阳/零类别和接种状态解释，缺类别的折不能计算分类AUC。2020产量测试只反映KAL单点；Phase I→II为品种/地点/设计共同变化。任何预处理、特征选择和BLUE均须训练折内完成。本次只审计支持，未构建或运行模型。

## 9. 是否有真正独立外部验证？

有界检查三个官方候选，详12。QLB Zenodo 3697514和Sabir 8248543属于BRIWECS既有试验/子集，NOT_INDEPENDENT。BASF 6521175小表已实际下载并解析，具有56个TrialId、处理、日期、病害评分和ERTRNE收获指标，但无品种/重复/样地ID。MCALCR未处理对照恒100，属于相对值；MWCGS的SampleUnit=HA仍未给出明确质量单位，不能凭缩写当dt/ha。病害分叶位/调查日、害物代码和输入交集也需统一。故为潜在独立试验级候选，**当前没有验证通过可直接用于正式外部测试的数据集**；未把内部留出包装成外部验证。检索到此停止，不为找外部数据无限扩展。

## 10. 最适合的科学问题与需要修改的设计

建议问题：在明确品种覆盖和管理支持域后，冬小麦产量及选定病害–产量关联在未见地点×年份中能否保持稳定，且这种稳定性是否因管理情景和接种制度而变化？目前仅证明数据支持评估，未证明预测精度或论文新颖性。建议主队列优先Phase I，非KAL的2018–2019作为阶段迁移扩展；KAL 2018–2020因编号/版本争议先暂缓，而不是把2020单点当完整未来年份验证。34,213为完整审计候选规模，最终主模型队列会因这些预设限制进一步缩小。

修改：放弃六病完整病例作为主线；以产量为主、条锈病/Septoria等各自可用队列为次要终点。区分原始非零和整数分值、AO和即时严重度；排除争议身份，保留缺失掩码。限定管理比较于真实共同环境/品种，QLB WF/NF及跨地点水分效应暂缓。整季天气仅用于回顾性分析；早期预测和独立外部验证先不作论文承诺。v1是正式冻结依据，GitHub更新另作版本敏感性核查，不能无记录替换。

## 完成状态与边界

COMPLETED：官方获取、版本/哈希、粒度/键/覆盖/管理/阶段/天气与泄漏审计、全部12张要求表、三个候选有界检索、裁定和复现入口。

PARTIAL：逐地随机化/空间设计重建、调查与接种精确日历、病害跨地点评分定义统一、版本差异完整历史解释、BASF质量单位及独立性确认。均未写成已证实结果。

BLOCKED：当前全数据无泄漏早期病害预报、六病完整队列的跨环境验证、正式独立外部验证；属于研究目标受限而非本次执行未获取核心文件。

NOT APPLICABLE：训练、调参、SHAP、因果估计、湿实验和Phase 2执行，按用户边界未开展。

证据：规模均来自本地表与脚本。来源语义参见[主论文](https://doi.org/10.1038/s41597-024-04332-7)、[QLB原研究方法](https://doi.org/10.1038/s41598-020-77200-0)、[QLB数据](https://zenodo.org/records/3697514)、[BASF候选](https://zenodo.org/records/6521175)、[Sabir数据](https://zenodo.org/records/8248543)。
'''
 (r/'BRIWECS_Phase1_Scientific_Feasibility_Report.md').write_text(science,encoding='utf-8')
 decisions=pd.DataFrame([
 ['A Disease prediction','CONDITIONAL GO','单病回顾性严重度可研究；发生分类需纠正整数截断，早期预报NO-GO；DTR/Fusarium限次要'],
 ['B Yield prediction','GO','34,213可靠产量，28地点×年份，6地点/年份；限定合法基线输入与严格分组验证，非精度保证'],
 ['C Disease–yield association','CONDITIONAL GO','12,993至少一病配对，22环境；需评分定义、接种和管理分层；完整六病跨环境NO-GO'],
 ['D Management comparison','CONDITIONAL GO','有真实共同环境/品种；QLB杀菌剂接种混杂、水分单地点、GGE版本标签争议需限制'],
 ['E External validation','NO-GO','当前正式外部验证未就绪；2个非独立子集，BASF候选仅CONDITIONAL，单位/输入尚不兼容']
 ],columns=['Track','Decision','Evidence / condition'])
 savepath=a/'track_decisions.csv';decisions.to_csv(savepath,index=False,encoding='utf-8-sig')
 decision=f'''# BRIWECS Phase 1 GO / NO-GO Decision

## Overall Decision

**CONDITIONAL GO**。核心可成立路线是严格环境留出的产量预测及限定病害的回顾性关联，必须缩减原始全目标方案。不是全六病早期预测、一般化管理因果分析或现成独立外部验证的GO。

## Track A–E Decisions

{md(decisions)}
## Evidence Supporting GO

官方v1的24性状有效单元格恰526,751，完整来源和哈希可核验；{s['reliable_units']:,}个唯一可追溯记录单位中，{s['coordinate_complete_units']:,}个坐标完整且无冲突，另499个仅具唯一记录键。可靠产量{s['strict_yield']:,}，至少一病+产量{s['any_pair']:,}，28总环境、22病害环境，可执行LOEO及LOYO/LOLO审计。

## Critical Limitations

病害NA→0、正小数→0及跨地点AO语义；142行身份/坐标争议；六病+产量416条仅HAN 2019；QLB接种与杀菌剂处理混杂；2020仅KAL无病害；GGE 2016/17水分标签有版本差异；调查日期不能完整关联；收获及产量推导变量泄漏。

## Unresolved Non-blocking Issues

GitHub比v1多34,971个性状值与210行已定位总量，完整历史解释未闭合；品种元数据文件改名、缺disease_record合并表；部分Excel日期和天气来源插值误差；BASF单位、身份独立性及输入交集。这些不能通过伪造或补零解决，也不阻断限定产量路线。

## Recommended Primary Scientific Question

产量与选定病害–产量关联能否在未见地点×年份上保持可迁移性，且管理支持域、品种组成与接种制度如何限制这种迁移？

## Recommended Scope Reduction

删除六病完整病例跨环境主线与全数据早期病害预报承诺；DTR/Fusarium仅次要或分层；QLB WF/NF不解释为单一杀菌剂效应；不推断一般化灌溉/遮雨效应；外部验证暂缓；GGE 2016/17处理比较先核查或排除。

## Proposed Phase 2 Dataset

以Figshare v1为冻结来源，`audit/audited_observations.csv`为审计索引，reliable_unit=True为合法记录候选；首选coordinates_complete=True且coordinate_collision=False子集，另499条无坐标单位仅作敏感性扩展。主分析优先Phase I，非KAL的2018–2019作阶段迁移扩展；KAL 2018–2020在品种编号和版本问题解决前暂缓。因此34,213是全部审计候选的产量规模，不是承诺主模型实际纳入规模。Seedyield主终点，六病各自可用掩码，TFI只在六项完整时描述。采用原始病害值重新定义标签前不得直接运行发生分类。管理/品种/环境索引由06及management_cultivar_support提供；病害接种状态由08按环境×病害关联，不依据地点一刀切。

## Phase 2 Readiness

**PARTIAL / READY WITH CONDITIONS**：产量数据可在剔除身份冲突、固定输入时间点和排除产量派生变量后进入开发；病害、管理和外部验证仍需上述范围收缩。11提供的是候选分组支持清单，不是已经训练或验证的模型性能。未自动进入Phase 2。

## Final Scientific Verdict

这套公开数据支持一条可靠、可复现的限定研究路线；原先“六病齐全、跨环境早期预报、全面管理比较、独立外部验证”的完整故事现阶段不成立。论文质量取决于正确处理观测单位、标签来源及外推边界，不能依靠总性状值规模证明。Phase 1至此停止。
'''
 (r/'BRIWECS_Phase1_GO_NO_GO_Decision.md').write_text(decision,encoding='utf-8')
 # Inspectable lightweight companion notebook delegates to the single, versioned audit entrypoint.
 nb=dict(cells=[dict(cell_type='markdown',metadata={},source=['# BRIWECS Phase 1\nRead-only audit; inspect scripts/run_phase1.py and reports. No model training.']),dict(cell_type='code',metadata={},execution_count=None,outputs=[],source=['from pathlib import Path\nimport subprocess, sys\nroot = Path.cwd()\nif root.name == "scripts": root = root.parent\nsubprocess.run([sys.executable, str(root / "scripts/run_phase1.py")], check=True)'])],metadata=dict(kernelspec=dict(display_name='Python 3',language='python',name='python3')),nbformat=4,nbformat_minor=5)
 (root/'Phase1_Audit.ipynb').write_text(json.dumps(nb,indent=2),encoding='utf-8')
 inv=[]
 for f in root.rglob('*'):
  if f.is_file() and not any(p in ['raw','.git'] for p in f.relative_to(root).parts) and f.name!='output_inventory.csv':inv.append(dict(path=f.relative_to(root).as_posix(),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 pd.DataFrame(inv).to_csv(a/'output_inventory.csv',index=False,encoding='utf-8-sig')
if __name__=='__main__':make_reports(Path(__file__).resolve().parents[1])
