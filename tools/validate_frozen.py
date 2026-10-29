"""Read-only replay of frozen evidence. No estimator fitting or inference."""
from pathlib import Path
import ast,csv,json,hashlib,sys,platform,importlib.metadata
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1];checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,pass_check=bool(ok),detail=str(detail)))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ih(x):return hashlib.sha256('\n'.join(sorted(map(str,x))).encode()).hexdigest()
for row in csv.DictReader((R/'provenance/Lean_Input_Receipt.csv').open(encoding='utf-8-sig')):
 p=R/row['repository_path'];check('authentic_byte_copy:'+row['repository_path'],p.exists() and sha(p)==row['sha256'])
for p in R.rglob('*.py'):
 try:ast.parse(p.read_text(encoding='utf-8-sig'));check('python_syntax:'+p.relative_to(R).as_posix(),True)
 except Exception as e:check('python_syntax:'+str(p),False,e)
versions={n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','scikit-learn','statsmodels','catboost','matplotlib','pyarrow','xgboost']}
recorded=json.loads((R/'Phase2B_Closure/config/runtime_versions.json').read_text())
for n,v in versions.items():
 expected=recorded.get(n,recorded.get('sklearn' if n=='scikit-learn' else n))
 if expected:check('runtime_version:'+n,v==expected,f'{v}; historical={expected}')
co=pd.read_parquet(R/'Phase2A/data/cohort/phaseI_yield_cohort.parquet');freeze=json.loads((R/'Phase2A/config/cohort_freeze.json').read_text())
check('yield_cohort_size',len(co)==24991);check('yield_sample_ids',co.sample_id.is_unique and ih(co.sample_id)==freeze['sample_id_set_sha256'])
check('yield_environments',co.site_year.nunique()==15)
s=pd.read_parquet(R/'Phase2B_Pilot/data/stripe_rust_paired_cohort.parquet')
check('stripe_cohort',len(s)==8539 and s.sample_id.is_unique);check('stripe_envs',s.site_year.nunique()==13)
check('natural_cohort',s.inoculation_status.eq('NATURAL').sum()==5899)
check('natural_envs',s.loc[s.inoculation_status.eq('NATURAL'),'site_year'].nunique()==10)
check('true_zero_scores',s.disease_severity.eq(0).sum()==453)
check('restored_fractional_scores',((s.published_disease_value==0)&(s.disease_severity>0)&(s.disease_severity<1)).sum()==199)
for entry in json.loads((R/'Phase2B_Pilot/config/cohort_freeze.json').read_text()):
 name='stripe_rust' if entry['disease']=='Stripe_rust' else 'septoria'
 p=R/f'Phase2B_Pilot/data/{name}_paired_cohort.parquet';d=pd.read_parquet(p)
 check(name+'_frozen_file',sha(p)==entry['cohort_file_sha256']);check(name+'_id_hash',ih(d.sample_id)==entry['sample_id_sha256'])
folds=pd.read_csv(R/'Phase2B_Pilot/splits/disease_outer_fold_manifest.csv')
for row in folds[folds.disease.eq('Stripe_rust')&folds.validation_strategy.eq('LOSYO')].itertuples():
 d=s[s.inoculation_status.eq('NATURAL')] if row.scope=='NATURAL' else s
 test=d[d.site_year.eq(row.test_group)];train=d[~d.site_year.eq(row.test_group)]
 check('LOSYO:'+row.scope+':'+row.fold_id,ih(test.sample_id)==row.test_id_sha256 and ih(train.sample_id)==row.train_id_sha256 and len(test)==row.test_n and len(train)==row.train_n and not set(test.sample_id)&set(train.sample_id))
pred=pd.read_parquet(R/'Phase2B_Pilot/results/predictions/disease_pilot_oof_predictions.parquet')
check('prediction_unique_keys',not pred.duplicated(['disease','scope','validation_strategy','model_name','context','sample_id']).any())
d=pred[pred.disease.eq('Stripe_rust')&pred.validation_strategy.eq('LOSYO')&pred.model_name.eq('CatBoost')&pred.scope.eq('FULL')]
j=d[d.context.eq('M0')].merge(d[d.context.eq('M1')],on='sample_id',suffixes=('_0','_1'),validate='one_to_one')
check('paired_full_predictions',len(j)==8539 and set(j.sample_id)==set(s.sample_id) and np.array_equal(j.observed_yield_0,j.observed_yield_1))
e0=j.predicted_yield_0-j.observed_yield_0;e1=j.predicted_yield_1-j.observed_yield_1
pooled=float(np.mean(abs(e0))-np.mean(abs(e1)))
points=[]
for sy,g in j.groupby('site_year_0'):
 e0=g.predicted_yield_0-g.observed_yield_0;e1=g.predicted_yield_1-g.observed_yield_1
 points.append(dict(site_year=sy,n=len(g),delta_MAE=np.mean(abs(e0))-np.mean(abs(e1)),delta_MSE=np.mean(e0**2)-np.mean(e1**2),delta_bias_squared=e0.mean()**2-e1.mean()**2,delta_centered_variance=np.var(e0,ddof=0)-np.var(e1,ddof=0)))
points=pd.DataFrame(points);agg=pd.read_csv(R/'Phase2B_Closure/results/tables/error_component_aggregation.csv')
for weighting in ['ENVIRONMENT_EQUAL','PLOT_WEIGHTED']:
 row=agg[agg.scope.eq('FULL')&agg.model_name.eq('CatBoost')&agg.comparison.eq('M1_vs_M0')&agg.weighting.eq(weighting)].iloc[0]
 for metric in ['delta_MAE','delta_MSE','delta_bias_squared','delta_centered_variance']:
  val=np.average(points[metric],weights=points.n if weighting=='PLOT_WEIGHTED' else None)
  check('frozen_accounting:'+weighting+':'+metric,np.isclose(val,row[metric],atol=1e-10),f'{val:.12g}')
check('pooled_MAE_increment',np.isclose(pooled,1.5810435434576888),pooled)
rank=pd.read_csv(R/'Phase2B_Closure/results/tables/paired_ranking_changes.csv');rank=rank[rank.scope.eq('FULL')&rank.model_name.eq('CatBoost')]
for level,count in [('PLOT',5),('CULTIVAR_MEAN',4)]:
 q=rank[rank.level.eq(level)];check('frozen_rank_count:'+level,len(q)==13 and q.delta_spearman.gt(0).sum()==count)
nat=pred[pred.disease.eq('Stripe_rust')&pred.validation_strategy.eq('LOSYO')&pred.model_name.eq('CatBoost')&pred.context.eq('M1')]
common=set(s.loc[s.inoculation_status.eq('NATURAL'),'sample_id'])
for scope,expected in [('NATURAL',12.9661065),('FULL',13.4457362)]:
 q=nat[nat.scope.eq(scope)&nat.sample_id.isin(common)];v=float(np.mean(abs(q.predicted_yield-q.observed_yield)))
 check('matched_natural_MAE:'+scope,len(q)==5899 and set(q.sample_id)==common and np.isclose(v,expected,atol=1e-6),v)
report=dict(status='PASS' if all(x['pass_check'] for x in checks) else 'FAIL',checks=len(checks),failures=[x for x in checks if not x['pass_check']],python=platform.python_version(),dependencies=versions,new_model_fits=0,new_inference_runs=0,validation_scope='Static compilation, byte hashes and deterministic replay of existing predictions; no clean-room fit reproduction',details=checks)
(R/'validation_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='details'},indent=2));sys.exit(0 if report['status']=='PASS' else 1)
