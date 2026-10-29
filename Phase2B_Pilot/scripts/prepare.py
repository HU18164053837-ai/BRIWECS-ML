from common import *
import warnings
warnings.filterwarnings('ignore',category=UserWarning,module='openpyxl')
def read_raw(f):
 for enc in ['utf-8-sig','cp1252']:
  try:return pd.read_csv(f,sep=';' if ';' in f.open(encoding=enc).readline() else ',',encoding=enc,dtype=str)
  except UnicodeDecodeError:pass
 raise ValueError(f)
def prepare():
 checks=[]
 for base,manifest in [(P1,P1/'audit/source_file_manifest.csv'),(A,A/'reports/output_manifest.csv')]:
  for row in pd.read_csv(manifest).itertuples():
   p=base/row.relative_path;ok=p.exists() and sha(p)==row.sha256;checks.append(dict(path=p.relative_to(P1).as_posix(),sha256=row.sha256,match=ok))
 csv(checks,'config/inherited_hash_check.csv');assert all(x['match'] for x in checks),'Inherited inputs changed; inspect verification'
 ac=json.loads((A/'config/training_config.json').read_text());config=dict(stage='Phase2B-Pilot',version='2.2',time_frozen_utc=now(),source='Official Figshare 27910269.v1 only',source_phase2a_config_sha256=sha(A/'config/training_config.json'),seed=42,models={m:ac['models'][m] for m in ['Ridge','CatBoost']},categorical=ac['F1']['categorical'],F1_numeric=ac['F1']['numeric'],weather=ac['F2']['numeric_additions'],F2_time_window='Target-year January 1 to March 31, unchanged from Phase2A',task='RETROSPECTIVE ASSOCIATION-ORIENTED ESTIMATION',contexts={'M0':dict(weather=False,disease=False),'M1':dict(weather=False,disease=True),'M2':dict(weather=True,disease=False),'M3':dict(weather=True,disease=True)},validation=['LOSYO','LOYO','LOLO'],scope_policy='Disease-specific natural infection only; UNKNOWN never natural; Septoria natural scope identical to FULL and reused without extra fits',endpoint='Single already reported raw location-table plot endpoint; retain original decimal, no max/date/correlation selection; AO/survey semantics vary and remain explicit',selection='Phase2A frozen cohort AND Phase1 supported disease mask AND uniquely recoverable original numeric value in [0,100]; no outcome/error filtering',delta_MAE_definition='MAE_without - MAE_with; positive improves',delta_RMSE_definition='RMSE_without - RMSE_with',delta_R2_definition='R2_with - R2_without',hyperparameter_search='NONE; exact Phase2A fixed configurations; no eval_set or early stopping',bootstrap_env=2000,bootstrap_cultivar_within_env=500,bootstrap_seed=42,association='Descriptive Spearman; environment and environment-management centering used only after prediction, never model inputs',resources=dict(cpu_workers=2,catboost_threads=4,blas_threads=2),phase2c_authorized=False)
 cp=ROOT/'config/model_config.json'
 if cp.exists():
  old=json.loads(cp.read_text());config['time_frozen_utc']=old['time_frozen_utc'];assert old==config,'Frozen configuration changed'
 else:js(config,'config/model_config.json')
 # Read and receipt the actual inherited artifacts, never use Phase2A scores as paired controls.
 receipt=[]
 for p in [*sorted((A/'reports').glob('BRIWECS_Phase2A*.md')),A/'data/cohort/phaseI_yield_cohort.parquet',A/'data/splits/outer_fold_manifest.csv',A/'results/metrics/model_performance_summary.csv',A/'results/predictions/out_of_fold_predictions.parquet',P1/'audit/05_disease_yield_overlap.csv',P1/'audit/08_inoculation_phase_audit.csv',P1/'audit/10_leakage_risk_registry.csv',P1/'scripts/run_phase1.py',RAW/'scripts/pre-processing/data_cleaning.R']:
  receipt.append(dict(path=p.relative_to(P1).as_posix(),sha256=sha(p),purpose='Read-only inheritance/provenance, no old predictions used as controls'))
 csv(receipt,'config/input_receipt.csv')
 oldpred=pd.read_parquet(A/'results/predictions/out_of_fold_predictions.parquet',columns=['sample_id','validation_strategy','feature_set','model_name'])
 js(dict(old_oof_rows=len(oldpred),purpose='Existence and schema check only; no predictions reused'),'config/phase2a_oof_receipt.json');del oldpred
 cohort=pd.read_parquet(A/'data/cohort/phaseI_yield_cohort.parquet');assert len(cohort)==24991 and cohort.sample_id.is_unique
 af=json.loads((A/'config/cohort_freeze.json').read_text());assert ih(cohort.sample_id)==af['sample_id_set_sha256']
 features=pd.read_parquet(A/'data/features/F2_features.parquet');assert features.sample_id.is_unique
 inoc=pd.read_csv(P1/'audit/08_inoculation_phase_audit.csv');inoc=inoc.rename(columns={'Location':'location','Year':'year','Treatment':'management','trait':'disease','inoculation':'inoculation_status'})
 inoc=inoc[['location','year','management','disease','inoculation_status']];assert not inoc.duplicated(['location','year','management','disease']).any()
 audit=pd.read_csv(P1/'audit/audited_observations.csv').set_index('published_row');qc=[];summary=[];envs=[];cache={};frozen=[]
 for disease in ['Stripe_rust','Septoria']:
  rows=[]
  for row in cohort.itertuples():
   supported=pd.notna(getattr(row,disease));reason='NOT_PHASE1_SUPPORTED' if not supported else 'CANDIDATE';raw_value=np.nan
   if supported:
    if row.source_file not in cache:cache[row.source_file]=read_raw(RAW/'data/locations'/row.source_file)
    raw=cache[row.source_file].iloc[int(row.source_row)-2]
    raw_value=pd.to_numeric(str(raw[disease]).replace(',','.'),errors='coerce')
    # Verify the original row is the same cultivar and trial year, not merely matching row index.
    assert str(raw.BRISONr).replace('BRISONR','BRISONr')==row.cultivar_id
    assert int(raw.Year)==row.year
    assert audit.loc[row.published_row,'unit_id']==row.unit_id
    reason='VALID_ORIGINAL_VALUE' if np.isfinite(raw_value) and 0<=raw_value<=100 else 'ORIGINAL_VALUE_NOT_RECOVERABLE_OR_INVALID'
    if reason=='VALID_ORIGINAL_VALUE':assert np.isclose(np.trunc(raw_value),getattr(row,disease)),'Unexpected original/published severity mismatch'
   qc.append(dict(sample_id=row.sample_id,published_row=row.published_row,disease=disease,phase1_supported=supported,published_value=getattr(row,disease),raw_original_value=raw_value,status=reason,positive_truncated_to_zero=bool(supported and np.isfinite(raw_value) and raw_value>0 and getattr(row,disease)==0),source_file=row.source_file,source_row=row.source_row))
   if reason=='VALID_ORIGINAL_VALUE':rows.append(dict(sample_id=row.sample_id,disease_severity=float(raw_value),disease=disease,published_disease_value=getattr(row,disease)))
  selected=pd.DataFrame(rows);c=cohort.drop(columns=['Stripe_rust','Septoria','Leaf_rust','Powdery_mildew','DTR','Fusarium','TFI_complete'],errors='ignore').merge(selected,on='sample_id',validate='one_to_one').merge(inoc,on=['location','year','management','disease'],how='left',validate='many_to_one')
  c['inoculation_status']=c.inoculation_status.fillna('UNKNOWN');assert c.inoculation_status.isin(['NATURAL','ARTIFICIAL','MIXED_OR_UNCERTAIN','UNKNOWN']).all()
  c['disease_source']='Figshare_v1/data/locations/'+c.source_file+':row='+c.source_row.astype(int).astype(str)
  c['disease_qc_status']='PHASE1_SUPPORTED_AND_ORIGINAL_DECIMAL_VERIFIED';c['disease_unit']='percent leaf area (reported endpoint; QLB stripe rust AO history)'
  c['endpoint_type']=np.where(c.location.eq('QLB')&c.disease.eq('Stripe_rust'),'REPORTED_PLOT_SCORE_WITH_AO_ANTECEDENT','REPORTED_PLOT_DISEASE_SCORE')
  c['survey_date_status']='NOT_COMPLETE_PLOT_ENDPOINT_TIMING; NO PROSPECTIVE CLAIM'
  c['other_pathogen_artificial_regime']=c.qlb_inoculation_regime & c.disease.eq('Septoria')
  assert c.sample_id.is_unique and c.unit_id.is_unique
  name='stripe_rust' if disease=='Stripe_rust' else 'septoria';c.to_parquet(ROOT/f'data/{name}_paired_cohort.parquet',index=False)
  ft=c[['sample_id','disease_severity']].merge(features,on='sample_id',validate='one_to_one');ft.to_parquet(ROOT/f'data/{name}_features.parquet',index=False)
  for scope,g in [('FULL',c),('NATURAL',c[c.inoculation_status.eq('NATURAL')])]:
   summary.append(dict(disease=disease,scope=scope,n=len(g),site_years=g.site_year.nunique(),locations=g.location.nunique(),years=g.year.nunique(),cultivars=g.cultivar_id.nunique(),managements=g.management.nunique(),zeros=int(g.disease_severity.eq(0).sum()),positive=int(g.disease_severity.gt(0).sum()),zero_proportion=g.disease_severity.eq(0).mean(),positive_proportion=g.disease_severity.gt(0).mean(),minimum=g.disease_severity.min(),median=g.disease_severity.median(),maximum=g.disease_severity.max(),artificial=int(g.inoculation_status.eq('ARTIFICIAL').sum()),unknown=int(g.inoculation_status.isin(['UNKNOWN','MIXED_OR_UNCERTAIN']).sum()),status='ESTIMABLE' if g.site_year.nunique()>1 and g.disease_severity.nunique()>1 else 'NOT ESTIMABLE',training_alias='FULL' if disease=='Septoria' and scope=='NATURAL' else ''))
   for (sy,mg,ist),h in g.groupby(['site_year','management','inoculation_status']):envs.append(dict(disease=disease,scope=scope,site_year=sy,management=mg,inoculation_status=ist,n=len(h),cultivars=h.cultivar_id.nunique(),severity_mean=h.disease_severity.mean(),severity_sd=h.disease_severity.std(),yield_mean=h.observed_yield.mean(),yield_sd=h.observed_yield.std(),zeros=int(h.disease_severity.eq(0).sum()),positive=int(h.disease_severity.gt(0).sum())))
  frozen.append(dict(disease=disease,n=len(c),sample_id_sha256=ih(c.sample_id),severity_vector_sha256=hashlib.sha256(c.disease_severity.to_numpy().tobytes()).hexdigest(),cohort_file_sha256=sha(ROOT/f'data/{name}_paired_cohort.parquet'),feature_file_sha256=sha(ROOT/f'data/{name}_features.parquet')))
 csv(qc,'data/disease_label_qc.csv');su=csv(summary,'data/disease_cohort_summary.csv');csv(envs,'data/disease_environment_coverage.csv')
 ff=ROOT/'config/cohort_freeze.json'
 if ff.exists():assert json.loads(ff.read_text())==frozen,'Frozen disease cohort changed'
 else:js(frozen,'config/cohort_freeze.json')
 # Exact inherited environment identities, including empty test groups retained as UNSCORABLE.
 original=pd.read_csv(A/'data/splits/outer_fold_manifest.csv');original=original[original.validation_strategy.isin(config['validation'])]
 foldrows=[];assignments=[]
 for disease,scope in scopes():
  c=load_cohort(disease,scope)
  for (strategy,fid),a in original.groupby(['validation_strategy','fold_id'],sort=False):
   te=c[c.sample_id.isin(a.sample_id)];tr=c[~c.sample_id.isin(a.sample_id)]
   assert not set(tr.sample_id)&set(te.sample_id)
   groupcol={'LOSYO':'site_year','LOYO':'year','LOLO':'location'}[strategy]
   assert not set(tr[groupcol])&set(te[groupcol])
   valid=len(te)>0 and len(tr)>1
   foldrows.append(dict(disease=disease,scope=scope,validation_strategy=strategy,fold_id=fid,test_group=a.test_group.iloc[0],train_n=len(tr),test_n=len(te),train_site_years=tr.site_year.nunique(),test_site_years=te.site_year.nunique(),test_cultivars=te.cultivar_id.nunique(),unseen_test_cultivars=len(set(te.cultivar_id)-set(tr.cultivar_id)),test_disease_positive_proportion=te.disease_severity.gt(0).mean(),test_disease_unique_values=te.disease_severity.nunique(),test_yield_mean=te.observed_yield.mean(),test_yield_sd=te.observed_yield.std(),test_yield_min=te.observed_yield.min(),test_yield_max=te.observed_yield.max(),train_id_sha256=ih(tr.sample_id),test_id_sha256=ih(te.sample_id),status='READY' if valid else 'UNSCORABLE',reason='' if valid else 'NO_PAIRED_TEST_OBSERVATIONS',evidence_strength='LOW_GROUP_COUNT: 3 harvest years' if strategy=='LOYO' else 'LIMITED: only 4-5 locations' if strategy=='LOLO' else 'LIMITED_SITE_YEAR_COUNT'))
   assignments.extend(dict(disease=disease,scope=scope,validation_strategy=strategy,fold_id=fid,sample_id=s,site_year=sy,role='TEST') for s,sy in zip(te.sample_id,te.site_year))
 folds=csv(foldrows,'splits/disease_outer_fold_manifest.csv');csv(assignments,'splits/disease_test_sample_assignments.csv')
 csv([dict(disease=d,source='Official Figshare v1 raw location tables',unit='percent leaf area',body_part='leaf',endpoint='Original reported plot endpoint; no date-specific endpoint chosen',repeated_measurement_policy='One sample_id per plot, no duplicated disease dates or raw time-series rows added',precision='Original decimals restored by exact source row; integer-cleaned values retained only for QC',cross_site_semantics='QLB stripe-rust antecedent uses AO; other sites may use survey endpoint; not assumed fully harmonized',date='Partial management calendars, no per-plot synchronized observation cutoff',inoculation='Phase1 disease-specific audit; QLB NF stripe rust ARTIFICIAL, Septoria disease-specific NATURAL with other-pathogen regime flagged') for d in ['Stripe_rust','Septoria']],'data/disease_data_dictionary.csv')
 gate=pd.read_csv(A/'data/features/F2_temporal_legality.csv');assert gate.window_supported.all();csv(gate,'config/inherited_F2_temporal_legality.csv')
 (ROOT/'reports/01_Disease_Cohort_Freeze.md').write_text(f'''# Disease cohort freeze

Frozen before fitting, official Figshare v1 and Phase2A primary cohort only. All inherited raw and delivered files were SHA256 checked. Every selected plot severity independently recovered from original source_file/source_row and checked against the audited identity and published truncated value. No missing-to-zero labels admitted, no plot duplicated, no outlier/error-based trimming, no time point selected for performance.

{table(su)}
Stripe-rust natural scope uses only explicit NATURAL; unknown never natural. Septoria FULL equals its disease-specific NATURAL scope, so the identical fits are reused rather than trained twice; QLB NF has other deliberately introduced pathogens, separately flagged. These labels do not prove an entirely unmanipulated microbial environment.

Endpoint is the previously reported raw plot score, with original precision restored. Measurement timing/summary semantics are incompletely synchronized. QLB stripe-rust AO antecedent is explicit; there is no invented common assessment date. All M0–M3 are RETROSPECTIVE ASSOCIATION-ORIENTED ESTIMATION, never pre-sowing or pre-disease forecasts.

F1/F2 are unchanged from Phase2A, but every control is retrained on exactly the same disease/scope training set as its augmented model. F2 keeps only Jan1–Mar31 weather; actual disease timing may be later/unsynchronized. No yield-derived covariates or environment/target means enter predictors.

{table(folds.groupby(['disease','scope','validation_strategy','status']).size().rename('folds').reset_index())}
All 15 primary environment identities remain in the manifest; disease-empty groups are UNSCORABLE with reasons, not bad-performance exclusions. Low group counts limit secondary validation evidence.
''',encoding='utf-8')
 print(su.to_string(index=False),flush=True);print(folds.groupby('status').size().to_string(),flush=True)
if __name__=='__main__':prepare()
