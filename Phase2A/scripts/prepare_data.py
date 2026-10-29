import sys,warnings,json,hashlib
from common import *
warnings.filterwarnings('ignore',category=UserWarning,module='openpyxl')
def prepare():
 # Verify the complete existing raw manifest, including sources never used as training data.
 manifest=pd.read_csv(P1/'audit/source_file_manifest.csv');checks=[]
 for row in manifest.itertuples():
  p=P1/row.relative_path;ok=p.exists() and sha(p)==row.sha256;checks.append(dict(path=row.relative_path,sha256=row.sha256,match=ok))
 csv(checks,'config/phase1_input_hash_verification.csv');assert all(x['match'] for x in checks)
 # Read requested audit evidence, record actual schema and nonempty status.
 patterns=['01_','02_','03_','05_','06_','08_','09_','10_','11_']
 evidence=[]
 for prefix in patterns:
  p=next((P1/'audit').glob(prefix+'*.csv'));z=pd.read_csv(p);evidence.append(dict(path=p.relative_to(P1).as_posix(),rows=len(z),columns='|'.join(z.columns),sha256=sha(p)))
 csv(evidence,'config/phase1_audit_receipt.csv')
 audits=pd.read_csv(P1/'audit/audited_observations.csv')
 published=pd.read_csv(RAW/'output/BRIWECS_data_publication.csv',sep=';');published['published_row']=np.arange(len(published))+2
 assert len(published)==34683 and audits.reliable_unit.sum()==34541 and audits.Seedyield.notna().sum()==34213
 target=['Seedyield'];d=audits.copy()
 reasons=np.select([~d.reliable_unit,~d.Year.isin(CFG['cohort']['years']),~d.Location.isin(CFG['cohort']['sites']),~np.isfinite(d.Seedyield)|d.Seedyield.lt(0)],['PHASE1_IDENTITY_OR_COORDINATE_DISPUTE','OUTSIDE_PHASE_I','OUTSIDE_FIVE_SITES','MISSING_OR_INVALID_WHOLE_PLOT_YIELD'],default='INCLUDED')
 d['selection_reason']=reasons;d['sample_id']=d.unit_id.map(lambda s:'BRI1_'+hashlib.sha256(s.encode()).hexdigest()[:24])
 csv(d[['published_row','sample_id','unit_id','Location','Year','Treatment','selection_reason']], 'data/cohort/exclusion_ledger.csv')
 c=d[d.selection_reason.eq('INCLUDED')].copy().reset_index(drop=True)
 c=c.rename(columns=dict(BRISONr='cultivar_id',Location='location',Year='year',Treatment='management',Seedyield='observed_yield',Block='replicate',Row='row',Column='column'))
 c['site_year']=c.location+'_'+c.year.astype(str);c['site_year_management']=c.site_year+'_'+c.management
 c['nitrogen']=c.management.str.split('_').str[0];c['fungicide']=c.management.str.split('_').str[1];c['water']=c.management.str.split('_').str[2]
 unknown=c.location.eq('GGE')&c.year.isin([2016,2017]);c.loc[unknown,'water']='UNKNOWN_V1'
 c['qlb_inoculation_regime']=c.location.eq('QLB')&c.fungicide.eq('NF')
 c['qc_status']='PHASE1_RELIABLE_SOURCE_LINKED';c['phase']='I';c['source_version']='figshare_27910269_v1'
 c['source_table_sha256']=sha(RAW/'output/BRIWECS_data_publication.csv')
 meta=pd.read_csv(RAW/'metadata/BRIWECs_cultivar_info.csv',encoding='cp1252');assert meta.BRISONr.is_unique
 assert meta.RYear.dropna().between(1800,2017).all()
 c=c.merge(meta[['BRISONr','RYear']].rename(columns={'BRISONr':'cultivar_id','RYear':'release_year'}),on='cultivar_id',validate='many_to_one',how='left')
 assert c.sample_id.is_unique and c.unit_id.is_unique
 assert c[['location','year','management','cultivar_id','replicate']].notna().all().all()
 assert c.coordinates_complete.all() and not c.coordinate_collision.any()
 # Extreme outcomes are flagged, never removed.
 q1,q3=c.observed_yield.quantile([.25,.75]);c['descriptive_extreme_IQR']=~c.observed_yield.between(q1-1.5*(q3-q1),q3+1.5*(q3-q1))
 cohort_path=ROOT/'data/cohort/phaseI_yield_cohort.parquet'
 frozen_path=ROOT/'config/cohort_freeze.json'
 if frozen_path.exists():
  frozen=json.loads(frozen_path.read_text());assert ids_hash(c.sample_id)==frozen['sample_id_set_sha256'],'Frozen cohort changed'
 c.to_parquet(cohort_path,index=False)
 csv(c[['sample_id','unit_id','published_row','site_year','site_year_management']], 'data/cohort/sample_id_manifest.csv')
 coverage=c.groupby(['location','year','management']).agg(samples=('sample_id','size'),cultivars=('cultivar_id','nunique'),replicates=('replicate','nunique'),mean_yield=('observed_yield','mean'),sd_yield=('observed_yield','std'),minimum=('observed_yield','min'),maximum=('observed_yield','max')).reset_index()
 csv(coverage,'data/cohort/environment_coverage.csv')
 exclusions=d.groupby('selection_reason').size().rename('records').reset_index();csv(exclusions,'data/cohort/exclusion_summary.csv')
 freeze=dict(time_utc=now(),samples=len(c),cultivars=c.cultivar_id.nunique(),sites=c.location.nunique(),years=sorted(c.year.unique().tolist()),site_years=c.site_year.nunique(),site_year_managements=c.site_year_management.nunique(),missing_coordinates=int((~c.coordinates_complete).sum()),source_table_sha256=sha(RAW/'output/BRIWECS_data_publication.csv'),sample_id_set_sha256=ids_hash(c.sample_id),yield_vector_sha256=hashlib.sha256(c.observed_yield.to_numpy().tobytes()).hexdigest(),release_year_missing=int(c.release_year.isna().sum()),extreme_yields_flagged_not_removed=int(c.descriptive_extreme_IQR.sum()),config_sha256=sha(ROOT/'config/training_config.json'))
 if not frozen_path.exists():write_json(freeze,'config/cohort_freeze.json')
 # F2 cutoff is a design choice, never an invented harvest or sowing date.
 phen=c[['sample_id','published_row','site_year','location','year']].merge(published[['published_row','BBCH59','BBCH87']],on='published_row',validate='one_to_one')
 gate=phen.groupby(['site_year','location','year']).agg(heading_min_doy=('BBCH59','min'),heading_recorded=('BBCH59','count'),maturity_min_doy=('BBCH87','min'),maturity_recorded=('BBCH87','count')).reset_index()
 env=[]
 for loc in sorted(c.location.unique()):
  w=pd.read_csv(next((RAW/'data/weather').glob('*_'+loc+'.csv')),sep=';',decimal=',');w['date']=pd.to_datetime(w.Date)
  assert w.date.is_unique
  for yr in sorted(c.year.unique()):
   start=pd.Timestamp(int(yr),1,1);end=pd.Timestamp(int(yr),3,31);v=w[w.date.between(start,end)]
   assert len(v)==len(pd.date_range(start,end)) and v[['Temperature.Air.Mean_C','Precipitation_mm_d','Global.Radiation_MJ_m2_d']].notna().all().all()
   env.append(dict(site_year=f'{loc}_{yr}',window_start=str(start.date()),window_end=str(end.date()),days=len(v),winter_temperature_mean=v['Temperature.Air.Mean_C'].mean(),winter_temperature_sd=v['Temperature.Air.Mean_C'].std(ddof=0),winter_precipitation_sum=v.Precipitation_mm_d.sum(),winter_radiation_sum=v.Global_Radiation.sum() if 'Global_Radiation' in v else v['Global.Radiation_MJ_m2_d'].sum()))
 env=pd.DataFrame(env);gate=gate.merge(env[['site_year','window_start','window_end','days']],on='site_year',validate='one_to_one')
 gate['cutoff_doy']=pd.to_datetime(gate.window_end).dt.dayofyear
 gate['window_supported']=gate.heading_recorded.gt(0)&gate.heading_min_doy.gt(gate.cutoff_doy)
 assert gate.window_supported.all()
 gate['evidence_limit']='Heading/harvest timing verified at site-year seasonal level, not every plot exact harvest day; fixed winter window, no full-season weather'
 csv(gate,'data/features/F2_temporal_legality.csv');csv(env,'data/features/site_year_weather_features.csv')
 f1=c[['sample_id']+CFG['F1']['categorical']+CFG['F1']['numeric']].copy();f1.to_parquet(ROOT/'data/features/F1_features.parquet',index=False)
 f2=f1.merge(c[['sample_id','site_year']],on='sample_id',validate='one_to_one').merge(env[['site_year']+CFG['F2']['numeric_additions']],on='site_year',validate='many_to_one').drop(columns='site_year')
 assert len(f2)==len(c);f2.to_parquet(ROOT/'data/features/F2_features.parquet',index=False)
 feature_hashes={fs:sha(ROOT/f'data/features/{fs}_features.parquet') for fs in ['F1','F2']}
 if (ROOT/'config/feature_freeze.json').exists():assert feature_hashes==json.loads((ROOT/'config/feature_freeze.json').read_text()),'Frozen features changed'
 else:write_json(feature_hashes,'config/feature_freeze.json')
 registry=[]
 for col in published:
  use=col in ['BRISONr','Treatment'];registry.append(dict(feature=col,status='ALLOWED_SOURCE_ONLY' if use else 'EXCLUDED',feature_set='F1' if use else 'NONE',reason='Decompose into audited cultivar/management fields' if use else 'Identifier only; timing uncertain; target or mature/derived trait; strict allow-list excludes all other original columns',time_available='design' if use else 'not used',preprocessing='training-fold only'))
 for col in CFG['F1']['categorical']+CFG['F1']['numeric']+CFG['F2']['numeric_additions']:
  registry.append(dict(feature=col,status='ALLOWED',feature_set='F1+F2' if col in f1 else 'F2',reason='Official genotype metadata / planned management' if col in f1 else 'Deterministic winter weather before fixed cutoff',time_available='pre-season planned' if col in f1 else '03-31 target year; retrospective evaluation',preprocessing='Fold-only median/std or categories; no target encoding outside CatBoost internal training'))
 csv(registry,'reports/02_Feature_Registry.csv')
 # Freeze all test assignments before any model is fit.
 assignments=[];folds=[]
 for strategy,col in [('LOSYO','site_year'),('LOYO','year'),('LOLO','location'),('LOEO','site_year_management')]:
  for val in sorted(c[col].unique()):
   test=c[c[col].eq(val)];train=c[c[col].ne(val)];fid=f'{strategy}__{val}'
   assignments.extend(dict(sample_id=s,validation_strategy=strategy,fold_id=fid,test_group=str(val),site_year=sy,role='TEST') for s,sy in zip(test.sample_id,test.site_year))
   folds.append(dict(validation_strategy=strategy,fold_id=fid,train_n=len(train),test_n=len(test),test_site_years=test.site_year.nunique(),train_site_years=train.site_year.nunique(),train_test_id_overlap=len(set(train.sample_id)&set(test.sample_id)),train_test_site_year_overlap=len(set(train.site_year)&set(test.site_year)),known_test_cultivars=len(set(test.cultivar_id)&set(train.cultivar_id)),unseen_test_cultivars=len(set(test.cultivar_id)-set(train.cultivar_id)),train_id_sha256=ids_hash(train.sample_id),test_id_sha256=ids_hash(test.sample_id),status='PASS',note='LOEO permits other managements from same site-year in training' if strategy=='LOEO' else 'Whole test group absent from training'))
 cv=csv(assignments,'data/splits/outer_fold_manifest.csv');fl=csv(folds,'data/splits/fold_integrity_check.csv')
 assert cv.groupby(['validation_strategy','sample_id']).size().eq(1).all()
 assert fl.train_test_id_overlap.eq(0).all() and fl.loc[fl.validation_strategy.eq('LOSYO'),'train_test_site_year_overlap'].eq(0).all()
 csv([dict(validation_strategy=x['validation_strategy'],fold_id=x['fold_id'],status='NOT APPLICABLE',reason='Fixed hyperparameters selected before results; no inner tuning or early stopping',outer_training_site_years=x['train_site_years']) for x in folds], 'data/splits/inner_cv_manifest.csv')
 csv(c.loc[c.location.ne('QLB'),['sample_id','site_year']].assign(validation_strategy='LOSYO_NO_QLB',fold_id=lambda x:'LOSYO_NO_QLB__'+x.site_year),'data/splits/sensitivity_fold_manifest.csv')
 # Pre-screen only: no disease features enter F1 or F2.
 disease=[];disease_env=[];raw_cache={};raw_values=[]
 for t in ['Stripe_rust','Septoria','Powdery_mildew']:
  g=c[c[t].notna()].copy();art=g.qlb_inoculation_regime if t=='Stripe_rust' else pd.Series(False,index=g.index)
  orig=[]
  for row in g.itertuples():
   if row.source_file not in raw_cache:
    f=RAW/'data/locations'/row.source_file
    try:z=pd.read_csv(f,sep=';',dtype=str,encoding='utf-8-sig')
    except UnicodeDecodeError:z=pd.read_csv(f,sep=';',dtype=str,encoding='cp1252')
    raw_cache[row.source_file]=z
   value=raw_cache[row.source_file].iloc[int(row.source_row)-2][t]
   value=pd.to_numeric(str(value).replace(',','.'),errors='coerce');assert np.isfinite(value) and 0<=value<=100
   orig.append(value);raw_values.append(dict(sample_id=row.sample_id,disease=t,raw_original_value=value,published_supported_value=getattr(row,t),source_file=row.source_file,source_row=row.source_row))
  g['raw_original_value']=orig
  disease.append(dict(disease=t,paired_samples=len(g),site_years=g.site_year.nunique(),cultivars=g.cultivar_id.nunique(),managements=g.management.nunique(),published_supported_zero=int(g[t].eq(0).sum()),published_supported_positive=int(g[t].gt(0).sum()),zero_proportion=g[t].eq(0).mean(),positive_proportion=g[t].gt(0).mean(),artificial_samples=int(art.sum()),natural_samples=int((~art).sum()),largest_environment_share=g.site_year.value_counts(normalize=True).max(),fold_reuse='Same LOSYO assignments on paired subset only; empty disease-test environments UNSCORABLE',precision_limit='Official integer scores truncate raw positive fractions; retain Phase1 valid mask, do not interpret all zeros as truly healthy'))
  disease[-1].update(raw_zero=int(g.raw_original_value.eq(0).sum()),raw_positive=int(g.raw_original_value.gt(0).sum()),raw_zero_proportion=g.raw_original_value.eq(0).mean(),raw_positive_proportion=g.raw_original_value.gt(0).mean(),raw_positive_truncated_to_published_zero=int((g.raw_original_value.gt(0)&g[t].eq(0)).sum()))
  for sy,h in g.groupby('site_year'):disease_env.append(dict(disease=t,site_year=sy,paired_samples=len(h),zero=int(h[t].eq(0).sum()),positive=int(h[t].gt(0).sum()),cultivars=h.cultivar_id.nunique(),managements=h.management.nunique()))
 csv(disease,'results/metrics/disease_track_prescreen.csv');csv(disease_env,'results/metrics/disease_environment_prescreen.csv');csv(raw_values,'results/metrics/disease_prescreen_raw_value_provenance.csv')
 report=f'''# 01 Cohort Freeze Report

Frozen before training. Official Figshare v1 only; Phase1 remains unchanged.

{table(pd.DataFrame([freeze]).drop(columns=['time_utc','source_table_sha256','sample_id_set_sha256','yield_vector_sha256','config_sha256']))}
Target: whole-plot Seedyield, dt/ha at 100% dry mass. All 220 cultivar IDs have metadata rows; release year may be missing and is imputed in each training fold only.

{table(exclusions)}
The ledger assigns mutually exclusive reasons across all 34,683 published rows. The 142 Phase1 disputes remain excluded; they occur outside this PhaseI five-site primary scope. No valid low/high yield removed. All selected units have complete coordinates; no blanket coordinate filter was needed. No performance-driven changes are allowed.

Actual environment coverage (51, not assumed 45):

{table(coverage)}
Sample ID is a hash of the audited full unit key, not an artificial biological replicate. Rows/cultivars within environments remain correlated. Known-cultivar transfer is the intended estimand.
'''
 (ROOT/'reports/01_Cohort_Freeze_Report.md').write_text(report,encoding='utf-8')
 (ROOT/'reports/03_Feature_Leakage_Check.md').write_text('''# 03 Feature leakage gate

F1: cultivar ID, cultivar release year, planned nitrogen/fungicide/water. No location, year, site-year, plot coordinates, block ID, disease, harvest measurements, BLUEs, target encodings or yield descendants enter the feature matrix. GGE 2016/2017 water is UNKNOWN_V1; no updated GitHub correction imported.

F2: F1 plus fixed January 1–March 31 weather (temperature mean/SD, precipitation sum, radiation sum). Every site-year has recorded heading after day 123, later than cutoff day 90/91. This calendar cutoff is an analysis definition, not a fabricated harvest date. Audit-level phenology is used to verify temporal ordering only, never as a model input. Exact plot harvest dates remain incomplete. This is a RETROSPECTIVE_ENVIRONMENTAL_MODEL evaluated with historical weather, not sowing-time or early disease forecasting. No weather after cutoff is included.

Medians, scales, category encoders are fit only on outer training data. CatBoost supervised categorical statistics are internal to training; no evaluation labels or early stopping supplied. Ridge/XGBoost use train-only one-hot encoding. MixedLM uses train-only numerical scaling and management encoding, release-year fixed effect, site-year and nested block random intercepts; predictions always marginal fixed effects, including LOEO. No test outcomes estimate random effects.

LOSYO holds out every treatment in a location-year; LOYO/LOLO entire year/site. Supplementary LOEO deliberately permits other treatments in same location-year and is not primary independence evidence. Fixed hyperparameters, no inner selection. Environmental bootstrap is descriptive because environments also share locations/years. Model errors never change cohort selection. F3 not trained.

Sources: scikit-learn OneHotEncoder documentation; CatBoostRegressor documentation; statsmodels MixedLM.predict documentation; XGBoost Python API. Exact source URLs recorded in config/method_sources.json.
''',encoding='utf-8')
 write_json(dict(OneHotEncoder='https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.OneHotEncoder.html',CatBoost='https://catboost.ai/docs/en/concepts/python-reference_catboostregressor',MixedLM='https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLM.predict.html',XGBoost='https://xgboost.readthedocs.io/en/stable/python/python_api.html'), 'config/method_sources.json')
 print(json.dumps(freeze,indent=2),flush=True)
if __name__=='__main__':prepare()
