from common import *
import joblib
from catboost import CatBoostRegressor
def verify():
 config=cfg();confighash=sha(ROOT/'config/model_config.json');frozen=json.loads((ROOT/'config/cohort_freeze.json').read_text());cohorts={}
 for row in frozen:
  name='stripe_rust' if row['disease']=='Stripe_rust' else 'septoria';c=pd.read_parquet(ROOT/f'data/{name}_paired_cohort.parquet');assert c.sample_id.is_unique and c.unit_id.is_unique
  assert ih(c.sample_id)==row['sample_id_sha256'] and hashlib.sha256(c.disease_severity.to_numpy().tobytes()).hexdigest()==row['severity_vector_sha256']
  assert sha(ROOT/f'data/{name}_paired_cohort.parquet')==row['cohort_file_sha256'] and sha(ROOT/f'data/{name}_features.parquet')==row['feature_file_sha256']
  assert c.disease_severity.between(0,100).all() and c.disease_qc_status.eq('PHASE1_SUPPORTED_AND_ORIGINAL_DECIMAL_VERIFIED').all();cohorts[row['disease']]=c
 for row in pd.read_csv(ROOT/'config/inherited_hash_check.csv').itertuples():assert sha(P1/row.path)==row.sha256,'Phase1/2A source changed'
 d=pd.read_parquet(ROOT/'results/predictions/disease_pilot_oof_predictions.parquet');status=pd.read_csv(ROOT/'results/metrics/model_fit_status.csv');folds=pd.read_csv(ROOT/'splits/disease_outer_fold_manifest.csv');assign=pd.read_csv(ROOT/'splits/disease_test_sample_assignments.csv')
 assert len(status)==len(folds)*8 and status.config_sha256.eq(confighash).all()
 assert not status.job.duplicated().any() and not d.duplicated(['disease','scope','validation_strategy','fold_id','model_name','context','sample_id']).any()
 assert np.isfinite(d.predicted_yield).all();count=0;reloaded=0
 for vals,g in d.groupby(['disease','scope','validation_strategy','fold_id','model_name','context']):
  disease,scope,strategy,fid,model,context=vals;c=load_cohort(disease,scope);a=assign[assign.disease.eq(disease)&assign.scope.eq(scope)&assign.fold_id.eq(fid)];expected=set(a.sample_id)
  assert set(g.sample_id)==expected
  tr=c[~c.sample_id.isin(expected)];col={'LOSYO':'site_year','LOYO':'year','LOLO':'location'}[strategy]
  assert not set(tr[col])&set(g[col]);st=status[status.job.eq(f'{disease}__{scope}__{fid}__{model}__{context}')].iloc[0]
  assert st.status.startswith('SUCCESS') and st.train_id_sha256==ih(tr.sample_id) and st.test_id_sha256==ih(g.sample_id)
  ref=c.set_index('sample_id').loc[g.sample_id];assert np.allclose(g.observed_yield,ref.observed_yield) and np.allclose(g.disease_severity,ref.disease_severity)
  if scope=='NATURAL':assert g.inoculation_status.eq('NATURAL').all()
  if strategy=='LOSYO':
   artifact=ROOT/'models'/f'{disease}__{scope}__{fid}__{model}__{context}';name='stripe_rust' if disease=='Stripe_rust' else 'septoria';f=pd.read_parquet(ROOT/f'data/{name}_features.parquet').set_index('sample_id');ct=config['contexts'][context];cats=config['categorical'];nums=config['F1_numeric']+(config['weather'] if ct['weather'] else [])+(['disease_severity'] if ct['disease'] else []);x=f.loc[g.sample_id,cats+nums].copy()
   if model=='Ridge':
    fitted=joblib.load(artifact.with_suffix('.joblib'));assert list(fitted.feature_names_in_)==cats+nums;pred=fitted.predict(x)
   else:
    imp=joblib.load(artifact.with_suffix('.imputer.joblib'));x[nums]=imp.transform(x[nums]);fitted=CatBoostRegressor();fitted.load_model(str(artifact.with_suffix('.cbm')));assert fitted.feature_names_==cats+nums;pred=fitted.predict(x,thread_count=4)
   assert np.allclose(pred,g.predicted_yield);reloaded+=1
  count+=1
 integrity=pd.read_csv(ROOT/'results/metrics/paired_prediction_integrity.csv');assert integrity.train_ids_identical.all();assert integrity.test_ids_identical.all(),'Pair coverage differs; must investigate failures'
 # Independent direct paired arithmetic validates delta sign and no Phase2A control reuse.
 perf=pd.read_csv(ROOT/'results/metrics/disease_incremental_performance.csv')
 for vals,g in d.groupby(['disease','scope','validation_strategy','model_name']):
  disease,scope,strategy,model=vals
  for without,withd in [('M0','M1'),('M2','M3')]:
   b=g[g.context.eq(without)];v=g[g.context.eq(withd)];q=b.merge(v[['sample_id','fold_id','predicted_yield']],on=['sample_id','fold_id'],suffixes=('_b','_v'),validate='one_to_one')
   delta=np.abs(q.predicted_yield_b-q.observed_yield).mean()-np.abs(q.predicted_yield_v-q.observed_yield).mean();row=perf[perf.disease.eq(disease)&perf.scope.eq(scope)&perf.validation_strategy.eq(strategy)&perf.model_name.eq(model)&perf.comparison.eq(f'{withd}_vs_{without}')&perf.aggregation.eq('POOLED')].iloc[0];assert np.isclose(delta,row.delta_MAE)
 result=dict(status='PASS',real_prediction_chunks=count,primary_models_reloaded=reloaded,task_states=len(status),successful_fits=int(status.status.str.startswith('SUCCESS').sum()),unscorable=int(status.status.eq('UNSCORABLE').sum()),failed_fits=int(status.status.eq('FAILED').sum()),source_artifacts_unchanged=True,paired_train_test_identity=True,original_decimal_severity_verified=True)
 js(result,'reports/verification_receipt.json');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':verify()
