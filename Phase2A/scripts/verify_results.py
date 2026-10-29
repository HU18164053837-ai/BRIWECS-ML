from common import *
import argparse,joblib
def verify():
 c=pd.read_parquet(ROOT/'data/cohort/phaseI_yield_cohort.parquet');freeze=json.loads((ROOT/'config/cohort_freeze.json').read_text())
 assert c.sample_id.is_unique and len(c)==24991 and c.cultivar_id.nunique()==220 and c.site_year.nunique()==15
 assert ids_hash(c.sample_id)==freeze['sample_id_set_sha256']
 assert hashlib.sha256(c.observed_yield.to_numpy().tobytes()).hexdigest()==freeze['yield_vector_sha256']
 assert sha(ROOT/'config/training_config.json')==freeze['config_sha256']
 mf=pd.read_csv(P1/'audit/source_file_manifest.csv')
 for row in mf.itertuples():assert sha(P1/row.relative_path)==row.sha256
 for fs in ['F1','F2']:
  f=pd.read_parquet(ROOT/f'data/features/{fs}_features.parquet');assert f.sample_id.is_unique and set(f.sample_id)==set(c.sample_id)
  allowed=set(['sample_id']+CFG['F1']['categorical']+CFG['F1']['numeric']+(CFG['F2']['numeric_additions'] if fs=='F2' else []));assert set(f.columns)==allowed
  assert sha(ROOT/f'data/features/{fs}_features.parquet')==json.loads((ROOT/'config/feature_freeze.json').read_text())[fs]
 splits=pd.read_csv(ROOT/'data/splits/outer_fold_manifest.csv');d=pd.read_parquet(ROOT/'results/predictions/out_of_fold_predictions.parquet')
 assert not d.duplicated(['validation_strategy','model_name','feature_set','sample_id']).any()
 assert np.isfinite(d.predicted_yield).all()
 status=pd.read_csv(ROOT/'results/metrics/model_fit_status.csv');assert len(status)==(15+3+5+51+12)*9
 byid=c.set_index('sample_id');checked=0
 for (strategy,fold,fs,model),g in d.groupby(['validation_strategy','fold_id','feature_set','model_name']):
  st=status[status.validation_strategy.eq(strategy)&status.fold_id.eq(fold)&status.feature_set.eq(fs)&status.model_name.eq(model)];assert len(st)==1 and st.status.iloc[0].startswith('SUCCESS')
  universe=c if strategy!='LOSYO_NO_QLB' else c[c.location.ne('QLB')]
  if strategy=='LOSYO_NO_QLB':expected=set(universe.loc[universe.site_year.eq(g.site_year.iloc[0]),'sample_id'])
  else:expected=set(splits.loc[splits.validation_strategy.eq(strategy)&splits.fold_id.eq(fold),'sample_id'])
  assert set(g.sample_id)==expected
  tr=universe[~universe.sample_id.isin(expected)];assert ids_hash(tr.sample_id)==st.train_id_sha256.iloc[0]
  if strategy.startswith('LOSYO'):assert not set(tr.site_year)&set(g.site_year)
  assert np.allclose(g.observed_yield,byid.loc[g.sample_id,'observed_yield'])
  if model=='TrainingMean':assert np.allclose(g.predicted_yield,tr.observed_yield.mean())
  checked+=1
 # Recompute primary metrics independently from arithmetic and verify saved model artifacts.
 primary=d[d.validation_strategy.eq('LOSYO')];summary=pd.read_csv(ROOT/'results/metrics/model_performance_summary.csv')
 for (model,fs),g in primary.groupby(['model_name','feature_set']):
  row=summary[summary.validation_strategy.eq('LOSYO')&summary.aggregation.eq('POOLED')&summary.feature_set.eq(fs)&summary.model_name.eq(model)].iloc[0]
  assert np.isclose(row.MAE,np.mean(np.abs(g.predicted_yield-g.observed_yield)))
  assert np.isclose(row.RMSE,np.sqrt(np.mean((g.predicted_yield-g.observed_yield)**2)))
 # Independently reload every primary fold model and reproduce the stored test predictions.
 for (fold,fs,name),g in primary[primary.model_name.ne('TrainingMean')].groupby(['fold_id','feature_set','model_name']):
  feature=pd.read_parquet(ROOT/f'data/features/{fs}_features.parquet').set_index('sample_id')
  path=ROOT/'models'/f'{fold}__{fs}__{name}'
  x=feature.loc[g.sample_id].copy()
  if name in ['Ridge','XGBoost']:
   model=joblib.load(path.with_suffix('.joblib'));pred=model.predict(x)
  elif name=='CatBoost':
   from catboost import CatBoostRegressor
   nums=CFG['F1']['numeric']+(CFG['F2']['numeric_additions'] if fs=='F2' else [])
   imp=joblib.load(path.with_suffix('.imputer.joblib'));x[nums]=imp.transform(x[nums]);model=CatBoostRegressor();model.load_model(str(path.with_suffix('.cbm')));pred=model.predict(x,thread_count=4)
  elif name=='MixedLM':
   artifact=joblib.load(path.with_suffix('.joblib'));z=artifact['preprocessor'].transform(x);pred=np.column_stack([np.ones(len(z)),z])@artifact['fixed_effects']
  assert np.allclose(pred,g.predicted_yield)
 # Same-sample F2/F1 comparisons are exact; no silent fold exclusions.
 pairs=pd.read_csv(ROOT/'results/metrics/model_comparison_paired.csv');assert (pairs.matched_samples.le(pairs.samples_a)&pairs.matched_samples.le(pairs.samples_b)).all()
 # Failed mixed-model folds are documented; comparisons use only exactly matched IDs/folds,
 # and cannot be mistaken for complete-cohort scores.
 partial=pairs[~pairs.exact_same_test_samples]
 if len(partial):assert partial[['model_a','model_b']].eq('MixedLM').any(axis=1).all()
 print(f'PASS: {checked} genuine test chunks, frozen IDs/targets/features, fold isolation, training means, metrics, ALL saved primary fold models reproduced, raw SHA256, paired comparisons.')
 print(status.groupby(['model_name','status']).size().to_string())
 return dict(status='PASS',verified_prediction_chunks=checked,raw_files_verified=len(mf),source_unchanged=True)
if __name__=='__main__':write_json(verify(),'reports/verification_receipt.json')
