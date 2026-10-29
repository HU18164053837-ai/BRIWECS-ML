import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('MKL_NUM_THREADS','2')
from common import *
import argparse,time,traceback,warnings,joblib
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from catboost import CatBoostRegressor
from xgboost import XGBRegressor
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

def preprocessing(cats,nums,sparse=True,drop=None):
 return ColumnTransformer([('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=sparse,drop=drop),cats),('num',Pipeline([('impute',SimpleImputer(strategy='median',keep_empty_features=True)),('scale',StandardScaler())]),nums)],sparse_threshold=1.0 if sparse else 0.0)

def fit_predict(name,fs,train,test,features,model_dir):
 cats=CFG['F1']['categorical'];nums=CFG['F1']['numeric']+(CFG['F2']['numeric_additions'] if fs=='F2' else [])
 if name=='TrainingMean':return np.full(len(test),train.observed_yield.mean()),dict(training_mean=float(train.observed_yield.mean()))
 Xtr=features.loc[train.sample_id,cats+nums].copy();Xte=features.loc[test.sample_id,cats+nums].copy();y=train.observed_yield.to_numpy()
 if name=='CatBoost':
  imp=SimpleImputer(strategy='median',keep_empty_features=True);Xtr[nums]=imp.fit_transform(Xtr[nums]);Xte[nums]=imp.transform(Xte[nums]);model=CatBoostRegressor(**CFG['models'][name]);model.fit(Xtr,y,cat_features=cats)
  pred=model.predict(Xte)
  if model_dir is not None:model.save_model(str(model_dir.with_suffix('.cbm')));joblib.dump(imp,model_dir.with_suffix('.imputer.joblib'))
  return pred,dict(iterations=model.tree_count_,numeric_medians=imp.statistics_.tolist())
 if name in ['Ridge','XGBoost']:
  prep=preprocessing(cats,nums);est=Ridge(**CFG['models'][name]) if name=='Ridge' else XGBRegressor(**CFG['models'][name]);model=Pipeline([('preprocess',prep),('model',est)]);model.fit(Xtr,y);pred=model.predict(Xte)
  if model_dir is not None:joblib.dump(model,model_dir.with_suffix('.joblib'),compress=3)
  return pred,dict(numeric_medians=prep.named_transformers_['num'].named_steps['impute'].statistics_.tolist(),encoded_dimensions=len(prep.get_feature_names_out()))
 if name=='MixedLM':
  # Lower-dimensional structural baseline, explicitly no cultivar ID fixed effect.
  mcats=['nitrogen','fungicide','water'];prep=preprocessing(mcats,nums,sparse=False,drop='first');A=prep.fit_transform(Xtr);B=prep.transform(Xte)
  xx=pd.DataFrame(A,columns=[f'x{i}' for i in range(A.shape[1])]);xx['outcome']=y;xx['site_year']=train.site_year.to_numpy();xx['management_block']=(train.management+'|'+train.replicate.astype(str)).to_numpy()
  formula='outcome ~ '+' + '.join(f'x{i}' for i in range(A.shape[1]));model=sm.MixedLM.from_formula(formula,groups='site_year',re_formula='1',vc_formula={'management_block':'0+C(management_block)'},data=xx,use_sparse=True)
  fit=model.fit(reml=True,method='lbfgs',maxiter=100,disp=False)
  if not fit.converged:raise RuntimeError('MixedLM did not converge under predefined lbfgs/maxiter=100; no test-informed fallback')
  pred=np.column_stack([np.ones(len(B)),B])@fit.fe_params.to_numpy()
  if not np.isfinite(pred).all():raise RuntimeError('MixedLM nonfinite marginal prediction')
  meta=dict(converged=bool(fit.converged),formula=formula,coefficients=fit.fe_params.to_dict(),site_year_random_variance=float(fit.cov_re.iloc[0,0]),block_variance=fit.vcomp.tolist(),prediction='fixed effects only',numeric_medians=prep.named_transformers_['num'].named_steps['impute'].statistics_.tolist())
  if model_dir is not None:joblib.dump(dict(preprocessor=prep,fixed_effects=fit.fe_params.to_numpy(),metadata=meta),model_dir.with_suffix('.joblib'),compress=3)
  return pred,meta
 raise ValueError(name)

def train(only=None,worker=False,collect_only=False):
 c=pd.read_parquet(ROOT/'data/cohort/phaseI_yield_cohort.parquet')
 freeze=json.loads((ROOT/'config/cohort_freeze.json').read_text());assert ids_hash(c.sample_id)==freeze['sample_id_set_sha256']
 assert sha(ROOT/'config/training_config.json')==freeze['config_sha256'],'Configuration changed after cohort freeze'
 features={fs:pd.read_parquet(ROOT/f'data/features/{fs}_features.parquet').set_index('sample_id') for fs in ['F1','F2']}
 split=pd.read_csv(ROOT/'data/splits/outer_fold_manifest.csv');sensitivity=pd.read_csv(ROOT/'data/splits/sensitivity_fold_manifest.csv');split=pd.concat([split,sensitivity],ignore_index=True)
 chunks=ROOT/'results/predictions/chunks';chunks.mkdir(exist_ok=True);statusdir=ROOT/'logs/fit_status';statusdir.mkdir(exist_ok=True)
 count=0
 with threadpool_limits(limits=2,user_api='blas'):
  for (strategy,fold),g in split.groupby(['validation_strategy','fold_id'],sort=False):
   if collect_only or (only and strategy not in only):continue
   universe=c if strategy!='LOSYO_NO_QLB' else c[c.location.ne('QLB')]
   test=universe[universe.sample_id.isin(g.sample_id)].copy();train=universe[~universe.sample_id.isin(g.sample_id)].copy()
   assert not set(train.sample_id)&set(test.sample_id)
   if strategy.startswith('LOSYO'):assert not set(train.site_year)&set(test.site_year)
   jobs=[('TrainingMean','F0')]+[(name,fs) for fs in ['F1','F2'] for name in ['Ridge','MixedLM','CatBoost','XGBoost']]
   for name,fs in jobs:
    job=f'{fold}__{fs}__{name}';statpath=statusdir/(job+'.json');predpath=chunks/(job+'.parquet');count+=1
    if statpath.exists():
     old=json.loads(statpath.read_text());assert old['config_sha256']==freeze['config_sha256'] and old['train_id_sha256']==ids_hash(train.sample_id)
     if old['status'].startswith('SUCCESS') and predpath.exists():continue
     if old['status']=='FAILED':continue # failure preserved; no performance-driven retries
    lock=statusdir/(job+'.lock')
    try:
     with lock.open('x') as claim:claim.write(str(os.getpid()))
    except FileExistsError:continue
    status=dict(job=job,model_name=name,feature_set=fs,validation_strategy=strategy,fold_id=fold,train_n=len(train),test_n=len(test),time_utc=now(),train_id_sha256=ids_hash(train.sample_id),test_id_sha256=ids_hash(test.sample_id),config_sha256=freeze['config_sha256'],cohort_id_sha256=freeze['sample_id_set_sha256'])
    t=time.monotonic();model_dir=ROOT/'models'/job if strategy=='LOSYO' else None
    print(f'START {job} train={len(train)} test={len(test)}',flush=True)
    try:
     with warnings.catch_warnings(record=True) as caught:
      warnings.simplefilter('always');pred,details=fit_predict(name,fs,train,test,features.get(fs),model_dir)
     assert len(pred)==len(test) and np.isfinite(pred).all()
     out=test[['sample_id','observed_yield','site_year','location','year','management','qlb_inoculation_regime']].copy();out['predicted_yield']=pred;out['model_name']=name;out['feature_set']=fs;out['validation_strategy']=strategy;out['fold_id']=fold
     tmp=predpath.with_suffix(f'.{os.getpid()}.tmp');out.to_parquet(tmp,index=False);tmp.replace(predpath)
     warning_text=' | '.join(dict.fromkeys(str(w.message) for w in caught));status.update(status='SUCCESS_WITH_WARNINGS' if warning_text else 'SUCCESS',warnings=warning_text,details=details)
    except Exception as exc:
     status.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc());(ROOT/'logs'/f'{job}.error.txt').write_text(status['traceback'],encoding='utf-8')
    status['seconds']=time.monotonic()-t;tempstat=statpath.with_suffix(f'.{os.getpid()}.tmp');tempstat.write_text(json.dumps(status,ensure_ascii=False,indent=2,default=str),encoding='utf-8');tempstat.replace(statpath);lock.unlink()
    print(f'END {job} {status["status"]} seconds={status["seconds"]:.2f}',flush=True)
 if worker:return
 allstatus=[json.loads(p.read_text(encoding='utf-8')) for p in statusdir.glob('*.json')]
 csv([{k:v for k,v in s.items() if k not in ['details','traceback']} for s in allstatus],'results/metrics/model_fit_status.csv')
 frames=[pd.read_parquet(p) for p in sorted(chunks.glob('*.parquet'))]
 if frames:pd.concat(frames,ignore_index=True).to_parquet(ROOT/'results/predictions/out_of_fold_predictions.parquet',index=False)
 print(f'Complete: {len(allstatus)} fit records; {len(frames)} actual prediction chunks',flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--strategies',nargs='*');ap.add_argument('--worker',action='store_true');ap.add_argument('--collect-only',action='store_true');args=ap.parse_args();train(args.strategies,args.worker,args.collect_only)
