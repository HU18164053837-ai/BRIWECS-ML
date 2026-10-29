import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('MKL_NUM_THREADS','2')
from common import *
import argparse,time,warnings,traceback,joblib
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from catboost import CatBoostRegressor
from threadpoolctl import threadpool_limits
def fit_predict(model_name,context,train,test,features,artifact):
 config=cfg();ct=config['contexts'][context];cats=config['categorical'];nums=config['F1_numeric']+(config['weather'] if ct['weather'] else [])+(['disease_severity'] if ct['disease'] else [])
 x=features.loc[train.sample_id,cats+nums].copy();z=features.loc[test.sample_id,cats+nums].copy();y=train.observed_yield.to_numpy()
 if model_name=='Ridge':
  prep=ColumnTransformer([('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=True),cats),('num',Pipeline([('impute',SimpleImputer(strategy='median',keep_empty_features=True)),('scale',StandardScaler())]),nums)],sparse_threshold=1.0)
  model=Pipeline([('preprocess',prep),('model',Ridge(**config['models']['Ridge']))]);model.fit(x,y);pred=model.predict(z)
  if artifact:joblib.dump(model,artifact.with_suffix('.joblib'),compress=3)
  details=dict(numeric_medians=prep.named_transformers_['num'].named_steps['impute'].statistics_.tolist(),feature_columns=cats+nums)
 else:
  imp=SimpleImputer(strategy='median',keep_empty_features=True);x[nums]=imp.fit_transform(x[nums]);z[nums]=imp.transform(z[nums]);model=CatBoostRegressor(**config['models']['CatBoost']);model.fit(x,y,cat_features=cats);pred=model.predict(z)
  if artifact:model.save_model(str(artifact.with_suffix('.cbm')));joblib.dump(imp,artifact.with_suffix('.imputer.joblib'))
  details=dict(numeric_medians=imp.statistics_.tolist(),feature_columns=cats+nums,iterations=model.tree_count_)
 return pred,details
def train(worker=False,collect_only=False):
 config=cfg();config_hash=sha(ROOT/'config/model_config.json');folds=pd.read_csv(ROOT/'splits/disease_outer_fold_manifest.csv');assign=pd.read_csv(ROOT/'splits/disease_test_sample_assignments.csv')
 chunks=ROOT/'results/predictions/chunks';chunks.mkdir(exist_ok=True);states=ROOT/'logs/fit_status';states.mkdir(exist_ok=True)
 with threadpool_limits(limits=2,user_api='blas'):
  for disease,scope in scopes():
   c=load_cohort(disease,scope);name='stripe_rust' if disease=='Stripe_rust' else 'septoria';features=pd.read_parquet(ROOT/f'data/{name}_features.parquet').set_index('sample_id');fr=folds[folds.disease.eq(disease)&folds.scope.eq(scope)]
   for row in fr.itertuples():
    if collect_only:continue
    sub=assign[assign.disease.eq(disease)&assign.scope.eq(scope)&assign.fold_id.eq(row.fold_id)];te=c[c.sample_id.isin(sub.sample_id)].copy();tr=c[~c.sample_id.isin(sub.sample_id)].copy()
    assert ih(tr.sample_id)==row.train_id_sha256 and ih(te.sample_id)==row.test_id_sha256
    for model in ['Ridge','CatBoost']:
     for context in ['M0','M1','M2','M3']:
      job=f'{disease}__{scope}__{row.fold_id}__{model}__{context}';sp=states/(job+'.json');pp=chunks/(job+'.parquet')
      if sp.exists():
       old=json.loads(sp.read_text());assert old['config_sha256']==config_hash and old['train_id_sha256']==row.train_id_sha256
       if old['status'] in ['FAILED','UNSCORABLE'] or (old['status'].startswith('SUCCESS') and pp.exists()):continue
      lock=states/(job+'.lock')
      try:
       with lock.open('x') as f:f.write(str(os.getpid()))
      except FileExistsError:continue
      status=dict(job=job,disease=disease,scope=scope,validation_strategy=row.validation_strategy,fold_id=row.fold_id,model_name=model,context=context,train_n=len(tr),test_n=len(te),train_id_sha256=row.train_id_sha256,test_id_sha256=row.test_id_sha256,config_sha256=config_hash,time_utc=now());t=time.monotonic()
      print('START '+job,flush=True)
      try:
       if row.status=='UNSCORABLE':status.update(status='UNSCORABLE',reason=row.reason)
       else:
        artifact=ROOT/'models'/job if row.validation_strategy=='LOSYO' else None
        with warnings.catch_warnings(record=True) as caught:
         warnings.simplefilter('always');pred,details=fit_predict(model,context,tr,te,features,artifact)
        assert len(pred)==len(te) and np.isfinite(pred).all()
        out=te[['sample_id','observed_yield','disease_severity','site_year','location','year','management','cultivar_id','inoculation_status','other_pathogen_artificial_regime']].copy();out['predicted_yield']=pred;out['disease']=disease;out['scope']=scope;out['validation_strategy']=row.validation_strategy;out['fold_id']=row.fold_id;out['model_name']=model;out['context']=context
        tmp=pp.with_suffix(f'.{os.getpid()}.tmp');out.to_parquet(tmp,index=False);tmp.replace(pp)
        warning=' | '.join(dict.fromkeys(str(w.message) for w in caught));status.update(status='SUCCESS_WITH_WARNINGS' if warning else 'SUCCESS',warnings=warning,details=details)
      except Exception as exc:
       status.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc());(ROOT/'logs'/f'{job}.error.txt').write_text(status['traceback'],encoding='utf-8')
      status['seconds']=time.monotonic()-t;tmp=sp.with_suffix(f'.{os.getpid()}.tmp');tmp.write_text(json.dumps(status,ensure_ascii=False,indent=2,default=str),encoding='utf-8');tmp.replace(sp);lock.unlink();print(f'END {job} {status["status"]} {status["seconds"]:.2f}s',flush=True)
 if worker:return
 states_all=[json.loads(f.read_text()) for f in states.glob('*.json')];csv([{k:v for k,v in s.items() if k not in ['details','traceback']} for s in states_all],'results/metrics/model_fit_status.csv')
 frames=[pd.read_parquet(f) for f in sorted(chunks.glob('*.parquet'))];pd.concat(frames,ignore_index=True).to_parquet(ROOT/'results/predictions/disease_pilot_oof_predictions.parquet',index=False)
 print(f'{len(states_all)} task states; {len(frames)} real test prediction chunks',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--worker',action='store_true');p.add_argument('--collect-only',action='store_true');args=p.parse_args();train(args.worker,args.collect_only)
