from common import *
from scipy.stats import spearmanr
def rho(x,y):return float(spearmanr(x,y).statistic) if pd.Series(x).nunique()>1 and pd.Series(y).nunique()>1 else np.nan
def run():
 checks=[]
 for base,manifest in [(P,P/'reports/SHA256_delivery_manifest.csv'),(A,A/'reports/output_manifest.csv')]:
  for r in pd.read_csv(manifest).itertuples():
   ok=sha(base/r.relative_path)==r.sha256;checks.append(dict(path=str((base/r.relative_path).relative_to(P1)),sha256=r.sha256,match=ok))
 csv(checks,'input_integrity.csv');assert all(r['match'] for r in checks)
 cfg=dict(seed=42,association_cluster_bootstrap=499,bootstrap_unit='site-year pairs, re-absorb fixed effects in every sample',new_ML_fits=0,transfer='Reuse Pilot FULL/NATURAL LOSYO M0/M1, same natural-domain test IDs and exact inherited training sets',error_variance_ddof=0,centering='Diagnostic only; never calibration',disease_scale='pooled: 10 nominal reported score units, not unconditionally harmonized percentage points',frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 if not (ROOT/'config/analysis.json').exists():js(cfg,'config/analysis.json')
 df=pd.read_parquet(P/'results/predictions/disease_pilot_oof_predictions.parquet');assert not df.duplicated(['disease','scope','validation_strategy','model_name','context','sample_id']).any()
 df=df[df.disease.eq('Stripe_rust')&df.validation_strategy.eq('LOSYO')].copy();co=pd.read_parquet(P/'data/stripe_rust_paired_cohort.parquet');assert len(co)==8539
 # Preserve source values and receipt prior reports/configuration.
 co.to_parquet(ROOT/'data/stripe_rust_inherited.parquet',index=False)
 receipt=[]
 for base in [A,P]:
  for f in list((base/'reports').glob('BRIWECS*.md'))+list((base/'config').glob('*.json'))+list((base/'scripts').glob('*.py')):receipt.append(dict(path=str(f.relative_to(P1)),sha256=sha(f),bytes=f.stat().st_size))
 csv(receipt,'inheritance_receipt.csv')
 # Decomposition is exact using population variance, not an MAE decomposition.
 dec=[];ranking=[];delta=[]
 for (scope,model,ctx,sy),g in df.groupby(['scope','model_name','context','site_year']):
  y=g.observed_yield.to_numpy();pred=g.predicted_yield.to_numpy();e=pred-y;m=met(y,pred);v=np.var(e,ddof=0);b=e.mean()**2
  assert np.isclose(m['MSE'],b+v,atol=1e-10)
  dec.append(dict(scope=scope,model_name=model,context=ctx,site_year=sy,observed_mean=y.mean(),predicted_mean=pred.mean(),yield_sd=np.std(y,ddof=0),bias_squared=b,error_variance=v,centered_MAE=np.abs(e-e.mean()).mean(),centered_RMSE=np.sqrt(v),**m,identity_error=m['MSE']-b-v,centering_use='POSTHOC_DIAGNOSTIC_NOT_DEPLOYABLE_CALIBRATION'))
  ranking.append(dict(scope=scope,model_name=model,context=ctx,site_year=sy,level='PLOT',n=len(g),cultivars=g.cultivar_id.nunique(),spearman=rho(y,pred),replicate_support='actual sample records; cultivar repeats correlated'))
  z=g.groupby('cultivar_id').agg(y=('observed_yield','mean'),pred=('predicted_yield','mean'),plots=('sample_id','size'))
  ranking.append(dict(scope=scope,model_name=model,context=ctx,site_year=sy,level='CULTIVAR_MEAN',n=len(z),cultivars=len(z),spearman=rho(z.y,z.pred),replicated_cultivars=int(z.plots.ge(2).sum()),low_quartile_MAE=abs(z.loc[z.y.le(z.y.quantile(.25)),'pred']-z.loc[z.y.le(z.y.quantile(.25)),'y']).mean(),high_quartile_MAE=abs(z.loc[z.y.ge(z.y.quantile(.75)),'pred']-z.loc[z.y.ge(z.y.quantile(.75)),'y']).mean(),replicate_support='Descriptive cultivar means; unequal management composition, not BLUEs'))
 dec=csv(dec,'01_environment_error_decomposition.csv');rank=csv(ranking,'02_within_environment_ranking.csv')
 for (scope,model),g in dec.groupby(['scope','model_name']):
  for b,v in [('M0','M1'),('M2','M3')]:
   j=g[g.context.eq(b)].merge(g[g.context.eq(v)],on='site_year',suffixes=('_b','_v'),validate='one_to_one')
   for r in j.itertuples():
    db=r.bias_squared_b-r.bias_squared_v;dv=r.error_variance_b-r.error_variance_v;dm=r.MSE_b-r.MSE_v;assert np.isclose(dm,db+dv)
    delta.append(dict(scope=scope,model_name=model,comparison=v+'_vs_'+b,site_year=r.site_year,n=r.n_b,delta_MSE=dm,delta_bias_squared=db,delta_centered_variance=dv,delta_MAE=r.MAE_b-r.MAE_v,delta_centered_MAE=r.centered_MAE_b-r.centered_MAE_v,delta_centered_RMSE=r.centered_RMSE_b-r.centered_RMSE_v))
 delta=csv(delta,'03_disease_delta_error_components.csv')
 # Transfer: full-trained and natural-trained predictions tested on identical natural records.
 states=pd.read_csv(P/'results/metrics/model_fit_status.csv');natids=set(co.loc[co.inoculation_status.eq('NATURAL'),'sample_id']);man=[];rows=[];summary=[];drop=[]
 for model in ['CatBoost','Ridge']:
  for sy in sorted(co.loc[co.inoculation_status.eq('NATURAL'),'site_year'].unique()):
   expected=set(co.loc[co.site_year.eq(sy)&co.inoculation_status.eq('NATURAL'),'sample_id']);predsets={}
   for strategy,scope in [('A_NATURAL_ONLY','NATURAL'),('B_FULL_ELIGIBLE','FULL')]:
    tr=co[co.site_year.ne(sy)];tr=tr[tr.inoculation_status.eq('NATURAL')] if scope=='NATURAL' else tr
    for context in ['M0','M1']:
     g=df[df.scope.eq(scope)&df.model_name.eq(model)&df.site_year.eq(sy)&df.context.eq(context)&df.sample_id.isin(natids)].copy();assert set(g.sample_id)==expected
     st=states[states.disease.eq('Stripe_rust')&states.scope.eq(scope)&states.validation_strategy.eq('LOSYO')&states.model_name.eq(model)&states.context.eq(context)&states.fold_id.eq('LOSYO__'+sy)].iloc[0];assert st.train_id_sha256==ih(tr.sample_id)
     man.append(dict(model_name=model,site_year=sy,strategy=strategy,context=context,test_n=len(g),test_id_sha256=ih(g.sample_id),train_n=len(tr),train_environments=tr.site_year.nunique(),artificial_training_fraction=tr.inoculation_status.eq('ARTIFICIAL').mean(),train_id_sha256=ih(tr.sample_id),reuse_status='EXACT_PILOT_MODEL_AND_OOF_REUSE',test_status='NATURAL_ONLY'))
     predsets[(strategy,context)]=g
    base=predsets[(strategy,'M0')];aug=predsets[(strategy,'M1')];j=base.merge(aug[['sample_id','predicted_yield']],on='sample_id',suffixes=('_b','_v'),validate='one_to_one');bm=met(j.observed_yield,j.predicted_yield_b);vm=met(j.observed_yield,j.predicted_yield_v)
    rows.append(dict(model_name=model,site_year=sy,strategy=strategy,n=len(j),MAE_without=bm['MAE'],MAE_with=vm['MAE'],delta_MAE=bm['MAE']-vm['MAE'],RMSE_without=bm['RMSE'],RMSE_with=vm['RMSE'],R2_without=bm['R2'],R2_with=vm['R2']))
 rows=pd.DataFrame(rows)
 comparisons=[]
 for model,g in rows.groupby('model_name'):
  for strategy,h in g.groupby('strategy'):
   lo,hi=ci(h.delta_MAE);pl,ph=ci(h.delta_MAE,h.n)
   summary.append(dict(model_name=model,strategy=strategy,n=int(h.n.sum()),environments=len(h),pooled_delta_MAE=np.average(h.delta_MAE,weights=h.n),macro_delta_MAE=h.delta_MAE.mean(),macro_CI95_lo=lo,macro_CI95_hi=hi,pooled_CI95_lo=pl,pooled_CI95_hi=ph,improved_environments=int(h.delta_MAE.gt(0).sum())))
  j=g[g.strategy.eq('A_NATURAL_ONLY')].merge(g[g.strategy.eq('B_FULL_ELIGIBLE')],on=['site_year','model_name'],suffixes=('_A','_B'));j['strategy_delta_difference_B_minus_A']=j.delta_MAE_B-j.delta_MAE_A
  comparisons.extend(j.to_dict('records'));lo,hi=ci(j.strategy_delta_difference_B_minus_A);pl,ph=ci(j.strategy_delta_difference_B_minus_A,j.n_A)
  summary.append(dict(model_name=model,strategy='B_MINUS_A_INCREMENT',n=int(j.n_A.sum()),environments=len(j),pooled_delta_MAE=np.average(j.strategy_delta_difference_B_minus_A,weights=j.n_A),macro_delta_MAE=j.strategy_delta_difference_B_minus_A.mean(),macro_CI95_lo=lo,macro_CI95_hi=hi,pooled_CI95_lo=pl,pooled_CI95_hi=ph,improved_environments=int(j.strategy_delta_difference_B_minus_A.gt(0).sum())))
  for removed in j.site_year:
   h=j[j.site_year.ne(removed)];drop.append(dict(model_name=model,removed_environment=removed,macro_A=h.delta_MAE_A.mean(),macro_B=h.delta_MAE_B.mean(),pooled_A=np.average(h.delta_MAE_A,weights=h.n_A),pooled_B=np.average(h.delta_MAE_B,weights=h.n_A),macro_B_minus_A=h.strategy_delta_difference_B_minus_A.mean(),pooled_B_minus_A=np.average(h.strategy_delta_difference_B_minus_A,weights=h.n_A)))
 csv(man,'07_matched_test_domain_manifest.csv');csv(summary,'08_matched_training_strategy_metrics.csv');csv(comparisons,'09_environment_level_strategy_comparison.csv');csv(rows,'matched_strategy_environment_metrics.csv');csv(drop,'matched_domain_delete_environment.csv')
 js(dict(actual_new_ML_fits=0,reused_primary_models=80,reason='Both training strategies and identical natural test samples exactly matched Pilot fingerprints'),'logs/retraining_status.json')
 js(dict(mse_identities_max_error=float(dec.identity_error.abs().max()),transfer_identical_test_IDs=True,new_ML_fits=0),'reports/prediction_verification.json')
 print(pd.DataFrame(summary).to_string(index=False),flush=True)
if __name__=='__main__':run()
