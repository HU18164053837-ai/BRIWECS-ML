from common import *
def met(y,p):
 y=np.asarray(y,float);p=np.asarray(p,float);e=p-y;t=np.square(y-y.mean()).sum();r2=1-np.square(e).sum()/t if len(y)>1 and t>np.finfo(float).eps*max(1,np.square(y).sum()) else np.nan
 return dict(n=len(y),MAE=np.abs(e).mean(),RMSE=np.sqrt(np.square(e).mean()),R2=r2,MBE=e.mean(),R2_status='DEFINED' if np.isfinite(r2) else 'UNDEFINED')
def interval(values,weights=None,B=2000):
 v=np.asarray(values,float);ok=np.isfinite(v);v=v[ok]
 if len(v)==0:return np.nan,np.nan
 rng=np.random.default_rng(42);ii=rng.integers(0,len(v),(B,len(v)))
 if weights is None:s=v[ii].mean(axis=1)
 else:
  w=np.asarray(weights,float)[ok];s=(v[ii]*w[ii]).sum(axis=1)/w[ii].sum(axis=1)
 return np.quantile(s,[.025,.975])
def evaluate():
 d=pd.read_parquet(ROOT/'results/predictions/disease_pilot_oof_predictions.parquet');assert not d.duplicated(['disease','scope','validation_strategy','model_name','context','sample_id']).any()
 statuses=pd.read_csv(ROOT/'results/metrics/model_fit_status.csv');absolute=[];foldmetrics=[];environment=[];performance=[];uncertainty=[];pairchecks=[];drop=[]
 group=['disease','scope','validation_strategy','model_name']
 for vals,g in d.groupby(group):
  keys=dict(zip(group,vals));absenv={}
  for context,x in g.groupby('context'):
   em=[]
   for sy,h in x.groupby('site_year'):
    m=met(h.observed_yield,h.predicted_yield);em.append(dict(site_year=sy,**m));absolute.append(dict(**keys,context=context,aggregation='SITE_YEAR',site_year=sy,**m))
   aa=pd.DataFrame(em);absolute.append(dict(**keys,context=context,aggregation='POOLED',site_year='ALL',**met(x.observed_yield,x.predicted_yield)));absolute.append(dict(**keys,context=context,aggregation='ENVIRONMENT_MACRO',site_year='ALL',n=len(x),MAE=aa.MAE.mean(),RMSE=aa.RMSE.mean(),R2=aa.R2.mean(),MBE=aa.MBE.mean(),R2_status='MEAN_OF_DEFINED_ENV_R2'))
   for fid,h in x.groupby('fold_id'):foldmetrics.append(dict(**keys,context=context,fold_id=fid,**met(h.observed_yield,h.predicted_yield)))
  for without,withd in [('M0','M1'),('M2','M3')]:
   ba=g[g.context.eq(without)];au=g[g.context.eq(withd)];comp=f'{withd}_vs_{without}'
   paired=ba.merge(au[['sample_id','fold_id','predicted_yield']],on=['sample_id','fold_id'],suffixes=('_without','_with'),validate='one_to_one')
   for fid in sorted(set(ba.fold_id)|set(au.fold_id)):
    u=ba[ba.fold_id.eq(fid)];v=au[au.fold_id.eq(fid)];st=statuses[(statuses.disease.eq(keys['disease']))&(statuses.scope.eq(keys['scope']))&(statuses.validation_strategy.eq(keys['validation_strategy']))&(statuses.model_name.eq(keys['model_name']))&statuses.fold_id.eq(fid)]
    st=st.set_index('context');same_train=without in st.index and withd in st.index and st.loc[without,'train_id_sha256']==st.loc[withd,'train_id_sha256']
    same_test=set(u.sample_id)==set(v.sample_id)
    pairchecks.append(dict(**keys,comparison=comp,fold_id=fid,train_ids_identical=same_train,test_ids_identical=same_test,n_without=len(u),n_with=len(v),status='PASS' if same_train and same_test and len(u) else 'PARTIAL_OR_UNSCORABLE'))
   if len(paired)==0:
    performance.append(dict(**keys,comparison=comp,aggregation='POOLED',status='NOT ESTIMABLE',n=0));continue
   assert np.isfinite(paired.observed_yield).all()
   em=[]
   for sy,h in paired.groupby('site_year'):
    base=met(h.observed_yield,h.predicted_yield_without);aug=met(h.observed_yield,h.predicted_yield_with)
    rr=dict(**keys,comparison=comp,site_year=sy,location=h.location.iloc[0],year=int(h.year.iloc[0]),n=len(h),MAE_without=base['MAE'],MAE_with=aug['MAE'],RMSE_without=base['RMSE'],RMSE_with=aug['RMSE'],R2_without=base['R2'],R2_with=aug['R2'],MBE_without=base['MBE'],MBE_with=aug['MBE'],delta_MAE=base['MAE']-aug['MAE'],delta_RMSE=base['RMSE']-aug['RMSE'],delta_R2=aug['R2']-base['R2'],relative_MAE_gain_percent=100*(base['MAE']-aug['MAE'])/base['MAE'],positive_severity_proportion=h.disease_severity.gt(0).mean(),inoculation_status='|'.join(sorted(h.inoculation_status.unique())),status='SCORABLE')
    environment.append(rr);em.append(rr)
   em=pd.DataFrame(em);base=met(paired.observed_yield,paired.predicted_yield_without);aug=met(paired.observed_yield,paired.predicted_yield_with)
   shared=dict(**keys,comparison=comp,n=len(paired),environments=len(em),outer_test_groups=paired.fold_id.nunique(),exact_same_prediction_coverage=len(paired)==len(ba)==len(au),improved_environments=int(em.delta_MAE.gt(0).sum()),worsened_environments=int(em.delta_MAE.lt(0).sum()),unchanged_environments=int(em.delta_MAE.eq(0).sum()),status='COMPLETE' if len(paired)==len(ba)==len(au) else 'PARTIAL',task='RETROSPECTIVE ASSOCIATION-ORIENTED ESTIMATION')
   performance.append(dict(**shared,aggregation='POOLED',MAE_without=base['MAE'],MAE_with=aug['MAE'],RMSE_without=base['RMSE'],RMSE_with=aug['RMSE'],R2_without=base['R2'],R2_with=aug['R2'],MBE_without=base['MBE'],MBE_with=aug['MBE'],delta_MAE=base['MAE']-aug['MAE'],delta_RMSE=base['RMSE']-aug['RMSE'],delta_R2=aug['R2']-base['R2'],relative_MAE_gain_percent=100*(base['MAE']-aug['MAE'])/base['MAE']))
   performance.append(dict(**shared,aggregation='ENVIRONMENT_MACRO',MAE_without=em.MAE_without.mean(),MAE_with=em.MAE_with.mean(),RMSE_without=em.RMSE_without.mean(),RMSE_with=em.RMSE_with.mean(),R2_without=em.R2_without.mean(),R2_with=em.R2_with.mean(),MBE_without=em.MBE_without.mean(),MBE_with=em.MBE_with.mean(),delta_MAE=em.delta_MAE.mean(),delta_RMSE=em.delta_RMSE.mean(),delta_R2=em.delta_R2.mean(),relative_MAE_gain_percent=100*em.delta_MAE.mean()/em.MAE_without.mean()))
   lo,hi=interval(em.delta_MAE);pl,ph=interval(em.delta_MAE,em.n);rl,rh=interval(em.delta_RMSE)
   gains=em.delta_MAE.clip(lower=0);positive_total=gains.sum();top2=gains.nlargest(2).sum()/positive_total if positive_total>0 else np.nan
   uncertainty.append(dict(**shared,macro_delta_MAE=em.delta_MAE.mean(),macro_CI95_lo=lo,macro_CI95_hi=hi,pooled_delta_MAE=base['MAE']-aug['MAE'],pooled_CI95_lo=pl,pooled_CI95_hi=ph,macro_delta_RMSE=em.delta_RMSE.mean(),macro_delta_RMSE_CI95_lo=rl,macro_delta_RMSE_CI95_hi=rh,two_largest_positive_gain_share=top2,uncertainty_unit='paired site-year clusters; 2000 percentile resamples; shared site/year limits independence',evidence_strength='LIMITED: 3 test-year groups' if keys['validation_strategy']=='LOYO' else 'LIMITED: 4-5 test-location groups' if keys['validation_strategy']=='LOLO' else 'LIMITED: 10-13 site-years'))
   if keys['validation_strategy']=='LOSYO':
    for removed in em.site_year:
     rest=em[em.site_year.ne(removed)];macro=rest.delta_MAE.mean();pooled=np.average(rest.delta_MAE,weights=rest.n)
     drop.append(dict(**keys,comparison=comp,removed_environment=removed,remaining_environments=len(rest),macro_delta_MAE=macro,pooled_delta_MAE=pooled,macro_sign_flip=bool(np.sign(macro)!=np.sign(em.delta_MAE.mean())),pooled_sign_flip=bool(np.sign(pooled)!=np.sign(base['MAE']-aug['MAE']))))
 csv(absolute,'results/metrics/absolute_model_performance.csv');csv(foldmetrics,'results/metrics/fold_model_results.csv');perf=csv(performance,'results/metrics/disease_incremental_performance.csv');environment=csv(environment,'results/metrics/environment_incremental_metrics.csv');uncertainty=csv(uncertainty,'results/metrics/paired_uncertainty_analysis.csv');csv(pairchecks,'results/metrics/paired_prediction_integrity.csv');csv(drop,'results/metrics/drop_one_environment_sensitivity.csv')
 # Natural-infection effects always compared with newly trained matched controls.
 ins=perf[perf.validation_strategy.eq('LOSYO')].merge(uncertainty[['disease','scope','validation_strategy','model_name','comparison','macro_CI95_lo','macro_CI95_hi','two_largest_positive_gain_share']],on=['disease','scope','validation_strategy','model_name','comparison'],validate='many_to_one')
 ins['interpretation']='CONFOUNDED: natural subset changes environment composition; no inoculation causal effect'
 ins['training_alias']=''
 alias=ins[ins.disease.eq('Septoria')&ins.scope.eq('FULL')].copy();alias['scope']='NATURAL';alias['training_alias']='FULL: identical samples, exact reused fits, not independent result';ins=pd.concat([ins,alias],ignore_index=True)
 csv(ins,'results/metrics/inoculation_sensitivity.csv')
 # Management outcome and errors remain descriptive, not treatment effects.
 management=[]
 for vals,g in d[d.validation_strategy.eq('LOSYO')].groupby(['disease','scope','model_name','context','management']):management.append(dict(zip(['disease','scope','model_name','context','management'],vals),n=len(g),environments=g.site_year.nunique(),severity_mean=g.disease_severity.mean(),positive_proportion=g.disease_severity.gt(0).mean(),**{k:v for k,v in met(g.observed_yield,g.predicted_yield).items() if k!='n'}))
 csv(management,'results/metrics/management_error_diagnostics.csv')
 print(perf[perf.validation_strategy.eq('LOSYO')][['disease','scope','model_name','comparison','aggregation','n','delta_MAE','delta_RMSE','delta_R2','improved_environments','environments']].to_string(index=False),flush=True)
if __name__=='__main__':evaluate()
