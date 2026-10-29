import os
os.environ.setdefault('MPLBACKEND','Agg')
from common import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size':9,'savefig.dpi':160})

def metrics(y,p):
 y=np.asarray(y);p=np.asarray(p);err=p-y;tss=np.square(y-y.mean()).sum()
 r2=1-np.square(err).sum()/tss if len(y)>1 and tss>np.finfo(float).eps*max(1,np.square(y).sum()) else np.nan
 return dict(n=len(y),MAE=float(np.abs(err).mean()),RMSE=float(np.sqrt(np.square(err).mean())),R2=float(r2),MBE=float(err.mean()),observed_sd=float(np.std(y)),R2_status='DEFINED' if np.isfinite(r2) else 'UNDEFINED_LOW_VARIANCE')
def bootstrap_mean(x,rep=2000):
 x=np.asarray(x,float);x=x[np.isfinite(x)]
 if not len(x):return (np.nan,np.nan)
 rng=np.random.default_rng(42);b=x[rng.integers(0,len(x),(rep,len(x)))].mean(axis=1)
 return np.quantile(b,[.025,.975])
def evaluate():
 d=pd.read_parquet(ROOT/'results/predictions/out_of_fold_predictions.parquet')
 assert not d.duplicated(['validation_strategy','feature_set','model_name','sample_id']).any()
 key=['validation_strategy','feature_set','model_name'];cohort=pd.read_parquet(ROOT/'data/cohort/phaseI_yield_cohort.parquet');expected=cohort.groupby('site_year').size()
 env=[];summary=[];folds=[]
 for vals,g in d.groupby(key,sort=True):
  row=dict(zip(key,vals));er=[]
  for sy,z in g.groupby('site_year'):
   stat=metrics(z.observed_yield,z.predicted_yield);x=dict(**row,site_year=sy,location=z.location.iloc[0],year=int(z.year.iloc[0]),test_fold_count=z.fold_id.nunique(),expected_n=int(expected[sy]),**stat);er.append(x);env.append(x)
  ed=pd.DataFrame(er);pooled=metrics(g.observed_yield,g.predicted_yield);lo,hi=bootstrap_mean(ed.MAE)
  summary.append(dict(**row,aggregation='POOLED',environments=len(ed),expected_total=len(cohort) if vals[0]!='LOSYO_NO_QLB' else int(cohort.location.ne('QLB').sum()),**pooled,MAE_env_bootstrap_lo=np.nan,MAE_env_bootstrap_hi=np.nan))
  summary.append(dict(**row,aggregation='ENVIRONMENT_MACRO',n=len(g),environments=len(ed),expected_total=len(cohort) if vals[0]!='LOSYO_NO_QLB' else int(cohort.location.ne('QLB').sum()),MAE=ed.MAE.mean(),RMSE=ed.RMSE.mean(),R2=ed.R2.mean(),MBE=ed.MBE.mean(),observed_sd=np.nan,R2_status='MEAN_OF_DEFINED_ENVIRONMENT_R2',MAE_env_bootstrap_lo=lo,MAE_env_bootstrap_hi=hi))
  fm=[]
  for fid,z in g.groupby('fold_id'):
   q=dict(**row,fold_id=fid,**metrics(z.observed_yield,z.predicted_yield));folds.append(q);fm.append(q)
  fm=pd.DataFrame(fm)
  summary.append(dict(**row,aggregation='TEST_GROUP_MACRO',n=len(g),environments=len(fm),expected_total=len(cohort) if vals[0]!='LOSYO_NO_QLB' else int(cohort.location.ne('QLB').sum()),MAE=fm.MAE.mean(),RMSE=fm.RMSE.mean(),R2=fm.R2.mean(),MBE=fm.MBE.mean(),observed_sd=np.nan,R2_status='MEAN_OF_DEFINED_TEST_GROUP_R2',MAE_env_bootstrap_lo=np.nan,MAE_env_bootstrap_hi=np.nan))
 env=csv(env,'results/metrics/environment_level_metrics.csv');summary=csv(summary,'results/metrics/model_performance_summary.csv');csv(folds,'results/metrics/fold_level_metrics.csv')
 # Matched comparison; negative delta = A improves upon B.
 comparison=[];differences=[]
 def compare(strategy,ma,fa,mb,fb,label):
  aa=d[d.validation_strategy.eq(strategy)&d.model_name.eq(ma)&d.feature_set.eq(fa)];bb=d[d.validation_strategy.eq(strategy)&d.model_name.eq(mb)&d.feature_set.eq(fb)]
  joined=aa[['sample_id','site_year','fold_id','observed_yield','predicted_yield']].merge(bb[['sample_id','fold_id','observed_yield','predicted_yield']],on=['sample_id','fold_id'],suffixes=('_a','_b'),validate='one_to_one')
  if not len(joined):return
  assert np.allclose(joined.observed_yield_a,joined.observed_yield_b)
  delta=[]
  for sy,g in joined.groupby('site_year'):
   am=metrics(g.observed_yield_a,g.predicted_yield_a);bm=metrics(g.observed_yield_b,g.predicted_yield_b)
   val=dict(comparison=label,validation_strategy=strategy,model_a=ma,feature_a=fa,model_b=mb,feature_b=fb,site_year=sy,n=len(g),delta_MAE=am['MAE']-bm['MAE'],delta_RMSE=am['RMSE']-bm['RMSE'],delta_R2=am['R2']-bm['R2'],MAE_a=am['MAE'],MAE_b=bm['MAE'],relative_MAE_improvement=(bm['MAE']-am['MAE'])/bm['MAE'])
   delta.append(val);differences.append(val)
  de=pd.DataFrame(delta);lo,hi=bootstrap_mean(de.delta_MAE)
  comparison.append(dict(comparison=label,validation_strategy=strategy,model_a=ma,feature_a=fa,model_b=mb,feature_b=fb,matched_samples=len(joined),samples_a=len(aa),samples_b=len(bb),exact_same_test_samples=len(joined)==len(aa)==len(bb),environments=len(de),macro_delta_MAE=de.delta_MAE.mean(),delta_MAE_CI95_lo=lo,delta_MAE_CI95_hi=hi,macro_delta_RMSE=de.delta_RMSE.mean(),macro_delta_R2=de.delta_R2.mean(),environments_improved=int(de.delta_MAE.lt(0).sum()),direction_consistency=de.delta_MAE.lt(0).mean(),uncertainty='Paired site-year bootstrap; shared years/sites limit independence; descriptive not confirmatory'))
 for strategy in sorted(d.validation_strategy.unique()):
  for fs in ['F1','F2']:
   for ma in ['CatBoost','XGBoost']:
    compare(strategy,ma,fs,'TrainingMean','F0',f'{ma}_{fs}_vs_TrainingMean')
    compare(strategy,ma,fs,'Ridge',fs,f'{ma}_{fs}_vs_Ridge_{fs}')
  for ma in ['Ridge','MixedLM','CatBoost','XGBoost']:compare(strategy,ma,'F2',ma,'F1',f'{ma}_F2_vs_F1')
 csv(comparison,'results/metrics/model_comparison_paired.csv');diff=csv(differences,'results/metrics/environment_paired_differences.csv')
 # Prespecified descriptive diagnostics, never used to refit.
 diag=[];primary=d[d.validation_strategy.eq('LOSYO')].copy();primary['yield_range']=pd.cut(primary.observed_yield,[-np.inf,cohort.observed_yield.quantile(.25),cohort.observed_yield.quantile(.75),np.inf],labels=['LOW_Q1','MIDDLE','HIGH_Q4'])
 for col in ['location','year','management','yield_range']:
  for vals,g in primary.groupby(['feature_set','model_name',col],observed=True):
   fs,model,group=vals;diag.append(dict(feature_set=fs,model_name=model,stratum_type=col,stratum=str(group),**metrics(g.observed_yield,g.predicted_yield)))
 csv(diag,'results/metrics/error_stratification.csv')
 weather=pd.read_csv(ROOT/'data/features/site_year_weather_features.csv');weatherflags=[]
 for row in weather.itertuples():
  other=weather[weather.site_year.ne(row.site_year)];outside=[]
  for variable in CFG['F2']['numeric_additions']:
   if not other[variable].min()<=getattr(row,variable)<=other[variable].max():outside.append(variable)
  weatherflags.append(dict(site_year=row.site_year,weather_outside_training_range=len(outside)>0,outside_variables='|'.join(outside),outside_variable_count=len(outside)))
 flags=csv(weatherflags,'results/metrics/weather_extrapolation_flags.csv')
 wd=env[env.validation_strategy.eq('LOSYO')].merge(weather,on='site_year',validate='many_to_one').merge(flags,on='site_year',validate='many_to_one');csv(wd,'results/metrics/environment_weather_diagnostics.csv')
 outliers=primary.merge(cohort[['sample_id','descriptive_extreme_IQR']],on='sample_id',validate='many_to_one');csv([dict(feature_set=fs,model_name=m,extreme_yield=bool(flag),**metrics(g.observed_yield,g.predicted_yield)) for (fs,m,flag),g in outliers.groupby(['feature_set','model_name','descriptive_extreme_IQR'])],'results/metrics/extreme_yield_diagnostics.csv')
 make_figures(primary,env,diff,cohort)
 print(summary[summary.aggregation.eq('POOLED')][key+['n','MAE','RMSE','R2','MBE']].to_string(index=False),flush=True)

def savefig(fig,name):
 fig.savefig(ROOT/'results/figures'/f'{name}.png',bbox_inches='tight');fig.savefig(ROOT/'results/figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def make_figures(primary,env,diff,c):
 coverage=c.groupby(['site_year','management']).size().unstack().sort_index();fig,ax=plt.subplots(figsize=(10,7));im=ax.imshow(coverage.fillna(0),aspect='auto',cmap='Blues');ax.set_xticks(range(len(coverage.columns)),coverage.columns,rotation=55,ha='right');ax.set_yticks(range(len(coverage)),coverage.index)
 for i in range(len(coverage)):
  for j in range(len(coverage.columns)):
   val=coverage.iloc[i,j];ax.text(j,i,'—' if pd.isna(val) else str(int(val)),ha='center',va='center',fontsize=8,color='black')
 fig.colorbar(im,ax=ax,label='Actual independent record units');ax.set_title('Figure 2A-1: frozen Phase I coverage (absent combinations not imputed)');savefig(fig,'Figure_2A_1_Cohort_Coverage')
 e=env[env.validation_strategy.eq('LOSYO')].copy();e['label']=e.model_name+' '+e.feature_set;labels=sorted(e.label.unique(),key=lambda s:('TrainingMean' not in s,s));fig,ax=plt.subplots(figsize=(11,5));values=[e.loc[e.label.eq(l),'MAE'].to_numpy() for l in labels];ax.boxplot(values,tick_labels=labels,showfliers=False)
 for i,val in enumerate(values):ax.scatter(np.full(len(val),i+1)+np.linspace(-.12,.12,len(val)),val,s=16,alpha=.65)
 ax.set_ylabel('Site-year MAE (dt/ha)');ax.tick_params(axis='x',rotation=50);ax.set_title('Figure 2A-2: LOSYO errors, all held-out environments');savefig(fig,'Figure_2A_2_Model_Errors')
 fig,axes=plt.subplots(3,3,figsize=(11,10));miny=min(primary.observed_yield.min(),primary.predicted_yield.min());maxy=max(primary.observed_yield.max(),primary.predicted_yield.max())
 for ax,label in zip(axes.flat,labels):
  model,fs=label.split(' ');g=primary[primary.model_name.eq(model)&primary.feature_set.eq(fs)];ax.hexbin(g.observed_yield,g.predicted_yield,gridsize=45,mincnt=1,bins='log',cmap='viridis');ax.plot([miny,maxy],[miny,maxy],'r--',lw=.8);met=metrics(g.observed_yield,g.predicted_yield);ax.set_title(f'{label}: R²={met["R2"]:.3f}, n={len(g)}');ax.set_xlabel('Observed (dt/ha)');ax.set_ylabel('OOF predicted (dt/ha)')
 fig.suptitle('Figure 2A-3: actual versus genuine LOSYO test predictions');fig.tight_layout();savefig(fig,'Figure_2A_3_Observed_Predicted')
 b=e.pivot(index='site_year',columns='label',values='MBE').reindex(columns=labels);lim=max(abs(b.min().min()),abs(b.max().max()));fig,ax=plt.subplots(figsize=(11,7));im=ax.imshow(b,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim);ax.set_xticks(range(len(b.columns)),b.columns,rotation=55,ha='right');ax.set_yticks(range(len(b)),b.index)
 for i in range(len(b)):
  for j in range(len(b.columns)):ax.text(j,i,f'{b.iloc[i,j]:.1f}',ha='center',va='center',fontsize=7)
 fig.colorbar(im,ax=ax,label='Predicted − observed (dt/ha)');ax.set_title('Figure 2A-4: all site-year systematic bias');savefig(fig,'Figure_2A_4_Bias_Heatmap')
 delta=diff[diff.validation_strategy.eq('LOSYO')&diff.feature_a.eq('F2')&diff.feature_b.eq('F1')];fig,ax=plt.subplots(figsize=(10,5))
 for i,(model,g) in enumerate(delta.groupby('model_a')):ax.scatter(np.full(len(g),i)+np.linspace(-.10,.10,len(g)),g.delta_MAE,s=24,alpha=.7);ax.scatter([i],[g.delta_MAE.mean()],s=80,marker='D',c='black')
 ax.set_xticks(range(delta.model_a.nunique()),sorted(delta.model_a.unique()));ax.axhline(0,color='red',ls='--');ax.set_ylabel('Paired site-year MAE(F2) − MAE(F1), dt/ha');ax.set_title('Figure 2A-5: fixed winter weather increment, same test samples');savefig(fig,'Figure_2A_5_F2_vs_F1')
 fig,ax=plt.subplots(figsize=(7,4));ax.hist(c.observed_yield,bins=55,color='#487c9b');ax.set_xlabel('Whole-plot grain yield (dt/ha)');ax.set_ylabel('Frozen records');ax.set_title('Outcome distribution; extreme records retained');savefig(fig,'Diagnostic_Yield_Distribution')
if __name__=='__main__':evaluate()
