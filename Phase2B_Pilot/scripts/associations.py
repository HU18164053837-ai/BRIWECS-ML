from common import *
from scipy.stats import spearmanr
import warnings
def rho(x,y):
 if pd.Series(x).nunique()<2 or pd.Series(y).nunique()<2:return np.nan
 return float(spearmanr(x,y).statistic)
def cluster_ci(g,B=500):
 r=rho(g.disease_severity,g.observed_yield)
 if not np.isfinite(r):return np.nan,np.nan
 ids=g.cultivar_id.unique();index={s:np.flatnonzero(g.cultivar_id.to_numpy()==s) for s in ids};x=g.disease_severity.to_numpy();y=g.observed_yield.to_numpy();rng=np.random.default_rng(42);res=[]
 with warnings.catch_warnings():
  warnings.simplefilter('ignore')
  for _ in range(B):
   sample=rng.choice(ids,size=len(ids),replace=True);ii=np.concatenate([index[s] for s in sample]);res.append(rho(x[ii],y[ii]))
 res=np.asarray(res);res=res[np.isfinite(res)]
 return np.quantile(res,[.025,.975]) if len(res) else (np.nan,np.nan)
def associations():
 env=[];overall=[];mgmt=[];art=[]
 for disease,scope in scopes():
  c=load_cohort(disease,scope);s=[]
  for sy,g in c.groupby('site_year'):
   r=rho(g.disease_severity,g.observed_yield);lo,hi=cluster_ci(g)
   row=dict(disease=disease,scope=scope,site_year=sy,n=len(g),cultivars=g.cultivar_id.nunique(),managements=g.management.nunique(),severity_mean=g.disease_severity.mean(),severity_sd=g.disease_severity.std(),severity_min=g.disease_severity.min(),severity_max=g.disease_severity.max(),zero_proportion=g.disease_severity.eq(0).mean(),positive_proportion=g.disease_severity.gt(0).mean(),yield_mean=g.observed_yield.mean(),yield_sd=g.observed_yield.std(),spearman=r,CI95_lo=lo,CI95_hi=hi,status='DEFINED' if np.isfinite(r) else 'UNDEFINED_NO_VARIATION',uncertainty='500 cultivar-cluster bootstrap samples within site-year; descriptive',inoculation_status='|'.join(sorted(g.inoculation_status.unique())))
   env.append(row);s.append(row)
  e=pd.DataFrame(s);xm=c.disease_severity-c.groupby('site_year').disease_severity.transform('mean');ym=c.observed_yield-c.groupby('site_year').observed_yield.transform('mean');xmm=c.disease_severity-c.groupby(['site_year','management']).disease_severity.transform('mean');ymm=c.observed_yield-c.groupby(['site_year','management']).observed_yield.transform('mean')
  overall.append(dict(disease=disease,scope=scope,n=len(c),environments=c.site_year.nunique(),pooled_spearman=rho(c.disease_severity,c.observed_yield),between_environment_spearman=rho(e.severity_mean,e.yield_mean),within_environment_centered_spearman=rho(xm,ym),within_environment_management_centered_spearman=rho(xmm,ymm),negative_environment_correlations=int(e.spearman.lt(0).sum()),positive_environment_correlations=int(e.spearman.gt(0).sum()),undefined_correlations=int(e.spearman.isna().sum()),interpretation='Descriptive, no causal or independent-genotype inference; centered true yield never used by prediction models'))
  for (sy,m),g in c.groupby(['site_year','management']):mgmt.append(dict(disease=disease,scope=scope,site_year=sy,management=m,n=len(g),cultivars=g.cultivar_id.nunique(),severity_mean=g.disease_severity.mean(),severity_sd=g.disease_severity.std(),yield_mean=g.observed_yield.mean(),spearman=rho(g.disease_severity,g.observed_yield),status='DEFINED' if np.isfinite(rho(g.disease_severity,g.observed_yield)) else 'UNDEFINED_NO_VARIATION'))
  if disease=='Stripe_rust' and scope=='FULL':
   for sy,g in c[c.inoculation_status.eq('ARTIFICIAL')].groupby('site_year'):art.append(dict(site_year=sy,n=len(g),cultivars=g.cultivar_id.nunique(),severity_mean=g.disease_severity.mean(),yield_mean=g.observed_yield.mean(),spearman=rho(g.disease_severity,g.observed_yield),status='CONFOUNDED: artificial regime restricted to QLB NF; descriptive only; no separate artificial model'))
 csv(env,'results/metrics/within_environment_association.csv');csv(overall,'results/metrics/between_within_association_summary.csv');csv(mgmt,'results/metrics/management_association.csv');csv(art,'results/metrics/artificial_regime_descriptive.csv')
 print(pd.DataFrame(overall).to_string(index=False),flush=True)
if __name__=='__main__':associations()
