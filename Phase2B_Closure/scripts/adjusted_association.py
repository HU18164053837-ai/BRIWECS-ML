from common import *
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
import itertools

def variance_ratio(e,d,w):
 a=e[(w>0)&(d.disease_severity.to_numpy()>=d.disease_severity.median())];b=e[(w>0)&(d.disease_severity.to_numpy()<d.disease_severity.median())]
 return float(np.var(a)/max(np.var(b),1e-12)) if len(a)>1 and len(b)>1 else np.nan
def codes(d,model):
 if model=='A0':return [np.zeros(len(d),dtype=int)]
 group=d.site_year if model=='A1' else d.site_year+'|'+d.management
 g=[pd.factorize(group)[0]]
 if model=='A3':g.append(pd.factorize(d.cultivar_id)[0])
 return g
def residualize(z,groups,w):
 v=np.asarray(z,float).copy();v[w==0]=0
 for iteration in range(500):
  before=v.copy()
  for g in groups:
   den=np.bincount(g,weights=w,minlength=int(g.max())+1)
   for j in range(v.shape[1]):
    num=np.bincount(g,weights=w*v[:,j],minlength=len(den));mean=np.divide(num,den,out=np.zeros_like(num),where=den>0);v[:,j]-=mean[g]
   v[w==0]=0
  if np.max(np.abs(v-before))<1e-9:return v,iteration+1,True
 return v,500,False
def estimate(d,model,w=None):
 w=np.ones(len(d)) if w is None else w;g=codes(d,model);v,it,conv=residualize(d[['disease_severity','observed_yield']].to_numpy(),g,w);x,y=v[:,0],v[:,1];ss=np.dot(w*x,x)
 if ss<1e-8:return np.nan,dict(status='NOT_IDENTIFIABLE_INSUFFICIENT_RESIDUAL_SEVERITY',iterations=it,converged=conv),x,y
 b=np.dot(w*x,y)/ss;e=y-b*x
 ng=[len(np.unique(a[w>0])) for a in g]
 if len(g)==1:rank=ng[0]+1
 else:
  n1=int(g[0].max())+1;n2=int(g[1].max())+1;i=g[0][w>0];j=g[1][w>0]+n1;adj=coo_matrix((np.ones(2*len(i)),(np.r_[i,j],np.r_[j,i])),shape=(n1+n2,n1+n2));nc,labels=connected_components(adj,directed=False);active=np.unique(np.r_[i,j]);components=len(np.unique(labels[active]));rank=sum(ng)-components+1
 return b,dict(status='IDENTIFIABLE' if conv else 'PARTIAL_ABSORPTION_CONVERGENCE',iterations=it,converged=conv,absorbed_design_rank_including_severity=rank,residual_df=int((w>0).sum())-rank,residual_severity_sd=float(np.sqrt(np.average(x*x,weights=w))),residual_RMSE=float(np.sqrt(np.average(e*e,weights=w))),residual_skew=float(pd.Series(e[w>0]).skew()),residual_variance_ratio_hi_vs_lo_severity=variance_ratio(e,d,w)),x,y
def run():
 co=pd.read_parquet(P/'data/stripe_rust_paired_cohort.parquet');assert len(co)==8539 and co.disease_severity.eq(0).sum()==453
 cfg=json.loads((ROOT/'config/analysis.json').read_text());B=cfg['association_cluster_bootstrap'];rows=[];sensitive=[];within=[];bs=[]
 for scope,d in [('FULL',co),('NATURAL',co[co.inoculation_status.eq('NATURAL')])]:
  d=d.reset_index(drop=True);clusters=pd.factorize(d.site_year)[0];G=d.site_year.nunique();rng=np.random.default_rng(42)
  for model in ['A0','A1','A2','A3']:
   b,diag,x,y=estimate(d,model);boot=[]
   for k in range(B):
    count=np.bincount(rng.integers(0,G,G),minlength=G);bb,dd,_,_=estimate(d,model,count[clusters].astype(float));boot.append(bb);bs.append(dict(scope=scope,model=model,replicate=k,coefficient=bb,status=dd['status']))
   valid=np.asarray(boot);valid=valid[np.isfinite(valid)];lo,hi=np.quantile(valid,[.025,.975]) if len(valid) else (np.nan,np.nan)
   rows.append(dict(scope=scope,model=model,n=len(d),environments=G,cultivars=d.cultivar_id.nunique(),coefficient_per_score_unit=b,CI95_lo=lo,CI95_hi=hi,association_per_10_score_units=10*b,per10_CI95_lo=10*lo,per10_CI95_hi=10*hi,bootstrap_valid=len(valid),bootstrap_requested=B,uncertainty='Site-year pairs bootstrap; fixed effects re-absorbed each replicate; percentile descriptive CI, only 10/13 dependent environments',unit_interpretation='10 nominal reported severity score units; heterogeneous temporal endpoints; pooled not unconditional percentage-point effect',**diag))
   for removed in sorted(d.site_year.unique()):
    z=d[d.site_year.ne(removed)].reset_index(drop=True);bb,dd,_,_=estimate(z,model);sensitive.append(dict(scope=scope,model=model,analysis='DELETE_ONE_SITE_YEAR',removed=removed,n=len(z),environments=z.site_year.nunique(),coefficient=bb,association_per10=bb*10,**dd))
   print(scope,model,'beta=',b,'CI=',lo,hi,flush=True)
  for sy,g in d.groupby('site_year'):
   for model in ['A0','A2','A3']:
    b,diag,_,_=estimate(g.reset_index(drop=True),model);within.append(dict(scope=scope,site_year=sy,model=model,n=len(g),cultivars=g.cultivar_id.nunique(),coefficient=b,association_per10=b*10,scale='within same site-year reported leaf-area score; timing remains source-specific',**diag))
 # Same-site estimates address source-score semantics; no new predictive models.
 for site,g in co.groupby('location'):
  b,diag,_,_=estimate(g.reset_index(drop=True),'A3');sensitive.append(dict(scope='SITE_SPECIFIC',model='A3',analysis='SAME_SITE_SCORE_SENSITIVITY',removed=site,n=len(g),environments=g.site_year.nunique(),coefficient=b,association_per10=b*10,scale='per10 reported leaf-area score units at same site; no causal loss, no cross-site harmonization assumption',**diag))
 # QLB has a documented AO endpoint; keep its unit-specific estimate separate.
 q=co[co.location.eq('QLB')].reset_index(drop=True);cl=pd.factorize(q.site_year)[0];qb,qd,_,_=estimate(q,'A3');draws=[]
 for sample in itertools.product(range(3),repeat=3):
  bb,_,_,_=estimate(q,'A3',np.bincount(sample,minlength=3)[cl].astype(float));draws.append(bb)
 lo,hi=np.nanquantile(draws,[.025,.975])
 csv([dict(scope='QLB_COMPARABLE_AO',n=len(q),environments=3,model='A3',per10_AO_percentage_points=qb*10,per10_CI95_lo=lo*10,per10_CI95_hi=hi*10,uncertainty='Exact enumeration of 27 ordered site-year pairs resamples; ONLY THREE environments, highly discrete descriptive interval',**qd)],'QLB_comparable_AO_association.csv')
 csv(rows,'04_adjusted_association_models.csv');csv(within,'05_within_environment_effects.csv');csv(sensitive,'06_cultivar_management_sensitivity.csv');csv(bs,'association_cluster_bootstrap_draws.csv')
 # Verify absorbed A3 point estimate using an independently constructed dense dummy projection.
 for scope,d in [('FULL',co),('NATURAL',co[co.inoculation_status.eq('NATURAL')])]:
  Z=pd.get_dummies(pd.DataFrame({'env_management':d.site_year+'|'+d.management,'cultivar':d.cultivar_id}),dtype=float).to_numpy();v=d[['disease_severity','observed_yield']].to_numpy();coef=np.linalg.lstsq(Z,v,rcond=None)[0];r=v-Z@coef;beta=np.dot(r[:,0],r[:,1])/np.dot(r[:,0],r[:,0]);ref=next(a['coefficient_per_score_unit'] for a in rows if a['scope']==scope and a['model']=='A3');assert np.isclose(beta,ref,atol=1e-7)
 js(dict(status='PASS',independent_dense_FWL_agrees=True,bootstrap_method='environment pairs, not ordinary row bootstrap',replicates=B),'reports/association_verification.json')
if __name__=='__main__':run()
