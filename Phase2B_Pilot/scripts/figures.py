from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
matplotlib.rcParams.update({'font.size':9,'savefig.dpi':170})
def save(fig,name):
 for suffix in ['png','pdf']:fig.savefig(ROOT/'results/figures'/f'{name}.{suffix}',bbox_inches='tight')
 plt.close(fig)
def figures():
 env=pd.read_csv(ROOT/'results/metrics/environment_incremental_metrics.csv');u=pd.read_csv(ROOT/'results/metrics/paired_uncertainty_analysis.csv');a=pd.read_csv(ROOT/'results/metrics/within_environment_association.csv')
 u=u[u.validation_strategy.eq('LOSYO')];e=env[env.validation_strategy.eq('LOSYO')&env.scope.eq('FULL')];cover=pd.read_csv(ROOT/'data/disease_environment_coverage.csv')
 fig,axes=plt.subplots(1,2,figsize=(12,5.5),sharey=True)
 order=sorted({f'{loc}_{year}' for loc in ['GGE','HAN','KAL','KIE','QLB'] for year in [2015,2016,2017]})
 colors={'NATURAL':'#3a9474','ARTIFICIAL':'#d38a36','UNKNOWN':'#999999','MIXED_OR_UNCERTAIN':'#c6a2cf'}
 for ax,dis in zip(axes,['Stripe_rust','Septoria']):
  z=cover[cover.disease.eq(dis)&cover.scope.eq('FULL')].groupby(['site_year','inoculation_status']).n.sum().unstack(fill_value=0).reindex(order,fill_value=0);left=np.zeros(len(order))
  for label,color in colors.items():
   n=z[label].to_numpy() if label in z else np.zeros(len(order));ax.barh(order,n,left=left,label=label,color=color);left+=n
  for i,n in enumerate(left):ax.text(n+7,i,str(int(n)) if n else '—',va='center',fontsize=8)
  ax.set_title(dis.replace('_',' '));ax.set_xlabel('Actual paired plots');ax.invert_yaxis()
 axes[1].legend(handles=[Patch(color=colors[s],label=s) for s in ['NATURAL','ARTIFICIAL']],loc='lower right',fontsize=8)
 fig.suptitle('B1 Disease cohort coverage (— = no paired observations; not zero disease)');fig.tight_layout();save(fig,'Figure_B1_Disease_Cohort_Coverage')
 full=u[u.scope.eq('FULL')].copy();full['label']=full.disease+' | '+full.model_name+' | '+full.comparison;full=full.sort_values(['disease','model_name','comparison']).reset_index(drop=True)
 fig,ax=plt.subplots(figsize=(11,5.7));yy=np.arange(len(full));err=np.row_stack([full.macro_delta_MAE-full.macro_CI95_lo,full.macro_CI95_hi-full.macro_delta_MAE]);ax.errorbar(full.macro_delta_MAE,yy,xerr=err,fmt='o',capsize=3,label='Site-year macro, 95% descriptive cluster interval');ax.scatter(full.pooled_delta_MAE,yy+.14,marker='D',s=22,color='#c7782e',label='Pooled');ax.axvline(0,color='black',ls='--',lw=.8);ax.set_yticks(yy,full.label);ax.invert_yaxis();ax.set_xlabel('ΔMAE = MAE(without disease) − MAE(with disease), dt/ha; positive improves');ax.set_title('B2 Genuine matched LOSYO incremental performance');ax.legend(fontsize=8);fig.tight_layout();save(fig,'Figure_B2_Incremental_MAE')
 limit=max(abs(e.delta_MAE.min()),abs(e.delta_MAE.max()));fig,axes=plt.subplots(1,2,figsize=(12,7),layout='constrained')
 for ax,dis in zip(axes,['Stripe_rust','Septoria']):
  z=e[e.disease.eq(dis)].copy();z['label']=z.model_name+' '+z.comparison;v=z.pivot(index='site_year',columns='label',values='delta_MAE').reindex(order);im=ax.imshow(np.ma.masked_invalid(v.to_numpy()),aspect='auto',cmap='RdBu',vmin=-limit,vmax=limit);ax.set_xticks(range(len(v.columns)),v.columns,rotation=50,ha='right');ax.set_yticks(range(len(v)),v.index);ax.set_title(dis.replace('_',' '))
  for i in range(len(v)):
   for j in range(len(v.columns)):ax.text(j,i,'—' if pd.isna(v.iloc[i,j]) else f'{v.iloc[i,j]:.1f}',ha='center',va='center',fontsize=8)
 fig.colorbar(im,ax=axes,label='ΔMAE (dt/ha), positive improves',shrink=.8);fig.suptitle('B3 Benefit heterogeneity across all planned site-years');save(fig,'Figure_B3_Environment_Benefit_Heterogeneity')
 fig,axes=plt.subplots(1,2,figsize=(12,6.5),layout='constrained')
 for ax,dis in zip(axes,['Stripe_rust','Septoria']):
  z=a[a.disease.eq(dis)&a.scope.eq('FULL')].sort_values('site_year').reset_index(drop=True);ok=z.spearman.notna();y=np.arange(len(z));err=np.row_stack([z.loc[ok,'spearman']-z.loc[ok,'CI95_lo'],z.loc[ok,'CI95_hi']-z.loc[ok,'spearman']]);err=np.maximum(err,0)
  ax.errorbar(z.loc[ok,'spearman'],y[ok],xerr=err,fmt='o',capsize=2);ax.set_yticks(y,[f'{r.site_year} (n={r.n})' for r in z.itertuples()]);ax.invert_yaxis();ax.axvline(0,color='black',ls='--',lw=.8);ax.set_xlim(-1.05,1.05);ax.set_title(dis.replace('_',' '));ax.set_xlabel('Within-site-year Spearman ρ; cultivar-cluster interval')
  for i,row in z.loc[~ok].iterrows():ax.text(0,i,'UNDEFINED',ha='center')
 fig.suptitle('B4 Disease–yield association: 500 cultivar-cluster resamples; not causal');save(fig,'Figure_B4_Within_Environment_Association')
 z=u[u.disease.eq('Stripe_rust')].copy();labels=sorted(set(z.model_name+' '+z.comparison));fig,ax=plt.subplots(figsize=(11,5))
 for scope,off,color in [('FULL',-.12,'#4377a1'),('NATURAL',.12,'#3a9474')]:
  zz=z[z.scope.eq(scope)].copy();zz['label']=zz.model_name+' '+zz.comparison;zz=zz.set_index('label').reindex(labels);y=np.arange(len(labels))+off;err=np.row_stack([zz.macro_delta_MAE-zz.macro_CI95_lo,zz.macro_CI95_hi-zz.macro_delta_MAE]);ax.errorbar(zz.macro_delta_MAE,y,xerr=err,fmt='o',capsize=3,color=color,label=f'{scope}: {int(zz.n.iloc[0])} plots, {int(zz.environments.iloc[0])} site-years')
 ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.axvline(0,color='black',ls='--',lw=.8);ax.set_xlabel('Matched macro ΔMAE (dt/ha); positive improves');ax.set_title('B5 Stripe-rust inoculation sensitivity: separately refitted matched controls');ax.legend(fontsize=8);fig.text(.5,.01,'Different environment coverage; FULL minus NATURAL is not an inoculation causal effect.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,1));save(fig,'Figure_B5_Inoculation_Sensitivity')
if __name__=='__main__':figures()
