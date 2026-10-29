from common import *
def diagnostics():
 d=pd.read_parquet(ROOT/'results/predictions/disease_pilot_oof_predictions.parquet');d=d[d.disease.eq('Stripe_rust')&d.scope.eq('FULL')&d.validation_strategy.eq('LOSYO')];rows=[]
 for model,g in d.groupby('model_name'):
  for b,v in [('M0','M1'),('M2','M3')]:
   x=g[g.context.eq(b)].merge(g[g.context.eq(v)][['sample_id','fold_id','predicted_yield']],on=['sample_id','fold_id'],suffixes=('_b','_v'),validate='one_to_one');x['absolute_error_reduction']=abs(x.predicted_yield_b-x.observed_yield)-abs(x.predicted_yield_v-x.observed_yield);total=x.absolute_error_reduction.sum()
   for regime,h in x.groupby('inoculation_status'):
    rows.append(dict(model_name=model,comparison=f'{v}_vs_{b}',regime=regime,n=len(h),environments=h.site_year.nunique(),pooled_delta_MAE=h.absolute_error_reduction.mean(),macro_delta_MAE=h.groupby('site_year').absolute_error_reduction.mean().mean(),total_absolute_error_reduction=h.absolute_error_reduction.sum(),signed_fraction_of_full_total=h.absolute_error_reduction.sum()/total if total!=0 else np.nan,interpretation='Partition of SAME full-cohort fitted OOF predictions; not separately fitted natural model; signed fractions may exceed 1 if another group worsens; CONFOUNDED'))
 csv(rows,'results/metrics/inoculation_regime_contribution.csv')
 e=pd.read_csv(ROOT/'results/metrics/environment_incremental_metrics.csv');x=e[e.validation_strategy.eq('LOSYO')].copy();x['absolute_R2_change']=abs(x.delta_R2);x=x.sort_values('absolute_R2_change',ascending=False);csv(x[['disease','scope','model_name','comparison','site_year','n','R2_without','R2_with','delta_R2','delta_MAE']],'results/metrics/R2_denominator_diagnostics.csv')
if __name__=='__main__':diagnostics()
