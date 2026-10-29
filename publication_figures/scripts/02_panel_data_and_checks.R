# This file extends the previously verified Figure 3 arithmetic; it fits nothing.
cohort <- read_source('S1_Phase2A_cohort_exclusions.csv')
yield_n <- sum(cohort$selection_reason=='INCLUDED')
coverage <- subset(read_source('S2_Disease_cohort_coverage.csv'),disease=='Stripe_rust' & scope=='FULL')
cover <- data.frame(site_year=env,natural=0,artificial=0,total=0)
for(i in seq_along(env)){
 z<-subset(coverage,site_year==env[i]);cover$natural[i]<-sum(z$n[z$inoculation_status=='NATURAL']);cover$artificial[i]<-sum(z$n[z$inoculation_status=='ARTIFICIAL']);cover$total[i]<-sum(z$n)
}
cohort_display<-data.frame(cohort=c('Yield plots','Stripe-rust paired','Natural infection','QLB artificial-related'),plots=c(yield_n,sum(cover$total),sum(cover$natural),sum(cover$artificial)),site_years=c(15,13,10,3))
label_qc<-subset(read_source('S3_Disease_label_QC.csv'),disease=='Stripe_rust' & status=='VALID_ORIGINAL_VALUE')
restored<-sum(label_qc$raw_original_value>0 & label_qc$published_value==0);zeros<-sum(label_qc$raw_original_value==0)
split<-subset(read_source('S9_Disease_split_manifest.csv'),disease=='Stripe_rust' & scope=='FULL' & validation_strategy=='LOSYO' & test_group %in% env)
split<-split[match(env,split$test_group),]
foldmatrix<-expand.grid(site_year=env,held_out_site_year=env,stringsAsFactors=FALSE)
foldmatrix$role<-ifelse(foldmatrix$site_year==foldmatrix$held_out_site_year,'TEST','TRAIN')
manifest<-read_source('S8_Matched_domain_manifest.csv')
train_ranges<-do.call(rbind,lapply(c('A_NATURAL_ONLY','B_FULL_ELIGIBLE'),function(s){z<-subset(manifest,strategy==s);data.frame(strategy=s,min_train=min(z$train_n),max_train=max(z$train_n),test_plots=5899,test_site_years=10)}))
write_data(cohort_display,'Figure1_A_Cohort_Counts.csv');write_data(cover,'Figure1_B_Environment_Coverage.csv')
write_data(data.frame(restored_positive_fractional=restored,genuine_zeros=zeros),'Figure1_A_Label_QC.csv')
write_data(split,'Figure1_C_Actual_Folds.csv');write_data(foldmatrix,'Figure1_C_Fold_Matrix.csv');write_data(train_ranges,'Figure1_C_Matched_Training_Ranges.csv')
fig2env<-subset(paired,scope=='FULL' & comparison=='M1_vs_M0' & model_name %in% c('CatBoost','Ridge'))
fig2env<-fig2env[order(match(fig2env$model_name,c('CatBoost','Ridge')),match(fig2env$site_year,env)),]
performance<-read_source('S4_All_disease_model_results.csv');uncertainty<-read_source('S4b_Paired_uncertainty.csv')
ag2<-list()
for(alg in c('CatBoost','Ridge'))for(weight in c('POOLED','ENVIRONMENT_MACRO')){
 z<-subset(performance,disease=='Stripe_rust' & scope=='FULL' & validation_strategy=='LOSYO' & comparison=='M1_vs_M0' & model_name==alg & aggregation==weight)
 u<-subset(uncertainty,disease=='Stripe_rust' & scope=='FULL' & validation_strategy=='LOSYO' & comparison=='M1_vs_M0' & model_name==alg)
 stopifnot(nrow(z)==1,nrow(u)==1)
 prefix<-if(weight=='POOLED')'pooled' else 'macro'
 ag2[[length(ag2)+1]]<-data.frame(model_name=alg,weighting=weight,n=z$n,environments=z$environments,delta_MAE=z$delta_MAE,lo=u[[paste0(prefix,'_CI95_lo')]],hi=u[[paste0(prefix,'_CI95_hi')]])
}
fig2agg<-do.call(rbind,ag2)
deletion<-subset(read_source('S10_Drop_environment_sensitivity.csv'),disease=='Stripe_rust' & scope=='FULL' & validation_strategy=='LOSYO' & model_name=='CatBoost' & comparison=='M1_vs_M0' & removed_environment=='QLB_2015')
stopifnot(nrow(deletion)==1)
sens<-data.frame(analysis=c('All 13 site-years','QLB 2015 omitted'),environments=c(13,12),macro_delta_MAE=c(fig2agg$delta_MAE[fig2agg$model_name=='CatBoost' & fig2agg$weighting=='ENVIRONMENT_MACRO'],deletion$macro_delta_MAE))
write_data(fig2env,'Figure2_A_Environment_Increments.csv');write_data(fig2agg,'Figure2_B_Aggregates_And_Intervals.csv');write_data(sens,'Figure2_C_QLB2015_Sensitivity.csv')
domain<-read_source('S8c_Matched_domain_environment_comparison.csv')
natural_env<-env[!grepl('^QLB',env)]
domain<-domain[order(match(domain$model_name,c('CatBoost','Ridge')),match(domain$site_year,natural_env)),]
abserr<-read_source('S14_Matched_domain_absolute_MAE.csv')
metrics<-read_source('S8b_Matched_domain_metrics.csv')
write_data(domain,'Figure4_A_Matched_Environment_Increments.csv');write_data(abserr,'Figure4_B_Absolute_MAE.csv');write_data(metrics,'Figure4_C_Macro_Increments_And_Contrasts.csv');write_data(manifest,'Figure4_Matched_Test_Fingerprints.csv')
# Audits against sources and the manuscript's frozen numeric claims.
audit('F1 yield plots',yield_n,24991,0);audit('F1 paired plots',sum(cover$total),8539,0);audit('F1 natural plots',sum(cover$natural),5899,0);audit('F1 artificial-related plots',sum(cover$artificial),2640,0)
audit('F1 restored fractions',restored,199,0);audit('F1 true zeros',zeros,453,0)
audit('F1 13 complete folds',nrow(split),13,0);audit('F1 fold-matrix rows',nrow(foldmatrix),169,0)
for(i in 1:13){audit(paste('F1',env[i],'test count'),split$test_n[i],cover$total[i],0);audit(paste('F1',env[i],'train count'),split$train_n[i],8539-cover$total[i],0);audit(paste('F1',env[i],'train site-years'),split$train_site_years[i],12,0)}
audit('F2 both algorithms all environments',nrow(fig2env),26,0)
for(alg in c('CatBoost','Ridge')){
 z<-subset(fig2env,model_name==alg);audit(paste('F2',alg,'environments'),nrow(z),13,0)
 for(w in c('POOLED','ENVIRONMENT_MACRO')){
  v<-if(w=='POOLED')weighted.mean(z$delta_MAE,z$n) else mean(z$delta_MAE)
  target<-subset(fig2agg,model_name==alg & weighting==w)$delta_MAE
  audit(paste('F2',alg,w,'MAE aggregation'),v,target)
 }
}
audit('F2 all-environment sensitivity reference',sens$macro_delta_MAE[1],.434,.0005);audit('F2 QLB deletion sensitivity',sens$macro_delta_MAE[2],-.143,.0005)
audit('F4 20 algorithm-environment rows',nrow(domain),20,0);audit('F4 eight absolute-error records',nrow(abserr),8,0);audit('F4 six macro/contrast records',nrow(metrics),6,0)
for(alg in c('CatBoost','Ridge')){
 z<-subset(domain,model_name==alg);stopifnot(setequal(z$site_year,natural_env),all(z$n_A==z$n_B))
 audit(paste('F4',alg,'natural tests'),sum(z$n_A),5899,0)
 for(i in 1:nrow(z)) for(s in c('A','B'))audit(paste('F4',alg,z$site_year[i],s,'absolute increment'),z[[paste0('MAE_without_',s)]][i]-z[[paste0('MAE_with_',s)]][i],z[[paste0('delta_MAE_',s)]][i])
 for(s in c('A','B'))for(w in c('POOLED','ENVIRONMENT_MACRO')){
  target_weight<-if(w=='POOLED')'Pooled plots' else 'Equal site-year'
  a<-abserr[abserr$algorithm==alg & abserr$strategy==s & abserr$weighting==target_weight,]
  stopifnot(nrow(a)==1)
  f<-if(w=='POOLED')function(v)weighted.mean(v,z$n_A) else mean
  audit(paste('F4',alg,s,w,'M0'),f(z[[paste0('MAE_without_',s)]]),a$M0_MAE)
  audit(paste('F4',alg,s,w,'M1'),f(z[[paste0('MAE_with_',s)]]),a$M1_MAE)
  audit(paste('F4',alg,s,w,'increment'),a$M0_MAE-a$M1_MAE,a$delta_MAE)
 }
 a<-subset(metrics,model_name==alg & strategy=='A_NATURAL_ONLY');b<-subset(metrics,model_name==alg & strategy=='B_FULL_ELIGIBLE');ba<-subset(metrics,model_name==alg & strategy=='B_MINUS_A_INCREMENT')
 audit(paste('F4',alg,'B-A macro contrast'),b$macro_delta_MAE-a$macro_delta_MAE,ba$macro_delta_MAE)
}
fingerprint_groups<-base::split(manifest,paste(manifest$model_name,manifest$site_year,sep='|'))
stopifnot(length(fingerprint_groups)==20)
for(key in names(fingerprint_groups)){
 z<-fingerprint_groups[[key]];stopifnot(nrow(z)==4,length(unique(z$test_id_sha256))==1,length(unique(z$test_n))==1)
 audit(paste('F4 identical A/B M0/M1 test fingerprints',key),length(unique(z$test_id_sha256)),1,0)
}
# The aggregates are within-environment component averages, not a pooled global bias decomposition.
agg$definition<-c('sum_h component_h / 13','sum_h n_h * component_h / 8539')
agg$interpretation<-'Average of within-site-year decompositions; not global pooled MBE squared / global error variance'
write_data(agg,'Figure3_B_Weighted_Components.csv')
write.csv(qa,file.path(root,'qa','R_Numerical_Consistency_Checks.csv'),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),'','Font files:',font_table$path[font_table$family==family]),file.path(root,'reports','R_Session_And_Fonts.txt'))
