#!/usr/bin/env Rscript
# Deterministic plotting of frozen results; no model fitting or inference.
suppressPackageStartupMessages({library(grid);library(svglite);library(ragg)})
args0 <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args0[grepl('^--file=',args0)][1])
root <- normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
src <- file.path(root,'..','frozen_SI')
options(digits=17)
read_source <- function(name) read.csv(file.path(src,name),check.names=FALSE,fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE)
write_data <- function(x,name) write.csv(x,file.path(root,'source_data','panels',name),row.names=FALSE,fileEncoding='UTF-8',na='')
font_table <- systemfonts::system_fonts()
stopifnot(any(font_table$family=='Times New Roman'))
family <- 'Times New Roman'
env <- c('GGE_2015','GGE_2016','GGE_2017','HAN_2015','HAN_2016','HAN_2017','KAL_2015','KAL_2016','KIE_2015','KIE_2016','QLB_2015','QLB_2016','QLB_2017')
labels <- sub('_',' ',env)
paired <- read_source('S6b_Paired_error_components.csv')
d <- subset(paired,scope=='FULL' & model_name=='CatBoost' & comparison=='M1_vs_M0')
stopifnot(nrow(d)==13, !anyDuplicated(d$site_year), setequal(d$site_year,env))
d <- d[match(env,d$site_year),]
absolute <- read_source('S6_Environment_error_decomposition.csv')
m0 <- subset(absolute,scope=='FULL' & model_name=='CatBoost' & context=='M0')
m1 <- subset(absolute,scope=='FULL' & model_name=='CatBoost' & context=='M1')
m0 <- m0[match(env,m0$site_year),];m1 <- m1[match(env,m1$site_year),]
stopifnot(nrow(m0)==13,nrow(m1)==13,all(m0$n==m1$n),all(m0$n==d$n),sum(d$n)==8539)
ranking <- read_source('S6c_Within_environment_ranking.csv')
r <- subset(ranking,scope=='FULL' & model_name=='CatBoost' & context %in% c('M0','M1'))
stopifnot(nrow(r)==52)
rho <- data.frame(site_year=env,plot_M0=NA_real_,plot_M1=NA_real_,cultivar_M0=NA_real_,cultivar_M1=NA_real_)
for (lv in c('PLOT','CULTIVAR_MEAN')) for (cx in c('M0','M1')) {
  z <- subset(r,level==lv & context==cx);stopifnot(nrow(z)==13,!anyDuplicated(z$site_year))
  rho[[paste(if(lv=='PLOT') 'plot' else 'cultivar',cx,sep='_')]] <- z$spearman[match(env,z$site_year)]
}
rho$delta_plot <- rho$plot_M1-rho$plot_M0;rho$delta_cultivar <- rho$cultivar_M1-rho$cultivar_M0
agg <- data.frame(weighting=c('Equal site-year','Plot-weighted'),environments=13,plots=8539,
  delta_bias_squared=c(mean(d$delta_bias_squared),weighted.mean(d$delta_bias_squared,d$n)),
  delta_centered_variance=c(mean(d$delta_centered_variance),weighted.mean(d$delta_centered_variance,d$n)),
  delta_MSE=c(mean(d$delta_MSE),weighted.mean(d$delta_MSE,d$n)))
write_data(d,'Figure3_A_Environment_Components.csv');write_data(agg,'Figure3_B_Weighted_Components.csv');write_data(rho,'Figure3_C_Ranking_Changes.csv')
write_data(r,'Figure3_C_Original_Rankings.csv');write_data(rbind(m0,m1),'Figure3_A_Absolute_Components.csv')
qa <- data.frame(check=character(),observed=numeric(),expected=numeric(),tolerance=numeric(),status=character())
audit <- function(check,obs,expected=0,tol=1e-10) {
  stopifnot(length(obs)==1,is.finite(obs));ok <- abs(obs-expected)<=tol
  qa <<- rbind(qa,data.frame(check=check,observed=obs,expected=expected,tolerance=tol,status=if(ok)'PASS' else 'ISSUE'))
  if(!ok)stop(check)
}
audit('13 complete primary environments',nrow(d),13,0);audit('8539 primary plots',sum(d$n),8539,0)
for(i in seq_along(env)) {
  audit(paste(env[i],'M0 exact identity'),m0$MSE[i]-m0$bias_squared[i]-m0$error_variance[i])
  audit(paste(env[i],'M1 exact identity'),m1$MSE[i]-m1$bias_squared[i]-m1$error_variance[i])
  audit(paste(env[i],'MSE subtraction'),m0$MSE[i]-m1$MSE[i]-d$delta_MSE[i])
  audit(paste(env[i],'bias subtraction'),m0$bias_squared[i]-m1$bias_squared[i]-d$delta_bias_squared[i])
  audit(paste(env[i],'variance subtraction'),m0$error_variance[i]-m1$error_variance[i]-d$delta_centered_variance[i])
  audit(paste(env[i],'paired identity'),d$delta_MSE[i]-d$delta_bias_squared[i]-d$delta_centered_variance[i])
  for(lv in c('plot','cultivar'))audit(paste(env[i],lv,'rank subtraction'),rho[[paste0(lv,'_M1')]][i]-rho[[paste0(lv,'_M0')]][i]-rho[[paste0('delta_',lv)]][i])
}
audit('Plot ranking improvements',sum(rho$delta_plot>0),5,0);audit('Cultivar ranking improvements',sum(rho$delta_cultivar>0),4,0)
audit('Bias component improvements',sum(d$delta_bias_squared>0),8,0);audit('Variance component improvements',sum(d$delta_centered_variance>0),6,0)
audit('Equal-site-year bias',agg$delta_bias_squared[1],45.227,0.0005)
audit('Equal-site-year variance',agg$delta_centered_variance[1],-17.304,0.0005)
audit('Equal-site-year net',agg$delta_MSE[1],27.923,0.0005)
audit('Plot-weighted bias',agg$delta_bias_squared[2],51.152,0.0005)
audit('Plot-weighted variance',agg$delta_centered_variance[2],11.236,0.0005)
audit('Plot-weighted net',agg$delta_MSE[2],62.388,0.0005)
for(i in 1:2) audit(paste(agg$weighting[i],'identity'),agg$delta_MSE[i]-agg$delta_bias_squared[i]-agg$delta_centered_variance[i])
write.csv(qa,file.path(root,'qa','Figure3_Data_Equivalence_Audit.csv'),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),'','Font files:',font_table$path[font_table$family==family]),file.path(root,'reports','R_Session_And_Fonts.txt'))

