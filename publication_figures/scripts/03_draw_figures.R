# Style inherited from the approved diverging-accounting candidate.
ink<-'#20262B';neutral<-'#6F7780';pale<-'#D8DDE0';family<-'Times New Roman';render_scale<-1;grayscale_render<-FALSE
paint<-function(col){if(!grayscale_render || is.na(col))return(col);grDevices::gray(sum((grDevices::col2rgb(col)/255)*c(.2126,.7152,.0722)))}
pal_color<-c(bias='#35576B',variance='#BC775C',net='#20262B')
pal_gray<-c(bias='#3E3E3E',variance='#777777',net='#111111')
txt<-function(s,x,y,size=7.6,just='left',bold=FALSE,col=ink)grid.text(s,x,y,default.units='native',just=just,gp=gpar(fontfamily=family,fontsize=size*render_scale,fontface=if(bold)'bold' else 'plain',col=paint(col)))
ln<-function(x0,y0,x1,y1,col=pale,lwd=.75,lty=1)grid.lines(c(x0,x1),c(y0,y1),default.units='native',gp=gpar(col=paint(col),lwd=lwd*render_scale,lty=lty,lineend='butt'))
pt<-function(x,y,pch=16,col=ink,size=1.55,fill=col)grid.points(x,y,default.units='native',pch=pch,size=unit(size*render_scale,'mm'),gp=gpar(col=paint(col),fill=paint(fill),lwd=.9*render_scale))
rect<-function(x0,x1,y,h,fill,col=NA)grid.rect(x=(x0+x1)/2,y=y,width=abs(x1-x0),height=h,default.units='native',gp=gpar(fill=paint(fill),col=paint(col),lwd=.6*render_scale))
mp<-function(v,lim,left,right)left+(v-lim[1])/diff(lim)*(right-left)
axis_x<-function(left,right,y,lim,ticks,label=NULL,label_y=y+8){
 ln(left,y,right,y,ink,.8)
 for(v in ticks){x<-mp(v,lim,left,right);ln(x,y,x,y+1.1,ink,.8);txt(format(v,digits=3,trim=TRUE,scientific=FALSE),x,y+3.4,7.4,'centre')}
 if(!is.null(label))txt(label,(left+right)/2,label_y,7.6,'centre')
}
panel<-function(letter,title,x,y){txt(letter,x,y,10,bold=TRUE);txt(title,x+5,y,8.6,bold=TRUE)}
zero<-function(left,right,lim,top,bottom){x<-mp(0,lim,left,right);ln(x,top,x,bottom,neutral,.8)}
separator<-function(left,right,ys,cuts=c(3,6,8,10))for(i in cuts)ln(left,mean(ys[c(i,i+1)]),right,mean(ys[c(i,i+1)]),pale,.5)
alg_colors<-function(gray)c(CatBoost=if(gray)'#333333' else '#35576B',Ridge=if(gray)'#777777' else '#BC775C')
hatch<-function(x0,x1,y,h,col=ink){
 rect(x0,x1,y,h,'white',col)
 if(x1>x0)for(x in seq(x0+.5,x1,by=1.8)){a<-max(x0,x-.6);b<-min(x1,x+.6);ln(a,y+h/2,b,y-h/2,col,.55)}
}
draw3<-function(gray=FALSE){
 grayscale_render<<-gray
 pal<-if(gray)pal_gray else pal_color
 pushViewport(viewport(xscale=c(0,178),yscale=c(174,0),clip='off'))
 panel('A','Within-site-year error components',2,6);panel('C','Rank changes',128,6)
 rect(26,30,14,1.25,pal['bias']);txt('ΔMBE²',31.5,14,7.6)
 rect(53,57,14,1.25,pal['variance']);txt('ΔVar(e)',58.5,14,7.6)
 pt(89,14,18,pal['net'],1.65);txt('ΔMSE',92,14,7.6)
 pt(130,13,16,ink,1.4);txt('Plot: 5/13 improved',133,13,7.3)
 pt(130,19,0,ink,1.5,fill='white');txt('Cultivar mean: 4/13',133,19,7.3)
 ys<-seq(29,95,length.out=13);left<-26;right<-119;lim<-c(-140,260)
 separator(left,right,ys);separator(133,174,ys);zero(left,right,lim,25,99)
 rz<-mp(0,c(-.42,.19),133,174);ln(rz,25,rz,99,neutral,.8)
 for(i in 1:13){
  txt(labels[i],23,ys[i],7.6,'right')
  xb<-mp(d$delta_bias_squared[i],lim,left,right);xv<-mp(d$delta_centered_variance[i],lim,left,right);x0<-mp(0,lim,left,right)
  rect(x0,xb,ys[i]-1.05,1.25,pal['bias']);rect(x0,xv,ys[i]+1.05,1.25,pal['variance'])
  pt(xb,ys[i]-1.05,16,pal['bias'],.85);pt(xv,ys[i]+1.05,17,pal['variance'],1)
  pt(mp(d$delta_MSE[i],lim,left,right),ys[i],18,pal['net'],1.6)
  if(abs(d$delta_bias_squared[i])<6)txt(sprintf('%+.3f',d$delta_bias_squared[i]),xb+2,ys[i]-1.3,7.2,col=pal['bias'])
  if(abs(d$delta_centered_variance[i])<6)txt(sprintf('%+.3f',d$delta_centered_variance[i]),xv+3,ys[i]+1.9,7.2,col=pal['variance'])
  xp<-mp(rho$delta_plot[i],c(-.42,.19),133,174);xc<-mp(rho$delta_cultivar[i],c(-.42,.19),133,174)
  ln(xp,ys[i],xc,ys[i],pale,.65);pt(xp,ys[i]-.4,16,ink,1.35);pt(xc,ys[i]+.4,0,ink,1.5,fill='white')
 }
 axis_x(left,right,102,lim,c(-100,0,100,200),'M0 − M1 change [(dt ha⁻¹)²]',111)
 axis_x(133,174,102,c(-.42,.19),c(-.4,-.2,0,.1),'Δ Spearman ρ',111)
 panel('B','Weighted means of within-site-year components',2,122)
 bl<-31;br<-164;lim2<-c(-25,70);zero(bl,br,lim2,129,154)
 for(i in 1:2){
  y<-c(136,150)[i];txt(if(i==1)'Equal site-year\nwₕ = 1/13' else 'Plot-count weights\nwₕ = nₕ/8,539',27,y,7.3,'right')
  vals<-c(agg$delta_bias_squared[i],agg$delta_centered_variance[i],agg$delta_MSE[i]);x0<-mp(0,lim2,bl,br)
  xb<-mp(vals[1],lim2,bl,br);xv<-mp(vals[2],lim2,bl,br);xn<-mp(vals[3],lim2,bl,br)
  rect(x0,xb,y-1.8,if(i==1)1.8 else 1.3,pal['bias']);rect(x0,xv,y+1.8,if(i==1)1.8 else 1.3,pal['variance']);pt(xn,y,18,pal['net'],if(i==1)1.9 else 1.5)
  txt(sprintf('%+.3f',vals[1]),xb+1.5,y-2,if(i==1)8 else 7.4,bold=i==1,col=pal['bias'])
  txt(sprintf('%+.3f',vals[2]),if(vals[2]<0)xv-1.4 else xv+1.5,y+2,if(i==1)8 else 7.4,if(vals[2]<0)'right' else 'left',bold=i==1,col=pal['variance'])
  txt(sprintf('%+.3f',vals[3]),xn+1.5,y+1.1,if(i==1)8 else 7.4,bold=i==1)
 }
 axis_x(bl,br,158,lim2,c(-20,0,20,40,60),'Weighted mean of within-site-year changes [(dt ha⁻¹)²]',168)
 popViewport()
}
draw1<-function(gray=FALSE){
 grayscale_render<<-gray
 natural_col<-if(gray)'#777777' else '#7D8E97'
 pushViewport(viewport(xscale=c(0,178),yscale=c(198,0),clip='off'))
 panel('A','Cohort sizes',2,6)
 al<-40;ar<-155;lim<-c(0,30000);yy<-c(17,28,39,50)
 for(i in 1:4){
  txt(cohort_display$cohort[i],36,yy[i],7.6,'right');x<-mp(cohort_display$plots[i],lim,al,ar)
  if(i==4)hatch(al,x,yy[i],3,ink)else rect(al,x,yy[i],3,if(i==1)'#D8DDE0' else natural_col)
  txt(format(cohort_display$plots[i],big.mark=',',trim=TRUE),x+1.7,yy[i],7.5,bold=TRUE)
  txt(paste0(cohort_display$site_years[i],' site-years'),174,yy[i],7.4,'right')
 }
 axis_x(al,ar,57,lim,c(0,10000,20000,30000),'Plots',65)
 txt('Score QC: 199 positive fractions restored; 453 genuine zero scores',40,73,7.4)
 panel('B','Environment coverage',2,84);panel('C','LOSYO holdouts',103,84)
 rect(27,31,92,1.5,natural_col);txt('Natural',33,92,7.3);hatch(56,60,92,1.5);txt('QLB artificial-related',62,92,7.3)
 ys<-seq(103,169,length.out=13);bl<-26;br<-89;lim2<-c(0,1400);separator(bl,br,ys)
 for(i in 1:13){
  txt(labels[i],23,ys[i],7.3,'right');xn<-mp(cover$natural[i],lim2,bl,br);xt<-mp(cover$total[i],lim2,bl,br)
  rect(bl,xn,ys[i],2.8,natural_col);if(cover$artificial[i]>0)hatch(xn,xt,ys[i],2.8)
  txt(as.character(cover$total[i]),92,ys[i],7.2)
 }
 axis_x(bl,br,176,lim2,c(0,700,1400),'Plots',185)
 # The cells represent complete eligible site-year groups in actual LOSYO folds.
 ml<-125;mt<-101;cell<-3.6
 for(i in 1:13){
  txt(labels[i],122,mt+(i-.5)*cell,7.1,'right')
  for(j in 1:13)rect(ml+(j-1)*cell+.15,ml+j*cell-.15,mt+(i-.5)*cell,cell-.3,if(i==j)ink else '#E7EAEC')
 }
 for(j in 1:13)txt(as.character(j),ml+(j-.5)*cell,151,7,'centre')
 txt('Held-out fold (row order)',148.5,157,7.2,'centre')
 rect(109,112,164,2.1,'#E7EAEC');txt('Train',114,164,7.2);rect(141,144,164,2.1,ink);txt('Test',146,164,7.2)
 txt('Same M0/M1 records in each fold',105,173,7.3)
 txt('Natural tests: 5,899 plots / 10 site-years',105,180,7.3)
 txt(sprintf('A train/fold: %s–%s',format(train_ranges$min_train[1],big.mark=',',trim=TRUE),format(train_ranges$max_train[1],big.mark=',',trim=TRUE)),105,186,7.3)
 txt(sprintf('B train/fold: %s–%s',format(train_ranges$min_train[2],big.mark=',',trim=TRUE),format(train_ranges$max_train[2],big.mark=',',trim=TRUE)),105,192,7.3)
 popViewport()
}
draw2<-function(gray=FALSE){
 grayscale_render<<-gray
 colors<-alg_colors(gray);pushViewport(viewport(xscale=c(0,178),yscale=c(174,0),clip='off'))
 panel('A','Environment-level ΔMAE',2,6);panel('B','Aggregate ΔMAE',112,6)
 pt(28,14,16,colors['CatBoost']);txt('CatBoost',31,14,7.5);pt(61,14,15,colors['Ridge']);txt('Ridge',64,14,7.5)
 ys<-seq(26,98,length.out=13);left<-26;right<-98;lim<-c(-16,8);zero(left,right,lim,22,102);separator(left,right,ys)
 cb<-subset(fig2env,model_name=='CatBoost');rg<-subset(fig2env,model_name=='Ridge')
 for(i in 1:13){
  txt(labels[i],23,ys[i],7.6,'right');a<-cb$delta_MAE[match(env[i],cb$site_year)];b<-rg$delta_MAE[match(env[i],rg$site_year)]
  xa<-mp(a,lim,left,right);xb<-mp(b,lim,left,right);ln(xa,ys[i],xb,ys[i],pale,.7)
  pt(xa,ys[i]-.5,16,colors['CatBoost']);pt(xb,ys[i]+.5,15,colors['Ridge'])
 }
 stopifnot(all(fig2env$delta_MAE>=lim[1] & fig2env$delta_MAE<=lim[2]))
 axis_x(left,right,106,lim,c(-15,-10,-5,0,5),'ΔMAE = MAE(M0) − MAE(M1) [dt ha⁻¹]',115)
 bl<-119;br<-173;lim2<-range(c(fig2agg$lo,fig2agg$hi));lim2<-c(floor(lim2[1]-.1),ceiling(lim2[2]+.1));zero(bl,br,lim2,22,102)
 by<-c(30,50,73,93)
 for(i in 1:4){
  z<-fig2agg[i,];co<-colors[z$model_name];txt(paste(z$model_name,if(z$weighting=='POOLED')'· pooled' else '· equal site-year'),114,by[i]-6,7.3)
  lo<-mp(z$lo,lim2,bl,br);hi<-mp(z$hi,lim2,bl,br);x<-mp(z$delta_MAE,lim2,bl,br)
  ln(lo,by[i],hi,by[i],co,.9);ln(lo,by[i]-1,lo,by[i]+1,co,.9);ln(hi,by[i]-1,hi,by[i]+1,co,.9)
  pt(x,by[i],if(z$model_name=='CatBoost')16 else 15,co);txt(sprintf('%+.3f',z$delta_MAE),br,by[i]+4,7.2,'right')
 }
 axis_x(bl,br,106,lim2,pretty(lim2,n=4),'ΔMAE [dt ha⁻¹]',115)
 panel('C','QLB 2015 deletion sensitivity',2,127)
 cl<-49;cr<-162;lc<-c(-.22,.52);zero(cl,cr,lc,134,153)
 for(i in 1:2){y<-c(139,150)[i];txt(if(i==1)'All site-years (13)' else 'QLB 2015 omitted (12)',44,y,7.3,'right');x<-mp(sens$macro_delta_MAE[i],lc,cl,cr);pt(x,y,16,colors['CatBoost'],1.7);txt(sprintf('%+.3f',sens$macro_delta_MAE[i]),x+2.3,y,7.5,bold=TRUE)}
 axis_x(cl,cr,158,lc,c(-.2,0,.2,.4),'Equal-site-year ΔMAE [dt ha⁻¹]',168)
 popViewport()
}
draw4<-function(gray=FALSE){
 grayscale_render<<-gray
 colors<-alg_colors(gray);pushViewport(viewport(xscale=c(0,178),yscale=c(201,0),clip='off'))
 panel('A','Same-domain environment increments',2,6)
 txt('CatBoost',58,21,8,'centre',TRUE,col=colors['CatBoost']);txt('Ridge',145,21,8,'centre',TRUE,col=colors['Ridge'])
 pt(27,14,0,ink,1.45,fill='white');txt('A: natural-only training',30,14,7.2);pt(77,14,15,ink,1.45);txt('B: all-eligible training',80,14,7.2)
 ys<-seq(28,82,length.out=10);lim<-c(-16,6)
 for(alg in c('CatBoost','Ridge')){
  left<-if(alg=='CatBoost')26 else 118;right<-if(alg=='CatBoost')89 else 174;labelx<-left-3
  zero(left,right,lim,24,86);separator(left,right,ys,c(3,6,8));z<-subset(domain,model_name==alg);z<-z[match(natural_env,z$site_year),]
  stopifnot(all(c(z$delta_MAE_A,z$delta_MAE_B)>=lim[1]),all(c(z$delta_MAE_A,z$delta_MAE_B)<=lim[2]))
  for(i in 1:10){
   txt(sub('_',' ',natural_env[i]),labelx,ys[i],7.3,'right');xa<-mp(z$delta_MAE_A[i],lim,left,right);xb<-mp(z$delta_MAE_B[i],lim,left,right)
   ln(xa,ys[i],xb,ys[i],pale,.8);pt(xa,ys[i]-.45,if(alg=='CatBoost')1 else 0,colors[alg],1.6,fill='white');pt(xb,ys[i]+.45,if(alg=='CatBoost')16 else 15,colors[alg],1.6)
  }
  axis_x(left,right,90,lim,c(-15,-10,-5,0,5),'ΔMAE [dt ha⁻¹]',99)
 }
 panel('B','Absolute MAE',2,111);panel('C','Macro increments and contrast',100,111)
 pt(43,119,1,neutral,1.6,fill='white');txt('M0',46,119,7.2);pt(60,119,17,ink,1.65);txt('M1',63,119,7.2)
 bl<-43;br<-91;lim2<-c(11,14.5);bby<-seq(129,178,length.out=8)
 for(i in 1:8){
  z<-abserr[i,];co<-colors[z$algorithm];weight<-if(z$weighting=='Pooled plots')'pooled' else 'macro';label<-paste(if(z$algorithm=='CatBoost')'CB' else 'Ridge',z$strategy,weight,sep=' · ')
  txt(label,39,bby[i],7.2,'right');x0<-mp(z$M0_MAE,lim2,bl,br);x1<-mp(z$M1_MAE,lim2,bl,br);ln(x0,bby[i],x1,bby[i],pale,.8);pt(x0,bby[i],1,neutral,1.65,fill='white');pt(x1,bby[i],17,co,1.65)
 }
 ln(bl,mean(bby[c(4,5)]),br,mean(bby[c(4,5)]),pale,.5)
 axis_x(bl,br,186,lim2,c(11,12,13,14),'MAE [dt ha⁻¹]',195)
 cl<-132;cr<-174;lim3<-c(-4.6,1.6);zero(cl,cr,lim3,123,182);cy<-c(129,139,149,159,169,179)
 for(i in 1:6){
  z<-metrics[i,];co<-colors[z$model_name];s<-if(z$strategy=='A_NATURAL_ONLY')'A' else if(z$strategy=='B_FULL_ELIGIBLE')'B' else 'B − A';label<-paste(if(z$model_name=='CatBoost')'CB' else 'Ridge',s,sep=' · ')
  txt(label,128,cy[i],7.3,'right');lo<-mp(z$macro_CI95_lo,lim3,cl,cr);hi<-mp(z$macro_CI95_hi,lim3,cl,cr);x<-mp(z$macro_delta_MAE,lim3,cl,cr)
  ln(lo,cy[i],hi,cy[i],co,.9);ln(lo,cy[i]-1,lo,cy[i]+1,co,.9);ln(hi,cy[i]-1,hi,cy[i]+1,co,.9)
  pt(x,cy[i],if(z$model_name=='CatBoost')16 else 15,co,1.7)
 }
 axis_x(cl,cr,186,lim3,c(-4,-2,0,1),'ΔMAE [dt ha⁻¹]',195)
 popViewport()
}
