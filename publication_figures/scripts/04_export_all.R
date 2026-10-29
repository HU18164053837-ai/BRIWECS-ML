#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(grid);library(svglite);library(ragg)})
script_arg<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))][1])
script_dir<-dirname(script_arg)
source(file.path(script_dir,'01_load_frozen_data.R'),encoding='UTF-8')
source(file.path(script_dir,'02_panel_data_and_checks.R'),encoding='UTF-8')
source(file.path(script_dir,'03_draw_figures.R'),encoding='UTF-8')
heights<-c(198,174,174,201);draws<-list(draw1,draw2,draw3,draw4)
# Minimal baseline-TIFF encoder keeps grayscale output truly 8 bit, one sample per pixel.
# Pixel generation remains entirely in R/ragg; rows are written without resampling.
gray_tiff<-function(draw,path,height){
 cap<-ragg::agg_capture(width=178,height=height,units='mm',res=800,background='white')
 grid.newpage();draw(TRUE);im<-cap(native=TRUE);dev.off()
 hh<-nrow(im);ww<-ncol(im);nentries<-12L;xres<-8L+2L+12L*nentries+4L;yres<-xres+8L;data_offset<-yres+8L
 con<-file(path,'wb');on.exit(close(con))
 w2<-function(v)writeBin(as.integer(v),con,size=2,endian='little');w4<-function(v)writeBin(as.integer(v),con,size=4,endian='little')
 writeBin(charToRaw('II'),con);w2(42);w4(8);w2(nentries)
 entry<-function(tag,type,value){w2(tag);w2(type);w4(1);if(type==3){w2(value);w2(0)}else w4(value)}
 entry(256,4,ww);entry(257,4,hh);entry(258,3,8);entry(259,3,1);entry(262,3,1)
 entry(273,4,data_offset);entry(277,3,1);entry(278,4,hh);entry(279,4,ww*hh)
 entry(282,5,xres);entry(283,5,yres);entry(296,3,2);w4(0);w4(800);w4(1);w4(800);w4(1)
 for(i in 1:hh)writeBin(as.raw(bitwAnd(as.integer(im[i,]),255L)),con)
 invisible(NULL)
}
export_one<-function(i){
 draw<-draws[[i]];height<-heights[i];stem<-paste0('Figure',i,'_v1.0')
 svglite::svglite(file.path(root,'figures',paste0(stem,'.svg')),width=178/25.4,height=height/25.4,bg='white',system_fonts=list(serif='Times New Roman'));grid.newpage();draw();dev.off()
 grDevices::cairo_pdf(file.path(root,'figures',paste0(stem,'.pdf')),width=178/25.4,height=height/25.4,family='Times New Roman',bg='white');grid.newpage();draw();dev.off()
 ragg::agg_png(file.path(root,'figures',paste0(stem,'.png')),width=178,height=height,units='mm',res=300,background='white');grid.newpage();draw();dev.off()
 # 1200 dpi exceeds the 600 dpi requested floor and preserves fine quantitative line work.
 ragg::agg_tiff(file.path(root,'figures',paste0(stem,'.tiff')),width=178,height=height,units='mm',res=1200,compression='lzw',background='white');grid.newpage();draw();dev.off()
 ragg::agg_png(file.path(root,'qa',paste0(stem,'_Grayscale.png')),width=178,height=height,units='mm',res=300,background='white');grid.newpage();draw(TRUE);dev.off()
 gray_tiff(draw,file.path(root,'qa',paste0(stem,'_Grayscale_800dpi.tiff')),height)
 svglite::svglite(file.path(root,'qa',paste0(stem,'_Grayscale.svg')),width=178/25.4,height=height/25.4,bg='white',system_fonts=list(serif='Times New Roman'));grid.newpage();draw(TRUE);dev.off()
 grDevices::cairo_pdf(file.path(root,'qa',paste0(stem,'_Grayscale.pdf')),width=178/25.4,height=height/25.4,family='Times New Roman',bg='white');grid.newpage();draw(TRUE);dev.off()
 cat('Exported ',stem,' SVG/PDF/PNG/TIFF plus grayscale checks\n',sep='')
}
if('--integration-editorial-figures' %in% commandArgs(TRUE)){
 export_one(3);export_one(4)
}else if('--figure4-editorial-only' %in% commandArgs(TRUE)){
 export_one(4)
}else if('--figure3-only' %in% commandArgs(TRUE)){
 export_one(3)
}else{
 # Figure 3 is the approved design baseline and is finalized before the other figures.
 for(i in c(3,1,2,4))export_one(i)
 grDevices::cairo_pdf(file.path(root,'qa','Figure_Visual_QA.pdf'),width=210/25.4,height=297/25.4,family='Times New Roman',bg='white',onefile=TRUE)
 for(gray in c(FALSE,TRUE))for(i in 1:4){
  grid.newpage();grid.text(paste0('Figure ',i,' | 178 mm / 100% | ',if(gray)'grayscale' else 'colour'),x=.075,y=.96,just='left',gp=gpar(fontfamily=family,fontsize=10))
  pushViewport(viewport(x=.5,y=.53,width=unit(178,'mm'),height=unit(heights[i],'mm')));draws[[i]](gray);popViewport()
 }
 for(i in 1:4){
  grid.newpage();grid.text(paste0('Figure ',i,' | 25% and 10% hierarchy checks'),x=.075,y=.96,just='left',gp=gpar(fontfamily=family,fontsize=10))
  for(sc in c(.25,.1)){
   pushViewport(viewport(x=if(sc==.25).3 else .75,y=.60,width=unit(178*sc,'mm'),height=unit(heights[i]*sc,'mm')))
   render_scale<-sc;draws[[i]]();render_scale<-1;popViewport()
  }
 }
 dev.off()
}
jsonlite::write_json(list(status='DATA_PASS',numerical_checks=nrow(qa),backend='R',font='Times New Roman',width_mm=178,height_mm=heights,rgb_tiff_dpi=1200,grayscale_tiff_dpi=800,new_models=0,new_inference=0,randomness='None'),file.path(root,'qa','R_Export_Status.json'),auto_unbox=TRUE,pretty=TRUE)
