# Redraw an existing manuscript figure from its complete frozen source tables.
args <- commandArgs(trailingOnly=FALSE)
script <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)]),fixed=TRUE)
root <- dirname(dirname(dirname(normalizePath(script))))
optional_library <- Sys.getenv('ENDOMETRIAL_R_LIBRARY',unset='')
if(nzchar(optional_library)).libPaths(c(optional_library,.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(jsonlite);library(grid)})
src <- file.path(root,'figure_source')
out <- file.path(root,'results/figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
options(device=function(...)grDevices::pdf(file.path(out,'diagnostic-default-device.pdf'),...))
scatter <- read.csv(file.path(src,'Source_clinical_a_all157_contrast_effects.csv'))
direct <- read.csv(file.path(src,'Source_clinical_b_all157_ENDO_UPP_intervals.csv'))
counts <- read.csv(file.path(src,'Source_clinical_c_fixed_family_sensitivity_counts.csv'))
selected <- direct[direct$BH_3U<.05,]
selected$symbol <- factor(selected$symbol,levels=selected$symbol[order(selected$effect)])
selected$evidence <- factor(selected$BY_3U<.05,levels=c(FALSE,TRUE),labels=c('BH471 only','Also BY471'))
stopifnot(nrow(scatter)==157L,nrow(direct)==157L,nrow(counts)==4L,nrow(selected)==7L)
theme_set(theme_classic(base_size=7,base_family='Helvetica')+
 theme(plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=6.7),
 axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 plot.caption=element_text(size=5.9,hjust=0),plot.margin=margin(12,8,7,8)))
rng <- range(c(scatter$effect,scatter$UPP_minus_NUPP),na.rm=TRUE);limits <- c(floor(rng[1]),ceiling(rng[2]))
a <- ggplot(scatter[scatter$estimable,],aes(effect,UPP_minus_NUPP,colour=direct_selected))+
 geom_abline(slope=1,intercept=0,linetype=2,linewidth=.35,colour='#777777')+
 geom_hline(yintercept=0,linewidth=.25,colour='#C1C5C8')+geom_vline(xintercept=0,linewidth=.25,colour='#C1C5C8')+
 geom_point(size=1.2,alpha=.85)+
 scale_colour_manual(values=c('FALSE'='#A4ABB2','TRUE'='#28788A'),
 labels=c('Other measured features','Direct ENDO-UPP: BH471 < 0.05'),name=NULL)+
 coord_fixed(xlim=limits,ylim=limits)+
 labs(title='Comparator-dependent effects',subtitle='144 samples; cycle-balanced contrasts',
 x='ENDO - NUPP (log2 signal)',y='UPP - NUPP (log2 signal)',
 caption='156 estimable features; 1 retained with P = 1.\nBoth axes share NUPP; no independent correlation test.')+
 guides(colour=guide_legend(nrow=2,byrow=TRUE))+
 theme(legend.position='bottom',legend.text=element_text(size=6),legend.key.height=grid::unit(2,'mm'))
b <- ggplot(selected,aes(effect,symbol))+
 geom_vline(xintercept=0,linetype=2,linewidth=.3,colour='#888888')+
 geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation='y',width=.18,linewidth=.4,colour='#577C91')+
 geom_point(aes(shape=evidence),size=2,stroke=.6,colour='#28788A')+
 scale_shape_manual(values=c(1,16),name=NULL)+
 labs(title='ENDO versus other pathology',subtitle='7 features selected by secondary BH471',
 x='ENDO - UPP (log2 signal)',y=NULL,
 caption='Pointwise 95% model intervals.\nTwo also pass BY471; absence of a hit is not equivalence.')+
 theme(legend.position='bottom',legend.text=element_text(size=6.5))
counts$label <- factor(paste0(counts$mode,' (n=',counts$N,')'),
                      levels=rev(paste0(counts$mode,' (n=',counts$N,')')))
long <- rbind(data.frame(label=counts$label,method='BH157',n=counts$BH),
              data.frame(label=counts$label,method='BY157',n=counts$BY))
c <- ggplot(counts,aes(y=label))+
 geom_vline(xintercept=157,linetype=3,linewidth=.3,colour='#AAAAAA')+
 geom_segment(aes(x=BY,xend=BH,yend=label),linewidth=.5,colour='#9EA6AD')+
 geom_point(data=long,aes(x=n,shape=method,colour=method),size=2,stroke=.7)+
 geom_text(data=long,aes(x=n,label=n),nudge_y=.23,size=2.2,family='Helvetica',colour='#333333')+
 scale_shape_manual(values=c('BH157'=16,'BY157'=1),name=NULL)+
 scale_colour_manual(values=c('BH157'='#28788A','BY157'='#636C75'),name=NULL)+
 scale_x_continuous(limits=c(0,165),breaks=c(0,40,80,120,157))+
 labs(title='Sensitivity of the primary ENDO-NUPP association',
 subtitle='Same fixed 157-gene family in each mode',x='Number with adjusted P < 0.05',y=NULL,
 caption='Models share participants. Global-QC sensitivity has identical primary samples and is not shown twice.')+
 theme(legend.position='right',legend.text=element_text(size=6.5))
top <- plot_grid(a,b,nrow=1,rel_widths=c(1.08,1),labels=c('a','b'),label_size=10,
                 label_fontfamily='Helvetica',label_fontface='bold')
figure <- plot_grid(top,c,ncol=1,rel_heights=c(1.55,1),labels=c('','c'),label_size=10,
                    label_fontfamily='Helvetica',label_fontface='bold')
name <- file.path(out,'Figure_clinical_comparator_robustness_v13');w <- 183;h <- 165
if(capabilities('cairo')) grDevices::cairo_pdf(filename=paste0(name,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica') else grDevices::pdf(file=paste0(name,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica');print(figure);dev.off()
svglite::svglite(paste0(name,'.svg'),width=w/25.4,height=h/25.4);print(figure);dev.off()
ragg::agg_png(paste0(name,'.png'),width=w,height=h,units='mm',res=600);print(figure);dev.off()
ragg::agg_tiff(paste0(name,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw');print(figure);dev.off()
ragg::agg_png(paste0(name,'_preview.png'),width=w,height=h,units='mm',res=300);print(figure);dev.off()
write_json(list(width_mm=w,height_mm=h,dpi=600,source_rows=c(a=157,b=157,c=4),backend='R',
 scope='Verified planned clinical gene comparisons; no diagnostic-model AUC or utility established; visualQA pending'),
 file.path(out,'clinical_figure_export_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_clinical_figure_session.txt'))
