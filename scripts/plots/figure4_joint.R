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
joint <- read.csv(file.path(src,'Source_second_a_all119_joint_associations.csv'))
joint$evidence <- factor(joint$evidence,levels=c('No joint rejection','BH119 only','Also BY119'))
strata <- read.csv(file.path(src,'Source_second_b_all119_model_derived_stratum_intervals.csv'))
counts <- read.csv(file.path(src,'Source_second_c_all6_fixed_mode_joint_counts.csv'))
genes <- sort(joint$symbol[joint$joint_BH_U2<.05])
display <- strata[strata$selected_primary_joint_BH,]
display$symbol <- factor(display$symbol,levels=rev(genes))
display$stratum <- factor(display$stratum,levels=c('Mean','Fertile','Infertile'))
stopifnot(nrow(joint)==119L,nrow(strata)==357L,nrow(counts)==6L,length(genes)==9L)
theme_set(theme_classic(base_size=7,base_family='Helvetica')+
 theme(plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=6.4),
 axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 plot.caption=element_text(size=5.8,hjust=0),plot.margin=margin(12,7,7,7)))
a <- ggplot(joint[joint$joint_estimable,],aes(primary144_effect,second_effect,colour=evidence,shape=evidence))+
 geom_hline(yintercept=0,linewidth=.3,colour='#B4BAC0')+geom_vline(xintercept=0,linewidth=.3,colour='#B4BAC0')+
 geom_point(size=1.45,alpha=.85,stroke=.6)+
 scale_colour_manual(values=c('#A4ABB2','#28788A','#9B633A'),name=NULL)+
 scale_shape_manual(values=c(16,1,16),name=NULL)+
 labs(title='Joint average associations',subtitle='119 fixed tests; 116 jointly estimable',
 x='GSE51981: ENDO - NUPP\n(processed log2 signal)',
 y='Raw GSE120103: ENDO - NONENDO\n(log2 signal)',
 caption='9 BH119; 3 BY119. No pooled effect or AUC.\nDifferent assays and clinical comparator mixtures.')+
 theme(legend.position='bottom',legend.text=element_text(size=5.8),legend.key.width=grid::unit(3,'mm'))+
 guides(colour=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))
b <- ggplot(display,aes(effect,symbol))+
 geom_vline(xintercept=0,linetype=2,linewidth=.3,colour='#9AA2A9')+
 geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation='y',width=.17,linewidth=.35,colour='#577C91')+
 geom_point(size=1.35,colour='#28788A')+
 facet_wrap(~stratum,nrow=1)+
 scale_x_continuous(breaks=function(x)pretty(x,n=3))+
 labs(title='Fertility-stratum context',subtitle='9 features selected by the primary joint BH119',
 x='Raw-cohort effect (log2 signal)',y=NULL,
 caption='Model-derived, pointwise 95% intervals; descriptive strata.\nThe mean is fertility-balanced, not population prevalence.')+
 theme(strip.background=element_blank(),strip.text=element_text(size=6.4),
       axis.text.y=element_text(size=6.2),panel.spacing.x=grid::unit(2,'mm'))
counts$label <- factor(c('Raw mean, primary (34)','Raw median (34)','Raw, source-QC subset (32)',
                         'Raw, FE-version adjusted (34)','Source processed (36)','Source, QC subset (32)'),
                       levels=rev(c('Raw mean, primary (34)','Raw median (34)','Raw, source-QC subset (32)',
                         'Raw, FE-version adjusted (34)','Source processed (36)','Source, QC subset (32)')))
long <- rbind(data.frame(label=counts$label,method='BH119',n=counts$BH),
              data.frame(label=counts$label,method='BY119',n=counts$BY))
c <- ggplot(counts,aes(y=label))+
 geom_segment(aes(x=BY,xend=BH,yend=label),linewidth=.5,colour='#A4ADB5')+
 geom_point(data=long,aes(x=n,colour=method,shape=method),size=2,stroke=.65)+
 geom_text(data=long,aes(x=n,label=n),nudge_y=.25,size=2.15,family='Helvetica')+
 scale_colour_manual(values=c('BH119'='#28788A','BY119'='#636C75'),name=NULL)+
 scale_shape_manual(values=c('BH119'=16,'BY119'=1),name=NULL)+
 scale_x_continuous(limits=c(0,15),breaks=c(0,3,6,9,12,15))+
 labs(title='All fixed processing and sample sensitivities',
 subtitle='Same-direction average-association test; family119 in every mode',
 x='Number meeting the joint adjustment',y=NULL,
 caption='Modes share profiles. Counts change with processing, sample availability and estimand; none replaces the raw primary.')+
 theme(legend.position='right',legend.text=element_text(size=6.5))
top <- plot_grid(a,b,nrow=1,rel_widths=c(1,1.2),labels=c('a','b'),label_size=10,
                 label_fontfamily='Helvetica',label_fontface='bold')
figure <- plot_grid(top,c,ncol=1,rel_heights=c(1.45,1.05),labels=c('','c'),label_size=10,
                    label_fontfamily='Helvetica',label_fontface='bold')
name <- file.path(out,'Figure_second_cohort_joint_context_v13');w <- 183;h <- 178
if(capabilities('cairo')) grDevices::cairo_pdf(filename=paste0(name,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica') else grDevices::pdf(file=paste0(name,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica');print(figure);dev.off()
svglite::svglite(paste0(name,'.svg'),width=w/25.4,height=h/25.4);print(figure);dev.off()
ragg::agg_png(paste0(name,'.png'),width=w,height=h,units='mm',res=600);print(figure);dev.off()
ragg::agg_tiff(paste0(name,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw');print(figure);dev.off()
ragg::agg_png(paste0(name,'_preview.png'),width=w,height=h,units='mm',res=300);print(figure);dev.off()
write_json(list(width_mm=w,height_mm=h,dpi=600,backend='R',source_rows=c(a=119,b=357,c=6),
 scope='Same-direction mean-association IUT, descriptive strata and fixed sensitivities; not original model validation; visualQA pending'),
 file.path(out,'second_cohort_figure_export_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_second_cohort_figure_session.txt'))
