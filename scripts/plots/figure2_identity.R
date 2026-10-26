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
groups <- read.csv(file.path(src,'Source_a_all42_measurement_groups.csv'))
pairs <- read.csv(file.path(src,'Source_b_all19_cross_series_pairs.csv'))
pairs <- pairs[order(pairs$sample_A),]
stopifnot(nrow(groups)==42L,nrow(pairs)==19L,sum(groups$records)==61L)
theme_set(theme_classic(base_size=7,base_family='Helvetica')+
 theme(plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=7),
 axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 plot.caption=element_text(size=6,hjust=0),plot.margin=margin(5,5,5,5)))
counts <- data.frame(multiplicity=factor(c('One GSM record','Two GSM records'),
 levels=c('Two GSM records','One GSM record')),groups=c(23,19),records=c(23,38))
a <- ggplot(counts,aes(groups,multiplicity))+
 geom_col(width=.48,fill='#577C91')+
 geom_text(aes(x=groups/2,label=paste0(groups,' groups\n',records,' GSM records')),
 colour='white',lineheight=1.1,size=2.25,family='Helvetica')+
 scale_x_continuous(limits=c(0,30),breaks=c(0,10,20,30),expand=expansion(mult=c(0,.02)))+
 labs(title='42 raw measurements in 61 records',subtitle='Complete raw-intensity vectors; fixed audit set',
 x='Number of distinct measurement groups',y=NULL,
 caption='20 records from GSE7305\n41 title-selected records from GSE7307\n\n19/20 GSE7305 measurements recur\nin the selected GSE7307 records.\n\nMeasurement groups are not patient counts.')+
 theme(axis.text.y=element_text(size=7),plot.margin=margin(15,7,20,7))
pairs$position <- rev(seq_len(nrow(pairs)))
pairs$tissue <- factor(pairs$group_A,levels=c('lesion','endometrial_comparator'),
 labels=c('Nominal lesion','Nominal endometrial comparator'))
b <- ggplot(pairs,aes(y=position,colour=tissue))+
 geom_segment(aes(x=.22,xend=.78,yend=position),linewidth=.5)+
 geom_point(aes(x=.22),size=1.5)+geom_point(aes(x=.78),size=1.5)+
 geom_text(aes(x=.18,label=sample_A),hjust=1,size=2.25,family='Helvetica',show.legend=FALSE)+
 geom_text(aes(x=.82,label=sample_B),hjust=0,size=2.25,family='Helvetica',show.legend=FALSE)+
 annotate('text',x=c(.10,.90),y=20.4,label=c('GSE7305','GSE7307'),
 size=2.6,fontface='bold',family='Helvetica')+
 scale_colour_manual(values=c('#28788A','#6B7077'),name=NULL)+
 coord_cartesian(xlim=c(-.2,1.2),ylim=c(.5,21),clip='off')+
 labs(title='19 identical CEL payload pairs',subtitle='All decompressed bytes equal; two intensity decoders agree',
 caption='10 lesion pairs and 9 endometrial comparator pairs.\nUnmatched in this fixed set: 1 GSE7305 record and 22 GSE7307 records.')+
 theme_void(base_size=7,base_family='Helvetica')+
 theme(plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=7),
 plot.caption=element_text(size=6,hjust=0),legend.position='bottom',legend.text=element_text(size=6.5),
 legend.key.width=grid::unit(4,'mm'),plot.margin=margin(15,7,5,7))
figure <- plot_grid(a,b,nrow=1,rel_widths=c(1,1.45),labels=c('a','b'),
 label_fontfamily='Helvetica',label_fontface='bold',label_size=10)
name <- file.path(out,'Figure_raw_measurement_identity_v13');w <- 183;h <- 122
if(capabilities('cairo')) grDevices::cairo_pdf(filename=paste0(name,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica') else grDevices::pdf(file=paste0(name,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica')
print(figure);dev.off()
svglite::svglite(paste0(name,'.svg'),width=w/25.4,height=h/25.4);print(figure);dev.off()
ragg::agg_png(paste0(name,'.png'),width=w,height=h,units='mm',res=600);print(figure);dev.off()
ragg::agg_tiff(paste0(name,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw');print(figure);dev.off()
ragg::agg_png(paste0(name,'_preview.png'),width=w,height=h,units='mm',res=300);print(figure);dev.off()
write_json(list(width_mm=w,height_mm=h,raster_dpi=600,plotting_backend='R',
 source_records=61,unique_raw_measurements=42,exact_cross_series_pairs=19,
 status='Raw measurement identity result; visual QA pending; manuscript not submission ready'),
 file.path(out,'export_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_figure_session.txt'))
