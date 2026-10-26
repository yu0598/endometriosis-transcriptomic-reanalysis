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
nodes <- read.csv(file.path(src,'Source_workflow_nodes.csv'))
edges <- read.csv(file.path(src,'Source_workflow_edges.csv'))
stopifnot(nrow(nodes)==12L,nrow(edges)==8L)
stopifnot(all(nodes$branch[match(edges$from,nodes$id)]==nodes$branch[match(edges$to,nodes$id)]))
width_mm <- 183; height_mm <- 170
palette <- c('Published features'='#577C91','Measurement archives'='#6B7077','Clinical estimands'='#28788A')
draw_text <- function(label,x,y,size=7,family='Helvetica',face='plain',colour='#262A2E',just='centre')
  grid.text(label,x=unit(x,'mm'),y=unit(y,'mm'),just=just,
    gp=gpar(fontfamily=family,fontface=face,fontsize=size,col=colour,lineheight=1.15))
draw_workflow <- function(){
  grid.newpage();pushViewport(viewport(width=unit(width_mm,'mm'),height=unit(height_mm,'mm')))
  draw_text('Study workflow: complementary provenance and clinical-context audits',91.5,165,8,face='bold')
  for(i in seq_len(3)){
    x <- c(32,92,153)[i];branch <- names(palette)[i]
    draw_text(letters[i],x-27,154,9,face='bold',just='left')
    draw_text(branch,x,154,8,face='bold',colour=palette[[branch]])
    grid.lines(unit(c(x-27,x+27),'mm'),unit(c(150.5,150.5),'mm'),gp=gpar(col=palette[[branch]],lwd=.7))
  }
  for(i in seq_len(nrow(edges))){
    xy <- do.call(rbind,strsplit(strsplit(edges$path[i],';',fixed=TRUE)[[1]],',',fixed=TRUE))
    grid.lines(unit(as.numeric(xy[,1]),'mm'),unit(as.numeric(xy[,2]),'mm'),
      arrow=arrow(length=unit(1.3,'mm'),type='closed'),
      gp=gpar(col='#78828B',fill='#78828B',lwd=.7,lty=if(edges$dashed[i])2 else 1))
  }
  for(i in seq_len(nrow(nodes))){
    n <- nodes[i,];colour <- palette[[n$branch]]
    grid.roundrect(unit(n$x,'mm'),unit(n$y,'mm'),unit(n$width,'mm'),unit(n$height,'mm'),
      r=unit(1.1,'mm'),gp=gpar(fill=if(n$id=='a5')'#F7F8F9' else '#FFFFFF',col=colour,lwd=.65))
    title_y <- n$y+n$height/2-3.8
    draw_text(n$title,n$x,title_y,7.3,face='bold',colour=colour)
    body_y <- if(n$id=='a5')n$y-2 else n$y-2.4
    draw_text(n$body,n$x,body_y,7)
  }
  draw_text('Measurement identity, not patient counts\nor paper-specific leakage claims',92,97,7,colour='#555F68')
  grid.lines(unit(c(65,119),'mm'),unit(c(87,87),'mm'),gp=gpar(col='#D3D8DC',lwd=.6,lty=3))
  draw_text('Source-processed 36 profiles retained\nRaw defects do not establish processed error',92,32,7,colour='#555F68')
  draw_text('Equal phase / fertility weighting\nAll fixed-family entries retained\nNo original-classifier or AUC test',153,36,7,colour='#555F68')
  draw_text('Parallel audits of provenance and estimands; no causal arrow from reuse to clinical performance.',91.5,10,7,colour='#414951')
  popViewport()
}
stem <- file.path(out,'Figure_study_workflow_v13')
if(capabilities('cairo')) grDevices::cairo_pdf(filename=paste0(stem,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,family='Helvetica') else grDevices::pdf(file=paste0(stem,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,family='Helvetica')
draw_workflow();dev.off()
svglite::svglite(paste0(stem,'.svg'),width=width_mm/25.4,height=height_mm/25.4);draw_workflow();dev.off()
ragg::agg_png(paste0(stem,'.png'),width=width_mm,height=height_mm,units='mm',res=600);draw_workflow();dev.off()
ragg::agg_tiff(paste0(stem,'.tiff'),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw');draw_workflow();dev.off()
ragg::agg_png(paste0(stem,'_preview.png'),width=width_mm,height=height_mm,units='mm',res=300);draw_workflow();dev.off()
