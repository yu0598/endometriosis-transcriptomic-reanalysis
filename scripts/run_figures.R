args <- commandArgs(trailingOnly=FALSE)
script <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)]),fixed=TRUE)
root <- dirname(dirname(normalizePath(script)))
for(name in c('figure1_workflow.R','figure2_identity.R','figure3_comparators.R','figure4_joint.R')){
  status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(file.path(root,'scripts/plots',name))))
  if(status!=0)stop('Figure redraw failed: ',name)
}
cat('FOUR_FIGURES_REDRAWN_FROM_FROZEN_SOURCE_TABLES\n')
