# Invoke the unchanged statistical bodies sequentially.
args <- commandArgs(trailingOnly=FALSE)
script <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)]),fixed=TRUE)
root <- dirname(dirname(normalizePath(script)))
if(normalizePath(getwd())!=root)stop('Run this command from the package root, with TMPDIR=runtime/tmp.')
for(name in c('run_primary.R','run_secondary.R')){
 extra <- if(name=='run_secondary.R')'--run' else character()
 status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(file.path('scripts',name)),extra))
 if(status!=0)stop('Rerun failed: ',name,' status ',status)
}
cat('ALL_FIXED_STATISTICAL_MODES_COMPLETE; run python3 scripts/compare_results.py next.\n')
