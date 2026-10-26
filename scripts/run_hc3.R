# Recompute the fixed HC3 sensitivities from fitted models and processed inputs.
# This does not fit new clinical contrasts or replace the moderated primary tests.
args <- commandArgs(trailingOnly=FALSE)
script <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)]),fixed=TRUE)
root <- dirname(dirname(normalizePath(script)))
source(file.path(root,'scripts/portable_io.R'))
verify_inputs()
out <- file.path(root,'results/hc3');dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_csv <- function(path)read.csv(path,colClasses='character',check.names=FALSE)
U <- read_csv(file.path(root,'data/clinical_preparation/U_all_identifiable_source_admissible_genes.csv'))
U2 <- read_csv(file.path(root,'data/clinical_preparation/U2_both_source_admissible_genes.csv'))

hc3 <- function(model,expression,universe,mode,column,secondary=FALSE){
  X <- model$design; C <- model$C[,column]; n <- nrow(X); df <- n-ncol(X)
  stopifnot(qr(X)$rank==ncol(X),df>0)
  inv <- solve(crossprod(X)); leverage <- rowSums((X%*%inv)*X)
  weights <- as.vector(X%*%inv%*%C)
  stopifnot(all(leverage<1),all(model$metadata$sample%in%colnames(expression)))
  idx <- match(universe$Entrez_ID,rownames(model$fit$coefficients));ok <- !is.na(idx)
  row <- universe
  row$mode <- mode
  if(secondary)row$contrast <- colnames(model$C)[column]
  row$N <- n
  if(secondary)row$ordinary_df <- df else row$df_residual <- df
  row$estimable <- ok
  for(field in c('effect','HC3_SE','HC3_t','HC3_lower95','HC3_upper95'))row[[field]] <- NA_real_
  row$P <- 1
  row$missing_reason <- ifelse(ok,'',if(secondary)'Absent/constant/nonfinite under fixed input and sample rules' else 'Absent unique-current-Entrez measurement or constant/nonfinite in fixed samples')
  if(secondary){
    row$family_size <- 119
    row$role <- 'Fixed HC3 sensitivity; does not replace moderated mean, interaction or joint IUT'
  }else row$inference <- 'HC3 sandwich; Student-t residual df approximation; same frozen contrast and U157 family'
  Y <- t(expression[match(universe$Entrez_ID[ok],rownames(expression)),model$metadata$sample,drop=FALSE])
  beta <- t(model$fit$coefficients[idx[ok],,drop=FALSE])
  residual <- Y-X%*%beta
  corrected <- residual/(1-leverage)
  effect <- as.vector(t(C)%*%beta)
  se <- sqrt(colSums((weights*corrected)^2))
  statistic <- effect/se
  row$effect[ok] <- effect;row$HC3_SE[ok] <- se;row$HC3_t[ok] <- statistic
  row$HC3_lower95[ok] <- effect-qt(.975,df)*se;row$HC3_upper95[ok] <- effect+qt(.975,df)*se
  row$P[ok] <- 2*pt(-abs(statistic),df)
  if(secondary){row$BH_U2 <- p.adjust(row$P,'BH');row$BY_U2 <- p.adjust(row$P,'BY')}
  else{row$BH_U <- p.adjust(row$P,'BH');row$BY_U <- p.adjust(row$P,'BY')}
  row
}

primary <- lapply(c('primary_mean_144','sensitivity_median_144','sensitivity_all148_additive'),function(mode){
  model <- readRDS(file.path(root,'results/primary',paste0(mode,'_model.rds')))
  aggregation <- if(mode=='sensitivity_median_144')'median' else 'mean'
  expression <- as.matrix(read.csv(gzfile(file.path(root,'data/clinical_preparation',paste0('GSE51981_',aggregation,'_probe_gene_matrix.csv.gz'))),row.names=1,check.names=FALSE))
  hc3(model,expression,U,mode,1)
})
primary <- do.call(rbind,primary);stopifnot(nrow(primary)==471L)
write.csv(primary,file.path(out,'primary_HC3_471.csv'),row.names=FALSE,na='')
modes <- c('rawmean34','rawmedian34','rawmean_noSourceQC32','rawmean_FEversionAdjusted34','sourceProcessed_mean36','sourceProcessed_mean_noQC32')
secondary <- do.call(rbind,lapply(modes,function(mode){
  model <- readRDS(file.path(root,'results/secondary',paste0(mode,'_model.rds')))
  rbind(hc3(model,model$expression,U2,mode,1,TRUE),hc3(model,model$expression,U2,mode,2,TRUE))
}))
stopifnot(nrow(secondary)==1428L)
write.csv(secondary,file.path(out,'second_HC3_1428.csv'),row.names=FALSE,na='')
cat('HC3_RECOMPUTED_471_PRIMARY_AND_1428_SECONDARY_ROWS; run tests/check_hc3.py\n')
