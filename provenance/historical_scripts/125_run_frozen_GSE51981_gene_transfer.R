args <- commandArgs(trailingOnly=FALSE)
script <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)]),fixed=TRUE)
root <- dirname(dirname(normalizePath(script)))
.libPaths(c(file.path(root,'runtime/R-library'),.libPaths()))
suppressPackageStartupMessages({library(limma);library(jsonlite);library(digest)})
prep <- file.path(root,'outputs/biomarker_transfer_v13/clinical_preparation')
out <- file.path(root,'outputs/biomarker_transfer_v13/clinical_analysis');dir.create(out,showWarnings=FALSE)
options(device=function(...)grDevices::pdf(file.path(out,'diagnostic-default-device.pdf'),...))
freeze <- fromJSON(file.path(prep,'primary_clinical_analysis_freeze.json'))
for(i in seq_len(nrow(freeze$files)))stopifnot(digest(file=file.path(root,freeze$files$path[i]),algo='sha256')==freeze$files$sha256[i])
qa <- fromJSON(file.path(prep,'independent_source_preparation_verification.json'))
stopifnot(qa$result=='PASS')
U <- read.csv(file.path(prep,'U_all_identifiable_source_admissible_genes.csv'),colClasses='character')
stopifnot(nrow(U)==157L,!anyDuplicated(U$Entrez_ID))
meta <- read.csv(file.path(prep,'GSE51981_fixed_clinical_design.csv'),colClasses='character')
flag <- function(x)tolower(x)=='true'
ids <- read.csv(file.path(prep,'orgHs_all_Entrez_to_symbol.csv'),colClasses='character')
valid_ids <- unique(ids$ENTREZID[!is.na(ids$SYMBOL)&nzchar(ids$SYMBOL)])
target_results <- list();secondary <- NULL;fit_summaries <- list()
fit_mode <- function(mode,aggregation,all148=FALSE){
  expr <- as.matrix(read.csv(gzfile(file.path(prep,paste0('GSE51981_',aggregation,'_probe_gene_matrix.csv.gz'))),
                            row.names=1,check.names=FALSE))
  m <- if(all148)meta else meta[flag(meta$include_primary),]
  stopifnot(nrow(m)==if(all148)148L else 144L)
  expr <- expr[,m$sample,drop=FALSE]
  eligible <- rownames(expr)%in%valid_ids&apply(expr,1,function(x)all(is.finite(x))&&sd(x)>0)
  expr <- expr[eligible,,drop=FALSE]
  if(all148){
    m$clinical_group <- factor(m$clinical_group,levels=c('NUPP','UPP','ENDO'))
    m$phase_collapsed <- factor(m$phase_collapsed,levels=c('P','Secretory','Unknown'))
    design <- model.matrix(~clinical_group+phase_collapsed,m)
    C <- matrix(0,nrow=ncol(design),ncol=3,dimnames=list(colnames(design),c('ENDO_minus_NUPP','ENDO_minus_UPP','UPP_minus_NUPP')))
    C['clinical_groupENDO',1:2] <- 1;C['clinical_groupUPP',2] <- -1;C['clinical_groupUPP',3] <- 1
  }else{
    cells <- c('NUPP_P','NUPP_ESE','NUPP_MSE','UPP_P','UPP_ESE','UPP_MSE','ENDO_P','ENDO_ESE','ENDO_MSE')
    m$primary_cell <- factor(m$primary_cell,levels=cells)
    stopifnot(all(table(m$primary_cell)>=2))
    design <- model.matrix(~0+primary_cell,m);colnames(design) <- cells
    C <- matrix(0,nrow=9,ncol=3,dimnames=list(cells,c('ENDO_minus_NUPP','ENDO_minus_UPP','UPP_minus_NUPP')))
    C[grep('^ENDO_',cells),1:2] <- 1/3;C[grep('^NUPP_',cells),c(1,3)] <- -1/3
    C[grep('^UPP_',cells),2] <- -1/3;C[grep('^UPP_',cells),3] <- 1/3
  }
  stopifnot(qr(design)$rank==ncol(design),all(is.finite(expr)))
  fit <- lmFit(expr,design)
  ef <- eBayes(contrasts.fit(fit,C),trend=TRUE,robust=TRUE)
  se <- sqrt(ef$s2.post)*ef$stdev.unscaled[,1]
  ordinary_se <- fit$sigma*ef$stdev.unscaled[,1]
  z <- data.frame(Entrez_ID=rownames(expr),effect=ef$coefficients[,1],SE=se,df=ef$df.total,
                  ordinary_residual_variance=fit$sigma^2,ordinary_df=fit$df.residual,
                  ordinary_SE=ordinary_se,ordinary_P=2*pt(-abs(ef$coefficients[,1]/ordinary_se),fit$df.residual),
                  ordinary_lower95=ef$coefficients[,1]-qt(.975,fit$df.residual)*ordinary_se,
                  ordinary_upper95=ef$coefficients[,1]+qt(.975,fit$df.residual)*ordinary_se,
                  posterior_variance=ef$s2.post,unscaled_SE=ef$stdev.unscaled[,1],t=ef$t[,1],P=ef$p.value[,1],
                  lower95=ef$coefficients[,1]-qt(.975,ef$df.total)*se,
                  upper95=ef$coefficients[,1]+qt(.975,ef$df.total)*se)
  z$equivalence_P <- pmax(pt((z$effect+.5)/se,z$df,lower.tail=FALSE),pt((z$effect-.5)/se,z$df))
  write.csv(z,file.path(out,paste0(mode,'_all_background_gene_primary_contrast.csv')),row.names=FALSE)
  target <- merge(U,z,by='Entrez_ID',all.x=TRUE,sort=FALSE)
  target <- target[match(U$Entrez_ID,target$Entrez_ID),]
  target$estimable <- !is.na(target$P);target$P[!target$estimable] <- 1
  target$equivalence_P[!target$estimable] <- 1;target$ordinary_P[!is.finite(target$ordinary_P)] <- 1
  target$BH_U <- p.adjust(target$P,'BH');target$BY_U <- p.adjust(target$P,'BY')
  target$ordinary_BH_U <- p.adjust(target$ordinary_P,'BH');target$ordinary_BY_U <- p.adjust(target$ordinary_P,'BY')
  target$equivalence_BH_U <- p.adjust(target$equivalence_P,'BH');target$equivalence_BY_U <- p.adjust(target$equivalence_P,'BY')
  target$mode <- mode;target$N <- nrow(m)
  target$missing_reason <- ifelse(target$estimable,'','Absent unique-current-Entrez measurement or constant/nonfinite in fixed samples')
  rownames(design) <- m$sample
  write.csv(design,file.path(out,paste0(mode,'_design_matrix.csv')))
  write.csv(C,file.path(out,paste0(mode,'_contrast_matrix.csv')))
  saveRDS(list(fit=fit,contrasts=ef,design=design,C=C,metadata=m,aggregation=aggregation),file.path(out,paste0(mode,'_model.rds')))
  if(mode=='primary_mean_144'){
    J <- matrix(0,9,4,dimnames=list(rownames(C),c('ENDO_ESEvsP','ENDO_MSEvsP','UPP_ESEvsP','UPP_MSEvsP')))
    for(j in 1:4){
      g <- if(j<=2)'ENDO' else 'UPP';p <- if(j%%2==1)'ESE' else 'MSE'
      J[paste0(g,'_',p),j] <- 1;J[paste0(g,'_P'),j] <- -1
      J[paste0('NUPP_',p),j] <- -1;J['NUPP_P',j] <- 1
    }
    interaction_fit <- eBayes(contrasts.fit(fit,J),trend=TRUE,robust=TRUE)
    sec <- list()
    for(j in 2:3)sec[[length(sec)+1]] <- data.frame(Entrez_ID=rownames(expr),test=colnames(C)[j],
        effect=ef$coefficients[,j],statistic=ef$t[,j],P=ef$p.value[,j])
    sec[[3]] <- data.frame(Entrez_ID=rownames(expr),test='group_by_phase_F4',effect=NA_real_,
                          statistic=interaction_fit$F,P=interaction_fit$F.p.value)
    sec <- do.call(rbind,sec)
    scaffold <- merge(U,data.frame(test=c(colnames(C)[2:3],'group_by_phase_F4')),by=NULL)
    secondary <<- merge(scaffold,sec,by=c('Entrez_ID','test'),all.x=TRUE,sort=FALSE)
    secondary$P[is.na(secondary$P)] <<- 1
    secondary$BH_3U <<- p.adjust(secondary$P,'BH');secondary$BY_3U <<- p.adjust(secondary$P,'BY')
    write.csv(J,file.path(out,'primary_group_by_phase_interaction_matrix.csv'))
    saveRDS(interaction_fit,file.path(out,'primary_group_by_phase_interaction_fit.rds'))
  }
  fit_summaries[[mode]] <<- list(N=nrow(m),design_rank=ncol(design),background_genes=nrow(expr),U=157,
    estimable=sum(target$estimable),BH_U_hits=sum(target$BH_U<.05),BY_U_hits=sum(target$BY_U<.05),
    bounded_mean_BH_U_hits=sum(target$equivalence_BH_U<.05),bounded_mean_BY_U_hits=sum(target$equivalence_BY_U<.05))
  target
}
target_results[[1]] <- fit_mode('primary_mean_144','mean')
target_results[[2]] <- fit_mode('sensitivity_median_144','median')
target_results[[3]] <- fit_mode('sensitivity_all148_additive','mean',all148=TRUE)
stopifnot(!any(flag(meta$global_distribution_flag)))
target_results[[4]] <- target_results[[1]];target_results[[4]]$mode <- 'global_QC_sensitivity_identical_primary'
write.csv(do.call(rbind,target_results),file.path(out,'all_fixed_U157_primary_and_sensitivity_results.csv'),row.names=FALSE)
write.csv(secondary,file.path(out,'all_secondary_3U471_results.csv'),row.names=FALSE)
write_json(list(completed_at_utc=format(Sys.time(),tz='UTC'),modes=fit_summaries,
 scope='Gene-level mean associations and fixed mean boundary; not diagnostic-model performance; independent numerical verification pending.'),
 file.path(out,'primary_51981_summary.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_primary_clinical_session.txt'))
print(fit_summaries)
