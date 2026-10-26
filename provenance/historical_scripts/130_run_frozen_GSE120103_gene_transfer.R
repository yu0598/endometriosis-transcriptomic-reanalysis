# Fixed second-cohort mean and interaction tests; no score, AUC or gene selection.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args %in% c('--preflight', '--freeze-only', '--run'))
root <- normalizePath(getwd())
stopifnot(basename(root) == 'IUA_TE_内膜修复研究_2026', Sys.getenv('TMPDIR') == 'runtime/tmp')
.libPaths(c(file.path(root, 'runtime/R-library'), .libPaths()))
suppressPackageStartupMessages({library(limma); library(jsonlite); library(digest)})
prep <- 'outputs/biomarker_transfer_v13/clinical_preparation'
rawprep <- 'outputs/biomarker_transfer_v13/GSE120103_raw_preparation'
primary <- 'outputs/biomarker_transfer_v13/clinical_analysis'
out <- file.path(primary, 'GSE120103_secondary')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(warn = 1, digits = 17, device = function(...) pdf(file.path(out, 'unused_default_device.pdf'), ...))
freeze_path <- file.path(out, 'secondary_statistical_freeze.json')
script_path <- 'scripts/130_run_frozen_GSE120103_gene_transfer.R'
sha <- function(path) digest(file = path, algo = 'sha256')
read_csv <- function(path) read.csv(path, colClasses = 'character', check.names = FALSE,
                                  fileEncoding = 'UTF-8-BOM', na.strings = c('', 'NA'))
flag <- function(x) !is.na(x) & tolower(x) == 'true'
write_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = '')
write_json_out <- function(x, name) write_json(x, file.path(out, name), pretty = TRUE,
                                              auto_unbox = TRUE, digits = NA, na = 'null')
cells <- c('NONENDO_Fertile', 'NONENDO_Infertile', 'ENDO_Fertile', 'ENDO_Infertile')
modes <- data.frame(mode = c('rawmean34', 'rawmedian34', 'rawmean_noSourceQC32',
                            'rawmean_FEversionAdjusted34', 'sourceProcessed_mean36', 'sourceProcessed_mean_noQC32'),
  input = c('raw', 'raw', 'raw', 'raw', 'source', 'source'),
  aggregation = c('mean', 'median', 'mean', 'mean', 'mean', 'mean'),
  drop_source_QC = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
  FE_adjust = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
  N = c(34L, 34L, 32L, 34L, 36L, 32L), stringsAsFactors = FALSE)
parameters <- list(U = 157L, U2 = 119L, cells = cells,
  mean_contrast = c(-.5, -.5, .5, .5), interaction_contrast = c(1, -1, -1, 1),
  FE_contrast_coefficient = 0, FE_indicator = 'FE9.5.3.1=1 vs FE10.7.3.1=0; fixed additive',
  background = 'Frozen current-valid Entrez with nonblank symbol; finite and sd>0 in each fixed sample subset',
  moderation = 'limma eBayes(trend=TRUE,robust=TRUE)',
  family = 'Complete U2=119 for each mean family and each separate SECONDARY interaction family, per mode',
  missing_P = 1, intervals = 'Two-sided 95% moderated t interval; ordinary OLS diagnostics separate',
  primary_component = 'Exact primary_mean_144_model.rds first contrast t/df, crosschecked against primary U157 CSV estimability',
  joint = 'jointUp=max(up1,up2); jointDown=max(down1,down2); jointP=min(1,2*min(jointUp,jointDown)); BH/BY119 per mode',
  one_sided_missing = 1, joint_primary_mode = 'rawmean34',
  other_modes = 'Five fixed sensitivities, never selectable replacements',
  no_primary_significance_filter = TRUE, no_AUC_or_scores = TRUE, no_new_GO = TRUE,
  interaction_cannot_rescue_mean_or_joint_primary = TRUE,
  raw_missing_targets = c('401546', '79605'), primary_constant_target = '4829',
  independent_HC3_sensitivity = list(
    contrasts = c('mean_ENDO_minus_NONENDO', 'SECONDARY_disease_by_fertility'),
    modes = 'Each of all six fixed modes, using its exact fixed design and samples',
    variance = 'HC3 heteroskedasticity-consistent covariance',
    reference_distribution = 'Student-t with ordinary residual degrees of freedom',
    multiplicity = 'Complete U2=119 BH/BY separately for each mode and contrast family; missing P=1',
    implementation = 'Computed later during independent verification, not by this moderated-fit script',
    interpretation = 'Fixed sensitivity; never replaces moderated primary or joint IUT'))
versions <- list(R = R.version.string, limma = as.character(packageVersion('limma')),
                 jsonlite = as.character(packageVersion('jsonlite')), digest = as.character(packageVersion('digest')))
U <- read_csv(file.path(prep, 'U_all_identifiable_source_admissible_genes.csv'))
U2 <- read_csv(file.path(prep, 'U2_both_source_admissible_genes.csv'))
stopifnot(nrow(U) == 157L, nrow(U2) == 119L, !anyDuplicated(U$Entrez_ID),
          !anyDuplicated(U2$Entrez_ID), all(U2$Entrez_ID %in% U$Entrez_ID))
raw_meta <- read_csv(file.path(rawprep, 'raw34_fixed_clinical_and_technical_metadata.csv'))
source_meta <- read_csv(file.path(prep, 'GSE120103_fixed_clinical_design.csv'))
stopifnot(nrow(raw_meta) == 34L, nrow(source_meta) == 36L)
raw_qa <- read_json(file.path(rawprep, 'independent_raw_preprocessing_verification.json'), simplifyVector = TRUE)
primary_qa <- read_json(file.path(primary, 'independent127_GSE51981_statistical_verification.json'), simplifyVector = TRUE)
stopifnot(raw_qa$result == 'PASS_WITH_EXTREME_TAIL_NUMERICAL_WARNING',
          primary_qa$result == 'PASS_WITH_BACKGROUND_OLS_NUMERICAL_WARNINGS')
raw_failed <- raw_qa$checks$check[!raw_qa$checks$passed]
stopifnot(identical(raw_failed, 'every_normexp_conditional_signal_plus_offset50'))

build_design <- function(config) {
  m <- if (config$input == 'raw') raw_meta else source_meta
  if (config$drop_source_QC) m <- m[!flag(m$global_distribution_flag), , drop = FALSE]
  stopifnot(nrow(m) == config$N, !anyDuplicated(m$sample), all(m$primary_cell %in% cells))
  m$primary_cell <- factor(m$primary_cell, levels = cells)
  counts <- as.integer(table(m$primary_cell))
  expected <- if (config$drop_source_QC) c(9L,9L,9L,5L) else if (config$input == 'source') rep(9L,4L) else c(9L,9L,9L,7L)
  stopifnot(identical(counts, expected), min(counts) >= 2L)
  design <- model.matrix(~0 + primary_cell, m); colnames(design) <- cells
  if (config$FE_adjust) {
    stopifnot(all(m$FE_version %in% c('9.5.3.1', '10.7.3.1')))
    design <- cbind(design, FE_9.5_vs_10.7 = as.integer(m$FE_version == '9.5.3.1'))
  }
  rownames(design) <- m$sample
  C <- cbind(mean_ENDO_minus_NONENDO = c(-.5,-.5,.5,.5),
             SECONDARY_disease_by_fertility = c(1,-1,-1,1))
  if (config$FE_adjust) C <- rbind(C, c(0,0))
  rownames(C) <- colnames(design)
  list(metadata = m, design = design, C = C, counts = counts,
       estimable = qr(design)$rank == ncol(design), rank = qr(design)$rank)
}
preflight <- lapply(seq_len(nrow(modes)), function(i) {
  d <- build_design(modes[i, ])
  data.frame(mode = modes$mode[i], N = nrow(d$metadata), design_columns = ncol(d$design),
             rank = d$rank, full_rank = d$estimable, NF = d$counts[1], NI = d$counts[2],
             EF = d$counts[3], EI = d$counts[4])
})
preflight <- do.call(rbind, preflight)
stopifnot(all(preflight$full_rank), identical(preflight$rank, c(4L,4L,4L,5L,4L,4L)))
if (args == '--preflight') {
  stopifnot(!file.exists(freeze_path))
  write_csv(preflight, 'secondary_design_preflight.csv')
  write_json_out(list(status = 'draft_parsed_metadata_design_only_no_group_fit', parameters = parameters,
    modes = modes, preflight = preflight, raw_QA = raw_qa$result, raw_QA_strict_warning = raw_failed,
    primary_QA = primary_qa$result, group_fits_started = FALSE), 'secondary_preflight.json')
  cat('PREFLIGHT_PASS: six fixed designs, counts 9/9/9/7 or 9/9/9/9 or 9/9/9/5; FE rank5; NO GROUP FIT.\n')
  quit(save = 'no', status = 0)
}
if (args == '--freeze-only') {
  stopifnot(!file.exists(freeze_path))
  rawfreeze <- read_json(file.path(rawprep, 'preprocessing_freeze.json'), simplifyVector = TRUE)
  capture.output(sessionInfo(), file = file.path(out, 'secondary_R_session_before_freeze.txt'))
  inputs <- unique(c('AGENTS.md',
    '33_临床基因迁移检验规则_v13_2.md', '34_临床检验执行补充_原始Agilent处理与新意边界_v13_3.md',
    '35_第二队列原始处理实施约定_v13_4.md', script_path,
    'scripts/125_run_frozen_GSE51981_gene_transfer.R', 'scripts/127_verify_GSE51981_models_independently.py',
    'scripts/128_reprocess_GSE120103_raw_Agilent.R', 'scripts/131_verify_raw_Agilent_preprocessing.py',
    rawfreeze$inputs$path, list.files(rawprep, recursive = TRUE, full.names = TRUE),
    list.files(prep, recursive = TRUE, full.names = TRUE),
    list.files(primary, recursive = FALSE, full.names = TRUE)[!dir.exists(list.files(primary, recursive = FALSE, full.names = TRUE))],
    file.path(out, 'secondary_R_session_before_freeze.txt'),
    file.path(out, c('secondary_preflight.json', 'secondary_design_preflight.csv')),
    'runtime/R-library/limma/DESCRIPTION', 'runtime/R-library/limma/R/limma.rdb',
    'runtime/R-library/limma/R/limma.rdx'))
  stopifnot(all(file.exists(inputs)), !any(dir.exists(inputs)))
  manifest <- data.frame(path = sort(inputs), stringsAsFactors = FALSE)
  manifest$bytes <- unname(file.info(manifest$path)$size)
  manifest$sha256 <- vapply(manifest$path, sha, character(1))
  write_json(list(frozen_at_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE), versions = versions,
    parameters = parameters, modes = modes, preflight = preflight,
    target_U_Entrez_ID = U$Entrez_ID, target_U2_Entrez_ID = U2$Entrez_ID,
    files = manifest, raw_preprocessing_QA = raw_qa$result, raw_preprocessing_strict_warning = raw_failed,
    primary_QA = primary_qa$result, secondary_group_fits_at_freeze = FALSE,
    primary_already_computed = TRUE, existing_primary_and_preparation_immutable = TRUE,
    scope = 'Second cohort association and same-direction IUT only, all fixed families; no diagnostic performance'),
    freeze_path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat('SECONDARY_FREEZE_WRITTEN', nrow(manifest), 'files; no secondary group fit.\n')
  quit(save = 'no', status = 0)
}
stopifnot(file.exists(freeze_path), !file.exists(file.path(out, 'secondary_run_started.json')))
freeze <- read_json(freeze_path, simplifyVector = TRUE)
stopifnot(identical(freeze$versions, versions), isTRUE(all.equal(freeze$parameters, parameters, check.attributes = FALSE)),
          isTRUE(all.equal(freeze$modes, modes, check.attributes = FALSE)),
          identical(freeze$target_U2_Entrez_ID, U2$Entrez_ID), identical(freeze$target_U_Entrez_ID, U$Entrez_ID))
verify_inputs <- function() {
  stopifnot(all(file.exists(freeze$files$path)), all(file.info(freeze$files$path)$size == freeze$files$bytes),
            identical(unname(vapply(freeze$files$path, sha, character(1))), freeze$files$sha256))
}
verify_inputs()
write_json_out(list(started_at_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
                    freeze_sha256 = sha(freeze_path)), 'secondary_run_started.json')
warnings_seen <- character(); current_mode <- 'primary_component_read_only'
main <- function() {
  # Exact first-cohort components; no significance or effect-size selection.
  primary_rds <- readRDS(file.path(primary, 'primary_mean_144_model.rds'))
  primary_csv <- read_csv(file.path(primary, 'all_fixed_U157_primary_and_sensitivity_results.csv'))
  primary_csv <- primary_csv[primary_csv$mode == 'primary_mean_144', ]
  stopifnot(nrow(primary_csv) == 157L, !anyDuplicated(primary_csv$Entrez_ID),
            setequal(primary_csv$Entrez_ID, U$Entrez_ID), nrow(primary_rds$metadata) == 144L,
            colnames(primary_rds$contrasts$coefficients)[1] == 'ENDO_minus_NUPP')
  ef1 <- primary_rds$contrasts
  primary_components <- U2
  j <- match(U2$Entrez_ID, rownames(ef1$coefficients)); k <- match(U2$Entrez_ID, primary_csv$Entrez_ID)
  primary_components$estimable <- !is.na(j)
  stopifnot(identical(primary_components$estimable, flag(primary_csv$estimable[k])))
  primary_components$effect <- ef1$coefficients[j, 1]
  primary_components$SE <- (sqrt(ef1$s2.post) * ef1$stdev.unscaled[, 1])[j]
  primary_components$t <- ef1$t[j, 1]; primary_components$df <- ef1$df.total[j]
  for (field in c('effect', 'SE', 't', 'df')) {
    a <- primary_components[[field]]; b <- as.numeric(primary_csv[[field]][k])
    ok <- primary_components$estimable
    stopifnot(isTRUE(all.equal(a[ok], b[ok], tolerance = 1e-12, check.attributes = FALSE)))
  }
  primary_components$P_up <- primary_components$P_down <- rep(1,119L)
  ok1 <- primary_components$estimable
  stopifnot(all(is.finite(primary_components$t[ok1])), all(primary_components$df[ok1] > 0))
  primary_components$P_up[ok1] <- pt(primary_components$t[ok1], primary_components$df[ok1], lower.tail = FALSE)
  primary_components$P_down[ok1] <- pt(primary_components$t[ok1], primary_components$df[ok1])
  primary_components$missing_reason <- ifelse(ok1, '', primary_csv$missing_reason[k])
  stopifnot(!primary_components$estimable[match('4829', primary_components$Entrez_ID)])
  write_csv(primary_components, 'all119_primary144_exact_RDS_components_and_CSV_estimability.csv')
  ids <- read_csv(file.path(prep, 'orgHs_all_Entrez_to_symbol.csv'))
  valid_ids <- unique(ids$ENTREZID[!is.na(ids$SYMBOL) & nzchar(ids$SYMBOL)])
  mean_targets <- interaction_targets <- joints <- backgrounds <- mode_summaries <- list()
  for (i in seq_len(nrow(modes))) {
    config <- modes[i, ]; current_mode <<- config$mode
    d <- build_design(config); m <- d$metadata; design <- d$design; C <- d$C
    if (config$input == 'raw') {
      expr <- readRDS(file.path(rawprep, paste0('Entrez_', config$aggregation, '_different_ProbeName_log2_matrix.rds')))
    } else {
      expr <- as.matrix(read.csv(gzfile(file.path(prep, 'GSE120103_mean_probe_gene_matrix.csv.gz')),
                                 row.names = 1, check.names = FALSE))
    }
    stopifnot(!anyDuplicated(rownames(expr)), !anyDuplicated(colnames(expr)), all(m$sample %in% colnames(expr)))
    expr <- expr[, m$sample, drop = FALSE]
    current <- rownames(expr) %in% valid_ids
    finite <- rowSums(!is.finite(expr)) == 0L
    nonconstant <- apply(expr, 1L, function(x) all(is.finite(x)) && sd(x) > 0)
    eligible <- current & finite & nonconstant
    admission <- data.frame(Entrez_ID = rownames(expr), current_valid_identifier = current,
                  all_values_finite = finite, nonconstant_in_fixed_samples = nonconstant, admitted = eligible,
                  reason = ifelse(!current, 'Not in frozen current-valid Entrez-with-symbol background',
                           ifelse(!finite, 'Nonfinite in fixed samples', ifelse(!nonconstant, 'Constant in fixed samples', 'Admitted'))))
    write_csv(admission, paste0(config$mode, '_all_input_gene_background_admission.csv'))
    expr <- expr[eligible, , drop = FALSE]
    stopifnot(nrow(expr) > 1L, all(is.finite(expr)))
    write_csv(data.frame(sample = m$sample, m[, setdiff(names(m), 'sample'), drop = FALSE]),
              paste0(config$mode, '_fixed_sample_metadata.csv'))
    write_csv(data.frame(sample = rownames(design), design, check.names = FALSE), paste0(config$mode, '_design_matrix.csv'))
    write_csv(data.frame(coefficient = rownames(C), C, check.names = FALSE), paste0(config$mode, '_contrast_matrix.csv'))
    if (!d$estimable) stop('Fixed design is not full rank; no ad hoc replacement allowed: ', config$mode)
    cat('GROUP_FIT_START', config$mode, ncol(expr), 'samples;', nrow(expr), 'background genes; rank', d$rank, '\n')
    fit <- lmFit(expr, design)
    ef <- eBayes(contrasts.fit(fit, C), trend = TRUE, robust = TRUE)
    saveRDS(list(fit = fit, contrasts = ef, design = design, C = C, metadata = m,
                 expression = expr, config = config, preprocessing_QA = raw_qa$result),
            file.path(out, paste0(config$mode, '_model.rds')), compress = 'gzip')
    write_csv(data.frame(Entrez_ID = rownames(expr), fit$coefficients, check.names = FALSE),
              paste0(config$mode, '_all_background_fitted_coefficients.csv'))
    targets <- list()
    for (column in 1:2) {
      se <- sqrt(ef$s2.post) * ef$stdev.unscaled[, column]
      ordinary_se <- fit$sigma * ef$stdev.unscaled[, column]
      effect <- ef$coefficients[, column]
      ordinary_t <- effect / ordinary_se
      z <- data.frame(Entrez_ID = rownames(expr), effect = effect, SE = se, df = ef$df.total,
        t = ef$t[, column], P = ef$p.value[, column], lower95 = effect - qt(.975, ef$df.total) * se,
        upper95 = effect + qt(.975, ef$df.total) * se,
        ordinary_residual_variance = fit$sigma^2, ordinary_df = fit$df.residual,
        ordinary_SE = ordinary_se, ordinary_t = ordinary_t,
        ordinary_P = 2 * pt(-abs(ordinary_t), fit$df.residual),
        ordinary_lower95 = effect - qt(.975, fit$df.residual) * ordinary_se,
        ordinary_upper95 = effect + qt(.975, fit$df.residual) * ordinary_se,
        posterior_variance = ef$s2.post, unscaled_SE = ef$stdev.unscaled[, column],
        average_expression = fit$Amean, mode = config$mode, contrast = colnames(C)[column], N = nrow(m))
      stopifnot(all(is.finite(z$effect)), all(is.finite(z$SE)), all(z$SE > 0),
                all(is.finite(z$t)), all(is.finite(z$P)), all(z$P >= 0 & z$P <= 1))
      write_csv(z, paste0(config$mode, if (column == 1L) '_all_background_mean_contrast.csv' else '_all_background_SECONDARY_interaction.csv'))
      target <- cbind(U2, z[match(U2$Entrez_ID, z$Entrez_ID), setdiff(names(z), 'Entrez_ID'), drop = FALSE])
      rownames(target) <- NULL
      target$estimable <- !is.na(target$P)
      missing_idx <- match(target$Entrez_ID, admission$Entrez_ID)
      target$missing_reason <- ifelse(target$estimable, '',
                             ifelse(is.na(missing_idx), 'Absent unique-current-Entrez measurement', admission$reason[missing_idx]))
      target$P[!target$estimable] <- 1
      target$ordinary_estimable <- target$estimable & is.finite(target$ordinary_P)
      target$ordinary_P[!target$ordinary_estimable] <- 1
      target$BH_U2 <- p.adjust(target$P, 'BH', n = 119L)
      target$BY_U2 <- p.adjust(target$P, 'BY', n = 119L)
      target$ordinary_BH_U2 <- p.adjust(target$ordinary_P, 'BH', n = 119L)
      target$ordinary_BY_U2 <- p.adjust(target$ordinary_P, 'BY', n = 119L)
      target$mode <- config$mode; target$contrast <- colnames(C)[column]; target$N <- nrow(m)
      target$family_size <- 119L
      target$family_scope <- if (column == 1L) 'Fixed mean-contrast U2 family' else 'Separate SECONDARY interaction U2 family'
      target$mode_role <- if (config$mode == 'rawmean34') 'primary_second_cohort_processing' else 'fixed_sensitivity_not_selectable'
      target$source_processing_limit <- if (config$input == 'raw') raw_qa$result else 'Source baseline/export unresolved; retained predefined sensitivity'
      stopifnot(nrow(target) == 119L, identical(target$Entrez_ID, U2$Entrez_ID))
      targets[[column]] <- target
    }
    t2 <- targets[[1]]
    if (config$input == 'raw') stopifnot(all(!t2$estimable[match(c('401546','79605'), t2$Entrez_ID)]))
    t2$P_up <- t2$P_down <- rep(1,119L)
    ok2 <- t2$estimable
    t2$P_up[ok2] <- pt(t2$t[ok2], t2$df[ok2], lower.tail = FALSE)
    t2$P_down[ok2] <- pt(t2$t[ok2], t2$df[ok2])
    joint <- U2
    joint$primary144_estimable <- primary_components$estimable
    joint$primary144_effect <- primary_components$effect; joint$primary144_t <- primary_components$t
    joint$primary144_df <- primary_components$df
    joint$primary144_up <- primary_components$P_up; joint$primary144_down <- primary_components$P_down
    joint$second_estimable <- t2$estimable; joint$second_effect <- t2$effect
    joint$second_t <- t2$t; joint$second_df <- t2$df
    joint$second_up <- t2$P_up; joint$second_down <- t2$P_down
    joint$joint_estimable <- joint$primary144_estimable & joint$second_estimable
    joint$jointUp <- pmax(joint$primary144_up, joint$second_up)
    joint$jointDown <- pmax(joint$primary144_down, joint$second_down)
    joint$jointP <- pmin(1, 2 * pmin(joint$jointUp, joint$jointDown))
    joint$joint_BH_U2 <- p.adjust(joint$jointP, 'BH', n = 119L)
    joint$joint_BY_U2 <- p.adjust(joint$jointP, 'BY', n = 119L)
    joint$direction_with_smaller_joint_one_sided_P <- ifelse(!joint$joint_estimable, 'not_estimable',
                          ifelse(joint$jointUp < joint$jointDown, 'up', ifelse(joint$jointDown < joint$jointUp, 'down', 'tie')))
    joint$same_observed_nonzero_sign <- joint$joint_estimable &
                          sign(joint$primary144_effect) == sign(joint$second_effect) & sign(joint$primary144_effect) != 0
    joint$primary_missing_reason <- primary_components$missing_reason; joint$second_missing_reason <- t2$missing_reason
    joint$mode <- config$mode; joint$N_second <- nrow(m); joint$family_size <- 119L
    joint$joint_role <- if (config$mode == 'rawmean34') 'PRIMARY_same_direction_joint_family' else 'fixed_joint_sensitivity_not_selectable'
    stopifnot(all(is.finite(joint$jointP)), all(joint$jointP[!joint$joint_estimable] == 1))
    mean_targets[[config$mode]] <- t2; interaction_targets[[config$mode]] <- targets[[2]]
    joints[[config$mode]] <- joint
    write_csv(t2, paste0(config$mode, '_all119_mean_contrast_results.csv'))
    write_csv(targets[[2]], paste0(config$mode, '_all119_SECONDARY_interaction_results.csv'))
    write_csv(joint, paste0(config$mode, '_all119_same_direction_IUT.csv'))
    mode_summaries[[config$mode]] <- list(N = nrow(m), cell_counts = setNames(d$counts, cells), design_rank = d$rank,
      background_genes = nrow(expr), fixed_family = 119L, mean_estimable = sum(t2$estimable),
      mean_BH_U2_hits = sum(t2$BH_U2 < .05), mean_BY_U2_hits = sum(t2$BY_U2 < .05),
      SECONDARY_interaction_BH_U2_hits = sum(targets[[2]]$BH_U2 < .05),
      SECONDARY_interaction_BY_U2_hits = sum(targets[[2]]$BY_U2 < .05),
      joint_estimable = sum(joint$joint_estimable), joint_BH_U2_hits = sum(joint$joint_BH_U2 < .05),
      joint_BY_U2_hits = sum(joint$joint_BY_U2 < .05), primary_mode = config$mode == 'rawmean34')
    cat('MODE_COMPLETE', config$mode, 'all 119 mean, interaction and joint rows retained\n')
  }
  write_csv(do.call(rbind, mean_targets), 'all_six_modes_714_U2_mean_results.csv')
  write_csv(do.call(rbind, interaction_targets), 'all_six_modes_714_U2_SECONDARY_interaction_results.csv')
  write_csv(do.call(rbind, joints), 'all_six_modes_714_U2_same_direction_IUT.csv')
  verify_inputs()
  capture.output(sessionInfo(), file = file.path(out, 'R_secondary_clinical_session.txt'))
  write_json_out(list(status = 'completed_pending_independent_statistical_QA',
    completed_at_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE), freeze_sha256 = sha(freeze_path),
    modes = mode_summaries, warnings = warnings_seen,
    preprocessing_QA = raw_qa$result, preprocessing_strict_warning = raw_failed,
    primary_QA = primary_qa$result, source_processed_baseline_export_unresolved = TRUE,
    all_original_primary_and_preparation_hashes_unchanged = TRUE,
    mean_family = '119 per fixed mode', SECONDARY_interaction_family = 'Separate 119 per fixed mode',
    joint_primary = 'rawmean34 versus exact primary_mean_144',
    five_other_joint_modes = 'Fixed sensitivities; never replace joint primary',
    no_diagnostic_performance_or_scores = TRUE), 'secondary_summary.json')
  cat('SECONDARY_ALL_FIXED_FAMILIES_COMPLETE; independent statistical QA pending.\n')
}
tryCatch(withCallingHandlers(main(), warning = function(w) {
  warnings_seen <<- c(warnings_seen, paste(current_mode, conditionMessage(w), sep = ': '))
}), error = function(e) {
  write_json_out(list(status = 'failed_no_automatic_retry', mode = current_mode,
    failed_at_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE), message = conditionMessage(e),
    warnings = warnings_seen), 'secondary_failure.json')
  stop(e)
})
