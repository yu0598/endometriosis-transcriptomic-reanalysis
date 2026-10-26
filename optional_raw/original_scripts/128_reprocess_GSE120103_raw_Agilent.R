# Frozen technical preprocessing only; no disease contrasts or candidate selection.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args %in% c('--freeze-only', '--run'))
root <- normalizePath(getwd())
stopifnot(basename(root) == 'IUA_TE_内膜修复研究_2026', Sys.getenv('TMPDIR') == 'runtime/tmp')
.libPaths(c(file.path(root, 'runtime/R-library'), .libPaths()))
suppressPackageStartupMessages({library(limma); library(jsonlite); library(digest)})
out <- 'outputs/biomarker_transfer_v13/GSE120103_raw_preparation'
rawdir <- 'data_raw/GSE120103_rawFE_v13'
sourcedir <- 'outputs/biomarker_transfer_v13/phase2_investigation/GSE120103_processing_sources'
clinical <- 'outputs/biomarker_transfer_v13/clinical_preparation'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(warn = 1, digits = 17, device = function(...) pdf(file.path(out, 'unused_default_device.pdf'), ...))
freeze_path <- file.path(out, 'preprocessing_freeze.json')
script_path <- 'scripts/128_reprocess_GSE120103_raw_Agilent.R'
sha <- function(path) digest(file = path, algo = 'sha256')
read_csv <- function(path) read.csv(path, colClasses = 'character', check.names = FALSE,
                                  fileEncoding = 'UTF-8-BOM', na.strings = NULL)
write_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = '')
write_json_out <- function(x, name) write_json(x, file.path(out, name), pretty = TRUE,
                                              auto_unbox = TRUE, digits = NA, na = 'null')
parameters <- list(
  cohort = 'GSE120103', source_files = 36L, admitted_files = 34L,
  raw_only_exclusions = c('GSM3393522', 'GSM3393523'),
  exclusion_rule = 'Incomplete public FEATURES records; source-processed 36-sample arm unchanged',
  spot_admission = 'ControlType=0 and finite gMedianSignal/gBGMedianSignal in all 34 admitted files',
  foreground = 'gMedianSignal', background = 'gBGMedianSignal',
  object_class = 'EListRaw', background_method = 'normexp', normexp_method = 'saddle', offset = 50,
  normalization = 'normalizeBetweenArrays(EListRaw,method=quantile): linear quantile then log2',
  duplicate_spots = 'limma::avereps on log2-normalized physical spots, grouped by ProbeName',
  annotation = 'Official frozen GPL6480.annot.gz; distinct trimmed Gene ID tokens split on ///',
  gene_admission = 'Exactly one positive numeric Entrez token, present in frozen orgHs Entrez-to-symbol export',
  gene_primary = 'Arithmetic mean across different ProbeName values',
  gene_sensitivity = 'Median across different ProbeName values',
  target_U2 = 119L, target_U2_action = 'Freeze and report availability only; never reduce family',
  no_intensity_cutoff = TRUE, no_per_feature_outlier_flag_filter = TRUE,
  no_biological_contrasts = TRUE, no_candidate_effect_selection = TRUE,
  no_cross_cohort_correction = TRUE, no_sample_or_value_imputation = TRUE,
  exact_baseline_of_source_processed_matrix = 'Unresolved; no linear/log mixture assumed')
versions <- list(R = R.version.string, limma = as.character(packageVersion('limma')),
                 jsonlite = as.character(packageVersion('jsonlite')),
                 digest = as.character(packageVersion('digest')))

if (args == '--freeze-only') {
  stopifnot(!file.exists(freeze_path))
  capture.output(sessionInfo(), file = file.path(out, 'preprocessing_R_session_before_freeze.txt'))
  capture.output({print(limma::backgroundCorrect); print(limma::backgroundCorrect.matrix)
                 print(limma::normalizeBetweenArrays); print(limma::normalizeQuantiles)
                 print(getS3method('avereps', 'default'))},
                 file = file.path(out, 'frozen_limma_implementation.txt'))
  inputs <- unique(c(
    'AGENTS.md', '34_临床检验执行补充_原始Agilent处理与新意边界_v13_3.md',
    '35_第二队列原始处理实施约定_v13_4.md', script_path,
    'scripts/GSE120103_rawFE_download_v13.py', 'scripts/GSE120103_rawFE_schema_inventory_v13.py',
    list.files(rawdir, recursive = TRUE, full.names = TRUE),
    list.files(sourcedir, recursive = TRUE, full.names = TRUE),
    'outputs/biomarker_transfer_v13/phase2_investigation/GSE120103_source_processing_followup.csv',
    'outputs/biomarker_transfer_v13/phase2_investigation/GSE120103_source_processing_followup.md',
    'sources/biomarker_transfer_v13/clinical_inputs/GPL6480.annot.gz',
    'sources/biomarker_transfer_v13/clinical_inputs/GSE120103_series_matrix.txt.gz',
    'sources/biomarker_transfer_v13/clinical_inputs/retrieval.json',
    'sources/endometriosis_resource_audit_v8/PMC6375207.xml',
    list.files(clinical, recursive = TRUE, full.names = TRUE),
    'runtime/R-library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite',
    file.path(out, 'preprocessing_R_session_before_freeze.txt'),
    file.path(out, 'frozen_limma_implementation.txt'),
    'runtime/R-library/limma/DESCRIPTION', 'runtime/R-library/limma/R/limma.rdb',
    'runtime/R-library/limma/R/limma.rdx'))
  stopifnot(all(file.exists(inputs)), !any(dir.exists(inputs)))
  manifest <- data.frame(path = sort(inputs), stringsAsFactors = FALSE)
  manifest$bytes <- unname(file.info(manifest$path)$size)
  manifest$sha256 <- vapply(manifest$path, sha, character(1))
  manifest$role <- ifelse(startsWith(manifest$path, paste0(clinical, '/')),
                          'immutable existing clinical_preparation input; no effect result read',
                          'frozen raw/source/implementation/identifier input')
  u2 <- read_csv(file.path(clinical, 'U2_both_source_admissible_genes.csv'))
  stopifnot(nrow(u2) == 119L, !anyDuplicated(u2$Entrez_ID))
  schema <- read_csv(file.path(sourcedir, 'raw_FE_schema_inventory.csv'))
  stopifnot(nrow(schema) == 36L, sum(schema$source_table_shape_complete == 'True') == 34L,
            setequal(schema$gsm[schema$source_table_shape_complete != 'True'], parameters$raw_only_exclusions))
  write_json(list(freeze_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
    scope = 'Technical normalization only, after primary GSE51981 effects; before secondary GSE120103 gene effects',
    parameters = parameters, versions = versions, target_U2_Entrez_ID = u2$Entrez_ID,
    inputs = manifest, clinical_preparation_immutable = TRUE,
    source_matrix_read_permission = 'Hash only; no numerical source matrix used in raw arm',
    independently_reviewed = FALSE), freeze_path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat('FREEZE_WRITTEN', nrow(manifest), 'inputs; 36 raw files; 119 unchanged U2 identifiers. No normalization run.\n')
  quit(save = 'no', status = 0)
}

stopifnot(file.exists(freeze_path), !file.exists(file.path(out, 'preprocessing_status.json')))
freeze <- read_json(freeze_path, simplifyVector = TRUE)
stopifnot(identical(freeze$versions, versions), isTRUE(all.equal(freeze$parameters, parameters, check.attributes = FALSE)))
input_manifest <- freeze$inputs
verify_inputs <- function() {
  stopifnot(all(file.exists(input_manifest$path)),
            all(unname(file.info(input_manifest$path)$size) == input_manifest$bytes),
            identical(unname(vapply(input_manifest$path, sha, character(1))), input_manifest$sha256))
}
verify_inputs()
write_json_out(list(started_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
                    freeze_sha256 = sha(freeze_path), scope = 'technical only'), 'preprocessing_started.json')
warnings_seen <- character()
main <- function() {
  schema <- read_csv(file.path(sourcedir, 'raw_FE_schema_inventory.csv'))
  retrieval <- read_csv(file.path(rawdir, 'rawFE_retrieval_manifest.csv'))
  source <- read_csv(file.path('outputs/biomarker_transfer_v13/phase2_investigation',
                               'GSE120103_source_processing_followup.csv'))
  schema <- schema[order(schema$gsm), ]; retrieval <- retrieval[match(schema$gsm, retrieval$gsm), ]
  source <- source[match(schema$gsm, source$gsm), ]
  stopifnot(identical(schema$gsm, retrieval$gsm), identical(schema$gsm, source$gsm))
  admitted <- schema$source_table_shape_complete == 'True' &
              schema$ordered_identity_equal_reference == 'True' & schema$malformed_feature_records == '0'
  stopifnot(sum(admitted) == 34L, setequal(schema$gsm[!admitted], parameters$raw_only_exclusions))
  sample_meta <- read_csv(file.path(clinical, 'GSE120103_fixed_clinical_design.csv'))
  sample_meta <- sample_meta[match(schema$gsm, sample_meta$sample), ]
  stopifnot(identical(sample_meta$sample, schema$gsm), all(sample_meta$series == 'GSE120103'))
  sample_meta$raw_arm_admitted <- admitted
  sample_meta$raw_arm_exclusion_reason <- ifelse(admitted, '', 'Public raw FEATURES table incomplete; no repair')
  sample_meta$source_processed_arm_retained <- TRUE
  sample_meta$raw_source_file <- schema$source_file
  sample_meta$raw_source_sha256 <- retrieval$compressed_sha256
  sample_meta$FE_version <- source$actual_FE_version
  sample_meta$FE_protocol <- source$actual_FE_protocol
  sample_meta$scanner <- source$actual_scanner
  sample_meta$grid_name <- source$actual_grid_name
  sample_meta$scan_date <- source$raw_Scan_Date
  sample_meta$microns_per_pixel_X <- source$raw_Scan_MicronsPerPixelX
  sample_meta$grid_genomic_build <- source$raw_Grid_GenomicBuild
  sample_meta$raw_arm_include_source_global_QC_sensitivity <- admitted & sample_meta$global_distribution_flag == 'False'
  write_csv(sample_meta, 'all36_raw_admission_and_source_metadata.csv')
  samples <- sample_meta[admitted, ]; rownames(samples) <- samples$sample
  write_csv(samples, 'raw34_fixed_clinical_and_technical_metadata.csv')
  counts <- as.data.frame(table(factor(samples$primary_cell,
               levels = c('NONENDO_Fertile', 'NONENDO_Infertile', 'ENDO_Fertile', 'ENDO_Infertile'))))
  names(counts) <- c('primary_cell', 'n')
  stopifnot(identical(counts$n, c(9L, 9L, 9L, 7L)), min(counts$n) >= 2L)
  write_csv(counts, 'raw34_clinical_cell_counts.csv')
  write_csv(as.data.frame(table(FE_version = samples$FE_version, clinical_cell = samples$primary_cell)),
            'raw34_FE_version_by_clinical_cell_counts.csv')
  readin <- list(); foreground <- background <- NULL; reference <- NULL
  for (i in which(admitted)) {
    item <- retrieval[i, ]; gsm <- item$gsm
    filename <- file.path(rawdir, item$file)
    stopifnot(sha(filename) == item$compressed_sha256,
              file.info(filename)$size == as.numeric(item$filelist_bytes))
    con <- gzfile(filename, 'rb')
    payload <- tryCatch(readBin(con, 'raw', n = as.numeric(item$uncompressed_bytes) + 1), finally = close(con))
    stopifnot(length(payload) == as.numeric(item$uncompressed_bytes), !any(payload == as.raw(0)),
              digest(payload, algo = 'sha256', serialize = FALSE) == item$uncompressed_sha256)
    lines <- strsplit(rawToChar(payload), '\n', fixed = TRUE)[[1]]
    lines <- sub('\r$', '', lines)
    start <- which(startsWith(lines, 'FEATURES\t'))
    stopifnot(length(start) == 1L)
    columns <- strsplit(lines[start], '\t', fixed = TRUE)[[1]][-1L]
    data_lines <- lines[seq.int(start + 1L, length(lines))]
    data_lines <- data_lines[startsWith(data_lines, 'DATA\t')]
    parts <- strsplit(data_lines, '\t', fixed = TRUE)
    stopifnot(length(parts) == 45015L, all(lengths(parts) == length(columns) + 1L))
    keep_columns <- c('FeatureNum', 'ProbeName', 'ControlType', 'gMedianSignal', 'gBGMedianSignal')
    stopifnot(all(keep_columns %in% columns))
    values <- do.call(rbind, lapply(parts, function(x) x[match(keep_columns, columns) + 1L]))
    identity <- data.frame(FeatureNum = as.integer(values[, 1]), ProbeName = values[, 2],
                            ControlType = values[, 3], stringsAsFactors = FALSE)
    stopifnot(!anyNA(identity), !anyDuplicated(identity$FeatureNum))
    json_probe <- encodeString(identity$ProbeName, quote = '"')
    identity_text <- paste0('[', identity$FeatureNum, ',', json_probe, ',"', identity$ControlType, '"]', collapse = '\n')
    identity_sha <- digest(charToRaw(identity_text), algo = 'sha256', serialize = FALSE)
    stopifnot(identity_sha == schema$ordered_identity_sha256[i])
    if (is.null(reference)) reference <- identity
    stopifnot(identical(reference, identity))
    fg <- as.numeric(values[, 4]); bg <- as.numeric(values[, 5])
    foreground <- cbind(foreground, fg); background <- cbind(background, bg)
    readin[[length(readin) + 1L]] <- data.frame(sample = gsm, raw_file = filename,
      compressed_bytes = file.info(filename)$size, compressed_sha256 = sha(filename),
      uncompressed_bytes = length(payload), uncompressed_sha256 = digest(payload, algo = 'sha256', serialize = FALSE),
      feature_rows = nrow(identity), feature_columns = length(columns),
      ordered_identity_sha256 = identity_sha, matches_frozen_identity = TRUE, matches_all34_order = TRUE,
      foreground_finite = sum(is.finite(fg)), background_finite = sum(is.finite(bg)))
    cat('READ_VERIFIED', gsm, nrow(identity), 'features\n'); flush.console()
  }
  colnames(foreground) <- colnames(background) <- samples$sample
  rownames(foreground) <- rownames(background) <- as.character(reference$FeatureNum)
  write_csv(do.call(rbind, readin), 'raw34_exact_readin_and_identity_checks.csv')
  ordinary <- reference$ControlType == '0'
  finite <- rowSums(!is.finite(foreground) | !is.finite(background)) == 0L
  retained <- ordinary & finite
  spot_admission <- reference
  spot_admission$foreground_nonfinite_samples <- rowSums(!is.finite(foreground))
  spot_admission$background_nonfinite_samples <- rowSums(!is.finite(background))
  spot_admission$retained <- retained
  spot_admission$reason <- ifelse(!ordinary, 'ControlType is not 0',
                                  ifelse(!finite, 'Nonfinite foreground/background in at least one admitted file', 'Retained'))
  write_csv(spot_admission, 'all45015_spot_admission_audit.csv')
  stopifnot(sum(retained) > 0L, !anyNA(reference$ProbeName[retained]),
            all(nzchar(reference$ProbeName[retained])))
  raw <- new('EListRaw', list(E = foreground[retained, , drop = FALSE],
               Eb = background[retained, , drop = FALSE], genes = reference[retained, ], targets = samples))
  stopifnot(all(is.finite(raw$E)), all(is.finite(raw$Eb)), is(raw, 'EListRaw'))
  saveRDS(list(foreground = foreground, background = background, identities = reference),
           file.path(out, 'all45015_raw_foreground_background.rds'), compress = 'gzip')
  saveRDS(raw, file.path(out, 'ordinary_raw_EListRaw.rds'), compress = 'gzip')
  export_matrix <- function(x, id_name, name) {
    tab <- data.frame(rownames(x), x, check.names = FALSE); names(tab)[1] <- id_name
    con <- gzfile(file.path(out, name), 'wt')
    tryCatch(write.csv(tab, con, row.names = FALSE, na = ''), finally = close(con))
  }
  export_matrix(foreground, 'FeatureNum', 'all45015_raw_foreground.csv.gz')
  export_matrix(background, 'FeatureNum', 'all45015_raw_background.csv.gz')
  export_matrix(raw$E, 'FeatureNum', 'ordinary_raw_foreground.csv.gz')
  export_matrix(raw$Eb, 'FeatureNum', 'ordinary_raw_background.csv.gz')
  cat('NORMALIZATION_START', nrow(raw$E), 'ordinary physical spots across', ncol(raw$E), 'samples\n')
  corrected <- backgroundCorrect(raw, method = 'normexp', normexp.method = 'saddle', offset = 50)
  stopifnot(is(corrected, 'EListRaw'), is.null(corrected$Eb), all(is.finite(corrected$E)), all(corrected$E > 0))
  saveRDS(corrected, file.path(out, 'ordinary_background_corrected_EListRaw.rds'), compress = 'gzip')
  export_matrix(corrected$E, 'FeatureNum', 'ordinary_background_corrected_linear.csv.gz')
  normalized <- normalizeBetweenArrays(corrected, method = 'quantile')
  stopifnot(is(normalized, 'EList'), all(is.finite(normalized$E)),
            identical(dimnames(normalized$E), dimnames(raw$E)))
  saveRDS(normalized, file.path(out, 'ordinary_normalized_log2_EList.rds'), compress = 'gzip')
  export_matrix(normalized$E, 'FeatureNum', 'ordinary_normalized_log2_spots.csv.gz')
  export_matrix(2^normalized$E, 'FeatureNum', 'ordinary_quantile_normalized_linear.csv.gz')
  probes <- avereps(normalized$E, ID = normalized$genes$ProbeName)
  stopifnot(!anyDuplicated(rownames(probes)), all(is.finite(probes)),
            identical(rownames(probes), unique(normalized$genes$ProbeName)))
  export_matrix(probes, 'ProbeName', 'ProbeName_mean_spot_log2_matrix.csv.gz')
  saveRDS(probes, file.path(out, 'ProbeName_mean_spot_log2_matrix.rds'), compress = 'gzip')
  spot_counts <- table(normalized$genes$ProbeName)
  write_csv(data.frame(ProbeName = rownames(probes),
                       physical_spots_averaged = as.integer(spot_counts[rownames(probes)])),
            'ProbeName_physical_spot_counts.csv')
  con <- gzfile('sources/biomarker_transfer_v13/clinical_inputs/GPL6480.annot.gz', 'rt')
  ann_lines <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
  a <- which(ann_lines == '!platform_table_begin'); b <- which(ann_lines == '!platform_table_end')
  stopifnot(length(a) == 1L, length(b) == 1L, b > a)
  ann <- read.delim(text = paste(ann_lines[seq.int(a + 1L, b - 1L)], collapse = '\n'),
                    sep = '\t', quote = '', comment.char = '', fill = FALSE,
                    colClasses = 'character', check.names = FALSE, na.strings = NULL)
  stopifnot(!anyDuplicated(ann$ID), all(c('ID', 'Gene symbol', 'Gene ID') %in% names(ann)))
  symbols <- read_csv(file.path(clinical, 'orgHs_all_Entrez_to_symbol.csv'))
  valid_ids <- unique(symbols$ENTREZID)
  unique_entrez <- function(x) {
    tokens <- unique(trimws(strsplit(x, '///', fixed = TRUE)[[1]])); tokens <- tokens[nzchar(tokens)]
    if (length(tokens) == 1L && grepl('^[1-9][0-9]*$', tokens)) tokens else ''
  }
  ann$unique_Entrez_ID <- vapply(ann[['Gene ID']], unique_entrez, character(1))
  ann$NCBI_identifier_status <- ifelse(!nzchar(ann$unique_Entrez_ID), 'blank_or_ambiguous_Entrez',
                     ifelse(ann$unique_Entrez_ID %in% valid_ids, 'present_in_frozen_current_Entrez_background',
                            'absent_from_frozen_current_Entrez_background'))
  write_csv(ann[, c('ID', 'Gene symbol', 'Gene ID', 'unique_Entrez_ID', 'NCBI_identifier_status')],
            'all_official_GPL6480_Entrez_mapping_audit.csv')
  idx <- match(rownames(probes), ann$ID)
  mapping <- data.frame(ProbeName = rownames(probes), physical_spots_averaged = as.integer(spot_counts[rownames(probes)]),
                        official_Gene_symbol = ifelse(is.na(idx), '', ann[['Gene symbol']][idx]),
                        official_Gene_ID = ifelse(is.na(idx), '', ann[['Gene ID']][idx]),
                        unique_Entrez_ID = ifelse(is.na(idx), '', ann$unique_Entrez_ID[idx]),
                        NCBI_identifier_status = ifelse(is.na(idx), 'ProbeName_absent_from_official_GPL6480',
                                                        ann$NCBI_identifier_status[idx]))
  mapping$retained_for_gene_aggregation <- mapping$NCBI_identifier_status == 'present_in_frozen_current_Entrez_background'
  write_csv(mapping, 'all_ProbeName_to_Entrez_mapping_and_admission.csv')
  valid <- mapping$retained_for_gene_aggregation
  entrez <- mapping$unique_Entrez_ID[valid]
  ids <- sort(unique(entrez))
  gene_mean <- avereps(probes[valid, , drop = FALSE], ID = entrez)[ids, , drop = FALSE]
  gene_median <- t(vapply(ids, function(id) apply(probes[valid, , drop = FALSE][entrez == id, , drop = FALSE],
                                                2L, median), numeric(ncol(probes))))
  dimnames(gene_median) <- list(ids, colnames(probes))
  stopifnot(identical(dimnames(gene_mean), dimnames(gene_median)), all(is.finite(gene_mean)), all(is.finite(gene_median)))
  for (mode in c('mean', 'median')) {
    matrix <- if (mode == 'mean') gene_mean else gene_median
    export_matrix(matrix, 'Entrez_ID', paste0('Entrez_', mode, '_different_ProbeName_log2_matrix.csv.gz'))
    saveRDS(matrix, file.path(out, paste0('Entrez_', mode, '_different_ProbeName_log2_matrix.rds')), compress = 'gzip')
  }
  gene_meta <- data.frame(Entrez_ID = ids,
               current_symbol = vapply(ids, function(id) paste(unique(symbols$SYMBOL[symbols$ENTREZID == id]), collapse = ';'), character(1)),
               different_ProbeNames = as.integer(table(entrez)[ids]))
  write_csv(gene_meta, 'gene_different_ProbeName_counts_and_identifiers.csv')
  u2 <- read_csv(file.path(clinical, 'U2_both_source_admissible_genes.csv'))
  stopifnot(nrow(u2) == 119L, identical(u2$Entrez_ID, freeze$target_U2_Entrez_ID))
  u2$measured_in_raw_mean_and_median <- u2$Entrez_ID %in% ids
  u2$family_size_unchanged <- 119L
  write_csv(u2, 'U2_119_frozen_gene_availability_only.csv')
  qc_matrix <- function(x, stage, scale) {
    do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
      y <- x[, j]; q <- quantile(y, c(0, .05, .25, .5, .75, .95, 1), names = FALSE)
      data.frame(sample = colnames(x)[j], stage = stage, scale = scale, features = length(y),
        finite = sum(is.finite(y)), min = q[1], p05 = q[2], p25 = q[3], median = q[4],
        p75 = q[5], p95 = q[6], max = q[7], mean = mean(y), SD = sd(y), IQR = q[5] - q[3])
    }))
  }
  qc <- rbind(qc_matrix(raw$E, 'ordinary_raw_foreground', 'linear'),
              qc_matrix(raw$Eb, 'ordinary_raw_background', 'linear'),
              qc_matrix(corrected$E, 'ordinary_background_corrected', 'linear'),
              qc_matrix(normalized$E, 'ordinary_quantile_normalized_spots', 'log2'),
              qc_matrix(probes, 'ProbeName_mean_spots', 'log2'),
              qc_matrix(gene_mean, 'Entrez_mean_different_probes', 'log2'),
              qc_matrix(gene_median, 'Entrez_median_different_probes', 'log2'))
  write_csv(qc, 'all34_sample_technical_distribution_QC_by_stage.csv')
  dimensions <- data.frame(stage = c('all_raw_features', 'ordinary_finite_spots', 'distinct_ProbeNames',
                                     'gene_admitted_ProbeNames', 'unique_current_Entrez_genes'),
                            rows = c(nrow(foreground), nrow(raw$E), nrow(probes), sum(valid), length(ids)), samples = 34L)
  write_csv(dimensions, 'preprocessing_stage_dimensions.csv')
  verify_inputs()
  capture.output(sessionInfo(), file = file.path(out, 'preprocessing_R_session_after_run.txt'))
  result <- list(status = 'completed_pending_independent_normalization_and_aggregation_QA',
    completed_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE), freeze_sha256 = sha(freeze_path),
    dimensions = dimensions, clinical_cell_counts = counts,
    all36_source_files_frozen = TRUE, raw_samples = 34L, excluded_raw_only = parameters$raw_only_exclusions,
    source_matrix_36_unchanged = TRUE, clinical_preparation_hashes_unchanged = TRUE,
    current_identifier_background = 'Frozen orgHs 3.23.1 Entrez-to-symbol export, EGSOURCEDATE 2026-Mar18; no guessed redirects',
    U2_family_size = 119L, U2_measured = sum(u2$measured_in_raw_mean_and_median),
    own_biological_contrasts = FALSE, candidate_effect_selection = FALSE, plots_created = FALSE,
    warnings = warnings_seen)
  write_json_out(result, 'preprocessing_status.json')
  files <- list.files(out, recursive = TRUE, full.names = TRUE)
  files <- files[!endsWith(files, 'preprocessing_output_manifest.csv')]
  write_csv(data.frame(path = files, bytes = unname(file.info(files)$size),
                       sha256 = vapply(files, sha, character(1))), 'preprocessing_output_manifest.csv')
  cat('COMPLETE_TECHNICAL_ONLY', nrow(raw$E), 'spots;', nrow(probes), 'ProbeNames;', length(ids),
      'current Entrez genes;', sum(u2$measured_in_raw_mean_and_median), 'of 119 U2 measured. Independent QA pending.\n')
}
tryCatch(withCallingHandlers(main(), warning = function(w) warnings_seen <<- c(warnings_seen, conditionMessage(w))),
  error = function(e) {
    write_json_out(list(status = 'failed_no_automatic_retry', failure_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
                        message = conditionMessage(e), warnings = warnings_seen,
                        own_biological_contrasts = FALSE, no_ad_hoc_changes = TRUE), 'preprocessing_failure.json')
    stop(e)
  })
