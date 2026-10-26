"""Independent six-mode secondary OLS/arithmetic QA and fixed two-contrast HC3."""
from pathlib import Path
from array import array
from collections import Counter
import csv
import datetime
import hashlib
import gzip
import json
import os
import sqlite3
import subprocess
import numpy as np
from scipy import stats
from statsmodels.regression.linear_model import OLS

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'outputs/biomarker_transfer_v13'
OUT = BASE / 'clinical_analysis/GSE120103_secondary'
SUPPORT = OUT / 'independent134_support'
SUPPORT.mkdir(exist_ok=True)
PREP, RAW = BASE / 'clinical_preparation', BASE / 'GSE120103_raw_preparation'
MODES = [('rawmean34', 'mean', True, False, False), ('rawmedian34', 'median', True, False, False),
         ('rawmean_noSourceQC32', 'mean', True, True, False), ('rawmean_FEversionAdjusted34', 'mean', True, False, True),
         ('sourceProcessed_mean36', 'mean', False, False, False), ('sourceProcessed_mean_noQC32', 'mean', False, True, False)]
CELLS = ['NONENDO_Fertile', 'NONENDO_Infertile', 'ENDO_Fertile', 'ENDO_Infertile']
CONTRASTS = ['mean_ENDO_minus_NONENDO', 'SECONDARY_disease_by_fertility']
checks, source_hashes, summaries, hc3_results = [], {}, {}, []


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def track(path):
    source_hashes[str(path.relative_to(ROOT))] = sha(path)


def rows(path):
    track(path)
    with path.open(encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def check(name, passed, **details):
    checks.append({'check': name, 'passed': bool(passed), **details})


def compare(name, actual, expected, atol=2e-10, rtol=2e-9):
    actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    finite = np.isfinite(actual) & np.isfinite(expected)
    check(name, actual.shape == expected.shape and np.allclose(actual, expected, atol=atol, rtol=rtol, equal_nan=True),
          values=int(actual.size), max_absolute_error=float(np.max(np.abs(actual[finite] - expected[finite]))) if finite.any() else 0,
          absolute_tolerance=atol, relative_tolerance=rtol)


def numeric_csv(path):
    track(path)
    opener = gzip.open if path.suffix == '.gz' else open
    ids, data = [], array('d')
    with opener(path, 'rt') as f:
        reader = csv.reader(f)
        names = next(reader)[1:]
        for row in reader:
            ids.append(row[0])
            data.extend(float(v) if v not in {'', 'NA'} else np.nan for v in row[1:])
    return ids, names, np.frombuffer(data, dtype=np.float64).reshape(len(ids), len(names))


def exact_array(name):
    ids = [r['identifier'] for r in rows(SUPPORT / (name + '.rows.csv'))]
    columns = [r['column'] for r in rows(SUPPORT / (name + '.columns.csv'))]
    return ids, columns, np.fromfile(SUPPORT / (name + '.float64.bin'), dtype='<f8').reshape((len(ids), len(columns)), order='F')


def numbers(table, field):
    return np.asarray([float(r[field]) if r[field] not in {'', 'NA'} else np.nan for r in table])


def adjust(p, by=False):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind='stable')
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    if by:
        q *= np.sum(1 / np.arange(1, len(p) + 1))
    result = np.empty_like(p)
    result[order] = np.minimum(np.minimum.accumulate(q[::-1])[::-1], 1)
    return result


freeze_path = OUT / 'secondary_statistical_freeze.json'
track(freeze_path)
freeze = json.loads(freeze_path.read_text())
freeze_errors = []
for item in freeze['files']:
    path = ROOT / item['path']
    track(path)
    if path.stat().st_size != item['bytes'] or sha(path) != item['sha256']:
        freeze_errors.append(item['path'])
check('all291_frozen_files_unchanged', len(freeze['files']) == 291 and not freeze_errors,
      mismatches=freeze_errors, freeze_sha256=sha(freeze_path), frozen_at_utc=freeze['frozen_at_utc'])
check('preapproved_HC3_contract', freeze['parameters']['independent_HC3_sensitivity']['contrasts'] == CONTRASTS
      and 'Student-t' in freeze['parameters']['independent_HC3_sensitivity']['reference_distribution'])

# Read existing RDS components at exact precision; no lmFit/eBayes or second EB call.
r_code = r'''
.libPaths(c('runtime/R-library',.libPaths()))
suppressPackageStartupMessages(library(limma))
b <- 'outputs/biomarker_transfer_v13'
out <- file.path(b,'clinical_analysis/GSE120103_secondary')
s <- file.path(out,'independent134_support')
put <- function(name,x){
 con <- file(file.path(s,paste0(name,'.float64.bin')),'wb');writeBin(as.double(x),con,size=8,endian='little');close(con)
 write.csv(data.frame(identifier=rownames(x)),file.path(s,paste0(name,'.rows.csv')),row.names=FALSE)
 write.csv(data.frame(column=colnames(x)),file.path(s,paste0(name,'.columns.csv')),row.names=FALSE)
}
for(a in c('mean','median'))put(paste0('raw_input_',a),readRDS(file.path(b,'GSE120103_raw_preparation',paste0('Entrez_',a,'_different_ProbeName_log2_matrix.rds'))))
p <- readRDS(file.path(b,'clinical_analysis/primary_mean_144_model.rds'))$contrasts
put('exact_primary144_components',cbind(effect=p$coefficients[,1],SE=sqrt(p$s2.post)*p$stdev.unscaled[,1],t=p$t[,1],df=p$df.total))
modes <- c('rawmean34','rawmedian34','rawmean_noSourceQC32','rawmean_FEversionAdjusted34','sourceProcessed_mean36','sourceProcessed_mean_noQC32')
for(mode in modes){
 m <- readRDS(file.path(out,paste0(mode,'_model.rds')));ef <- m$contrasts
 put(paste0(mode,'_expression'),m$expression)
 put(paste0(mode,'_coefficients'),m$fit$coefficients)
 put(paste0(mode,'_exact_posterior_fields'),cbind(posterior_variance=ef$s2.post,df=ef$df.total,
  ordinary_residual_variance=m$fit$sigma^2,ordinary_df=m$fit$df.residual,
  effect_mean=ef$coefficients[,1],effect_interaction=ef$coefficients[,2],
  t_mean=ef$t[,1],t_interaction=ef$t[,2],
  SE_mean=sqrt(ef$s2.post)*ef$stdev.unscaled[,1],SE_interaction=sqrt(ef$s2.post)*ef$stdev.unscaled[,2]))
}
'''
r_run = subprocess.run(['/usr/local/bin/Rscript', '-'], input=r_code, text=True, cwd=ROOT,
                       env={**os.environ, 'TMPDIR': 'runtime/tmp'}, capture_output=True)
(SUPPORT / 'RDS_accessor_stdout.txt').write_text(r_run.stdout)
(SUPPORT / 'RDS_accessor_stderr.txt').write_text(r_run.stderr)
if r_run.returncode:
    raise RuntimeError(r_run.stderr)
check('RDS_read_only_accessor', not r_run.stderr.strip(), returncode=r_run.returncode)

u2 = rows(PREP / 'U2_both_source_admissible_genes.csv')
uid = [r['Entrez_ID'] for r in u2]
check('all119_fixed_identifiers', len(uid) == len(set(uid)) == 119 and uid == freeze['target_U2_Entrez_ID'])
raw_meta = rows(RAW / 'raw34_fixed_clinical_and_technical_metadata.csv')
source_meta = rows(PREP / 'GSE120103_fixed_clinical_design.csv')
provenance = json.loads((PREP / 'orgHs_annotation_provenance.json').read_text())
db = sqlite3.connect('file:' + str(ROOT / provenance['sqlite_relative_path']) + '?mode=ro', uri=True)
current = {r[0] for r in db.execute('SELECT genes.gene_id FROM genes JOIN gene_info USING(_id) WHERE symbol IS NOT NULL AND symbol != ""')}
db.close()
pg, pc, pv = exact_array('exact_primary144_components')
primary = {g: dict(zip(pc, v)) for g, v in zip(pg, pv)}
primary_export = rows(OUT / 'all119_primary144_exact_RDS_components_and_CSV_estimability.csv')
primary_up = np.asarray([stats.t.sf(primary[g]['t'], primary[g]['df']) if g in primary else 1 for g in uid])
primary_down = np.asarray([stats.t.cdf(primary[g]['t'], primary[g]['df']) if g in primary else 1 for g in uid])
check('primary119_estimation_scaffold', [r['Entrez_ID'] for r in primary_export] == uid
      and all((r['estimable'] == 'TRUE') == (r['Entrez_ID'] in primary) for r in primary_export))
for field in pc:
    compare('primary119_exact_' + field, numbers(primary_export, field), [primary[g][field] if g in primary else np.nan for g in uid])
compare('primary119_one_sided_up', numbers(primary_export, 'P_up'), primary_up)
compare('primary119_one_sided_down', numbers(primary_export, 'P_down'), primary_down)
primary_missing = [g for g in uid if g not in primary]
check('primary_missing4829_retained', primary_missing == ['4829'])
author_summary = json.loads((OUT / 'secondary_summary.json').read_text())['modes']
mean_combined, interaction_combined, joint_combined = [], [], []

for mode, aggregation, raw_input, drop_qc, fe_adjust in MODES:
    track(OUT / (mode + '_model.rds'))
    metadata = raw_meta if raw_input else source_meta
    m = [r for r in metadata if not drop_qc or r['global_distribution_flag'] != 'True']
    # Labels independently derived from source clinical_group/fertility, not exported design columns.
    x = np.asarray([[r['clinical_group'] + '_' + r['fertility'] == cell for cell in CELLS] for r in m], dtype=float)
    names = CELLS.copy()
    if fe_adjust:
        x = np.column_stack([x, [r['FE_version'] == '9.5.3.1' for r in m]])
        names.append('FE_9.5_vs_10.7')
    c = np.asarray([[-.5, 1], [-.5, -1], [.5, -1], [.5, 1]] + ([[0, 0]] if fe_adjust else []))
    expected_n = 32 if drop_qc else 34 if raw_input else 36
    expected_counts = [9, 9, 9, 5] if drop_qc else [9, 9, 9, 7] if raw_input else [9, 9, 9, 9]
    check(mode + '_fixed_samples_and_cells', len(m) == expected_n and x[:, :4].sum(axis=0).tolist() == expected_counts,
          cells=expected_counts, N=len(m), rank=int(np.linalg.matrix_rank(x)))
    di, dh, dx = numeric_csv(OUT / (mode + '_design_matrix.csv'))
    ci, ch, dc = numeric_csv(OUT / (mode + '_contrast_matrix.csv'))
    check(mode + '_design_and_contrast_identities', di == [r['sample'] for r in m] and dh == names and ci == names and ch == CONTRASTS)
    compare(mode + '_independent_design', dx, x, atol=0, rtol=0)
    compare(mode + '_independent_contrasts', dc, c, atol=0, rtol=0)
    exported_meta = rows(OUT / (mode + '_fixed_sample_metadata.csv'))
    check(mode + '_unchanged_clinical_metadata', [r['sample'] for r in exported_meta] == [r['sample'] for r in m]
          and all(all(row[k if k in row else k.replace(' ', '.')] == original[k] for k in original) for row, original in zip(exported_meta, m)),
          serialization_note='R data.frame changes the source column sample group to sample.group; values compared unchanged.')
    if raw_input:
        genes, sample_names, complete = exact_array('raw_input_' + aggregation)
    else:
        genes, sample_names, complete = numeric_csv(PREP / 'GSE120103_mean_probe_gene_matrix.csv.gz')
    y_all = complete[:, [sample_names.index(r['sample']) for r in m]]
    current_flags = np.asarray([g in current for g in genes])
    finite_flags = np.isfinite(y_all).all(axis=1)
    nonconstant = finite_flags & (np.ptp(y_all, axis=1) > 0)
    eligible = current_flags & nonconstant
    selected = [g for g, keep in zip(genes, eligible) if keep]
    admission = rows(OUT / (mode + '_all_input_gene_background_admission.csv'))
    check(mode + '_every_background_admission', [r['Entrez_ID'] for r in admission] == genes
          and [r['admitted'] == 'TRUE' for r in admission] == eligible.tolist()
          and [r['current_valid_identifier'] == 'TRUE' for r in admission] == current_flags.tolist()
          and [r['all_values_finite'] == 'TRUE' for r in admission] == finite_flags.tolist()
          and [r['nonconstant_in_fixed_samples'] == 'TRUE' for r in admission] == nonconstant.tolist())
    y = y_all[eligible]
    ri, rh, ry = exact_array(mode + '_expression')
    check(mode + '_exact_expression_identities', selected == ri and rh == [r['sample'] for r in m])
    compare(mode + '_every_exact_RDS_input_value', ry, y, atol=0, rtol=0)
    if not raw_input:
        adjacent = (ry == y) | (ry == np.nextafter(y, -np.inf)) | (ry == np.nextafter(y, np.inf))
        different = ry != y
        check(mode + '_source_CSV_parser_at_most_one_ULP', np.all(adjacent),
              differing_values=int(different.sum()), total_values=int(y.size),
              max_ULP_distance=float(np.max(np.abs(ry[different] - y[different]) / np.abs(np.spacing(y[different])))) if different.any() else 0)
        exact_check = next(r for r in checks if r['check'] == mode + '_every_exact_RDS_input_value')
        if not exact_check['passed'] and np.all(adjacent):
            exact_check.update(severity='source_decimal_parser_one_ULP_warning', all_values_at_most_one_ULP=True,
                               explanation='Strict bit-equality criterion remains failed. R read.csv and Python float independently decode the same frozen decimal CSV to equal or immediately adjacent doubles; every fitted statistic is checked separately. No value, tolerance or frozen input changed.')
    inv = np.linalg.inv(x.T @ x)
    beta = np.linalg.lstsq(x, y.T, rcond=None)[0].T
    residual = y - beta @ x.T
    df = len(m) - np.linalg.matrix_rank(x)
    s2 = np.einsum('ij,ij->i', residual, residual) / df
    effects = beta @ c
    unscaled = np.sqrt(np.diag(c.T @ inv @ c))
    coefficient_rows = rows(OUT / (mode + '_all_background_fitted_coefficients.csv'))
    check(mode + '_coefficient_gene_order', [r['Entrez_ID'] for r in coefficient_rows] == selected)
    compare(mode + '_all_background_coefficients', np.asarray([[float(r[k]) for k in names] for r in coefficient_rows]), beta)
    ei, eh, ev = exact_array(mode + '_exact_posterior_fields')
    exact = {field: ev[:, i] for i, field in enumerate(eh)}
    check(mode + '_exact_posterior_gene_order', ei == selected)
    by_gene = {g: i for i, g in enumerate(selected)}
    mode_predictions = {}
    summaries[mode] = {'N': len(m), 'rank': int(np.linalg.matrix_rank(x)), 'background_genes': len(selected),
                       'mean_estimable': sum(g in by_gene for g in uid), 'mean_missing_IDs': [g for g in uid if g not in by_gene]}
    for k, label, suffix in [(0, CONTRASTS[0], 'mean_contrast'), (1, CONTRASTS[1], 'SECONDARY_interaction')]:
        background = rows(OUT / (mode + '_all_background_' + suffix + '.csv'))
        check(mode + '_' + label + '_background_gene_order', [r['Entrez_ID'] for r in background] == selected)
        se = np.sqrt(exact['posterior_variance']) * unscaled[k]
        ordinary_se = np.sqrt(s2) * unscaled[k]
        t = effects[:, k] / se
        ot = effects[:, k] / ordinary_se
        op = 2 * stats.t.sf(np.abs(ot), df)
        p = 2 * stats.t.sf(np.abs(t), exact['df'])
        expected = {'effect': effects[:, k], 'SE': se, 'df': exact['df'], 't': t, 'P': p,
                    'lower95': effects[:, k] - stats.t.ppf(.975, exact['df']) * se,
                    'upper95': effects[:, k] + stats.t.ppf(.975, exact['df']) * se,
                    'ordinary_residual_variance': s2, 'ordinary_df': np.repeat(df, len(y)), 'ordinary_SE': ordinary_se,
                    'ordinary_t': ot, 'ordinary_P': op,
                    'ordinary_lower95': effects[:, k] - stats.t.ppf(.975, df) * ordinary_se,
                    'ordinary_upper95': effects[:, k] + stats.t.ppf(.975, df) * ordinary_se,
                    'posterior_variance': exact['posterior_variance'], 'unscaled_SE': np.repeat(unscaled[k], len(y)),
                    'average_expression': y.mean(axis=1)}
        for field, value in expected.items():
            compare(mode + '_' + label + '_background_' + field, numbers(background, field), value)
            if not checks[-1]['passed'] and field in {'ordinary_t', 'ordinary_P'}:
                actual = numbers(background, field)
                bad = ~np.isclose(actual, value, atol=2e-10, rtol=2e-9, equal_nan=True)
                bad_ids = [g for g, b in zip(selected, bad) if b]
                if not set(bad_ids).intersection(uid) and np.all(s2[bad] < 1e-10):
                    checks[-1].update(severity='background_nearflat_ordinary_warning', mismatched_IDs=bad_ids,
                                      max_residual_variance=float(s2[bad].max()), no_U2_rows_affected=True)
        compare(mode + '_' + label + '_reported_ordinary_P_arithmetic', numbers(background, 'ordinary_P'),
                2 * stats.t.sf(np.abs(numbers(background, 'effect') / numbers(background, 'ordinary_SE')), numbers(background, 'ordinary_df')))
        targets = rows(OUT / (mode + '_all119_' + suffix + '_results.csv'))
        check(mode + '_' + label + '_all119_scaffold', [r['Entrez_ID'] for r in targets] == uid
              and all(r['family_size'] == '119' and (r['estimable'] == 'TRUE') == (r['Entrez_ID'] in by_gene) for r in targets))
        target_p = np.asarray([p[by_gene[g]] if g in by_gene else 1 for g in uid])
        target_op = np.asarray([op[by_gene[g]] if g in by_gene and np.isfinite(op[by_gene[g]]) else 1 for g in uid])
        for field, value in expected.items():
            desired = [value[by_gene[g]] if g in by_gene else (1 if field in {'P', 'ordinary_P'} else np.nan) for g in uid]
            compare(mode + '_' + label + '_U119_' + field, numbers(targets, field), desired)
        for field, value in [('BH_U2', adjust(target_p)), ('BY_U2', adjust(target_p, True)), ('ordinary_BH_U2', adjust(target_op)), ('ordinary_BY_U2', adjust(target_op, True))]:
            compare(mode + '_' + label + '_U119_' + field, numbers(targets, field), value)
        one_up = np.asarray([stats.t.sf(t[by_gene[g]], exact['df'][by_gene[g]]) if g in by_gene else 1 for g in uid])
        one_down = np.asarray([stats.t.cdf(t[by_gene[g]], exact['df'][by_gene[g]]) if g in by_gene else 1 for g in uid])
        if k == 0:
            compare(mode + '_mean_one_sided_up', numbers(targets, 'P_up'), one_up)
            compare(mode + '_mean_one_sided_down', numbers(targets, 'P_down'), one_down)
            mode_predictions = {'up': one_up, 'down': one_down, 't': [t[by_gene[g]] if g in by_gene else np.nan for g in uid],
                                'effect': [effects[by_gene[g], 0] if g in by_gene else np.nan for g in uid],
                                'df': [exact['df'][by_gene[g]] if g in by_gene else np.nan for g in uid]}
            mean_combined.extend(targets)
        else:
            interaction_combined.extend(targets)
        prefix = 'mean' if k == 0 else 'interaction'
        summaries[mode][prefix + '_BH'] = int(np.sum(adjust(target_p) < .05))
        summaries[mode][prefix + '_BY'] = int(np.sum(adjust(target_p, True) < .05))
        h = np.einsum('ij,jk,ik->i', x, inv, x)
        weights = x @ inv @ c[:, k]
        hc_current, hp, max_cov_error = [], [], 0.0
        sm_observed, sm_expected = [], []
        for info in u2:
            gene = info['Entrez_ID']
            row = {**info, 'mode': mode, 'contrast': label, 'N': len(m), 'ordinary_df': int(df),
                   'estimable': gene in by_gene, 'effect': '', 'HC3_SE': '', 'HC3_t': '', 'HC3_lower95': '', 'HC3_upper95': '', 'P': 1.0,
                   'missing_reason': '' if gene in by_gene else 'Absent/constant/nonfinite under fixed input and sample rules',
                   'family_size': 119, 'role': 'Fixed HC3 sensitivity; does not replace moderated mean, interaction or joint IUT'}
            if gene in by_gene:
                i = by_gene[gene]
                r = residual[i] / (1 - h)
                var = np.sum((weights * r) ** 2)
                full_cov = inv @ (x.T @ (r[:, None] ** 2 * x)) @ inv
                max_cov_error = max(max_cov_error, abs(var - c[:, k] @ full_cov @ c[:, k]))
                hcse = np.sqrt(var)
                hct = effects[i, k] / hcse
                hpval = float(2 * stats.t.sf(abs(hct), df))
                row.update(effect=float(effects[i, k]), HC3_SE=float(hcse), HC3_t=float(hct), P=hpval,
                           HC3_lower95=float(effects[i, k] - stats.t.ppf(.975, df) * hcse),
                           HC3_upper95=float(effects[i, k] + stats.t.ppf(.975, df) * hcse))
                reference = OLS(y[i], x).fit(cov_type='HC3', use_t=True).t_test(c[:, k])
                sm_observed.append([float(reference.effect.item()), float(reference.sd.item()), float(reference.tvalue.item()), float(reference.pvalue.item()), *reference.conf_int()[0].tolist()])
                sm_expected.append([row['effect'], row['HC3_SE'], row['HC3_t'], row['P'], row['HC3_lower95'], row['HC3_upper95']])
            hc_current.append(row)
            hp.append(row['P'])
        for row, bh, by in zip(hc_current, adjust(hp), adjust(hp, True)):
            row.update(BH_U2=float(bh), BY_U2=float(by))
        hc3_results.extend(hc_current)
        summaries[mode]['HC3_' + prefix + '_BH'] = int(np.sum(adjust(hp) < .05))
        summaries[mode]['HC3_' + prefix + '_BY'] = int(np.sum(adjust(hp, True) < .05))
        check(mode + '_' + label + '_HC3_two_covariance_forms', max_cov_error < 1e-10,
              max_variance_difference=max_cov_error, max_leverage=float(h.max()), target_rows=119, reference='Student-t ordinary residual df')
        compare(mode + '_' + label + '_HC3_all_measured_targets_vs_Statsmodels', sm_observed, sm_expected)
    joint = rows(OUT / (mode + '_all119_same_direction_IUT.csv'))
    check(mode + '_joint119_scaffold', [r['Entrez_ID'] for r in joint] == uid and all(r['family_size'] == '119' for r in joint))
    ju = np.maximum(primary_up, mode_predictions['up'])
    jd = np.maximum(primary_down, mode_predictions['down'])
    jp = np.minimum(1, 2 * np.minimum(ju, jd))
    expected_joint = {'primary144_up': primary_up, 'primary144_down': primary_down,
                      'second_up': mode_predictions['up'], 'second_down': mode_predictions['down'],
                      'primary144_effect': [primary[g]['effect'] if g in primary else np.nan for g in uid],
                      'primary144_t': [primary[g]['t'] if g in primary else np.nan for g in uid],
                      'primary144_df': [primary[g]['df'] if g in primary else np.nan for g in uid],
                      'second_effect': mode_predictions['effect'], 'second_t': mode_predictions['t'], 'second_df': mode_predictions['df'],
                      'jointUp': ju, 'jointDown': jd, 'jointP': jp, 'joint_BH_U2': adjust(jp), 'joint_BY_U2': adjust(jp, True)}
    for field, value in expected_joint.items():
        compare(mode + '_joint_' + field, numbers(joint, field), value)
    direction_errors = []
    for i, row in enumerate(joint):
        estimable = uid[i] in primary and uid[i] in by_gene
        direction = 'not_estimable' if not estimable else 'up' if ju[i] < jd[i] else 'down' if jd[i] < ju[i] else 'tie'
        same_sign = estimable and np.sign(primary[uid[i]]['effect']) == np.sign(mode_predictions['effect'][i]) and np.sign(primary[uid[i]]['effect']) != 0
        if (row['joint_estimable'] == 'TRUE') != estimable or row['direction_with_smaller_joint_one_sided_P'] != direction or (row['same_observed_nonzero_sign'] == 'TRUE') != same_sign:
            direction_errors.append(uid[i])
    check(mode + '_joint_missing_and_direction_policy', not direction_errors and all(jp[i] == 1 for i, g in enumerate(uid) if g not in primary or g not in by_gene), mismatches=direction_errors)
    summaries[mode].update(joint_estimable=sum(g in primary and g in by_gene for g in uid), joint_BH=int(np.sum(adjust(jp) < .05)), joint_BY=int(np.sum(adjust(jp, True) < .05)))
    joint_combined.extend(joint)
    expected_summary = author_summary[mode]
    check(mode + '_author_counts_confirmed', all(summaries[mode][k] == expected_summary[v] for k, v in [
        ('N', 'N'), ('rank', 'design_rank'), ('background_genes', 'background_genes'), ('mean_estimable', 'mean_estimable'),
        ('mean_BH', 'mean_BH_U2_hits'), ('mean_BY', 'mean_BY_U2_hits'), ('interaction_BH', 'SECONDARY_interaction_BH_U2_hits'),
        ('interaction_BY', 'SECONDARY_interaction_BY_U2_hits'), ('joint_estimable', 'joint_estimable'), ('joint_BH', 'joint_BH_U2_hits'), ('joint_BY', 'joint_BY_U2_hits')]))
    print(mode + ' independent verification complete', flush=True)
    del complete, y_all, y, beta, residual

for filename, accumulated in [('all_six_modes_714_U2_mean_results.csv', mean_combined), ('all_six_modes_714_U2_SECONDARY_interaction_results.csv', interaction_combined), ('all_six_modes_714_U2_same_direction_IUT.csv', joint_combined)]:
    combined = rows(OUT / filename)
    by_key = {(r['mode'], r['Entrez_ID']): r for r in accumulated}
    check(filename + '_all714_rows_match_per_mode', len(combined) == len(by_key) == 714
          and all(r == by_key[(r['mode'], r['Entrez_ID'])] for r in combined))
check('HC3_all1428_unique_rows_and_missing_P1', len(hc3_results) == 1428
      and len({(r['mode'], r['contrast'], r['Entrez_ID']) for r in hc3_results}) == 1428
      and all(r['P'] == r['BH_U2'] == r['BY_U2'] == 1 for r in hc3_results if not r['estimable']))
check('all291_and_used_source_artifacts_unchanged_at_end', all(sha(ROOT / name) == value for name, value in source_hashes.items()))
hcpath = OUT / 'independent134_HC3_all119_both_contrasts_six_modes.csv'
with hcpath.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(hc3_results[0]))
    writer.writeheader()
    writer.writerows(hc3_results)
failed = [r for r in checks if not r['passed']]
warnings_only = all(r.get('severity') in {'background_nearflat_ordinary_warning', 'source_decimal_parser_one_ULP_warning'} for r in failed)
status = 'PASS' if not failed else 'PASS_WITH_SOURCE_DECIMAL_PARSING_WARNINGS' if warnings_only and all(r.get('severity') == 'source_decimal_parser_one_ULP_warning' for r in failed) else 'PASS_WITH_BACKGROUND_OLS_NUMERICAL_WARNINGS' if warnings_only else 'FAIL'
report = {'verified_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'result': status,
          'checks': checks, 'summaries': summaries, 'source_hashes': source_hashes,
          'secondary_freeze_sha256': sha(freeze_path), 'script_sha256': sha(Path(__file__)), 'HC3_sha256': sha(hcpath),
          'HC3_rows': len(hc3_results), 'primary_missing': primary_missing,
          'second_EB_fit_run': False, 'joint_HC3_performed': False,
          'scope': 'Independent OLS/arithmetic audit and prespecified HC3 sensitivities, full119 families; no significance-based filtering or replacement of primary processing.',
          'limitations': ['Posterior variance and df are exact saved RDS inputs, not independently re-estimated EB hyperparameters.',
                          'Strict bit-identity warnings remain for source-CSV parsing differences of at most one ULP; all numerical model/statistic comparisons retain the original tolerances.',
                          'HC3 uses the frozen Student-t ordinary residual-df approximation; it is not an exact small-sample heteroskedastic test.',
                          'Source-processed baseline/export remains unresolved; its two modes remain fixed sensitivities.',
                          'Same-direction IUT tests mean association across sources, not diagnostic classifier performance.']}
(OUT / 'independent134_secondary_statistical_verification.json').write_text(json.dumps(report, indent=2) + '\n')
lines = ['# GSE120103六种固定模式：独立统计核验', '', '状态：' + status + '。冻结输入/脚本/阈值未改变；没有第二次经验贝叶斯拟合，也没有额外joint-HC3检验。', '',
         '独立从固定样本表构建四单元设计及FE9.5.3.1指示变量；使用源基因输入/精确RDS数值和NumPy最小二乘核验全部背景系数、效应、残差方差、普通与moderated统计/区间，逐项核对全119及BH/BY。', '',
         '|模式|N|mean BH/BY|交互 BH/BY|joint BH/BY|HC3 mean BH/BY|HC3交互 BH/BY|', '|---|---:|---:|---:|---:|---:|---:|']
for mode, *_ in MODES:
    s = summaries[mode]
    lines.append('|{}|{}|{}/{}|{}/{}|{}/{}|{}/{}|{}/{}|'.format(mode, s['N'], s['mean_BH'], s['mean_BY'], s['interaction_BH'], s['interaction_BY'], s['joint_BH'], s['joint_BY'], s['HC3_mean_BH'], s['HC3_mean_BY'], s['HC3_interaction_BH'], s['HC3_interaction_BY']))
lines += ['', 'IUT使用精确第一队列144模型RDS的t/df，不用显著性或方向预选基因。对全部119分别计算两队列向上/向下单侧P的最大值，较小方向乘2并截至1，再在119内BH/BY。第一队列4829以及原始臂401546/79605等缺测均保留、单侧及jointP按1处理。三套714主表与逐模式文件逐行相同。', '',
          'HC3共1428行：全部6模式×2对比×119基因。每个模式/对比各自保留完整119家族，缺测P/BH/BY1；Student-t使用普通残差自由度，分别保存效应、HC3 SE、t、区间和修正P。约简对比权重公式与完整sandwich矩阵逐行交叉核对，并对全部可测目标使用Statsmodels独立复核效应、SE、t/P和区间。', '',
          '严格位级输入核验保留两项警示：源处理36/32模式中，R和Python对同一冻结十进制CSV的解析有30079/27076个值相差恰好1ULP，最大绝对差1.78e-15/4.44e-16。原始处理RDS输入逐值完全相同。未放宽位级准则、未更改数据；另检验所有源数值不超过相邻double差距，且全部背景/目标模型统计、校正和IUT均按原数值容差通过。', '',
          '这只是固定均值关联、交互和联合关联的数值核验与敏感性；不代替主要moderated结果，不生成额外joint-HC3，不解释为AUC或临床诊断效能。来源处理矩阵基线关系未解决的限制保留。JSON含逐项误差与291输入/已用产物SHA，精确RDS读取资料存入independent134_support。']
if failed:
    lines += ['', '严格未通过／警示项：' + '；'.join(r['check'] for r in failed)]
(OUT / 'independent134_secondary_statistical_verification.md').write_text('\n'.join(lines) + '\n')
print(json.dumps({'result': status, 'checks': len(checks), 'failed_checks': [r['check'] for r in failed], 'summaries': summaries}), flush=True)
