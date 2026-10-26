"""Frozen GSE51981 OLS/arithmetic audit and fixed HC3 sensitivity; no second EB fit."""
from pathlib import Path
from array import array
from collections import Counter
from decimal import Decimal, localcontext
import csv
import datetime
import gzip
import hashlib
import json
import os
import sqlite3
import subprocess
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / 'outputs/biomarker_transfer_v13/clinical_preparation'
OUT = ROOT / 'outputs/biomarker_transfer_v13/clinical_analysis'
checks, source_hashes, summaries, high_precision_rows = [], {}, {}, []


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def track(path):
    source_hashes[str(path.relative_to(ROOT))] = sha(path)


def record(name, passed, **details):
    checks.append({'check': name, 'passed': bool(passed), **details})


def compare(name, actual, expected, atol=2e-10, rtol=2e-9):
    actual, expected = np.asarray(actual), np.asarray(expected)
    finite = np.isfinite(actual) & np.isfinite(expected)
    error = float(np.max(np.abs(actual[finite] - expected[finite]))) if finite.any() else 0.0
    record(name, np.allclose(actual, expected, atol=atol, rtol=rtol, equal_nan=True),
           values=int(actual.size), max_absolute_error=error, absolute_tolerance=atol, relative_tolerance=rtol)


def rows(path):
    track(path)
    with path.open() as f:
        return list(csv.DictReader(f))


def numeric_table(path):
    track(path)
    opener = gzip.open if path.suffix == '.gz' else open
    identifiers, data = [], array('d')
    with opener(path, 'rt') as f:
        reader = csv.reader(f)
        header = next(reader)[1:]
        for row in reader:
            identifiers.append(row[0])
            data.extend(float(x) if x not in {'NA', ''} else np.nan for x in row[1:])
    return identifiers, header, np.frombuffer(data, dtype=np.float64).reshape(len(identifiers), len(header))


def padjust(p, by=False):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind='stable')
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    if by:
        adjusted *= np.sum(1 / np.arange(1, len(p) + 1))
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1)
    return result


def decimal_ols_one(y, x, contrast):
    """50-digit normal-equation solution for audit-discrepant rows only."""
    with localcontext() as ctx:
        ctx.prec = 50
        n, k = x.shape
        dy = [Decimal.from_float(float(v)) for v in y]
        normal = [[Decimal(int(v)) for v in row] for row in x.T @ x]
        rhs = [sum((dy[i] for i in range(n) if x[i, j]), Decimal(0)) for j in range(k)]
        augmented = [normal[j] + [rhs[j]] for j in range(k)]
        for j in range(k):
            pivot = augmented[j][j]
            augmented[j] = [v / pivot for v in augmented[j]]
            for i in range(k):
                if i != j:
                    multiple = augmented[i][j]
                    augmented[i] = [v - multiple * w for v, w in zip(augmented[i], augmented[j])]
        beta = [row[-1] for row in augmented]
        dc = [Decimal(1) / 3 if v > 0 else -Decimal(1) / 3 if v < 0 else Decimal(0) for v in contrast] if k == 9 else [Decimal(int(v)) for v in contrast]
        effect = sum((a * b for a, b in zip(beta, dc)), Decimal(0))
        residual = [dy[i] - sum((beta[j] for j in range(k) if x[i, j]), Decimal(0)) for i in range(n)]
        s2 = sum((v * v for v in residual), Decimal(0)) / Decimal(n - k)
        return float(effect), float(s2)


freeze_path = PREP / 'primary_clinical_analysis_freeze.json'
track(freeze_path)
freeze = json.loads(freeze_path.read_text())
changed_frozen = []
for item in freeze['files']:
    p = ROOT / item['path']
    track(p)
    if sha(p) != item['sha256'] or p.stat().st_size != item['bytes']:
        changed_frozen.append(item['path'])
record('all59_frozen_files_unchanged', not changed_frozen, frozen_file_count=len(freeze['files']), mismatches=changed_frozen)

# Extract existing RDS fields and installed implementation only; never call a fit or EB estimator.
r_source = r'''
.libPaths(c('runtime/R-library',.libPaths()))
suppressPackageStartupMessages(library(limma))
out <- 'outputs/biomarker_transfer_v13/clinical_analysis'
x <- readRDS(file.path(out,'primary_group_by_phase_interaction_fit.rds'))
write.csv(data.frame(Entrez_ID=rownames(x$coefficients),posterior_variance=x$s2.post,
 df_prior=x$df.prior,ordinary_df=x$df.residual,df_total=x$df.total,F=x$F,P=x$F.p.value),
 file.path(out,'independent127_existing_interaction_RDS_fields.csv'),row.names=FALSE)
capture.output({print(packageVersion('limma'));print(limma::eBayes);print(limma:::classifyTestsF)},
 file=file.path(out,'independent127_installed_limma_F_implementation.txt'))
'''
env = os.environ.copy()
env['TMPDIR'] = 'runtime/tmp'
subprocess.run(['/usr/local/bin/Rscript', '-'], input=r_source, text=True, cwd=ROOT, env=env, check=True, capture_output=True)
for p in OUT.glob('*_model.rds'):
    track(p)
track(OUT / 'primary_group_by_phase_interaction_fit.rds')

metadata = rows(PREP / 'GSE51981_fixed_clinical_design.csv')
universe = rows(PREP / 'U_all_identifiable_source_admissible_genes.csv')
uid = [r['Entrez_ID'] for r in universe]
master = rows(OUT / 'all_fixed_U157_primary_and_sensitivity_results.csv')
secondary = rows(OUT / 'all_secondary_3U471_results.csv')
provenance = json.loads((PREP / 'orgHs_annotation_provenance.json').read_text())
db = sqlite3.connect('file:' + str(ROOT / provenance['sqlite_relative_path']) + '?mode=ro', uri=True)
valid_ids = {r[0] for r in db.execute('SELECT genes.gene_id FROM genes JOIN gene_info USING(_id) WHERE symbol IS NOT NULL AND symbol != ""')}
db.close()
hc3_rows, mode_predictions = [], {}
primary_fit = None
for mode, aggregation, all148 in [('primary_mean_144', 'mean', False), ('sensitivity_median_144', 'median', False), ('sensitivity_all148_additive', 'mean', True)]:
    genes, samples, complete = numeric_table(PREP / ('GSE51981_' + aggregation + '_probe_gene_matrix.csv.gz'))
    m = metadata if all148 else [r for r in metadata if r['include_primary'] == 'True']
    positions = [samples.index(r['sample']) for r in m]
    y_all = complete[:, positions]
    eligible = np.asarray([g in valid_ids for g in genes]) & np.isfinite(y_all).all(axis=1) & (np.ptp(y_all, axis=1) > 0)
    selected_genes = [g for g, keep in zip(genes, eligible) if keep]
    y = y_all[eligible]
    if all148:
        names = ['(Intercept)', 'clinical_groupUPP', 'clinical_groupENDO', 'phase_collapsedSecretory', 'phase_collapsedUnknown']
        x = np.asarray([[1, r['clinical_group'] == 'UPP', r['clinical_group'] == 'ENDO', r['phase_collapsed'] == 'Secretory', r['phase_collapsed'] == 'Unknown'] for r in m], dtype=float)
        c = np.asarray([[0, 0, 0], [0, -1, 1], [1, 1, 0], [0, 0, 0], [0, 0, 0]], dtype=float)
    else:
        names = [g + '_' + p for g in ['NUPP', 'UPP', 'ENDO'] for p in ['P', 'ESE', 'MSE']]
        x = np.asarray([[r['clinical_group'] + '_' + r['phase'] == cell for cell in names] for r in m], dtype=float)
        c = np.zeros((9, 3))
        c[:3, 0] = -1 / 3
        c[6:, :2] = 1 / 3
        c[3:6, 1] = -1 / 3
        c[:3, 2] = -1 / 3
        c[3:6, 2] = 1 / 3
    di, dh, dx = numeric_table(OUT / (mode + '_design_matrix.csv'))
    ci, ch, dc = numeric_table(OUT / (mode + '_contrast_matrix.csv'))
    record(mode + '_independent_design_identity', di == [r['sample'] for r in m] and dh == names and ci == names
           and ch == ['ENDO_minus_NUPP', 'ENDO_minus_UPP', 'UPP_minus_NUPP'])
    compare(mode + '_independent_design_values', dx, x, atol=1e-14, rtol=0)
    compare(mode + '_independent_contrast_values', dc, c, atol=1e-14, rtol=0)
    xtx_inverse = np.linalg.inv(x.T @ x)
    beta = np.linalg.lstsq(x, y.T, rcond=None)[0].T
    residual = y - beta @ x.T
    df = len(m) - np.linalg.matrix_rank(x)
    s2 = np.einsum('ij,ij->i', residual, residual) / df
    effect_all = beta @ c
    unscaled = np.sqrt(np.diag(c.T @ xtx_inverse @ c))
    bgids, bh, background = numeric_table(OUT / (mode + '_all_background_gene_primary_contrast.csv'))
    record(mode + '_background_eligibility', selected_genes == bgids, independent_count=len(selected_genes), exported_count=len(bgids))
    if selected_genes != bgids:
        raise ValueError('Background disagreement requires adjudication before arithmetic: ' + mode)
    bg = {name: background[:, i] for i, name in enumerate(bh)}
    effect = effect_all[:, 0]
    se = np.sqrt(bg['posterior_variance']) * unscaled[0]
    ordinary_se = np.sqrt(s2) * unscaled[0]
    ordinary_t = effect / ordinary_se
    moderated_t = effect / se
    p = 2 * stats.t.sf(np.abs(moderated_t), bg['df'])
    eqp = np.maximum(stats.t.sf((effect + .5) / se, bg['df']), stats.t.cdf((effect - .5) / se, bg['df']))
    independent = {'effect': effect, 'ordinary_residual_variance': s2, 'ordinary_df': np.repeat(df, len(y)),
                   'ordinary_SE': ordinary_se, 'ordinary_P': 2 * stats.t.sf(np.abs(ordinary_t), df),
                   'ordinary_lower95': effect - stats.t.ppf(.975, df) * ordinary_se,
                   'ordinary_upper95': effect + stats.t.ppf(.975, df) * ordinary_se,
                   'unscaled_SE': np.repeat(unscaled[0], len(y)), 'SE': se, 't': moderated_t, 'P': p,
                   'lower95': effect - stats.t.ppf(.975, bg['df']) * se,
                   'upper95': effect + stats.t.ppf(.975, bg['df']) * se, 'equivalence_P': eqp}
    for field, expected in independent.items():
        compare(mode + '_all_background_' + field, bg[field], expected)
    compare(mode + '_reported_ordinary_P_arithmetic', bg['ordinary_P'],
            2 * stats.t.sf(np.abs(bg['effect'] / bg['ordinary_SE']), bg['ordinary_df']))
    ordinary_bad = ~np.isclose(bg['ordinary_P'], independent['ordinary_P'], atol=2e-10, rtol=2e-9)
    for i in np.flatnonzero(ordinary_bad):
        deffect, ds2 = decimal_ols_one(y[i], x, c[:, 0])
        dse = np.sqrt(ds2) * unscaled[0]
        high_precision_rows.append({'mode': mode, 'Entrez_ID': bgids[i], 'in_U': bgids[i] in uid,
                                    'source_effect': float(bg['effect'][i]), 'numpy_effect': float(effect[i]), 'decimal50_effect': deffect,
                                    'source_variance': float(bg['ordinary_residual_variance'][i]), 'numpy_variance': float(s2[i]), 'decimal50_variance': ds2,
                                    'source_ordinary_P': float(bg['ordinary_P'][i]), 'numpy_ordinary_P': float(independent['ordinary_P'][i]),
                                    'decimal50_effect_variance_Scipy_P': float(2 * stats.t.sf(abs(deffect / dse), df)),
                                    'source_minus_decimal_effect': float(bg['effect'][i] - deffect),
                                    'numpy_minus_decimal_effect': float(effect[i] - deffect), 'sample_range': float(np.ptp(y[i]))})
    if ordinary_bad.any():
        check = next(z for z in checks if z['check'] == mode + '_all_background_ordinary_P')
        check.update(severity='background_OLS_floating_point_warning',
                     mismatched_IDs=[bgids[i] for i in np.flatnonzero(ordinary_bad)],
                     any_mismatch_in_U=any(bgids[i] in uid for i in np.flatnonzero(ordinary_bad)),
                     explanation='Strict comparison retained; near-flat background rows amplify machine-scale coefficient cancellation. Source arithmetic and every U ordinary P checked separately. No tolerance widened.')
    by_gene = {g: i for i, g in enumerate(bgids)}
    target = [next(r for r in master if r['mode'] == mode and r['Entrez_ID'] == gene) for gene in uid]
    for field in [k for k in independent if k not in {'P', 'ordinary_P', 'equivalence_P'}] + ['posterior_variance', 'df']:
        expected_field = independent[field] if field in independent else bg[field]
        expected_target = [expected_field[by_gene[g]] if g in by_gene else np.nan for g in uid]
        actual_target = [float(r[field]) if r[field] not in {'NA', ''} else np.nan for r in target]
        compare(mode + '_U157_' + field, actual_target, expected_target)
    target_p = np.asarray([p[by_gene[g]] if g in by_gene else 1 for g in uid])
    target_eq = np.asarray([eqp[by_gene[g]] if g in by_gene else 1 for g in uid])
    target_op = np.asarray([independent['ordinary_P'][by_gene[g]] if g in by_gene else 1 for g in uid])
    for pfield, values, bhfield, byfield in [('P', target_p, 'BH_U', 'BY_U'), ('equivalence_P', target_eq, 'equivalence_BH_U', 'equivalence_BY_U'), ('ordinary_P', target_op, 'ordinary_BH_U', 'ordinary_BY_U')]:
        compare(mode + '_U157_' + pfield, [float(r[pfield]) for r in target], values)
        compare(mode + '_U157_' + bhfield, [float(r[bhfield]) for r in target], padjust(values))
        compare(mode + '_U157_' + byfield, [float(r[byfield]) for r in target], padjust(values, True))
    missing = [g for g in uid if g not in by_gene]
    record(mode + '_fixed157_missing_and_estimability', len(target) == 157 and all((r['estimable'] == 'TRUE') == (r['Entrez_ID'] in by_gene) for r in target), missing_IDs=missing)
    leverage = np.einsum('ij,jk,ik->i', x, xtx_inverse, x)
    contrast_weights = x @ xtx_inverse @ c[:, 0]
    hc3_p, hc3_current, covariance_error = [], [], 0.0
    for info in universe:
        gene = info['Entrez_ID']
        row = {**info, 'mode': mode, 'N': len(m), 'df_residual': int(df), 'estimable': gene in by_gene,
               'effect': '', 'HC3_SE': '', 'HC3_t': '', 'HC3_lower95': '', 'HC3_upper95': '', 'P': 1.0,
               'missing_reason': '' if gene in by_gene else 'Absent unique-current-Entrez measurement or constant/nonfinite in fixed samples',
               'inference': 'HC3 sandwich; Student-t residual df approximation; same frozen contrast and U157 family'}
        if gene in by_gene:
            i = by_gene[gene]
            corrected = residual[i] / (1 - leverage)
            var = np.sum((contrast_weights * corrected) ** 2)
            full_cov = xtx_inverse @ (x.T @ (corrected[:, None] ** 2 * x)) @ xtx_inverse
            covariance_error = max(covariance_error, abs(var - c[:, 0] @ full_cov @ c[:, 0]))
            hcse = np.sqrt(var)
            hct = effect[i] / hcse
            hcp = float(2 * stats.t.sf(abs(hct), df))
            row.update(effect=float(effect[i]), HC3_SE=float(hcse), HC3_t=float(hct), P=hcp,
                       HC3_lower95=float(effect[i] - stats.t.ppf(.975, df) * hcse),
                       HC3_upper95=float(effect[i] + stats.t.ppf(.975, df) * hcse))
        hc3_p.append(row['P'])
        hc3_current.append(row)
    for row, q, qby in zip(hc3_current, padjust(hc3_p), padjust(hc3_p, True)):
        row['BH_U'], row['BY_U'] = float(q), float(qby)
    record(mode + '_HC3_reduced_and_full_sandwich_identity', covariance_error < 1e-10,
           max_variance_difference=covariance_error, family_size=157, max_leverage=float(leverage.max()))
    hc3_rows.extend(hc3_current)
    summaries[mode] = {'N': len(m), 'design_rank': int(np.linalg.matrix_rank(x)), 'ordinary_df': int(df),
                       'background_genes': len(y), 'U157_estimable': 157 - len(missing),
                       'BH_U_hits': int(np.sum(padjust(target_p) < .05)), 'BY_U_hits': int(np.sum(padjust(target_p, True) < .05)),
                       'TOST_BH_U_hits': int(np.sum(padjust(target_eq) < .05)), 'TOST_BY_U_hits': int(np.sum(padjust(target_eq, True) < .05)),
                       'HC3_BH_U_hits': int(np.sum(padjust(hc3_p) < .05)), 'HC3_BY_U_hits': int(np.sum(padjust(hc3_p, True) < .05)),
                       'missing_IDs': missing}
    mode_predictions[mode] = {r['Entrez_ID']: r for r in target}
    if mode == 'primary_mean_144':
        primary_fit = (x, xtx_inverse, beta, bg, by_gene, effect_all, unscaled)
    del complete, y_all, y, residual
    print(mode + ' verification complete', flush=True)

x, xtx_inverse, beta, bg, by_gene, effect_all, unscaled = primary_fit
j = np.zeros((9, 4))
for k, (group, phase) in enumerate([(2, 1), (2, 2), (1, 1), (1, 2)]):
    j[group * 3 + phase, k] = 1
    j[group * 3, k] = -1
    j[phase, k] = -1
    j[0, k] = 1
ji, jh, actual_j = numeric_table(OUT / 'primary_group_by_phase_interaction_matrix.csv')
compare('independent_F4_interaction_contrast_matrix', actual_j, j, atol=1e-14, rtol=0)
g = j.T @ xtx_inverse @ j
d = beta @ j
f = np.einsum('ij,jk,ik->i', d, np.linalg.inv(g), d) / (4 * bg['posterior_variance'])
rid, rh, rdata = numeric_table(OUT / 'independent127_existing_interaction_RDS_fields.csv')
rd = {field: rdata[:, i] for i, field in enumerate(rh)}
record('interaction_background_ID_order', rid == list(by_gene))
compare('interaction_posterior_matches_main', rd['posterior_variance'], bg['posterior_variance'])
compare('interaction_df_total_matches_main', rd['df_total'], bg['df'])
compare('independent_all_background_F4_quadratic_form', rd['F'], f)
f_df2 = rd['ordinary_df'] + rd['df_prior']
f_p = stats.f.sf(f, 4, f_df2)
compare('independent_all_background_F4_P', rd['P'], f_p)
record('installed_limma_F_df_cap_nuance', True,
       F_uses='ordinary_df+df_prior, from classifyTestsF, not capped df.total',
       main_t_uses='min(ordinary_df+df_prior, pooled residual df)',
       genes_where_df_differs=int(np.sum(np.abs(f_df2 - bg['df']) > 1e-9)),
       F_df2_min=float(f_df2.min()), F_df2_max=float(f_df2.max()))
sec_p, sec_statistics, sec_effects = [], [], []
for row in secondary:
    gene, test = row['Entrez_ID'], row['test']
    if gene not in by_gene:
        p, statistic, effect = 1.0, np.nan, np.nan
    else:
        i = by_gene[gene]
        if test == 'group_by_phase_F4':
            p, statistic, effect = f_p[i], f[i], np.nan
        else:
            k = 1 if test == 'ENDO_minus_UPP' else 2
            effect = effect_all[i, k]
            statistic = effect / (np.sqrt(bg['posterior_variance'][i]) * unscaled[k])
            p = 2 * stats.t.sf(abs(statistic), bg['df'][i])
    sec_p.append(p)
    sec_statistics.append(statistic)
    sec_effects.append(effect)
record('secondary_fixed_3U_scaffold', len(secondary) == 471 and Counter(r['test'] for r in secondary) == {'ENDO_minus_UPP': 157, 'UPP_minus_NUPP': 157, 'group_by_phase_F4': 157})
for field, expected in [('P', sec_p), ('statistic', sec_statistics), ('effect', sec_effects), ('BH_3U', padjust(sec_p)), ('BY_3U', padjust(sec_p, True))]:
    actual = [float(r[field]) if r[field] not in {'', 'NA'} else np.nan for r in secondary]
    compare('secondary471_' + field, actual, expected)
summaries['secondary_BH_471_family'] = {test: sum(float(r['BH_3U']) < .05 for r in secondary if r['test'] == test) for test in ['ENDO_minus_UPP', 'UPP_minus_NUPP', 'group_by_phase_F4']}
qc_rows = [r for r in master if r['mode'] == 'global_QC_sensitivity_identical_primary']
primary_rows = mode_predictions['primary_mean_144']
record('global_QC_sensitivity_is_exact_primary_copy', len(qc_rows) == 157
       and all(all(v == primary_rows[r['Entrez_ID']][k] for k, v in r.items() if k != 'mode') for r in qc_rows))
record('all_frozen_and_audited_sources_unchanged_at_end', all(sha(ROOT / name) == value for name, value in source_hashes.items()))
published_summary = json.loads((OUT / 'primary_51981_summary.json').read_text())['modes']
for mode in ['primary_mean_144', 'sensitivity_median_144', 'sensitivity_all148_additive']:
    expected_summary = published_summary[mode]
    our_summary = summaries[mode]
    record(mode + '_reported_counts_independently_confirmed',
           all(our_summary[k] == expected_summary[k] for k in ['N', 'design_rank', 'background_genes', 'BH_U_hits', 'BY_U_hits'])
           and our_summary['U157_estimable'] == expected_summary['estimable']
           and our_summary['TOST_BH_U_hits'] == expected_summary['bounded_mean_BH_U_hits']
           and our_summary['TOST_BY_U_hits'] == expected_summary['bounded_mean_BY_U_hits'])

hcpath = OUT / 'independent127_HC3_all_U157_three_modes.csv'
with hcpath.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(hc3_rows[0]))
    writer.writeheader()
    writer.writerows(hc3_rows)
precision_path = OUT / 'independent127_background_OLS_high_precision_diagnostics.csv'
if high_precision_rows:
    with precision_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(high_precision_rows[0]))
        writer.writeheader()
        writer.writerows(high_precision_rows)
failed = [c for c in checks if not c['passed']]
warnings_only = all(c.get('severity') == 'background_OLS_floating_point_warning' and not c['any_mismatch_in_U'] for c in failed)
status = 'PASS' if not failed else 'PASS_WITH_BACKGROUND_OLS_NUMERICAL_WARNINGS' if warnings_only else 'FAIL'
report = {'verified_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'result': status,
          'checks': checks, 'summaries': summaries, 'source_hashes': source_hashes,
          'script_sha256': sha(Path(__file__)), 'freeze_sha256': sha(freeze_path),
          'HC3_rows': len(hc3_rows), 'HC3_sha256': sha(hcpath),
          'strict_background_warning_count': len(failed) if warnings_only else 0,
          'high_precision_audit_rows': len(high_precision_rows),
          'second_empirical_Bayes_fit_run': False,
          'scope': 'Independent OLS/arithmetic verification and fixed157 HC3 sensitivity; no selection by effects, no diagnostic classifier/AUC claim.',
          'limits': ['Posterior variance/prior degrees of freedom are checked for arithmetic consistency, not independently re-estimated.',
                     'HC3 Student-t residual-df inference is a finite-sample approximation, not an exact heteroskedastic test.',
                     'The three sensitivity modes retain the same U157 family, with unavailable measurements assigned P1.',
                     'Independent numerical agreement is not independent biological validation.']}
(OUT / 'independent127_GSE51981_statistical_verification.json').write_text(json.dumps(report, indent=2) + '\n')
lines = ['# GSE51981冻结模型：独立数值核验与HC3敏感性', '', '核验状态：' + report['result'] + '。不改变冻结脚本、输入、阈值或基因家族。没有第二次经验贝叶斯拟合。', '',
         '从固定元数据独立构建9单元设计和148样本加性设计，以NumPy最小二乘恢复全部背景基因效应、残差方差、普通SE/df/区间/P；使用已导出的后验方差/df复算moderated t/P/区间和±0.5尺度界值TOST，并核对U157全部BH/BY。', '',
         '|模式|N|背景基因|可估计U|moderated BH/BY|TOST BH/BY|HC3 BH/BY|', '|---|---:|---:|---:|---:|---:|---:|']
for mode in ['primary_mean_144', 'sensitivity_median_144', 'sensitivity_all148_additive']:
    s = summaries[mode]
    lines.append('|{}|{}|{}|{}|{}/{}|{}/{}|{}/{}|'.format(mode, s['N'], s['background_genes'], s['U157_estimable'], s['BH_U_hits'], s['BY_U_hits'], s['TOST_BH_U_hits'], s['TOST_BY_U_hits'], s['HC3_BH_U_hits'], s['HC3_BY_U_hits']))
lines += ['', '交互F按独立J与(X′X)⁻¹的4维二次型及同一后验方差计算；显式核对已安装limma的F自由度规则（残差df+先验df）与t的截顶df.total是否不同。次要ENDO-UPP、UPP-NUPP、交互三类共同在471项家族校正，结果逐项比对；全局QC敏感性157行与主要结果逐项一致。', '',
          'HC3对三种模式的全157基因保存效应、sandwich SE、t、残差df、95%区间、P及BH/BY。用残差/(1-h)构造HC3，直接完整sandwich矩阵与对比权重约简公式独立交叉核验。采用Student-t残差自由度近似；它并非异方差下精确有限样本检验。缺测行保留P1，家族不缩减。', '',
          '本核验不重新估计经验贝叶斯超参数，也不把基因均值关联或HC3敏感性解释为分类器AUC、诊断效能或临床无效。JSON含逐项误差、固定输入和冻结文件SHA；HC3 CSV保存全部471行。', '',
          '保留严格背景普通OLS数值警示：144模式各4行、148模式3行（均不在U）近乎恒定，残差方差约2.9e-20至3.2e-12，1e-15级系数消去误差使普通P最大相差3.03e-5。没有放宽原核验容差；原严格项仍标未通过并单列warning。50位Decimal诊断逐行保存系数/方差/P，证实两种浮点实现对同一极小差值的舍入敏感性。原导出effect/SE/df的普通P算术全部复现；全部U普通P及校正、全部背景moderated算术均通过。这不影响已固定主要/次要家族的计数。']
failures = [r['check'] for r in checks if not r['passed']]
if failures:
    lines += ['', '未通过项：' + '；'.join(failures)]
(OUT / 'independent127_GSE51981_statistical_verification.md').write_text('\n'.join(lines) + '\n')
print(json.dumps({'result': report['result'], 'failed_checks': failures, 'checks': len(checks), 'summaries': summaries}), flush=True)
