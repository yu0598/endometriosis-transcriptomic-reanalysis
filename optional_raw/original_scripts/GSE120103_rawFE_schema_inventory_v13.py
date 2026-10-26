"""Inventory raw FE schemas and technical feature identities; no normalization."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

root = Path.cwd()
raw = root / 'data_raw/GSE120103_rawFE_v13'
out = root / 'outputs/biomarker_transfer_v13/phase2_investigation/GSE120103_processing_sources'
headers = out / 'raw_FE_headers'
headers.mkdir(exist_ok=True)
manifest = list(csv.DictReader(open(raw / 'rawFE_retrieval_manifest.csv')))
manifest.sort(key=lambda row: row['gsm'])
numeric_names = ['gMedianSignal', 'gMeanSignal', 'gBGMedianSignal', 'gProcessedSignal']
rows, metadata_rows, anomalies, reference = [], [], [], None

def fingerprint(items):
    data = '\n'.join(json.dumps(x, separators=(',', ':'), ensure_ascii=False) for x in items)
    return hashlib.sha256(data.encode()).hexdigest()

for item in manifest:
    gsm, filename = item['gsm'], item['file']
    preamble, feature_columns, parameter_columns = [], None, None
    section, params, identities = None, {}, []
    nul_bytes, malformed_rows, data_records = 0, 0, 0
    finite = {name: 0 for name in numeric_names}
    not_numeric = {name: 0 for name in numeric_names}
    with gzip.open(raw / filename, 'rt', encoding='utf-8-sig', errors='strict') as handle:
        for lineno, line in enumerate(handle, 1):
            nul_bytes += line.count('\x00')
            cells = line.rstrip('\r\n').split('\t')
            kind, values = cells[0], cells[1:]
            if feature_columns is None:
                preamble.append(line)
                if kind == 'FEPARAMS':
                    section, parameter_columns = kind, values
                elif kind == 'STATS':
                    section = kind
                elif kind == 'FEATURES':
                    section, feature_columns = kind, values
                    positions = {name: i for i, name in enumerate(feature_columns)}
                    assert all(name in positions for name in ['FeatureNum', 'ProbeName', 'ControlType'])
                elif kind == 'DATA' and section == 'FEPARAMS':
                    assert len(values) == len(parameter_columns)
                    params = dict(zip(parameter_columns, values))
                continue
            if kind != 'DATA':
                continue
            data_records += 1
            if len(values) != len(feature_columns) or '\x00' in line:
                malformed_rows += 1
                anomalies.append(dict(gsm=gsm, physical_line=lineno, expected_fields=len(feature_columns),
                                      observed_fields=len(values), nul_bytes=line.count('\x00'),
                                      observed_FeatureNum_prefix=values[positions['FeatureNum']].split('\x00', 1)[0][:80] if len(values)>positions['FeatureNum'] else '',
                                      observed_ProbeName_prefix=values[positions['ProbeName']].split('\x00', 1)[0][:80] if len(values)>positions['ProbeName'] else '',
                                      interpretation='Malformed source record; not repaired or supplied as numeric input'))
                continue
            identity = (int(values[positions['FeatureNum']]), values[positions['ProbeName']],
                        values[positions['ControlType']])
            identities.append(identity)
            for name in numeric_names:
                if name not in positions:
                    continue
                try:
                    finite[name] += math.isfinite(float(values[positions[name]]))
                except ValueError:
                    not_numeric[name] += 1
    assert feature_columns and identities
    header_file = headers / (gsm + '_raw_header.txt')
    header_file.write_text(''.join(preamble))
    for key, value in params.items():
        metadata_rows.append(dict(gsm=gsm, parameter=key, value=value))
    identity_map = {a: (b, c) for a, b, c in identities}
    if reference is None:
        reference = identities
        refmap = identity_map
        with gzip.open(out / 'GSM3393491_feature_identity_reference.tsv.gz', 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['FeatureNum', 'ProbeName', 'ControlType'])
            writer.writerows(reference)
    common = sorted(set(identity_map) & set(refmap))
    proto = {k: v for k, v in params.items() if any(x in k.lower() for x in
             ['version', 'protocol', 'grid', 'design', 'scan', 'date', 'barcode'])}
    controls = {}
    for _, _, control in identities:
        controls[control] = controls.get(control, 0) + 1
    row = dict(gsm=gsm, source_file=filename, source_compressed_sha256=item['compressed_sha256'],
               feature_rows=len(identities), feature_data_records_seen=data_records,
               malformed_feature_records=malformed_rows, nul_bytes=nul_bytes,
               source_table_shape_complete=(malformed_rows == 0 and nul_bytes == 0),
               unique_FeatureNum=len(identity_map),
               duplicate_FeatureNum=len(identities)-len(identity_map),
               unique_ProbeName=len({b for _, b, _ in identities}),
               ControlType_counts=json.dumps(controls, sort_keys=True),
               columns_count=len(feature_columns), columns_json=json.dumps(feature_columns),
               feature_schema_sha256=fingerprint(feature_columns),
               ordered_identity_sha256=fingerprint(identities),
               unordered_identity_sha256=fingerprint(sorted(identities)),
               ordered_identity_equal_reference=identities == reference,
               unordered_identity_equal_reference=sorted(identities) == sorted(reference),
               common_FeatureNum_with_reference=len(common),
               ProbeName_matches_at_common_FeatureNum=sum(identity_map[k][0] == refmap[k][0] for k in common),
               ControlType_matches_at_common_FeatureNum=sum(identity_map[k][1] == refmap[k][1] for k in common),
               full_identity_matches_at_common_FeatureNum=sum(identity_map[k] == refmap[k] for k in common),
               missing_reference_FeatureNum=len(set(refmap)-set(identity_map)),
               extra_relative_to_reference_FeatureNum=len(set(identity_map)-set(refmap)),
               protocol_metadata_json=json.dumps(proto),
               filename_protocol_token='GE1-v5_95_Feb07' if 'GE1-v5_95_Feb07' in filename else 'GE1_107_Sep09',
               header_path=str(header_file.relative_to(root)),
               header_sha256=hashlib.sha256(header_file.read_bytes()).hexdigest(),
               own_normalization=False, biological_group_statistics=False)
    for name in numeric_names:
        row[name+'_available'] = name in positions
        row[name+'_finite_rows'] = finite[name]
        row[name+'_not_numeric_rows'] = not_numeric[name]
    rows.append(row)
    print('Inventoried', gsm, len(identities), 'features', len(feature_columns), 'columns', flush=True)
for filename, content in [('raw_FE_schema_inventory.csv', rows), ('raw_FE_FEPARAMS_inventory.csv', metadata_rows), ('raw_FE_structural_anomalies.csv', anomalies)]:
    with open(out / filename, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(content[0]))
        writer.writeheader()
        writer.writerows(content)
print('COMPLETE', len(rows), 'schemas; no normalization or biological group statistics', flush=True)
