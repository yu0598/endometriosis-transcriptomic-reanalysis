"""Compare every fixed HC3 target against the archived independent implementation."""
from pathlib import Path
import csv
import gzip
import json
import math

ROOT = Path(__file__).resolve().parents[1]
def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as stream:
        return list(csv.DictReader(stream))

checks = []
for name, filename, keys, n in [
    ('primary', 'primary_HC3_471.csv', ['mode', 'Entrez_ID'], 471),
    ('secondary', 'second_HC3_1428.csv', ['mode', 'contrast', 'Entrez_ID'], 1428),
]:
    ref = read(ROOT / 'reference/HC3' / (filename + '.gz'))
    actual = read(ROOT / 'results/hc3' / filename)
    lookup = {tuple(r[k] for k in keys): r for r in actual}
    assert len(ref) == len(actual) == len(lookup) == n
    columns = ['N', 'effect', 'HC3_SE', 'HC3_t', 'HC3_lower95', 'HC3_upper95', 'P']
    columns += ['df_residual', 'BH_U', 'BY_U'] if name == 'primary' else ['ordinary_df', 'BH_U2', 'BY_U2']
    largest = 0.0
    values = 0
    for expected in ref:
        observed = lookup[tuple(expected[k] for k in keys)]
        assert expected['estimable'].lower() == observed['estimable'].lower()
        assert expected['symbol'] == observed['symbol']
        if observed['estimable'].lower() == 'false':
            assert float(observed['P']) == 1.0
            assert all(observed[col] in ('', 'NA') for col in
                       ['effect', 'HC3_SE', 'HC3_t', 'HC3_lower95', 'HC3_upper95'])
        for col in columns:
            a, b = expected[col], observed[col]
            if a in ('', 'NA'):
                assert b in ('', 'NA'), (name, expected['Entrez_ID'], col)
            else:
                a, b = float(a), float(b)
                largest = max(largest, abs(a-b)); values += 1
                assert math.isclose(a, b, abs_tol=2e-10, rel_tol=2e-9), (name, expected['mode'], expected['Entrez_ID'], col, a, b)
        assert (float(expected['P']) < .05) == (float(observed['P']) < .05)
    checks.append({'cohort': name, 'rows': n, 'numeric_values': values, 'max_absolute_error': largest, 'passed': True})

report = {'result': 'PASS', 'checks': checks, 'absolute_tolerance': 2e-10, 'relative_tolerance': 2e-9}
(ROOT / 'qa').mkdir(exist_ok=True)
(ROOT / 'qa/public_hc3_comparison.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
