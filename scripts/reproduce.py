"""Run the fixed processed-input analyses, comparisons, HC3 checks and figures."""
from pathlib import Path
import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--rscript', default=os.environ.get('RSCRIPT', 'Rscript'))
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
(root / 'runtime/tmp').mkdir(parents=True, exist_ok=True)
env = {**os.environ, 'TMPDIR': 'runtime/tmp', 'RSCRIPT': args.rscript}
steps = [
    [args.rscript, '--vanilla', 'scripts/run_all.R'],
    [sys.executable, 'scripts/compare_results.py'],
    [args.rscript, '--vanilla', 'scripts/run_hc3.R'],
    [sys.executable, 'tests/check_hc3.py'],
    [args.rscript, '--vanilla', 'scripts/run_figures.R'],
]
for step in steps:
    print('Running:', ' '.join(step), flush=True)
    subprocess.run(step, cwd=root, env=env, check=True)
print('Processed-input models, HC3 checks and four figure redraws completed.')
