"""Compare all canonical target rows and background coefficients; Python standard library only."""
from pathlib import Path
import csv,gzip,hashlib,json,math,os,subprocess,datetime
ROOT=Path(__file__).resolve().parents[1]
def load(path):
 with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return list(csv.DictReader(f))
checks=[]
def compare(name,refpath,actualpath,keys,required_rows=None):
 ref,actual=load(refpath),load(actualpath)
 lookup={tuple(r[k] for k in keys):r for r in actual}
 failures=[];nnum=0;largest=0.;maxrel=0.
 if len(actual)!=len(lookup) or len(actual)!=len(ref) or (required_rows is not None and len(ref)!=required_rows):failures.append('row count/uniqueness')
 for r in ref:
  key=tuple(r[k] for k in keys);a=lookup.get(key)
  if a is None:failures.append('missing key:'+str(key));continue
  if set(a)!=set(r):failures.append('column identity');continue
  for field,value in r.items():
   observed=a[field]
   if value in ('','NA'):
    if observed not in ('','NA'):failures.append(str(key)+':'+field+':missingness')
    continue
   if field in keys:
    if observed!=value:failures.append(str(key)+':'+field)
    continue
   try:x,y=float(observed),float(value)
   except ValueError:
    if observed!=value:failures.append(str(key)+':'+field+':metadata')
    continue
   nnum+=1;delta=abs(x-y);largest=max(largest,delta);maxrel=max(maxrel,delta/max(abs(y),1e-300))
   if not math.isclose(x,y,abs_tol=2e-10,rel_tol=2e-9):failures.append(str(key)+':'+field+':numeric')
   if field in ('P','BH_U','BY_U','BH_U2','BY_U2','jointP','joint_BH_U2','joint_BY_U2') and value=='1' and observed!='1':failures.append(str(key)+':'+field+':P1')
 checks.append(dict(table=name,reference_rows=len(ref),actual_rows=len(actual),numeric_cells=nnum,max_absolute_error=largest,max_relative_error=maxrel,passed=not failures,failures=failures[:30]))
targets=[('primary628','reference/targets/primary628.csv.gz','results/primary/all_fixed_U157_primary_and_sensitivity_results.csv',['mode','Entrez_ID'],628),
 ('primary_secondary471','reference/targets/primary_secondary471.csv.gz','results/primary/all_secondary_3U471_results.csv',['test','Entrez_ID'],471),
 ('second_mean714','reference/targets/second_mean714.csv.gz','results/secondary/all_six_modes_714_U2_mean_results.csv',['mode','Entrez_ID'],714),
 ('second_interaction714','reference/targets/second_interaction714.csv.gz','results/secondary/all_six_modes_714_U2_SECONDARY_interaction_results.csv',['mode','Entrez_ID'],714),
 ('joint714','reference/targets/joint714.csv.gz','results/secondary/all_six_modes_714_U2_same_direction_IUT.csv',['mode','Entrez_ID'],714)]
for name,ref,actual,keys,n in targets:compare(name,ROOT/ref,ROOT/actual,keys,n)
qa=ROOT/'qa';qa.mkdir(exist_ok=True);coef=qa/'rerun_coefficients';coef.mkdir(exist_ok=True)
r="""root<-normalizePath(getwd());source('scripts/portable_io.R');suppressPackageStartupMessages(library(limma))
for(mode in c('primary_mean_144','sensitivity_median_144','sensitivity_all148_additive')){
 x<-readRDS(file.path('results/primary',paste0(mode,'_model.rds')))$fit$coefficients
 write.csv(data.frame(Entrez_ID=rownames(x),x,check.names=FALSE),file.path('qa/rerun_coefficients',paste0(mode,'.csv')),row.names=FALSE)
}"""
run=subprocess.run([os.environ.get('RSCRIPT','Rscript'),'--vanilla','-'],input=r,text=True,cwd=ROOT,env={**os.environ,'TMPDIR':'runtime/tmp'},capture_output=True)
if run.returncode:raise RuntimeError(run.stderr)
for mode in ['primary_mean_144','sensitivity_median_144','sensitivity_all148_additive']:
 compare(mode+'_all_background_coefficients',ROOT/'reference/coefficients'/f'{mode}.csv.gz',coef/f'{mode}.csv',['Entrez_ID'])
for mode in ['rawmean34','rawmedian34','rawmean_noSourceQC32','rawmean_FEversionAdjusted34','sourceProcessed_mean36','sourceProcessed_mean_noQC32']:
 compare(mode+'_all_background_coefficients',ROOT/'reference/coefficients'/f'{mode}.csv.gz',ROOT/'results/secondary'/f'{mode}_all_background_fitted_coefficients.csv',['Entrez_ID'])
report=dict(verified_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),result='PASS' if all(c['passed'] for c in checks) else 'FAIL',checks=checks,target_rows=sum(c['reference_rows'] for c in checks[:5]),tolerances=dict(absolute=2e-10,relative=2e-9),source_parser_warning_policy='Preserve previously documented at-most-one-ULP source CSV parsing differences; do not change numerical thresholds on failure.',scope='All fixed target and fitted-coefficient rows; no result-based selection; HC3 reference tables supplied separately, not recomputed by this standard-library comparator.')
(qa/'rerun_numeric_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
if report['result']!='PASS':raise SystemExit(1)
