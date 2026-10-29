from common import *
import shutil

def run():
 checks=[]
 def check(name,condition,detail=''):
  checks.append(dict(check=name,status='PASS' if condition else 'FAIL',detail=detail))
  assert condition,name
 # No source or earlier delivery changes are accepted.
 integrity=pd.read_csv(ROOT/'results/tables/input_integrity.csv')
 check('Prior stage immutable files',all(sha(P1/r.path)==r.sha256 for r in integrity.itertuples()),f'{len(integrity)} inherited files')
 raw=pd.read_csv(P1/'audit/source_file_manifest.csv');ok=[]
 for r in raw.itertuples():
  p=P1/r.relative_path
  ok.append(p.exists() and sha(p)==r.sha256)
 check('Phase1 original source SHA manifest',all(ok),f'{len(ok)} source files')
 co=pd.read_parquet(ROOT/'data/stripe_rust_inherited.parquet');check('Frozen cohort unique IDs',len(co)==8539 and co.sample_id.nunique()==8539)
 restored=int((co.published_disease_value.eq(0)&co.disease_severity.gt(0)&co.disease_severity.lt(1)).sum())
 check('Restored positive fractional scores',restored==199,str(restored))
 check('Final true zeros',co.disease_severity.eq(0).sum()==453)
 dec=pd.read_csv(ROOT/'results/tables/01_environment_error_decomposition.csv');check('Exact population variance identity',np.allclose(dec.MSE,dec.bias_squared+dec.error_variance,atol=1e-10))
 delta=pd.read_csv(ROOT/'results/tables/03_disease_delta_error_components.csv');check('Paired MSE decomposition',np.allclose(delta.delta_MSE,delta.delta_bias_squared+delta.delta_centered_variance,atol=1e-10))
 man=pd.read_csv(ROOT/'results/tables/07_matched_test_domain_manifest.csv');check('80 reused model tasks',len(man)==80 and man.reuse_status.eq('EXACT_PILOT_MODEL_AND_OOF_REUSE').all())
 check('Identical natural tests',man.groupby('site_year').test_id_sha256.nunique().eq(1).all() and man.drop_duplicates('site_year').test_n.sum()==5899)
 assoc=pd.read_csv(ROOT/'results/tables/04_adjusted_association_models.csv');check('Association identified and converged',assoc.status.eq('IDENTIFIABLE').all() and assoc.converged.all() and assoc.bootstrap_valid.eq(499).all())
 check('Independent A3 verification',json.loads((ROOT/'reports/association_verification.json').read_text())['independent_dense_FWL_agrees'])
 matrix=pd.read_csv(ROOT/'results/tables/10_claim_evidence_novelty_matrix.csv');check('Seven original studies and three claims',matrix.DOI.nunique()==7 and len(matrix)==21)
 for i in range(1,11):check(f'Required table {i:02}',len(list((ROOT/'results/tables').glob(f'{i:02}_*.csv')))==1)
 for prefix in ['C1','C2','C3','C4']:
  for ext in ['png','pdf']:check(f'{prefix} {ext}',len(list((ROOT/'results/figures').glob(prefix+'*.'+ext)))==1)
 check('Six required scientific reports',len(list((ROOT/'reports').glob('*Findings.md')))==3 and (ROOT/'reports/04_Literature_Novelty_Assessment.md').exists() and len(list((ROOT/'reports').glob('BRIWECS*.md')))==2)
 shutil.copyfile(P/'config/model_config.json',ROOT/'config/inherited_pilot_model_config.json')
 shutil.copyfile(P/'config/cohort_freeze.json',ROOT/'config/inherited_pilot_cohort_freeze.json')
 shutil.copyfile(P/'config/disease_unit_evidence.csv',ROOT/'config/inherited_disease_unit_evidence.csv')
 shutil.copyfile(ROOT/'results/tables/07_matched_test_domain_manifest.csv',ROOT/'logs/model_fit_status.csv')
 csv(checks,'delivery_verification_checks.csv')
 js(dict(status='PASS',checks=len(checks),source_files_unchanged=len(raw),prior_stage_files_unchanged=len(integrity),new_ML_fits=0,actual_reused_model_tasks=80,restored_positive_fractional_scores=restored,true_zero_count=453,scientific_decision='Manuscript CONDITIONAL GO - narrow evaluation-boundary case only',phase2c_started=False,literature_access_limits=['P5 publisher-indexed Methods/Results only; direct HTML/PDF denied','P6 subscription preview, no complete Methods'],figure_visual_inspection='All four PNGs inspected; C2 label shortened after layout check'),'reports/delivery_verification.json')
 excluded={'reports/SHA256_delivery_manifest.csv','reports/output_inventory.csv','logs/run_verify_delivery.log','logs/verify_delivery.log'}
 files=[p for p in ROOT.rglob('*') if p.is_file() and not any(x in p.parts for x in ['__pycache__','matplotlib']) and p.relative_to(ROOT).as_posix() not in excluded]
 entries=[dict(relative_path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
 pd.DataFrame(entries).to_csv(ROOT/'reports/SHA256_delivery_manifest.csv',index=False,encoding='utf-8-sig')
 pd.DataFrame([dict(relative_path=e['relative_path'],bytes=e['bytes'],category=e['relative_path'].split('/')[0]) for e in entries]).to_csv(ROOT/'reports/output_inventory.csv',index=False,encoding='utf-8-sig')
 print(f'PASS: {len(checks)} checks; {len(entries)} files hashed; source manifests unchanged. No new ML training.')
if __name__=='__main__':run()
