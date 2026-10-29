from pathlib import Path
import os,subprocess,sys
root=Path(__file__).resolve().parents[1]
os.chdir(root)
os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('MKL_NUM_THREADS','2')
os.environ.setdefault('PYTHONIOENCODING','utf-8')
for script in ['prediction_analysis.py','adjusted_association.py','closure_delivery.py','verify_delivery.py']:
 with (root/'logs'/('run_'+script.replace('.py','.log'))).open('w',encoding='utf-8') as log:
  result=subprocess.run([sys.executable,str(root/'scripts'/script)],stdout=log,stderr=subprocess.STDOUT)
 if result.returncode:raise SystemExit(f'{script} failed; inspect logs')
 print(script,'PASS',flush=True)
