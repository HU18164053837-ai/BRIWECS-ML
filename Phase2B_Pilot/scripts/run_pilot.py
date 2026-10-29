from pathlib import Path
import subprocess,sys
r=Path(__file__).resolve().parent
for name in ['prepare.py','measurement_evidence.py','train.py','evaluate.py','associations.py','posthoc_diagnostics.py','figures.py','make_reports.py','verify.py','finalize.py']:
 subprocess.run([sys.executable,str(r/name)],check=True)
