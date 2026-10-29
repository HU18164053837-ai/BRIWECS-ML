import subprocess,sys
from pathlib import Path
r=Path(__file__).resolve().parent
for script in ['prepare_data.py','train_models.py','evaluate_models.py','make_reports.py','verify_results.py','finalize.py']:
 result=subprocess.run([sys.executable,str(r/script)],check=True)
