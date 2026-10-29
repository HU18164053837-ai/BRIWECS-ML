"""Acquire the immutable outcome source and optional pinned provenance repositories.
Never fits a model or executes original-author preprocessing.
"""
from pathlib import Path
import argparse,subprocess,sys,json,urllib.request,hashlib,shutil
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--with-provenance-repositories',action='store_true');args=ap.parse_args()
subprocess.run([sys.executable,str(R/'tools/acquire_figshare_v1.py')],check=True)
if args.with_provenance_repositories:
 records=json.loads((R/'provenance/Third_Party_Repositories.json').read_text())
 records.append(dict(repository='BRIWECS_Data_Publication_shallow',url=records[0]['url'],commit=records[0]['commit']))
 for row in records:
  dest=R/'raw/github'/row['repository'];dest.parent.mkdir(parents=True,exist_ok=True)
  newly_created=not dest.exists()
  if newly_created:subprocess.run(['git','clone','--no-checkout',row['url'],str(dest)],check=True)
  current=subprocess.run(['git','-C',str(dest),'rev-parse','HEAD'],capture_output=True,text=True)
  if current.returncode==0 and current.stdout.strip()==row['commit']:continue
  status=subprocess.run(['git','-C',str(dest),'status','--porcelain'],capture_output=True,text=True,check=True).stdout
  if not newly_created and status.strip():raise RuntimeError('Refusing to change a nonempty or dirty existing checkout: '+str(dest))
  subprocess.run(['git','-C',str(dest),'checkout','--detach',row['commit']],check=True)
  actual=subprocess.check_output(['git','-C',str(dest),'rev-parse','HEAD'],text=True).strip();assert actual==row['commit']
 print('Pinned provenance repositories materialised; official Figshare remains the outcome source.')
