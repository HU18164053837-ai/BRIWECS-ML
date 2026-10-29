from common import *
import sys,platform,importlib.metadata as im
def finalize():
 packages=sorted([(p.metadata['Name'],p.version) for p in im.distributions() if p.metadata['Name']])
 (ROOT/'requirements.txt').write_text('\n'.join(f'{n}=={v}' for n,v in packages)+'\n',encoding='utf-8')
 write_json(dict(python=sys.version,executable=sys.executable,platform=platform.platform(),packages=dict(packages),time_utc=now(),model_threads=4,blas_threads=2),'config/environment.json')
 files=[]
 for f in ROOT.rglob('*'):
  rel=f.relative_to(ROOT)
  if not f.is_file() or any(x in ['.venv','.pip_cache','__pycache__'] for x in rel.parts) or f.name=='output_manifest.csv':continue
  files.append(dict(relative_path=rel.as_posix(),file_size_bytes=f.stat().st_size,sha256=sha(f)))
 csv(files,'reports/output_manifest.csv')
if __name__=='__main__':finalize()
