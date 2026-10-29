from common import *
import sys,platform,importlib.metadata as im
def finalize():
 packages={d.metadata['Name']:d.version for d in im.distributions() if d.metadata['Name']};(ROOT/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in sorted(packages.items()))+'\n',encoding='utf-8');js(dict(python=sys.version,executable=sys.executable,packages=packages,platform=platform.platform(),time_utc=now(),dependency_policy='Existing Phase2A isolated environment reused without package installation or updates'),'config/environment.json')
 files=[]
 for f in ROOT.rglob('*'):
  rel=f.relative_to(ROOT)
  if not f.is_file() or '__pycache__' in rel.parts or 'matplotlib_cache' in rel.parts or f.name=='SHA256_delivery_manifest.csv':continue
  files.append(dict(relative_path=rel.as_posix(),file_size_bytes=f.stat().st_size,sha256=sha(f)))
 csv(files,'reports/SHA256_delivery_manifest.csv')
if __name__=='__main__':finalize()
