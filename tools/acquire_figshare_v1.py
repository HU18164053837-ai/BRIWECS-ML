"""Download only the immutable official archive, validate before safe extraction."""
from pathlib import Path
import hashlib,json,urllib.request,zipfile,shutil
R=Path(__file__).resolve().parents[1];lock=json.loads((R/'provenance/Figshare_v1_Lock.json').read_text())
with urllib.request.urlopen(lock['metadata_url'],timeout=30) as r:m=json.load(r)
assert m['version']==1 and m['doi']==lock['doi']
f=next(f for f in m['files'] if f['id']==lock['file_id']);assert f['computed_md5']==lock['md5'] and f['size']==lock['bytes']
dest=R/'raw/figshare'/lock['filename'];dest.parent.mkdir(parents=True,exist_ok=True)
if not dest.exists():
 tmp=dest.with_suffix('.download')
 with urllib.request.urlopen(lock['download_url'],timeout=120) as response,tmp.open('wb') as out:shutil.copyfileobj(response,out)
 assert tmp.stat().st_size==lock['bytes'] and hashlib.sha256(tmp.read_bytes()).hexdigest()==lock['sha256'];tmp.replace(dest)
assert dest.stat().st_size==lock['bytes'] and hashlib.md5(dest.read_bytes()).hexdigest()==lock['md5'] and hashlib.sha256(dest.read_bytes()).hexdigest()==lock['sha256']
target=(dest.parent/'extracted').resolve();target.mkdir(exist_ok=True)
with zipfile.ZipFile(dest) as z:
 assert z.testzip() is None
 for member in z.infolist():
  p=(target/member.filename).resolve();assert p.is_relative_to(target),'Unsafe ZIP member'
 z.extractall(target)
print('PASS: Figshare version 1 archive hashes, CRC and extraction paths')
