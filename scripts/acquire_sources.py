"""Optional bounded source acquisition, skips existing data; keeps official archive."""
from pathlib import Path
import datetime as dt, hashlib, json, subprocess, zipfile
ROOT=Path(__file__).resolve().parents[1]
def acquire():
 log=[]
 def download(url,target):
  target.parent.mkdir(parents=True,exist_ok=True)
  if target.exists():log.append(dict(url=url,path=str(target.relative_to(ROOT)),status='EXISTING_NOT_OVERWRITTEN'));return
  q=subprocess.run(['curl.exe','-L','--fail','--retry','2','--connect-timeout','20','--max-time','180','-C','-',url,'-o',str(target)],capture_output=True,text=True)
  log.append(dict(url=url,path=str(target.relative_to(ROOT)),status='COMPLETED' if q.returncode==0 else 'BLOCKED',returncode=q.returncode,message=q.stderr[-2500:],time_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
  if q.returncode:raise RuntimeError('Download failed; see acquisition log: '+url)
 try:
  download('https://api.figshare.com/v2/articles/27910269',ROOT/'sources/figshare_article.json')
  m=json.loads((ROOT/'sources/figshare_article.json').read_text());assert m['version']==1,'Version changed; select explicitly before proceeding'
  for f in m['files']:
   target=ROOT/'raw/figshare'/f['name'];download(f['download_url'],target)
   assert target.stat().st_size==f['size']
   assert hashlib.md5(target.read_bytes()).hexdigest()==f['computed_md5']
   if zipfile.is_zipfile(target) and not (target.parent/'extracted').exists():
    with zipfile.ZipFile(target) as z:
     for name in z.namelist():assert (target.parent/'extracted'/name).resolve().is_relative_to((target.parent/'extracted').resolve())
     z.extractall(target.parent/'extracted')
  for owner,repo in [('tillrose','BRIWECS_Data_Publication'),('Illustratien','Scientific_Data_Analyis')]:
   p=ROOT/'raw/github'/repo
   if p.exists():log.append(dict(repository=repo,status='EXISTING_NOT_MODIFIED'));continue
   q=subprocess.run(['git','-c','http.lowSpeedLimit=1000','-c','http.lowSpeedTime=30','clone','--depth','1',f'https://github.com/{owner}/{repo}.git',str(p)],capture_output=True,text=True,timeout=180)
   log.append(dict(repository=repo,status='COMPLETED' if q.returncode==0 else 'BLOCKED',message=q.stderr))
  download('https://zenodo.org/records/6521175/files/basf-wheat-diseases.txt',ROOT/'sources/basf-wheat-diseases.txt')
 finally:(ROOT/'sources/acquisition_rerun_log.json').write_text(json.dumps(log,ensure_ascii=False,indent=2),encoding='utf-8')
if __name__=='__main__':acquire()
