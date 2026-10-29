from common import *
import urllib.request, concurrent.futures
SOURCES=[
 ('P1','10.1038/s41598-020-77200-0','https://www.nature.com/articles/s41598-020-77200-0'),
 ('P2','10.1038/s41597-024-04332-7','https://www.nature.com/articles/s41597-024-04332-7'),
 ('P3','10.1038/s41467-026-72330-x','https://www.nature.com/articles/s41467-026-72330-x'),
 ('P4','10.1007/s11119-026-10462-2','https://link.springer.com/article/10.1007/s11119-026-10462-2'),
 ('P5','10.3390/rs18162694','https://www.mdpi.com/2072-4292/18/16/2694'),
 ('P6','10.1038/s41477-023-01516-8','https://www.nature.com/articles/s41477-023-01516-8'),
 ('P7','10.1002/tpg2.20532','https://pmc.ncbi.nlm.nih.gov/articles/PMC11726420/')]
def acquire(s):
 ident,doi,url=s;out=[];dest=ROOT/'data/literature';dest.mkdir(exist_ok=True)
 for kind,u in [('crossref','https://api.crossref.org/works/'+doi),('publisher',url)]:
  r=dict(id=ident,doi=doi,kind=kind,url=u,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
  try:
   req=urllib.request.Request(u,headers={'User-Agent':'BRIWECS-ML scientific literature audit/1.0'})
   with urllib.request.urlopen(req,timeout=35) as f:b=f.read();r.update(status=f.status,final_url=f.url)
   p=dest/(ident+('_crossref.json' if kind=='crossref' else '.html'));p.write_bytes(b);r.update(file=str(p.relative_to(ROOT)),sha256=sha(p),bytes=len(b))
  except Exception as e:r.update(status='UNAVAILABLE',error=str(e))
  out.append(r)
 return out
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=[r for batch in ex.map(acquire,SOURCES) for r in batch]
 csv(rows,'literature_source_receipt.csv');print(pd.DataFrame(rows)[['id','kind','status']].to_string(index=False))
