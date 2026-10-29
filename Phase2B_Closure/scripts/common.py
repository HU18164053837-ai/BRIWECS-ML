from pathlib import Path
import os,json,hashlib,datetime
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1];P1=ROOT.parent;P=P1/'Phase2B_Pilot';A=P1/'Phase2A'
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'config/matplotlib'));(ROOT/'config/matplotlib').mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ih(x):return hashlib.sha256('\n'.join(sorted(map(str,x))).encode()).hexdigest()
def csv(d,n):
 if not isinstance(d,pd.DataFrame):d=pd.DataFrame(d)
 d.to_csv(ROOT/'results/tables'/n,index=False,encoding='utf-8-sig');return d
def js(d,n):(ROOT/n).write_text(json.dumps(d,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def table(d):
 d=d.fillna('');f=lambda x:f'{x:.4f}' if isinstance(x,float) else str(x).replace('|',' / ').replace('\n',' ')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(map(f,r))+' |' for r in d.itertuples(index=False,name=None))+'\n'
def met(y,p):
 y=np.asarray(y);p=np.asarray(p);e=p-y;return dict(n=len(y),MAE=abs(e).mean(),MSE=(e*e).mean(),RMSE=np.sqrt((e*e).mean()),MBE=e.mean(),R2=1-(e*e).sum()/((y-y.mean())**2).sum() if np.var(y)>1e-12 else np.nan)
def ci(v,w=None):
 v=np.asarray(v);rng=np.random.default_rng(42);i=rng.integers(0,len(v),(2000,len(v)))
 b=v[i].mean(axis=1) if w is None else (v[i]*np.asarray(w)[i]).sum(axis=1)/np.asarray(w)[i].sum(axis=1)
 return np.quantile(b,[.025,.975])
