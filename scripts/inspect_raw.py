from pathlib import Path
import sys
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd
r=Path(__file__).resolve().parents[1]/'raw/figshare/extracted/Briwecs_data'
frames=[]
for f in sorted((r/'data/locations').glob('*.csv')):
 try:d=pd.read_csv(f,sep=';',dtype=str)
 except UnicodeDecodeError:d=pd.read_csv(f,sep=';',dtype=str,encoding='cp1252')
 if len(d.columns)==1:d=pd.read_csv(f,dtype=str)
 d['source_file']=f.name;d['source_row']=range(2,len(d)+2);frames.append(d)
d=pd.concat(frames,ignore_index=True)
print('raw',d.shape, d.columns.tolist())
print(d[['Location','Treatment']].drop_duplicates().to_string(index=False))
for f in sorted((r/'data/management').glob('*.xlsx')):
 x=pd.read_excel(f,header=None).fillna('')
 for i,row in x.iterrows():
  text=' | '.join(str(v) for v in row if str(v))
  if any(t in text.lower() for t in ['infection','inocul','harvest','sowing','irrigat','rust','fusarium']):print(f.name,i+1,text)
