from pathlib import Path
import pandas as pd
import sys
sys.stdout.reconfigure(encoding='utf-8')
r=Path(__file__).resolve().parents[2];base=r/'raw/figshare/extracted/Briwecs_data'
d=pd.read_csv(base/'output/BRIWECS_data_publication.csv',sep=';');d=d[d.Year.between(2015,2017)&d.Location.isin(['GGE','HAN','KAL','KIE','QLB'])]
print(d.groupby(['Location','Year']).agg(rows=('Seedyield','size'),yield_n=('Seedyield','count'),sowing_min=('Sowing_date','min'),sowing_max=('Sowing_date','max'),heading_min=('BBCH59','min'),heading_n=('BBCH59','count'),maturity_min=('BBCH87','min'),maturity_n=('BBCH87','count')).to_string())
for loc in ['GGE','HAN','KAL','KIE','QLB','DKI']:
 for year in [2015,2016,2017]:
  f=base/f'data/management/Management_information_{loc}_{year}.xlsx'
  if not f.exists():continue
  x=pd.read_excel(f,header=None).fillna('');print('\n',f.name)
  for i,row in x.iterrows():
   vals=[str(v) for v in row];text=' | '.join(v for v in vals if v)
   if any(q in text.lower() for q in ['harvest','sowing date','sowing:','ernte','irrigation','bewässer']):print(i+1,text)
print(pd.read_csv(base/'metadata/BRIWECs_cultivar_info.csv').describe(include='all').to_string())
