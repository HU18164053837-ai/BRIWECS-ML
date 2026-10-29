from common import *
import re,warnings
warnings.filterwarnings('ignore',category=UserWarning,module='openpyxl')
def measurement_evidence():
 units=pd.read_excel(RAW/'metadata/Unit.xlsx');csv(units[units.trait.isin(['Stripe_rust','Septoria'])],'config/disease_unit_evidence.csv')
 rows=[]
 for loc in ['GGE','HAN','KAL','KIE','QLB','DKI']:
  for year in [2015,2016,2017]:
   f=RAW/f'data/management/Management_information_{loc}_{year}.xlsx'
   if not f.exists():continue
   x=pd.read_excel(f,header=None).fillna('');header=None;treatment=''
   for i in range(len(x)):
    vals=x.iloc[i].tolist()
    if any(str(v)=='Infection_data' for v in vals):header=[str(v) for v in vals];continue
    m=re.search(r'\(([^)]+)\)',str(vals[0]))
    if m and re.search('HN|LN',m.group(1)):treatment=m.group(1)
    if header is None:continue
    for j,value in enumerate(vals):
     if value not in ['Stripe_rust','Septoria'] or header[j]!='Infection_data':continue
     evidence=[]
     for k in range(j+1,len(vals)):
      label=header[k]
      if 'Infection_data' in label or 'On-plot' in label or 'Irrigation' in label:break
      if any(word in label.lower() for word in ['date','measure']):evidence.append(dict(header=label,raw_value=str(vals[k]),column=k+1))
     rows.append(dict(file=f.name,excel_row=i+1,disease=str(value),treatment_raw=treatment,date_fields_json=json.dumps(evidence,ensure_ascii=False),has_inoculation_header=any('inoculation' in e['header'].lower() for e in evidence),nonempty_measure_date_fields=sum(bool(e['raw_value']) and ('measure' in e['header'].lower()) for e in evidence),interpretation='Treatment-calendar evidence only; not plot-specific endpoint; never choose max/most correlated date'))
 csv(rows,'data/measurement_calendar_evidence.csv')
 print(f'{len(rows)} source calendar entries; inoculation and measure-date columns explicitly separated; no endpoint re-selection.')
if __name__=='__main__':measurement_evidence()
