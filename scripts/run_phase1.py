"""Read-only Phase 1 audit. No model fitting, imputation, or raw-file mutation."""
from pathlib import Path
import argparse, datetime as dt, hashlib, json, platform, re, subprocess, sys, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')
ROOT=Path(__file__).resolve().parents[1]
DS=['Stripe_rust','Leaf_rust','Powdery_mildew','Septoria','DTR','Fusarium']
K=['Year','Location','Treatment','Block','Row','Column','BRISONr']
E=['Location','Year','Treatment']
NOW=dt.datetime.now(dt.timezone.utc).isoformat()
def save(df,name):
 if not isinstance(df,pd.DataFrame):df=pd.DataFrame(df)
 df.to_csv(ROOT/'audit'/name,index=False,encoding='utf-8-sig')
 return df
def readcsv(f):
 for enc in ['utf-8-sig','cp1252']:
  try:
   sep=';' if ';' in f.open(encoding=enc).readline() else ','
   return pd.read_csv(f,sep=sep,encoding=enc,dtype=str)
  except UnicodeDecodeError:continue
 raise ValueError(f)
def normalize(x):
 x=x.copy()
 for c in ['Year','Row','Column']:x[c]=pd.to_numeric(x[c],errors='coerce')
 x['BRISONr']=x.BRISONr.str.replace('BRISONR','BRISONr',regex=False)
 x['Block']=x.Block.replace({'1':'B1','2':'B2'})
 t=x.Treatment.str.replace(r'D{1,2}','D',regex=True)
 irrig=x.Location.eq('GGE')&x.Year.isin([2015,2018,2019])&~t.str.contains('D',na=False)
 t=t.where(~irrig,t+'_IR');t=t.where(~x.Location.eq('DKI'),t+'_RO')
 t=t.str.replace('_D','',regex=False);t=t.where(t.str.contains(r'(?:IR|RO)$',na=False),t+'_RF')
 x['Treatment']=t;x['Location']=x.Location.replace({'DKI':'KIE','KAD':'KAL'})
 return x
def flag_inoc(loc,year,t,trait):
 if loc=='QLB' and t in ['HN_NF_RF','LN_NF_RF'] and trait in ['Stripe_rust','Leaf_rust','Fusarium']:
  return 'ARTIFICIAL' if year in [2015,2016,2017] else 'MIXED_OR_UNCERTAIN'
 return 'NATURAL' # paper explicitly identifies the exception, disease-specific
def main():
 global ROOT
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);args=ap.parse_args();ROOT=args.root.resolve()
 for n in ['audit','reports','logs','sources','environment']: (ROOT/n).mkdir(exist_ok=True)
 r=ROOT/'raw/figshare/extracted/Briwecs_data'
 d=pd.read_csv(r/'output/BRIWECS_data_publication.csv',sep=';');d['published_row']=np.arange(len(d))+2
 units=pd.read_excel(r/'metadata/Unit.xlsx').set_index('trait')
 raw=[];parse=[]
 for f in sorted((r/'data/locations').glob('*.csv')):
  x=readcsv(f);parse.append(dict(file=f.name,rows=len(x),columns=len(x.columns),delimiter=';' if ';' in f.read_text(encoding='cp1252').splitlines()[0] else ','))
  x['source_file']=f.name;x['source_row']=np.arange(len(x))+2;raw.append(x)
 raw=pd.concat(raw,ignore_index=True);raw=normalize(raw)
 raw=raw[raw.BRISONr.notna()&~raw.BRISONr.isin(['BRISONr_?','BRISONr_NA'])&raw.Treatment.ne('LLN_WF_RF')].copy()
 save(parse,'raw_parse_inventory.csv')
 # Never merge two nonunique keys. Treat all ambiguous records as excluded.
 rawdup=raw.duplicated(K,keep=False);pubdup=d.duplicated(K,keep=False)
 ru=raw.loc[~rawdup,K+DS+['source_file','source_row']].copy()
 ru=ru.rename(columns={c:'raw_'+c for c in DS})
 linked=d.merge(ru,on=K,how='left',validate='many_to_one',indicator=True)
 assert len(linked)==len(d)
 linked['coordinates_complete']=linked[['Row','Column']].notna().all(axis=1)
 spatial_key=['Year','Location','Treatment','Block','Row','Column']
 linked['coordinate_collision']=linked.coordinates_complete&linked.duplicated(spatial_key,keep=False)
 linked['reliable_unit']=~pubdup.to_numpy()&linked._merge.eq('both')&~linked.coordinate_collision
 strict=linked.copy();maskcounts=[]
 for t in DS:
  v=pd.to_numeric(linked['raw_'+t].astype(str).str.replace(',','.',regex=False),errors='coerce')
  good=v.between(0,100)&linked.reliable_unit
  maskcounts.append(dict(trait=t,published_nonmissing=int(d[t].notna().sum()),raw_supported_nonmissing=int((good&d[t].notna()).sum()),published_zero_raw_missing_or_negative=int((d[t].eq(0)&~v.ge(0)).sum()),raw_noninteger_positive=int((v.gt(0)&v.mod(1).ne(0)).sum()),published_integer_truncation=int((v.gt(0)&d[t].notna()&v.ne(d[t])).sum()),raw_positive_truncated_to_zero=int((good&v.gt(0)&d[t].eq(0)).sum())))
  strict[t]=d[t].where(good) # retain official value; do not silently reintroduce raw precision
 strict['Seedyield']=strict.Seedyield.where(strict.reliable_unit)
 for t in units.index:
  if t in strict and t not in DS+['Seedyield']:strict[t]=strict[t].where(strict.reliable_unit)
 strict['TFI_complete']=strict[DS].sum(axis=1,min_count=6)
 strict['phase']=np.where(strict.Year.le(2017),'I','II')
 strict['unit_id']=strict[K].astype(str).fillna('NA').agg('|'.join,axis=1)
 save(maskcounts,'disease_missingness_provenance.csv')
 save(linked.loc[~linked.reliable_unit,K+['published_row','_merge','coordinate_collision']], 'unresolved_plot_links.csv')
 save(linked.loc[linked.coordinate_collision,K+['published_row','Seedyield']], 'coordinate_collisions.csv')
 save(d.loc[pubdup], 'ambiguous_published_records.csv')
 save(raw.loc[rawdup], 'ambiguous_raw_records.csv')
 save(strict[K+['published_row','unit_id','reliable_unit','coordinates_complete','coordinate_collision','source_file','source_row','phase']+DS+['Seedyield','TFI_complete']], 'audited_observations.csv')
 grain=[dict(metric='published_rows',value=len(d)),dict(metric='raw_rows_before_filter',value=sum(a['rows'] for a in parse)),dict(metric='cultivars',value=d.BRISONr.nunique()),dict(metric='locations',value=d.Location.nunique()),dict(metric='years',value=d.Year.nunique()),dict(metric='location_year',value=len(d[['Location','Year']].drop_duplicates())),dict(metric='location_year_management',value=len(d[E].drop_duplicates())),dict(metric='published_unique_full_keys',value=len(d[K].drop_duplicates())),dict(metric='published_ambiguous_rows',value=int(pubdup.sum())),dict(metric='reliable_source_linked_units',value=int(linked.reliable_unit.sum())),dict(metric='reliable_with_coordinates',value=int((linked.reliable_unit&linked.coordinates_complete).sum())),dict(metric='unresolved_rows',value=int((~linked.reliable_unit).sum())),dict(metric='exact_duplicate_rows',value=int(d.drop(columns='published_row').duplicated().sum())),dict(metric='repeated_dated_disease_assessments',value='NOT_IDENTIFIABLE: no observation-date column')]
 save(grain,'01_record_grain_audit.csv')
 save(grain+[dict(metric='coordinate_collision_rows',value=int(linked.coordinate_collision.sum()))],'01_record_grain_audit.csv')
 keys=[]
 for key in [E+['BRISONr'],E+['BRISONr','Block'],K,['Year','Location','Treatment','Block','Row','Column'],['Year','Location','Treatment','Row','Column']]:
  keys.append(dict(key='|'.join(key),rows=len(d),distinct_keys=len(d[key].drop_duplicates()),duplicate_rows=int(d.duplicated(key,keep=False).sum()),missing_key_rows=int(d[key].isna().any(axis=1).sum()),join_policy='Only unique full key to unique raw key; no many-to-many joins'))
 save(keys,'02_join_key_integrity.csv')
 cov=[]
 for basis,frame in [('published',d),('source_supported',strict)]:
  for t in units.index:
   if t not in frame:continue
   v=frame[t];s=frame[v.notna()]
   cov.append(dict(basis=basis,trait=t,unit=units.loc[t,'unit'],sample_source=units.loc[t,'sample.source'],raw_record_count=int(pd.to_numeric(raw[t].astype(str).str.replace(',','.',regex=False),errors='coerce').notna().sum()) if t in raw else np.nan,total_rows=len(frame),valid=int(v.notna().sum()),missing_rate=float(v.isna().mean()),zero_rate_valid=float(v.eq(0).sum()/v.notna().sum()) if v.notna().any() else np.nan,positive=int(v.gt(0).sum()),minimum=v.min(),q25=v.quantile(.25),median=v.median(),q75=v.quantile(.75),maximum=v.max(),years=s.Year.nunique(),sites=s.Location.nunique(),managements=s.Treatment.nunique(),environments=len(s[['Location','Year']].drop_duplicates()),out_of_range=int((~v.between(0,100)&v.notna()).sum()) if t in DS else int((v.lt(0)&v.notna()).sum()),repeated_dates='UNKNOWN: date absent'))
 coverage=save(cov,'03_trait_coverage.csv')
 em=[];over=[];support=[];phase=[]
 for env,g in strict.groupby(E):
  loc,year,treat=env
  row=dict(zip(E,env));row.update(rows=len(g),reliable_units=int(g.reliable_unit.sum()),cultivars=g.BRISONr.nunique(),replicates=g.Block.nunique(),yield_valid=int(g.Seedyield.notna().sum()),any_disease_yield=int((g.Seedyield.notna()&g[DS].notna().any(axis=1)).sum()),all_six_yield=int((g.Seedyield.notna()&g[DS].notna().all(axis=1)).sum()),phase='I' if year<=2017 else 'II',design_exception=('QLB/RHH 2017' if year==2017 and loc in ['QLB','RHH'] else 'GGE/KAL/RHH 2018' if year==2018 and loc in ['GGE','KAL','RHH'] else ''),coordinate_missing=int(g[['Row','Column']].isna().any(axis=1).sum()))
  for t in DS:
   valid=g[t].notna();paired=valid&g.Seedyield.notna();row[t+'_valid']=int(valid.sum());row[t+'_yield']=int(paired.sum())
   em.append(dict(zip(E,env),trait=t,valid=int(valid.sum()),positive=int(g[t].gt(0).sum()),zeros=int(g[t].eq(0).sum()),missing_rate=float((~valid).mean()),inoculation=flag_inoc(loc,year,treat,t)))
   phase.append(dict(zip(E,env),trait=t,phase=row['phase'],design_exception=row['design_exception'],cultivars=g.BRISONr.nunique(),reliable_units=int(g.reliable_unit.sum()),disease_valid=int(valid.sum()),yield_valid=row['yield_valid'],paired=int(paired.sum()),inoculation=flag_inoc(loc,year,treat,t),evidence='Scientific Data 2025 Methods; Scientific Reports 2020 Methods: 2015-2017 T1/T3',date_precision='season/phenological window only'))
   over.append(dict(zip(E,env),basis='source_supported',trait=t,paired_units=int(paired.sum())))
  row['disease_missing_rate']=float(g[DS].isna().to_numpy().mean());row['inoculation']='DISEASE_SPECIFIC: see 08'
  support.append(row)
 for basis,g in [('published',d),('source_supported',strict)]:
  for t in DS+['ALL_SIX','ANY_DISEASE']:
   ok=g[DS].notna().all(axis=1) if t=='ALL_SIX' else g[DS].notna().any(axis=1) if t=='ANY_DISEASE' else g[t].notna()
   s=g[ok&g.Seedyield.notna()];over.append(dict(Location='ALL',Year='ALL',Treatment='ALL',basis=basis,trait=t,paired_units=len(s),environments=len(s[['Location','Year']].drop_duplicates())))
 save(em,'04_disease_environment_matrix.csv');save(over,'05_disease_yield_overlap.csv');sup=save(support,'06_management_support_matrix.csv');save(phase,'08_inoculation_phase_audit.csv')
 cultivar_support=strict.groupby(E+['BRISONr'],dropna=False).agg(rows=('published_row','size'),reliable_units=('reliable_unit','sum'),yield_valid=('Seedyield','count'),**{t:(t,'count') for t in DS}).reset_index()
 save(cultivar_support,'management_cultivar_support.csv')
 comparisons=[]
 pairs=[('HN_WF_RF','HN_NF_RF'),('HN_NF_RF','LN_NF_RF'),('HN_WF_RF','LN_WF_RF'),('HN_WF_IR','HN_WF_RF'),('HN_WF_RO','HN_WF_RF')]
 for a,b in pairs:
  for (loc,year),g in strict.groupby(['Location','Year']):
   ga=g[g.Treatment.eq(a)&g.reliable_unit];gb=g[g.Treatment.eq(b)&g.reliable_unit];common=set(ga.BRISONr)&set(gb.BRISONr)
   if not len(ga) and not len(gb):continue
   rr=dict(Location=loc,Year=year,treatment_a=a,treatment_b=b,common_cultivars=len(common),units_a=len(ga),units_b=len(gb),shared_blocks='|'.join(sorted(set(ga.Block)&set(gb.Block))),inoculation_confounding=loc=='QLB' and (('_WF_' in a)!=('_WF_' in b)),status='COMPARABLE' if common else 'NOT_IDENTIFIABLE',interpretation='within-environment cultivar overlap; blocks are not plot pairs; not causal')
   if rr['inoculation_confounding']:rr['status']='LIMITED'
   for t in ['Seedyield']+DS:
    ca=set(ga.loc[ga[t].notna(),'BRISONr']);cb=set(gb.loc[gb[t].notna(),'BRISONr']);rr[t+'_common_cultivars']=len(ca&cb);rr[t+'_a']=int(ga[t].notna().sum());rr[t+'_b']=int(gb[t].notna().sum())
    if t in DS:
     pa=ga[t].notna()&ga.Seedyield.notna();pb=gb[t].notna()&gb.Seedyield.notna()
     rr[t+'_yield_units_a']=int(pa.sum());rr[t+'_yield_units_b']=int(pb.sum());rr[t+'_yield_common_cultivars']=len(set(ga.loc[pa,'BRISONr'])&set(gb.loc[pb,'BRISONr']))
   comparisons.append(rr)
 comp=save(comparisons,'07_management_comparison_feasibility.csv')
 # Daily weather: EU decimal comma, semicolon delimiter; calendar coverage is not forecast legality.
 weather=[]
 for f in sorted((r/'data/weather').glob('*')):
  w=pd.read_csv(f,sep=';',decimal=',');dates=pd.to_datetime(w.Date,errors='coerce');unique=dates.dropna().drop_duplicates();missing=pd.date_range(unique.min(),unique.max()).difference(unique)
  for c in w.columns:
   if c in ['Date','Location']:continue
   weather.append(dict(file=f.name,Location='|'.join(w.Location.unique()),variable=c,rows=len(w),valid=int(w[c].notna().sum()),earliest=str(dates.min().date()),latest=str(dates.max().date()),years='|'.join(map(str,sorted(dates.dt.year.unique()))),missing_dates=len(missing),duplicate_dates=int(dates.duplicated().sum()),minimum=w[c].min(),maximum=w[c].max(),prediction_legality='must truncate before explicit prediction cutoff; includes post-outcome years'))
 save(weather,'09_weather_temporal_coverage.csv')
 # Inspect original workbooks and keep exact cell evidence; no invented dates.
 cells=[];disease_dates=[];events=[]
 for f in sorted((r/'data/management').glob('*.xlsx')):
  x=pd.read_excel(f,header=None).fillna('');treatment=''
  for i in range(len(x)):
   vals=x.iloc[i].tolist()
   m=re.search(r'\(([^)]+)\)',str(vals[0]))
   if m and re.search(r'(HN|LN)',m.group(1)):treatment=m.group(1)
   for j,v in enumerate(vals):
    s=str(v)
    if any(word in s.lower() for word in ['harvest','sowing','inocul','irrigat','infection','measure date']):cells.append(dict(file=f.name,sheet=0,excel_row=i+1,excel_col=j+1,value=s,adjacent_values=' | '.join(str(a) for a in vals[j+1:j+4])))
    if s in DS:
     disease_dates.append(dict(file=f.name,excel_row=i+1,excel_col=j+1,treatment_raw=treatment,trait=s,date_raw=str(vals[j+1]) if j+1<len(vals) else '',date_status='PRESENT_UNPARSED' if j+1<len(vals) and str(vals[j+1]) else 'MISSING'))
 save(cells,'management_time_cell_evidence.csv');dd=save(disease_dates,'disease_date_evidence.csv')
 dd['date_parsed']=pd.to_datetime(dd.date_raw.where(dd.date_raw.str.fullmatch(r'\d{4}-\d{2}-\d{2} 00:00:00',na=False)),errors='coerce',format='%Y-%m-%d %H:%M:%S')
 dd['duplicate_template_entry']=dd.duplicated(['file','treatment_raw','trait','date_raw'],keep=False)
 save(dd,'disease_date_evidence.csv')
 save(dd.groupby('file').agg(template_trait_cells=('trait','size'),nonempty_date_cells=('date_raw',lambda s:s.ne('').sum()),parsed_date_cells=('date_parsed','count'),unique_dates=('date_parsed','nunique'),duplicate_template_cells=('duplicate_template_entry','sum')).reset_index(),'disease_date_file_summary.csv')
 for name in ['fertilizer','plant_protection','soil']:
  z=pd.read_excel(r/f'metadata/{name}.xlsx');events.append(dict(table=name,rows=len(z),date_nonmissing=int(z.Date.notna().sum()) if 'Date' in z else None,date_note='author extraction uses origin 1900-01-01; check original Excel dates before temporal features',columns='|'.join(z.columns)))
 save(events,'management_metadata_profile.csv')
 risks=[('R01','disease labels','HIGH','author replaces NA/negative disease values by zero within Y/L/M','mask with original raw source; strict tables exclude unsupported zeros'),('R02','early disease forecast','CRITICAL','plot table lacks disease assessment dates; workbook dates incomplete','NO-GO for proven leakage-free early forecast until dates and outcome definition resolved'),('R03','whole-season weather','HIGH','weather extends after target and through 2024','truncate by prespecified cutoff; retrospective weather only for association'),('R04','yield derived traits','CRITICAL','Grain=Seedyield*1000/TKW; Protein_yield=Seedyield*Crude_protein/100; Biomass/Straw derived via HI','exclude all yield descendants from yield predictors'),('R05','50 cm vs whole plot','HIGH','Seedyield_bio is mature 50 cm sample, different unit and sample source','do not substitute or count as independent plot'),('R06','random record split','HIGH','cultivar, blocks, year and site recur','split complete environments; fit preprocessing within training folds'),('R07','artificial inoculation','HIGH','QLB NF artificially inoculated; WF non-inoculated','disease-specific flag; QLB stratification or exclusion for fungicide comparison'),('R08','ambiguous identity','HIGH','duplicate or unmatched source keys','exclude ambiguous rows; row numbers are lineage only, not biological plot IDs'),('R09','summary outcome semantics','HIGH','QLB antecedent study reports AO from repeated scores; cleaning casts to integer','avoid pooled instantaneous severity claims; raw precision and temporal definition require Phase 2 review'),('R10','management actual versus planned','HIGH','observed fertilizer/fungicide events occur after sowing','pre-season model uses planned treatment only; actual events only before cutoff'),('R11','phase and year','HIGH','2020 only KAL; II panel/site composition differs','report composition shift; 2020 not representative future-year evidence'),('R12','date conversion','MEDIUM','Excel origin 1900-01-01 in author script may offset true workbook dates','use native workbook datetime; no automatic inherited-date legality'),('R13','prior BLUEs','HIGH','fitting adjusted means across full data could share test information','use plot outcomes or calculate summaries entirely within each fold')]
 save([dict(id=a,scope=b,severity=c,evidence=e,mitigation=f) for a,b,c,e,f in risks],'10_leakage_risk_registry.csv')
 folds=[]
 eligible=strict[strict.reliable_unit].copy()
 def fold(strategy,label,train,test):
  for target in ['Seedyield']+DS:
   tr=train[train[target].notna()];te=test[test[target].notna()]
   tc=set(tr.BRISONr);ec=set(te.BRISONr);tm=set(tr.Treatment);em=set(te.Treatment)
   folds.append(dict(strategy=strategy,test_group=label,target=target,train_units=len(tr),test_units=len(te),train_environments=len(tr[['Location','Year']].drop_duplicates()),test_environments=len(te[['Location','Year']].drop_duplicates()),test_cultivars=len(ec),seen_test_cultivars=len(tc&ec),unseen_test_cultivars=len(ec-tc),test_positive=int(te[target].gt(0).sum()),test_zero=int(te[target].eq(0).sum()),train_positive=int(tr[target].gt(0).sum()),train_zero=int(tr[target].eq(0).sum()),management_overlap='|'.join(sorted(tm&em)),unseen_management='|'.join(sorted(em-tm)),status='FEASIBLE_WITH_LIMITATIONS' if len(tr) and len(te) else 'NOT_IDENTIFIABLE',leakage_risk='environment split; no global imputation/BLUE fitting; AO/date caveat',independence='INTERNAL_HELD_OUT',limitation='2020 only KAL' if str(label)=='2020' else 'small environment count; shared cultivars; QLB intervention and phase composition shift'))
 for strategy,col in [('LOYO','Year'),('LOLO','Location'),('MANAGEMENT_TRANSFER','Treatment')]:
  for value in eligible[col].unique():fold(strategy,str(value),eligible[eligible[col].ne(value)],eligible[eligible[col].eq(value)])
 for (loc,y),g in eligible.groupby(['Location','Year']):fold('LOEO',f'{loc}_{y}',eligible[~(eligible.Location.eq(loc)&eligible.Year.eq(y))],g)
 fold('PHASE_I_TO_II','II',eligible[eligible.Year.le(2017)],eligible[eligible.Year.gt(2017)])
 validation=save(folds,'11_validation_feasibility.csv')
 # External scope deliberately bounded to three official candidates.
 ext=[]
 for name,doi,years,sites,diseases,yield_av,mgt,level,ind,compat,reason in [
  ('QLB BRIWECS adjusted means','10.5281/zenodo.3697514','2015-2017','1','stripe rust;leaf rust;powdery mildew;FHB','YES','YES','cultivar x treatment adjusted mean','NOT_INDEPENDENT','REJECT','Same BRIWECS QLB seasons; 178 German cultivars; repeated-score AO; not external'),
  ('56 European BASF fungicide trials','10.5281/zenodo.6521175','2017-2019','56 trials','SEPTTR;PUCCST;PUCCRE;PUCCRT;PUCCSI;PYRNTR;ERYSGA;ERYSGR;ERYSGT','YES_VERIFY_UNIT','YES','trial x treatment x date x leaf layer; plot IDs absent','LIKELY_INDEPENDENT_NOT_PROVEN','CONDITIONAL','Public table actually downloaded; no cultivar/replicate IDs; ERTRNE includes MCALCR relative yield and MWCGS grain weight per HA with no explicit mass unit; input intersection and units require reconciliation; trial-level test only'),
  ('Stage-specific GxE yield components','10.5281/zenodo.8248543','2014-2017','BRIWECS MET','not documented in deposited description','YES','PARTIAL','field yield components / climate','NOT_INDEPENDENT','REJECT','Paper explicitly identifies these climatic/trait data as earlier BRIWECS subset')]:
  ext.append(dict(dataset_name=name,DOI=doi,source='https://zenodo.org/records/'+doi.split('.')[-1],country='Europe' if 'BASF' in name else 'Germany',years=years,sites=sites,diseases=diseases,yield_available=yield_av,management_available=mgt,record_level=level,public_download='YES',independence_status=ind,compatibility=compat,reason=reason))
 save(ext,'12_external_dataset_candidates.csv')
 ef=ROOT/'sources/basf-wheat-diseases.txt'
 if ef.exists():
  ex=pd.read_csv(ef,sep='\t');save(ex.groupby(['Organism','Parameter','SampleUnit','Method','PlantPart'],dropna=False).agg(rows=('Value','size'),minimum=('Value','min'),maximum=('Value','max')).reset_index(),'external_basf_profile.csv')
  save(ex.groupby(['TrialId','Treatment','Parameter'],dropna=False).size().rename('rows').reset_index(),'external_basf_trial_support.csv')
 # Freeze both git sources, all actual files, exact byte comparisons.
 repos=[]
 for p in sorted((ROOT/'raw/github').iterdir()):
  if not p.is_dir():continue
  def git(*cmd):
   q=subprocess.run(['git','-C',str(p),*cmd],capture_output=True,text=True);return q.stdout.strip() if q.returncode==0 else 'UNAVAILABLE'
  repos.append(dict(repository=p.name,commit=git('rev-parse','HEAD'),branch=git('branch','--show-current'),remote=git('remote','get-url','origin'),status=git('status','--porcelain')))
 (ROOT/'sources/git_versions.json').write_text(json.dumps(repos,indent=2),encoding='utf-8')
 manifest=[]
 for f in sorted((ROOT/'raw').rglob('*')):
  if not f.is_file() or '.git' in f.parts:continue
  h=hashlib.sha256(f.read_bytes()).hexdigest();rel=f.relative_to(ROOT).as_posix();repo='figshare';version='v1';commit='';url='https://doi.org/10.6084/m9.figshare.27910269.v1'
  if 'github' in f.parts:
   repo=f.relative_to(ROOT/'raw/github').parts[0];rv=next(a for a in repos if a['repository']==repo);commit=rv['commit'];version=rv['branch'];url=rv['remote'].removesuffix('.git')+'/blob/'+commit+'/'+f.relative_to(ROOT/'raw/github'/repo).as_posix()
  manifest.append(dict(relative_path=rel,source_repository=repo,source_url=url,source_version=version,git_commit=commit,file_size_bytes=f.stat().st_size,sha256=h,retrieval_time_utc=dt.datetime.fromtimestamp(f.stat().st_mtime,dt.timezone.utc).isoformat(),validation_status='SHA256_COMPUTED',timestamp_note='local materialization mtime; ZIP members retain archive timestamps, not network retrieval time'))
 mf=save(manifest,'source_file_manifest.csv')
 archive=ROOT/'raw/figshare/Briwecs_data.zip';meta=json.loads((ROOT/'sources/figshare_article.json').read_text(encoding='utf-8'));expected=meta['files'][0]['computed_md5'];actual=hashlib.md5(archive.read_bytes()).hexdigest()
 assert actual==expected, 'Official archive MD5 mismatch'
 mf.loc[mf.relative_path.eq('raw/figshare/Briwecs_data.zip'),'validation_status']='OFFICIAL_MD5_MATCH_AND_SHA256_COMPUTED'
 mf['same_sha256_as_prior_file']=mf.sha256.duplicated()
 # Retrieval for archive members derives from archive download, not stored archive member times.
 archive_time=dt.datetime.fromtimestamp((ROOT/'raw/figshare/Briwecs_data.zip').stat().st_mtime,dt.timezone.utc).isoformat()
 mf.loc[mf.source_repository.eq('figshare'),'retrieval_time_utc']=archive_time
 mf.loc[mf.source_repository.eq('figshare'),'timestamp_note']='retrieval of containing official archive'
 mf.to_csv(ROOT/'audit/source_file_manifest.csv',index=False,encoding='utf-8-sig')
 matches=[];gr=ROOT/'raw/github/BRIWECS_Data_Publication'
 for f in r.rglob('*'):
  if not f.is_file():continue
  g=gr/f.relative_to(r)
  if g.exists():matches.append(dict(relative_path=f.relative_to(r).as_posix(),figshare_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),github_sha256=hashlib.sha256(g.read_bytes()).hexdigest(),status='IDENTICAL' if f.read_bytes()==g.read_bytes() else 'VERSION_DIFFERENCE'))
 save(matches,'figshare_github_comparison.csv')
 known=['BRIWECS_data_publication.csv','disease_record.xlsx','fertilizer.xlsx','plant_protection.xlsx','soil.xlsx','BRIWECS_BRISONr_information.csv','Unit.xlsx','scripts/run.R','scripts/pre-processing/data_cleaning.R','scripts/pre-processing/extract_management.R']
 save([dict(requested=n,actual_paths='|'.join(str(p.relative_to(r)) for p in r.rglob(Path(n).name)),status='PRESENT' if list(r.rglob(Path(n).name)) else 'MISSING_OR_RENAMED') for n in known],'expected_file_check.csv')
 # Informative descriptive figure, no paper figure production.
 c=coverage[(coverage.basis=='source_supported')&coverage.trait.isin(DS)]
 svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="300"><rect width="100%" height="100%" fill="white"/><text x="15" y="22" font-size="17">Source-supported disease coverage (green positive, blue zero, gray missing)</text>']
 for i,row in enumerate(c.itertuples()):
  y=45+i*39;svg.append(f'<text x="10" y="{y+17}" font-size="14">{row.trait}</text>');x=175
  for n,color in [(row.positive,'#24845a'),(row.valid-row.positive,'#4680b3'),(len(d)-row.valid,'#dadada')]:
   width=n/len(d)*640;svg.append(f'<rect x="{x}" y="{y}" width="{width}" height="24" fill="{color}"/>');x+=width
  svg.append(f'<text x="825" y="{y+17}" font-size="13">n={row.valid}, +={row.positive}</text>')
 svg.append('</svg>');(ROOT/'reports/disease_coverage.svg').write_text(''.join(svg),encoding='utf-8')
 # Version population counts, based on metadata trait list rather than guessed count.
 versions=[]
 for source,f in [('figshare_v1',r/'output/BRIWECS_data_publication.csv'),('github_current',gr/'output/BRIWECS_data_publication.csv')]:
  if f.exists():
   z=pd.read_csv(f,sep=';');traits=[t for t in units.index if t in z]
   versions.append(dict(source=source,rows=len(z),trait_columns=len(traits),nonnull_trait_cells=int(z[traits].notna().sum().sum()),nonnull_numeric_cells_excluding_keys=int(z.drop(columns=[a for a in K if a in z]).notna().sum().sum()),years='|'.join(map(str,sorted(z.Year.unique())))))
 save(versions,'version_population_counts.csv')
 z=pd.read_csv(gr/'output/BRIWECS_data_publication.csv',sep=';')
 ve=[]
 for source,frame in [('figshare_v1',d),('github_current',z)]:
  for env,g in frame.groupby(E):ve.append(dict(source=source,**dict(zip(E,env)),rows=len(g),nonnull_trait_cells=int(g[[t for t in units.index if t in g]].notna().sum().sum())))
 save(ve,'version_environment_counts.csv')
 vs=pd.DataFrame(ve).groupby(['source','Location','Year']).agg(rows=('rows','sum'),trait_cells=('nonnull_trait_cells','sum')).reset_index()
 vp=vs.pivot(index=['Location','Year'],columns='source',values='trait_cells').fillna(0)
 vp['delta']=vp.github_current-vp.figshare_v1
 save(vp.reset_index(),'version_delta_reconciliation.csv')
 code=[]
 for source,base in [('figshare_v1',r),('github_current',gr)]:
  for name in ['scripts/pre-processing/data_cleaning.R','scripts/pre-processing/functions.R','scripts/pre-processing/extract_management.R']:
   for line,text in enumerate((base/name).read_text(encoding='utf-8-sig').splitlines(),1):
    if re.search(r'NA|as.integer|write_delim|RHH_2018|BRISONr_|2015|2016|2017|KAL_table|write.xlsx\(disease',text):code.append(dict(source=source,path=name,line=line,text=text))
 save(code,'author_code_evidence.csv')
 # Actual per-phase panel sizes including design exceptions.
 ph=[]
 for (loc,y),g in strict.groupby(['Location','Year']):ph.append(dict(Location=loc,Year=y,phase='I' if y<=2017 else 'II',rows=len(g),cultivars=g.BRISONr.nunique(),valid_yield=int(g.Seedyield.notna().sum()),**{t:int(g[t].notna().sum()) for t in DS}))
 save(ph,'phase_actual_samples.csv')
 panels=[]
 for scope,frame in [('all_sites',d),('excluding_KAL',d[d.Location.ne('KAL')])]:
  p1=set(frame.loc[frame.Year.le(2017),'BRISONr']);p2=set(frame.loc[frame.Year.gt(2017),'BRISONr'])
  panels.append(dict(scope=scope,phase_I_unique_IDs=len(p1),phase_II_unique_IDs=len(p2),shared_IDs=len(p1&p2),interpretation='observed identifiers only; not proof of biological identity; KAL remapping in current author script requires sensitivity check'))
 save(panels,'phase_panel_overlap.csv')
 common_phase=len(set(d.loc[d.Year.le(2017),'BRISONr'])&set(d.loc[d.Year.gt(2017),'BRISONr']))
 summary=dict(published_rows=len(d),reliable_units=int(linked.reliable_unit.sum()),coordinate_complete_units=int((linked.reliable_unit&linked.coordinates_complete).sum()),coordinate_collision_rows=int(linked.coordinate_collision.sum()),unresolved_rows=int((~linked.reliable_unit).sum()),published_duplicate_rows=int(pubdup.sum()),cultivars=d.BRISONr.nunique(),sites=d.Location.nunique(),years=d.Year.nunique(),location_years=len(d[['Location','Year']].drop_duplicates()),management_environments=len(d[E].drop_duplicates()),managements=d.Treatment.nunique(),strict_yield=int(strict.Seedyield.notna().sum()),published_yield=int(d.Seedyield.notna().sum()),any_pair=int((strict[DS].notna().any(axis=1)&strict.Seedyield.notna()).sum()),six_pair=int((strict[DS].notna().all(axis=1)&strict.Seedyield.notna()).sum()),tfi_complete=int(strict.TFI_complete.notna().sum()),tfi_max=None if strict.TFI_complete.notna().sum()==0 else float(strict.TFI_complete.max()),cross_phase_cultivars=common_phase,disease_dates_present=int(dd.date_raw.ne('').sum()),source_files=len(mf),versions=versions)
 (ROOT/'audit/summary.json').write_text(json.dumps(summary,indent=2,default=int),encoding='utf-8')
 (ROOT/'environment/runtime.json').write_text(json.dumps(dict(python=sys.version,executable=sys.executable,pandas=pd.__version__,numpy=np.__version__,platform=platform.platform(),run_time_utc=NOW),indent=2),encoding='utf-8')
 print(json.dumps(summary,indent=2,default=int))
 from make_reports import make_reports
 make_reports(ROOT)
if __name__=='__main__':main()
