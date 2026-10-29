"""Independent consistency checks on delivered audit results, no model training."""
from pathlib import Path
import hashlib,json,subprocess
import pandas as pd
r=Path(__file__).resolve().parents[1];a=r/'audit';s=json.loads((a/'summary.json').read_text())
d=pd.read_csv(a/'audited_observations.csv');ds=['Stripe_rust','Leaf_rust','Powdery_mildew','Septoria','DTR','Fusarium']
assert len(d)==s['published_rows']
assert d.reliable_unit.sum()==s['reliable_units']
assert not d.loc[d.reliable_unit,'unit_id'].duplicated().any()
assert d.loc[~d.reliable_unit,ds+['Seedyield']].isna().all().all()
assert d.Seedyield.notna().sum()==s['strict_yield']
assert (d.Seedyield.notna()&d[ds].notna().any(axis=1)).sum()==s['any_pair']
assert (d.Seedyield.notna()&d[ds].notna().all(axis=1)).sum()==s['six_pair']
assert d.loc[d[ds].notna().all(axis=1),['Location','Year']].drop_duplicates().shape[0]==1
assert d.TFI_complete.notna().equals(d[ds].notna().all(axis=1))
assert d.loc[d.TFI_complete.notna(),'TFI_complete'].equals(d[ds].sum(axis=1,min_count=6).dropna())
mf=pd.read_csv(a/'source_file_manifest.csv')
for row in mf.itertuples():assert hashlib.sha256((r/row.relative_path).read_bytes()).hexdigest()==row.sha256
for repo in json.loads((r/'sources/git_versions.json').read_text()):
 q=subprocess.run(['git','-C',str(r/'raw/github'/repo['repository']),'status','--porcelain'],capture_output=True,text=True);assert q.returncode==0 and not q.stdout.strip()
for number in range(1,13):assert len(list(a.glob(f'{number:02d}_*.csv')))==1
for filename in ['BRIWECS_Phase1_Source_Audit.md','BRIWECS_Phase1_Scientific_Feasibility_Report.md','BRIWECS_Phase1_GO_NO_GO_Decision.md']:assert (r/'reports'/filename).stat().st_size>1000
ver=pd.read_csv(a/'version_population_counts.csv');assert int(ver.loc[ver.source.eq('figshare_v1'),'nonnull_trait_cells'].iloc[0])==526751
delta=pd.read_csv(a/'version_delta_reconciliation.csv');assert delta.delta.sum()==34971
print('PASS: unit identity, masking, overlaps, TFI, all raw SHA256, clean author repositories, 12 audit tables, 3 reports, version counts.')
print(delta[delta.delta.ne(0)].to_string(index=False))
g=pd.read_csv(r/'raw/github/BRIWECS_Data_Publication/output/BRIWECS_data_publication.csv',sep=';')
print('GitHub apparent shared phase cultivar IDs:',len(set(g.loc[g.Year.le(2017),'BRISONr'])&set(g.loc[g.Year.gt(2017),'BRISONr'])))
