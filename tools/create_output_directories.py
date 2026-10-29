from pathlib import Path
R=Path(__file__).resolve().parents[1]
for stage in ['Phase2A','Phase2B_Pilot','Phase2B_Closure']:
 for folder in ['logs','models','reports','results/figures','results/metrics','results/tables','results/predictions']:(R/stage/folder).mkdir(parents=True,exist_ok=True)
for folder in ['figures','qa','reports','source_data/panels']:(R/'publication_figures'/folder).mkdir(parents=True,exist_ok=True)
print('Output directories created; no analysis run')
