# BRIWECS-ML: Stripe-rust information and winter wheat yield estimation

This public repository accompanies the manuscript "Context-dependent incremental value of stripe-rust scores for winter wheat yield estimation: a multi-environment reanalysis". It contains historical analysis scripts, fixed model configurations, reconstructed disease-label records, cohort and validation assignments, frozen out-of-fold predictions, numerical results and publication-figure generation code.

The repository supports inspection and verification of reported analyses. Frozen-result verification has been completed; a complete clean-environment, raw-data-to-model retraining exercise has not been independently demonstrated. Public accessibility does not confer an open-source licence for original code. No version-archive DOI is claimed.

## Quick verification without fitting

Use Python 3.12.14 and the exact recorded dependencies in Phase2A/requirements.txt or Phase2B_Closure/requirements.txt. In an environment with those packages installed, run from the repository root:

    python tools/validate_frozen.py

The verifier checks byte receipts, syntax, recorded dependency versions, cohort IDs, LOSYO train/test hashes and deterministic errors from frozen predictions. It does not invoke any historical runner or fit method. See REPRODUCIBILITY_STATUS.md for VERIFIED versus AVAILABLE_NOT_RERUN versus BLOCKED. A fresh environment installation and a complete raw-data-to-model reproduction were not performed in this stage.

## Official source acquisition

Original outcome data are not re-distributed. Download the immutable official archive with:

    python tools/prepare_sources.py

It requests https://api.figshare.com/v2/articles/27910269/versions/1 and verifies file ID 50832927, byte size, MD5 and SHA256 before safe extraction under raw/figshare. The source DOI is https://doi.org/10.6084/m9.figshare.27910269.v1. The lock is provenance/Figshare_v1_Lock.json. The prepared outcome table is used directly; original-author preprocessing need not be rerun. Figure inputs and the historical frozen cohorts/OOF outputs are already retained, so the verification command needs no network.

Phase 1's complete historical source audit additionally inspects original-author Git history. To materialise its exact pinned repositories in a separate working copy, use:

    python tools/prepare_sources.py --with-provenance-repositories

This clones and checks out fixed commits; it never uses current GitHub values as model outcomes. The third duplicate checkout exists only for the historical acquisition comparison. Do not run scripts/acquire_sources.py for a new acquisition: it is retained unchanged as historical source, and used a then-current v1 endpoint and current branches. The locked tool above replaces only the acquisition entry point.

## Source layout and execution order

| Path | Scientific function |
| --- | --- |
| scripts/run_phase1.py; audit/ | Plot identity, raw-source linkage, disease masks, provenance and QC |
| Phase2A/scripts/prepare_data.py | Yield cohort, F1/F2 features, weather legality and original validation assignments |
| Phase2A/scripts/train_models.py | Actual fitted pipelines including Ridge and CatBoost; fixed config and train-fold preprocessing |
| Phase2B_Pilot/scripts/prepare.py | Disease-specific cohorts and original decimal-label restoration, preserving real zeros |
| Phase2B_Pilot/scripts/train.py; evaluate.py | Actual paired models, combined OOF outputs and evaluation |
| Phase2B_Closure/scripts/prediction_analysis.py | Environment MSE accounting, rankings and common natural test-domain comparison |
| Phase2B_Closure/scripts/adjusted_association.py | Actual adjusted fixed-effect association analysis and historical bootstrap code |
| publication_figures/scripts/01–04 | Accepted Figure 1–4 data accounting, drawing and vector/raster export |
| frozen_SI/ | Single canonical copy of all 32 frozen scientific tables, configs and registries used by plots |

All stage paths derive from __file__ rather than a user directory. The figure loader's only change is an input-path adapter to frozen_SI, removing its duplicate data directory; mathematical and drawing code are unchanged. R 4.4.1, grid, svglite, ragg, systemfonts, jsonlite and Times New Roman are required. Figure export directories can be created with `python tools/create_output_directories.py`; it creates directories only and never draws. The recorded R session is retained. The historical code uses Windows curl.exe, Git and a Times New Roman font dependency; operating-system portability and a clean installation remain uncertified. Python source imports and relative-path replay were tested in this stage.

## Future full reproduction

Use a separate writable working copy. Never overwrite the distributed frozen outputs. A future authorised full reconstruction starts with source acquisition and scripts/run_phase1.py, followed by Phase2A/scripts/run_phase2a.py, Phase2B_Pilot/scripts/run_pilot.py and Phase2B_Closure/scripts/run_closure.py. These historical runners train models and/or calculate inference and must not be used for the read-only verification requested here. Their historical delivery manifests include omitted original model binaries and logs; simply invoking a downstream prepare script against an old delivery manifest will fail. A fresh sequential run must regenerate each stage's outputs and manifests before the next stage. Timestamps and production manifests will differ; numerical identities are validated separately.

For a focused future refit from supplied frozen cohorts, the genuine train_models.py and train.py entry points and exact model configs are present; they are AVAILABLE_NOT_RERUN, not a claim of current training success. Frozen historical configuration includes a dry-mass label whose applicability to every original record was not independently verified; interpretation follows the manuscript's explicit moisture-basis limitation, without changing the historical config.

## Statistical scope and attribution

LOSYO holds out complete site-years within a previously observed cultivar network; it is not independent external or new-genotype validation. Natural and artificial-inoculation-related records are explicitly distinguished. Disease scores have heterogeneous endpoints and are not unconditionally harmonised percentage points. Error variance uses ddof=0. Equal-site-year and sample-weighted within-environment components are different summaries; the latter are not the global pooled bias–variance decomposition. Cultivar-mean rankings are descriptive, with unequal management composition. Fixed seeds and model parameters remain in the historical JSON configurations; no hyperparameter search was added.

Original Figshare data and data-derived results require CC BY 4.0 attribution to Wang and colleagues and indication of reanalysis changes. The small Holst feasibility-audit dataset is separately attributed and did not enter model cohorts. Third-party repository source trees are downloaded on demand and their original licence notices remain under third_party_licences. No MIT/GPL licence is granted to original reanalysis code: see LICENSE_PENDING.md and FILE_OWNERSHIP.csv. Original-code rights/licensing and any future DOI remain pending. The public URL is stated above. Author and funding statements belong to the manuscript.
