# Reproducibility state definitions

VERIFIED means an operation actually executed and passed in this closure stage: syntax/guarded imports, dependency/config checks, frozen cohort/fold/OOF verification, deterministic metric replay and reconstruction of figure numerical panels in a disposable scratch copy. It does not mean full training reproduction.

AVAILABLE_NOT_RERUN means genuine sources/configs and frozen inputs or fixed publicly archived acquisition targets are supplied, but the historical operation was not rerun: raw acquisition, cohort and label re-creation, all fits, bootstrap intervals, adjusted association estimation and figure rendering. Raw-data reconstruction first requires the documented version-locked download. Output-manifest regeneration requires a fresh sequential working run. Do not invoke full runners for verification.

BLOCKED applies to author/licence confirmation and public release where approval is absent. Clean-environment installation and cross-platform certification are AVAILABLE_NOT_RERUN, not demonstrated failures. No missing-input blocker affects the tested frozen replay path. No public-release licence has been selected, no repository has been published and no DOI has been obtained. Detailed stage states are in Reproducibility_Status.csv.
