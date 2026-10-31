# Reproducibility status

Public access, code licensing and computational reproducibility are separate.

## VERIFIED
Previously executed checks covered syntax/guarded imports, recorded dependency/configuration checks, frozen file SHA256 receipts, cohort and LOSYO assignments/test identities, deterministic metrics from frozen OOF predictions and figure numerical-panel reconstruction. VERIFIED describes these recorded operations, not a fresh end-to-end retraining exercise. Evidence: validation logs/results, Reproducibility_Status.csv and tools/validate_frozen.py.

## AVAILABLE_NOT_RERUN
Genuine historical code/configurations are supplied for original-data acquisition and cohort/label reconstruction, all Ridge/CatBoost fitting, bootstrap intervals and adjusted association estimation. These were not fully rerun in the frozen-output verification. A fresh installation, cross-platform certification and full raw-data-to-model reproduction remain undemonstrated. Run full pipelines only in a separate working copy; they overwrite outputs and perform fitting/inference. See README for fixed acquisition and sequential manifest reconstruction.

## PENDING
Original reanalysis-code rights/licence require rights-holder authorisation. A version-archive DOI has not been established; no DOI is invented. Third-party materials remain under their respective notices; any unresolved redistribution authority needs author review. Public access has already been authorised and the repository is public.

Reproducibility_Status.csv and historical verification reports describe the earlier local verification snapshot. Historical publication-approval status in those records must not be read as the current repository visibility. No historical scientific evidence has been rewritten.
