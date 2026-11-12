# Revision 5 (release v1.3.0), 12 November 2026

## Why
An external-style review of the manuscripts and of release v1.2.0 reproduced all 59 tables but raised scientific points that the code did not yet address. This revision adds the missing analyses. No primary result, event selection or numerical method changed.

## Added
- scripts/revision_checks.py: laboratory comparison, low-concentration screen, influence of large event errors, exposure-adjusted peak contrast, flow-resolution audit, and a check that the 2300 L/s sensitivity equals the flagged-window exclusion.
- results/archives/revision_results.zip and four summary tables in results/tables/.
- Two unit tests (laboratory pairing and agreement; recovery of a known effect by the exposure adjustment). The suite now has eight tests.
- scripts/run_pipeline.py and scripts/verify_outputs.py include the new step and archive.

## Findings carried into the manuscripts
1. Sensor-laboratory differences (median absolute 21% for local COD, 28% for global TSS) exceed the temporal reconstruction errors. The manuscripts now state that the errors are relative to the sensor record.
2. Eight windows contain local COD below 20 mg/L. Excluding them lowers the COD peak-gap mean from 8.19% to 7.18%; conclusions are unchanged in direction.
3. Means are driven by a few windows; medians are now reported beside means in the abstract and conclusions.
4. After linear adjustment for reference mass share, no separate peak-timing effect is detected. The manuscripts now describe the peak-gap excess as an exposure effect.
5. The 2300 L/s sensitivity and the flagged-window exclusion are one 65-window test.
6. Results at 3-6 minutes are presented as a consistency check only.

## Verification
A full rerun from the checksum-verified source archive passed eight tests and matched 66 numeric tables at rtol=atol=1e-10 and 40 release-file hashes (Python 3.13.18, NumPy 2.5.3, pandas 3.0.5, matplotlib 3.11.2). The log is verification/revision5_rerun.log.txt.

## Limits remaining
- The screen threshold is post hoc. The physical cause of low concentrations is unresolved.
- Only four laboratory samples fall inside retained windows, so event mass cannot be validated against laboratory data.
- The exposure adjustment is linear and descriptive, not causal.
- Live Zenodo download was not tested; the checksum-verified archive was supplied locally.
- No Zenodo DOI has been minted for this repository.
