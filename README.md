# Sampling intervals and sensor outages in COD and TSS mass reconstruction

Research supplement for the Graz-West R05 case study by Viet Nam Tran, Van Lang University. This repository is not evidence of conference acceptance or journal publication.

## Version and scope

Version 1.2.0 responds to the November 2026 pre-submission review. The repository keeps a flat script layout: run commands from the folder containing this README. Release v1.1.0 and its commit remain unchanged. Do not combine commands from this version with an uncommitted alternative src/data/results layout.

The primary analysis contains 66 rainfall-response windows and 228096 reconstructions. Flow stays available while concentration observations are reduced or masked. Mean absolute errors at 15-minute sampling with linear reconstruction are 2.45%/2.73% (COD/TSS) without outages, 3.86%/3.90% for random 60-minute gaps, and 8.19%/6.20% near peak flow. All means give equal weight to events.

## Reference quality and the flagged window

The reference was screened for data coverage and flow range; it is not independently validated pollutant mass. In event 20110604T1250, local COD was 0.1316 mg/L and global TSS was 10.7 mg/L at the flow peak. No observation-level quality flags, dated maintenance record or same-day laboratory sample was found in the supplied archive. A probe artefact is possible, but the cause cannot be established from these records. No unsupported universal lower concentration threshold is used.

The 66-event analysis is retained for transparency. Excluding this window gives no-gap errors of 2.33%/2.50%, random-gap errors of 3.73%/3.66%, and peak-gap errors of 7.53%/6.19%. Paired peak penalties become 5.09/3.64 percentage points. Full Table 2, Table 3, penalty and contrast sensitivities are in supplemental_results.zip. These checks do not establish that the other windows are free of sensor artefacts.

Both pollutant equivalents come from the same spectrometer and are not independent sensor replications. A separate nearest-original-observation test gives 3-, 6-, and 15-minute no-gap errors of 0.35/0.51%, 0.85/1.12%, and 2.64/2.99%. It selects the nearest native point within one minute, removes duplicate selections and drops scheduled points with no eligible observation. The reference itself remains interpolated, so this is not independent accuracy validation.

## Data source and licences

Original observations: Pichler, Hofer, Gamerith and Gruber, Zenodo version 0.1, https://doi.org/10.5281/zenodo.15783832. Dataset description: https://doi.org/10.3390/data11030045. Raw archive MD5: `307525c96633e977e144357e014acfc8`. The original data, including the small source-record extract in the supplement, retain CC BY 4.0 attribution requirements. Analysis code is MIT licensed; see LICENSE. AI tools assisted code development, analysis and writing.

## Reproduce from a clean clone

Use Python 3.12. Run in the repository root:

```text
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_outputs.py --rerun
```

The verifier runs prepare_data.py, test_analysis.py, analyse_graz.py, make_figures.py, additional_analysis.py, native_knot_analysis.py, uncertainty_analysis.py, make_supplemental_figures.py and qc_sensitivity.py. It then compares regenerated numeric tables with the archived release (rtol=atol=1e-10), checks release-file SHA256 hashes, and writes verification_report.json. Reproduction is numerical consistency, not independent validation of field measurements.

prepare_data.py downloads the original Zenodo archive. If a download is blocked, obtain the unchanged archive from the cited Zenodo record and place it at `graz/data-set.zip`; the script verifies its MD5 before extraction. A run using this already downloaded, checksum-verified archive is recorded in reproduction_log.txt. The source download endpoint is not guaranteed to be continuously available.

For separate steps, use `python run_pipeline.py`, followed by `python verify_outputs.py`. Analysis outputs go to analysis/ and supplemental_analysis/, never to the script directory. Importing a module does not run its analysis.

## Read saved results without running the analysis

```text
python -m zipfile -e analysis_results.zip .
python -m zipfile -e supplemental_results.zip .
```

The first archive contains analysis/ (primary data tables, replicate metrics, sensitivity tables and Figures 1–4). The second contains supplemental_analysis/ (additional analyses, Figure 5, exploratory correlation figure, native-sample and source-QC/exclusion outputs). Key tables are also supplied at the root. SHA256_manifest.json uses paths relative to the repository root; the manifest, verification report and reproduction log are excluded from its own list to avoid circular hashes and allow local verification to regenerate the report.

![Equal-event mean absolute error](figure3_outage_matrix.png)

## Interpretation and reproducibility notes

- The one-minute reference permits interpolation across at most six minutes. Original-point fractions have the same rounded median/range for COD and TSS (0.590; 0.323–0.896), although five windows differ. All source Q/COD/TSS timestamps are aligned to whole minutes; at_minutes nevertheless retains sub-minute timing for other inputs.
- The concentration shift test rebuilds its reference from shifted concentrations. It measures sampling-error stability, not clock-offset mass bias against an unchanged reference.
- Spearman results now describe correlation time tau versus error at each interval. Interval/tau ranks are the reverse of tau ranks within a fixed interval; the five intervals are not independent tests. Tau is computed on an interpolated series and searched only up to one quarter of each window.
- Primary means and all 228096 primary reconstruction rows are unchanged. Some bootstrap CI endpoints change because statistical streams are separated from plotting jitter and group iteration order. Correlation signs reverse with the stated predictor. See revision_changes.json.
- Peak/onset gaps retain their full 60-minute duration; peak gaps are shifted within bounds. A flow-weighted gap is not a mass-matched causal timing comparison. Native-knot quadrature changes integration support, not the virtual sampling source.
- At two decimals, TSS random-120 and peak-60 means both display 6.20%, from 6.2023879494 and 6.2019061963%. The TSS peak paired penalty (3.3896575463 points) and previous-value no-gap error (3.3938250241%) both display 3.39 by coincidence. Paired penalties use matching phase draws.
- No window exceeds 90% native points for both pollutants. The empty high_native_subset.csv is a documented negative diagnostic, not a failed or omitted calculation.
- Pichler and Muschalla (Water 2026, 18(20), 2486; https://doi.org/10.3390/w18202486) is an existing publication, confirmed on the publisher page on 3 November 2026.

The versioned GitHub release identifies the code and results. It has no Zenodo DOI. The unresolved physical cause of the flagged event should be reviewed with the supervisor or source-data authors before treating the manuscript as ready for submission.
