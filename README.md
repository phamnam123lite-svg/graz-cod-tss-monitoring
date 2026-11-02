# Effects of sampling interval and sensor outage timing on COD and TSS mass reconstruction in a combined sewer

Analysis code and results accompanying the manuscript by Viet Nam Tran, Van Lang University, Vietnam. This repository is a research supplement; it does not indicate conference acceptance or journal publication.

The study reuses public observations from Graz-West R05, Austria (2008–2011). It compares sampling schedules and simulated sensor outages in 66 retained rainfall-response windows. These are masks applied to existing observations, not new field experiments.

## Data source and attribution

Original data: Pichler, Hofer, Gamerith and Gruber (2026), version 0.1, [Zenodo record 15783832](https://doi.org/10.5281/zenodo.15783832). Dataset description: [10.3390/data11030045](https://doi.org/10.3390/data11030045). The source metadata specifies CC BY 4.0. Please cite the original creators when reusing the data. The manuscript author did not collect these field observations.

Raw measurements are downloaded from Zenodo by `prepare_data.py` and are not duplicated in this repository. The archive MD5 is `307525c96633e977e144357e014acfc8`.

## Reproduce the analysis

Use Python 3.12 in a virtual environment. From the repository folder:

```sh
python -m pip install -r requirements.txt
python prepare_data.py
python analyse_graz.py
python make_figures.py
python additional_analysis.py
python native_knot_analysis.py
python uncertainty_analysis.py
python make_supplemental_figures.py
```

The download requires internet access. If you already have the unchanged source archive, place it at `graz/data-set.zip`; the preparation script verifies and extracts it. Output files are written to `analysis/`. The analysis runs the main experiment and the sensitivity analyses with seed 20261101. Bootstrap samples resample complete events, not individual reconstruction replicates. Numerical results can differ slightly across platforms or software versions.

To read the saved outputs without rerunning:

```sh
python -m zipfile -e analysis_results.zip .
python -m zipfile -e supplemental_results.zip .
```

This extracts `analysis/`, including all result tables, event selection audits, replicate metrics and figures. A few key files are also provided at the repository root for easy inspection. `SHA256_manifest.json` records the checksums of the distributed files.

## Main findings and scope

There are 206 candidate windows, 66 retained windows and 228096 primary reconstructions. With linear reconstruction at 15-minute intervals, mean absolute mass errors are 2.45% for COD and 2.73% for TSS. A simulated 60-minute outage near peak flow raises these errors to 8.19% and 6.20%; randomly located outages of the same duration give 3.86% and 3.90%.

![Mean absolute mass errors by interval and outage scenario](figure3_outage_matrix.png)

COD and TSS are optical equivalents. Flow remains intact in the simulated concentration-outage experiments. The one-minute computational grid interpolates native observations only across gaps of at most six minutes; it does not represent independent one-minute field measurements. Analogue-flow windows reaching 2400 L/s, including buffers, are excluded. Digital flow is not used because its documented clock drift could distort matching with concentration.

The estimates describe sensor-derived, inlet-associated equivalent mass, including dry-weather background. They do not establish true pollutant discharge to receiving waters, universal monitoring intervals, or the accuracy of the source sensor calibration.

## Files

- `analyse_graz.py`: event screening, reconstruction experiments, event bootstrap summaries and sensitivity analyses. Includes checks of time conversion, gap handling and mass integration.
- `make_figures.py`: main figures, paired outage penalties and diagnostic traces.
- `main_summary.csv`: event-weighted main results.
- `events.csv`: retained windows.
- `paired_penalty_table.csv`: outage penalties relative to matched sampling phases.
- `analysis_results.zip`: full saved outputs, including selection audits, sensitivity tables and replicate metrics.
- `provenance.json`: parameters, source archive checksum and analysis software versions.

## Research transparency

AI assistance was used in developing the analysis code and drafting the accompanying research text. Sharing code and checks does not establish independent validation of the field measurements.

The original dataset retains its own CC BY 4.0 terms. No additional licence for the analysis code is assigned in this initial repository version.

## Manuscript version 3 supplement

The supplemental scripts and `supplemental_results.zip` now include every additional analysis reported in the full paper and ICUD extended abstract: event quantiles, paired contrasts, gap volume/mass exposure, flow-weighted random outages, flow-only regression, native-knot quadrature, seasonal/rainfall/flow subgroups and event-bootstrap Spearman intervals at all five sampling intervals. Files are stored under `supplemental_analysis/`; key CSV tables are also readable at the repository root.

Both COD and TSS native fractions have median 0.5898385565 across the 66 retained windows (range 0.3234873749–0.8961748634); individual event fractions can differ between pollutants. These fractions count exact original observations among the one-minute reference points, not a new set of independent measurements.

Rounding checks: the TSS peak-gap paired penalty is 3.3896575463 percentage points, while the previous-observation no-gap error is 3.3938250241%. TSS random-120 and peak-60 errors are 6.2023879494% and 6.2019061963%. These pairs happen to agree when rounded to two decimals. Paired penalties use matching phase draws rather than the all-phase no-outage mean.

Peak gaps are shifted within a window if necessary; their length is not shortened. All 66 peak and onset scenarios retain the full 60 minutes. Flow-weighted random placement is not mass-matched and does not isolate a causal timing effect. The native-knot comparison still uses linear assumptions and is not independent field validation. No second catchment or dry-weather-excess experiment is claimed.

See `verification_report.json` for numerical consistency checks. The versioned GitHub release and commit identify this snapshot; a GitHub release is not a Zenodo DOI. Original raw data remain attributed to their Zenodo creators. No manuscript, private contact details or portrait are included here.
