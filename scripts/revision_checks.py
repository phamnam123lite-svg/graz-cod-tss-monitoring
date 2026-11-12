"""Revision 5 checks requested in review.

Adds five analyses without changing any primary result:

1. Laboratory grab samples against the nearest original sensor observation.
2. A disclosed low-concentration screen applied to every retained window.
3. Influence of the largest event errors on equal-event means.
4. Peak-gap contrast adjusted for the reference mass share inside the gap.
5. Effective resolution of the analogue flow record.

It also records that the 2300 L/s saturation sensitivity and the flagged-window
exclusion remove the same single window, so they are one test, not two.

Outputs go to revision_analysis/. Nothing here validates field measurements.
"""
from pathlib import Path
from importlib.metadata import version
import json
import sys

import numpy as np
import pandas as pd

import analyse_graz as base

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "analysis"
SUPPLEMENT = ROOT / "supplemental_analysis"
OUT = ROOT / "revision_analysis"

FLAGGED = "20110604T1250"
LAB_TOLERANCE_MIN = 3          # largest native sensor spacing
SCREEN_THRESHOLDS_MGL = [10, 20, 30]
PRIMARY_SCREEN_MGL = 20
BOOTSTRAP_DRAWS = 2000
SCENARIOS = ["none", "random60", "peakQ60"]


def mean_ci(values):
    """Mean with a percentile 95% whole-event bootstrap interval."""
    values = np.asarray(values, float)
    rng = np.random.default_rng(base.SEED)
    draws = values[rng.integers(0, len(values), (BOOTSTRAP_DRAWS, len(values)))].mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(values.mean()), float(low), float(high)


def nearest_observation(series, when, tolerance_min):
    """Value and time offset (s) of the nearest valid original observation."""
    window = series.loc[when - pd.Timedelta(minutes=tolerance_min):
                        when + pd.Timedelta(minutes=tolerance_min)].dropna()
    if window.empty:
        return np.nan, np.nan
    offsets = (window.index - when).total_seconds()
    pick = int(np.argmin(np.abs(offsets)))
    return float(window.iloc[pick]), float(offsets[pick])


def agreement(lab, sensor):
    """Summary of sensor-minus-laboratory differences relative to laboratory."""
    lab = np.asarray(lab, float)
    sensor = np.asarray(sensor, float)
    relative = 100 * (sensor - lab) / lab
    ranks = pd.Series(lab).rank().corr(pd.Series(sensor).rank())
    return {
        "n_pairs": int(len(lab)),
        "median_signed_diff_pct": float(np.median(relative)),
        "mean_abs_diff_pct": float(np.abs(relative).mean()),
        "median_abs_diff_pct": float(np.median(np.abs(relative))),
        "pearson_r": float(np.corrcoef(lab, sensor)[0, 1]),
        "spearman_rho": float(ranks),
    }


def laboratory_comparison(events):
    """Pair each grab sample with the nearest original sensor reading."""
    sensors = {
        "COD_local": base.read_series("pollutant data/calibrated.csv", "COD"),
        "COD_global": base.read_series("pollutant data/spectrometer.csv", "COD"),
        "TSS_global": base.read_series("pollutant data/spectrometer.csv", "TSS"),
    }
    lab = pd.read_csv(base.RAW / "lab data/samples.csv")
    lab["time_UTC"] = pd.to_datetime(lab.time, utc=True)
    lab["lab_COD"] = lab[["COD1", "COD2", "COD3"]].mean(axis=1)
    lab["lab_TSS"] = lab[["TSS1", "TSS2", "TSS3"]].mean(axis=1)
    rows = []
    for _, sample in lab.iterrows():
        row = {"time_UTC": sample.time_UTC, "lab_COD_mgL": sample.lab_COD, "lab_TSS_mgL": sample.lab_TSS,
               "inside_retained_window": any(e["start"] <= sample.time_UTC <= e["end"] for e in events)}
        for label, series in sensors.items():
            row[label + "_mgL"], row[label + "_offset_s"] = nearest_observation(series, sample.time_UTC, LAB_TOLERANCE_MIN)
        rows.append(row)
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT / "lab_sensor_pairs.csv", index=False)

    summary = []
    for label, lab_column in [("COD_local", "lab_COD_mgL"), ("COD_global", "lab_COD_mgL"), ("TSS_global", "lab_TSS_mgL")]:
        ok = pairs.dropna(subset=[label + "_mgL", lab_column])
        summary.append({"sensor_series": label, **agreement(ok[lab_column], ok[label + "_mgL"])})
    summary = pd.DataFrame(summary)
    summary.to_csv(OUT / "lab_sensor_summary.csv", index=False)
    info = {
        "lab_samples": int(len(lab)),
        "sampling_days": int(lab.time_UTC.dt.date.nunique()),
        "first_sample_UTC": str(lab.time_UTC.min()),
        "last_sample_UTC": str(lab.time_UTC.max()),
        "samples_inside_retained_windows": int(pairs.inside_retained_window.sum()),
        "pairing_tolerance_min": LAB_TOLERANCE_MIN,
        "note": "The local COD calibration was derived from laboratory analyses of this station, so local-COD agreement is not out-of-sample. Global COD and TSS use the manufacturer calibration.",
    }
    return summary, info


def low_concentration_screen(events):
    """Count original local-COD observations below each threshold in every window."""
    cod = base.read_series("pollutant data/calibrated.csv", "COD")
    flow = base.read_series("hydraulic data/inflow-analogue.csv", "flow_rate")
    all_cod = cod.dropna()
    rows = []
    for event in events:
        inside = cod.loc[event["start"]:event["end"]].dropna()
        row = {"event_id": event["event_id"], "COD_min_mgL": float(inside.min()), "Q_peak_Ls": event["Q_peak_Ls"]}
        for threshold in SCREEN_THRESHOLDS_MGL:
            low = inside[inside < threshold]
            row[f"n_below_{threshold}"] = int(len(low))
        low = inside[inside < PRIMARY_SCREEN_MGL]
        row["flow_during_low_max_Ls"] = float(flow.reindex(low.index).max()) if len(low) else np.nan
        rows.append(row)
    windows = pd.DataFrame(rows)
    windows.to_csv(OUT / "signal_screen_windows.csv", index=False)
    record_share = {str(t): float(100 * (all_cod < t).mean()) for t in SCREEN_THRESHOLDS_MGL}
    return windows, record_share


def screen_sensitivity(windows, event_metrics):
    """Headline 15-minute linear errors after removing screened windows."""
    view = event_metrics[(event_metrics.interval_min == 15) & (event_metrics.method == "linear")
                         & event_metrics.scenario.isin(SCENARIOS)]
    samples = {"all_66": set()}
    for threshold in SCREEN_THRESHOLDS_MGL:
        samples[f"exclude_COD_below_{threshold}"] = set(windows.loc[windows[f"n_below_{threshold}"] > 0, "event_id"])
    rows = []
    for name, removed in samples.items():
        kept = view[~view.event_id.isin(removed)]
        for (pollutant, scenario), group in kept.groupby(["pollutant", "scenario"]):
            errors = group.mean_abs_pct.to_numpy()
            mean, low, high = mean_ci(errors)
            penalty = mean_ci(group.paired_outage_penalty_pp)
            rows.append({"sample": name, "pollutant": pollutant, "scenario": scenario, "n_events": len(group),
                         "n_removed": len(removed), "mean_abs_pct": mean, "ci_low": low, "ci_high": high,
                         "median_abs_pct": float(np.median(errors)),
                         "paired_penalty_pp": penalty[0], "penalty_ci_low": penalty[1], "penalty_ci_high": penalty[2]})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "signal_screen_sensitivity.csv", index=False)
    return table


def influence(event_metrics):
    """How much the largest event errors drive the equal-event mean."""
    view = event_metrics[(event_metrics.interval_min == 15) & (event_metrics.method == "linear")]
    rows = []
    for (pollutant, scenario), group in view.groupby(["pollutant", "scenario"]):
        errors = np.sort(group.mean_abs_pct.to_numpy())[::-1]
        trim = int(np.floor(0.1 * len(errors)))
        rows.append({"pollutant": pollutant, "scenario": scenario, "n_events": len(errors),
                     "mean_abs_pct": float(errors.mean()), "median_abs_pct": float(np.median(errors)),
                     "trimmed_mean_10pct": float(errors[trim:len(errors) - trim].mean()),
                     "mean_without_top5": float(errors[5:].mean()),
                     "top5_share_of_sum_pct": float(100 * errors[:5].sum() / errors.sum()),
                     "largest_event_error_pct": float(errors[0])})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "influence_summary.csv", index=False)
    return table


def within_event_fit(frame):
    """OLS with event fixed effects: error ~ mass share + flow-weighted + peak."""
    columns = ["gap_mass_share_pct", "is_flowweighted", "is_peak"]
    centred = frame[columns + ["mean_abs_pct"]] - frame.groupby("event_id")[columns + ["mean_abs_pct"]].transform("mean")
    return np.linalg.lstsq(centred[columns].to_numpy(), centred.mean_abs_pct.to_numpy(), rcond=None)[0]


def exposure_adjusted_contrast():
    """Peak-minus-uniform error difference at equal mean reference mass share.

    Each event contributes three placements (uniform, flow-weighted, peak). Event
    fixed effects remove between-event differences; events are resampled whole.
    The adjustment is linear and descriptive, not a causal timing estimate.
    """
    metrics = pd.read_csv(SUPPLEMENT / "additional_event_metrics.csv")
    metrics = metrics[(metrics.method == "linear") & (metrics.scenario != "none")].copy()
    metrics["is_flowweighted"] = (metrics.scenario == "flowweighted60").astype(float)
    metrics["is_peak"] = (metrics.scenario == "peakQ60").astype(float)
    rows = []
    for pollutant, frame in metrics.groupby("pollutant"):
        frame = frame.sort_values(["event_id", "scenario"]).reset_index(drop=True)
        estimate = within_event_fit(frame)
        ids = frame.event_id.unique()
        blocks = {event_id: block for event_id, block in frame.groupby("event_id")}
        rng = np.random.default_rng(base.SEED)
        draws = []
        for _ in range(BOOTSTRAP_DRAWS):
            chosen = rng.integers(0, len(ids), len(ids))
            sample = pd.concat([blocks[ids[k]].assign(event_id=n) for n, k in enumerate(chosen)], ignore_index=True)
            draws.append(within_event_fit(sample))
        low, high = np.quantile(np.array(draws), [0.025, 0.975], axis=0)
        unadjusted = frame.pivot(index="event_id", columns="scenario", values="mean_abs_pct")
        rows.append({"pollutant": pollutant, "n_events": len(ids),
                     "unadjusted_peak_minus_uniform_pp": float((unadjusted.peakQ60 - unadjusted.random60).mean()),
                     "slope_pp_per_mass_share_pct": float(estimate[0]), "slope_ci_low": float(low[0]), "slope_ci_high": float(high[0]),
                     "adjusted_peak_minus_uniform_pp": float(estimate[2]), "adjusted_ci_low": float(low[2]), "adjusted_ci_high": float(high[2])})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "exposure_adjusted_contrast.csv", index=False)
    return table


def flow_resolution(events):
    """Share of consecutive identical analogue flow readings inside each window."""
    flow = base.read_series("hydraulic data/inflow-analogue.csv", "flow_rate")
    rows = []
    for event in events:
        inside = flow.loc[event["start"]:event["end"]].dropna()
        rows.append({"event_id": event["event_id"], "flow_readings": int(len(inside)),
                     "identical_to_previous_share": float((inside.diff() == 0).mean())})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "flow_resolution.csv", index=False)
    return table


def duplicate_sensitivity_check():
    """The 2300 L/s test and the flagged-window exclusion remove the same window."""
    main = pd.read_csv(ANALYSIS / "event_audit_gap6_sep6_tail3_sat2400.csv")
    lower = pd.read_csv(ANALYSIS / "event_audit_gap6_sep6_tail3_sat2300.csv")
    removed = sorted(set(main.loc[main.retained, "event_id"]) - set(lower.loc[lower.retained, "event_id"]))
    return {"windows_removed_by_2300_Ls_limit": removed, "identical_to_flagged_exclusion": removed == [FLAGGED]}


def main():
    OUT.mkdir(exist_ok=True)
    events = pd.read_csv(ANALYSIS / "events.csv", parse_dates=["start", "end"]).to_dict("records")
    event_metrics = pd.read_csv(ANALYSIS / "main_event_metrics.csv")

    lab_summary, lab_info = laboratory_comparison(events)
    windows, record_share = low_concentration_screen(events)
    screen = screen_sensitivity(windows, event_metrics)
    influence_table = influence(event_metrics)
    exposure = exposure_adjusted_contrast()
    resolution = flow_resolution(events)
    duplicate = duplicate_sensitivity_check()

    report = {
        "laboratory_comparison": lab_info,
        "low_concentration_screen": {
            "series": "locally calibrated COD equivalent, original observations inside each window",
            "thresholds_mgL": SCREEN_THRESHOLDS_MGL,
            "primary_threshold_mgL": PRIMARY_SCREEN_MGL,
            "share_of_whole_record_below_threshold_pct": record_share,
            "windows_flagged": {str(t): int((windows[f"n_below_{t}"] > 0).sum()) for t in SCREEN_THRESHOLDS_MGL},
            "flagged_event_ids_primary": windows.loc[windows[f"n_below_{PRIMARY_SCREEN_MGL}"] > 0, "event_id"].tolist(),
            "note": "Post-hoc sensitivity screen. Low values can be real dilution; the screen does not classify them as sensor failure.",
        },
        "flow_resolution": {
            "median_identical_share": float(resolution.identical_to_previous_share.median()),
            "min_identical_share": float(resolution.identical_to_previous_share.min()),
            "max_identical_share": float(resolution.identical_to_previous_share.max()),
        },
        "saturation_2300_vs_flagged_exclusion": duplicate,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "seed": base.SEED,
        "python": sys.version,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "matplotlib": version("matplotlib"),
    }
    (OUT / "revision_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(lab_summary.to_string(index=False), flush=True)
    print(screen[screen["sample"].isin(["all_66", f"exclude_COD_below_{PRIMARY_SCREEN_MGL}"])].to_string(index=False), flush=True)
    print(influence_table.to_string(index=False), flush=True)
    print(exposure.to_string(index=False), flush=True)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
