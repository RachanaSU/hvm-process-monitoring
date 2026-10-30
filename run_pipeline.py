"""Run the full pipeline end to end.

    python run_pipeline.py

1. Simulate 60 days of etch-module telemetry and load it into SQLite.
2. Extract and clean it with SQL + pandas (data-quality report).
3. SPC: X-bar/R charts, Western Electric rules, Cp/Cpk, ANOVA.
4. Reliability: fit the pump wear ODE and forecast failure.
5. DOE: optimise the recipe and confirm it.
6. FMEA: rank failure modes by RPN, before and after the actions.

Everything is written to outputs/ (CSV tables, summary.json, figures/).
"""
import json

import numpy as np
import pandas as pd

from src import doe, etl, reliability, simulate, spc
from src.config import OUT_DIR, ROOT, TOOLS
from src.db import driver_name

INJECTED = {"ETCH02": "pressure step after PM on day 32",
            "ETCH03": "chuck-temperature drift from day 20",
            "ETCH04": "pump wear (trips ~day 51)"}
EVENT_DAY = {"ETCH02": 32, "ETCH03": 20}


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else round(float(obj), 4)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    return obj


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = {}

    print(f"[1/6] Simulating telemetry and loading SQL database ({driver_name()})")
    summary["data"] = simulate.main()

    print("[2/6] Extracting and cleaning")
    metro, pumps, subgroups, dq, sensors, n_runs = etl.run()
    dq.to_csv(OUT_DIR / "data_quality_report.csv", index=False)
    print(f"  {n_runs:,} wafer sensor records, {len(metro):,} clean metrology wafers, {len(subgroups):,} complete lot subgroups, "
          f"{int(dq['rows_affected'].sum())} raw rows removed or repaired")
    summary["clean"] = {"metrology_wafers": len(metro), "lot_subgroups": len(subgroups),
                        "pump_readings": len(pumps), "rows_fixed": int(dq["rows_affected"].sum())}

    print("[3/6] Statistical process control")
    charts, signals, cap = spc.analyse(metro, subgroups)
    signals.to_csv(OUT_DIR / "spc_signals.csv", index=False)
    cap.to_csv(OUT_DIR / "capability.csv", index=False)
    detect = spc.first_detection(signals, EVENT_DAY)
    healthy = charts["ETCH01"][0]
    healthy_mon = healthy[healthy["day"] >= spc.BASELINE_DAYS]
    false_alarms = {"three_way_chart": int((healthy_mon["rule"] != "").sum()),
                    "textbook_A2_limits": spc.classic_false_alarms(healthy_mon, charts["ETCH01"][1]),
                    "lots_monitored": len(healthy_mon)}
    an_base, _ = spc.anova(metro, 0, 15)
    an_late, tukey = spc.anova(metro, 50, 60)
    tukey.to_csv(OUT_DIR / "anova_tukey_days50-60.csv", index=False)
    spc.figures(charts, metro, INJECTED)
    spc.capability_figure(cap)
    spc.anova_figure(metro, 50, 60, an_late)
    spc.root_cause_figure(sensors)
    for tool, d in detect.items():
        print(f"  {tool}: event day {d['event_day']} -> first signal day {d['first_signal_day']} "
              f"(lag {d['detection_lag_days']} d, rule {d['first_rules']})")
    print(f"  ETCH01 (healthy) false alarms over {false_alarms['lots_monitored']} lots: "
          f"{false_alarms['three_way_chart']} (textbook limits would give {false_alarms['textbook_A2_limits']})")
    print(f"  ANOVA days 0-15: p = {an_base['p_value']:.3f} | days 50-60: F = {an_late['F']}, p = {an_late['p_value']:.1e}")
    summary["spc"] = {"detection": detect, "false_alarms_ETCH01": false_alarms,
                      "capability": cap.to_dict(orient="records"),
                      "anova_baseline": an_base, "anova_days50_60": an_late}

    print("[4/6] Pump reliability model")
    rel, detail, rolling = reliability.analyse(pumps)
    rel.to_csv(OUT_DIR / "pump_forecast_summary.csv", index=False)
    rolling.to_csv(OUT_DIR / "pump_forecast_rolling.csv", index=False)
    life, _ = reliability.first_life(pumps, "ETCH04")
    life[["day", "pump_vibration_mm_s"]].to_csv(OUT_DIR / "pump_ETCH04_first_life.csv", index=False)
    reliability.figures(pumps, rel, detail, rolling)
    r4 = rel.set_index("tool_id").loc["ETCH04"]
    print(f"  ETCH04 forecast on day {reliability.FORECAST_DAY}: failure day {r4['predicted_trip_day']:.2f} "
          f"(90% {r4['trip_p05']:.2f}-{r4['trip_p95']:.2f}); actual {r4['actual_trip_day']:.2f}; "
          f"straight-line trend said {r4['linear_trend_trip_day']:.2f}; planned PM day {r4['recommended_pm_day']:.0f}")
    summary["reliability"] = rel.replace({np.nan: None}).to_dict(orient="records")

    print("[5/6] Design of experiments")
    d = doe.run()
    d["runs"].to_csv(OUT_DIR / "doe_runs.csv", index=False)
    d["effects_er"].to_csv(OUT_DIR / "doe_effects_etch_rate.csv", index=False)
    d["effects_nu"].to_csv(OUT_DIR / "doe_effects_nonuniformity.csv", index=False)
    doe.figures(d)
    b, a = d["before"], d["after"]
    print(f"  new recipe {d['new_recipe']}")
    print(f"  non-uniformity {b['nu_mean']:.2f}% -> {a['nu_mean']:.2f}% | Ppk(NU) {b['ppk_nu']:.2f} -> {a['ppk_nu']:.2f} | "
          f"Ppk(ER) {b['ppk_er']:.2f} -> {a['ppk_er']:.2f} | yield {b['yield_pct']:.1f}% -> {a['yield_pct']:.1f}%")
    summary["doe"] = {k: d[k] for k in ("doe_recipe", "new_recipe", "pred_er", "pred_nu", "curvature",
                                        "before", "first_confirmation", "after")}

    print("[6/6] FMEA")
    fmea = pd.read_csv(ROOT / "reports" / "fmea.csv")
    fmea["RPN"] = fmea["S"] * fmea["O"] * fmea["D"]
    fmea["RPN_after"] = fmea["S_after"] * fmea["O_after"] * fmea["D_after"]
    fmea = fmea.sort_values("RPN", ascending=False)
    fmea.to_csv(OUT_DIR / "fmea_ranked.csv", index=False)
    print(f"  total RPN {fmea['RPN'].sum()} -> {fmea['RPN_after'].sum()} after actions; "
          f"top risk: {fmea.iloc[0]['failure_mode']} (RPN {fmea.iloc[0]['RPN']})")
    summary["fmea"] = fmea[["id", "failure_mode", "RPN", "RPN_after"]].to_dict(orient="records")

    with open(OUT_DIR / "summary.json", "w") as f:
        json.dump(to_jsonable(summary), f, indent=2, default=str)
    print(f"Done. Results in {OUT_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
