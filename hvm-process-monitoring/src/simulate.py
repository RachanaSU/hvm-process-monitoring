"""Generate 60 days of simulated etch-module telemetry and load it into SQL.

Four chambers run the same recipe. Three problems are injected so the
analysis has something real to find:

* ETCH02 - step shift: after a preventive-maintenance (PM) visit on day 32 the
  throttle valve is mis-calibrated and chamber pressure sits 1.6 mTorr high.
* ETCH03 - slow drift: from day 20 the chiller degrades and chuck temperature
  creeps up by 0.25 C per day.
* ETCH04 - pump wear: the dry vacuum pump degrades (an accelerating wear
  process). Vibration rises until the pump trips, the chamber goes down for
  36 hours and the pump is replaced.

Raw data also carries the defects real logs have: missing values, duplicated
rows from a re-sent message, and sentinel values from sensor dropouts.
"""
import numpy as np
import pandas as pd

from . import process_model as pm
from .config import (LOTS_PER_DAY, N_DAYS, RECIPE, SAMPLED_PER_LOT, SEED,
                     TOOLS, VIB_TRIP, WAFERS_PER_LOT)
from .db import get_connection

START = pd.Timestamp("2026-07-01 06:00")
HOURS_PER_LOT = 24 / LOTS_PER_DAY

# Pump wear model (used to generate ETCH04 vibration): dD/dt = k*u(t)*(1 + D/Dc)
PUMP_WEAR = {
    "ETCH01": dict(v0=1.7, k=0.0020, dc=40.0),
    "ETCH02": dict(v0=1.9, k=0.0030, dc=40.0),
    "ETCH03": dict(v0=1.6, k=0.0025, dc=40.0),
    "ETCH04": dict(v0=1.8, k=0.0510, dc=1.5),
}
PUMP_DOWN_HOURS = 36


def pump_load(t_days):
    """Pump duty: busier during the day shift, never idle."""
    return 0.85 + 0.15 * np.sin(2 * np.pi * (t_days % 1.0))


def simulate_pumps(rng):
    """Hourly pump telemetry. Returns (dataframe, {tool: trip_time_days})."""
    hours = np.arange(0, N_DAYS * 24)
    rows, trips = [], {}
    for tool in TOOLS:
        p = PUMP_WEAR[tool]
        damage, down_until, tripped = 0.0, -1.0, False
        for h in hours:
            t = h / 24.0
            if t < down_until:      # chamber down, pump being replaced
                continue
            damage += p["k"] * pump_load(t) * (1 + damage / p["dc"]) / 24.0
            vib = p["v0"] + damage + rng.normal(0, 0.12)
            current = 11.0 + 0.55 * damage + 1.2 * pump_load(t) + rng.normal(0, 0.15)
            if vib >= VIB_TRIP and not tripped:
                trips[tool] = t
                tripped = True
                down_until = t + PUMP_DOWN_HOURS / 24.0
                damage = 0.0                 # new pump fitted
                p = dict(p, k=0.0020, dc=40.0)
            rows.append((tool, START + pd.Timedelta(hours=int(h)),
                         round(vib, 3), round(current, 3)))
    df = pd.DataFrame(rows, columns=["tool_id", "timestamp", "pump_vibration_mm_s", "pump_current_A"])
    return df, trips


def chamber_state(tool, t_days, pump_vib):
    """Actual chamber conditions for a wafer processed at time t (days)."""
    rf, pres, temp = RECIPE["rf_power_W"], RECIPE["pressure_mTorr"], RECIPE["chuck_temp_C"]
    pres_noise = 0.15
    if tool == "ETCH02" and t_days >= 32:
        pres += 1.6
    if tool == "ETCH03" and t_days >= 20:
        temp += 0.25 * (t_days - 20)
    if tool == "ETCH04":
        # A worn pump struggles to hold pressure: control gets noisier
        pres_noise = 0.15 + 0.12 * max(pump_vib - 2.0, 0)
    return rf, pres, temp, pres_noise


def simulate_wafers(rng, pumps, trips):
    run_rows, metro_rows = [], []
    pump_by_tool = {t: g.set_index("timestamp")["pump_vibration_mm_s"] for t, g in pumps.groupby("tool_id")}
    for tool in TOOLS:
        vib_series = pump_by_tool[tool]
        for day in range(N_DAYS):
            for lot_idx in range(LOTS_PER_DAY):
                t_lot = day + lot_idx * HOURS_PER_LOT / 24.0
                if tool in trips and trips[tool] <= t_lot < trips[tool] + PUMP_DOWN_HOURS / 24.0:
                    continue  # chamber down
                ts_lot = START + pd.Timedelta(days=t_lot)
                vib = float(vib_series.asof(ts_lot)) if ts_lot >= vib_series.index[0] else 1.8
                rf, pres, temp, pres_noise = chamber_state(tool, t_lot, vib)
                lot_id = f"{tool[-2:]}L{day:02d}{lot_idx}"
                lot_offset = rng.normal(0, pm.ER_LOT_SD)
                sampled = set(rng.choice(WAFERS_PER_LOT, SAMPLED_PER_LOT, replace=False))
                for w in range(WAFERS_PER_LOT):
                    ts = ts_lot + pd.Timedelta(minutes=6.5 * w)
                    m_rf = rf + rng.normal(0, 4.0)
                    m_p = pres + rng.normal(0, pres_noise)
                    m_t = temp + rng.normal(0, 0.2)
                    gas = 120 + rng.normal(0, 0.8)
                    wafer_id = f"{lot_id}W{w + 1:02d}"
                    er = pm.mean_etch_rate(m_rf, m_p, m_t) + lot_offset + rng.normal(0, pm.ER_WAFER_SD)
                    endpoint_s = 60 * 400 / er + rng.normal(0, 0.8)   # time to clear 400 nm
                    run_rows.append((wafer_id, lot_id, tool, ts, round(m_rf, 1), round(m_p, 3),
                                     round(m_t, 2), round(gas, 2), round(endpoint_s, 2)))
                    if w in sampled:
                        nu = pm.mean_nonuniformity(m_rf, m_p, m_t) + rng.normal(0, pm.NU_SD)
                        metro_rows.append((wafer_id, lot_id, tool, ts + pd.Timedelta(minutes=40),
                                           round(er, 2), round(nu, 3)))
    runs = pd.DataFrame(run_rows, columns=["wafer_id", "lot_id", "tool_id", "timestamp", "rf_power_W",
                                           "pressure_mTorr", "chuck_temp_C", "gas_flow_sccm", "endpoint_time_s"])
    metro = pd.DataFrame(metro_rows, columns=["wafer_id", "lot_id", "tool_id", "measured_at",
                                              "etch_rate_nm_min", "nonuniformity_pct"])
    return runs, metro


def add_log_defects(rng, df, value_cols, missing=0.005, glitch=0.001, dup=0.003):
    """Make raw logs look like real ones: gaps, sensor dropouts (-999) and re-sent rows."""
    df = df.copy()
    for col in value_cols:
        mask = rng.random(len(df)) < missing
        df.loc[mask, col] = np.nan
        mask = rng.random(len(df)) < glitch
        df.loc[mask, col] = -999.0
    dups = df.sample(frac=dup, random_state=int(rng.integers(1e9)))
    time_col = "timestamp" if "timestamp" in df.columns else "measured_at"
    return pd.concat([df, dups]).sort_values(time_col, kind="stable").reset_index(drop=True)


def main():
    rng = np.random.default_rng(SEED)
    pumps, trips = simulate_pumps(rng)
    runs, metro = simulate_wafers(rng, pumps, trips)

    runs_raw = add_log_defects(rng, runs, ["rf_power_W", "pressure_mTorr", "chuck_temp_C", "gas_flow_sccm"])
    metro_raw = add_log_defects(rng, metro, ["etch_rate_nm_min", "nonuniformity_pct"], missing=0.004, glitch=0.001)
    pumps_raw = add_log_defects(rng, pumps, ["pump_vibration_mm_s", "pump_current_A"], missing=0.004, glitch=0.001)

    con = get_connection()
    runs_raw.to_sql("wafer_runs", con, if_exists="replace", index=False)
    metro_raw.to_sql("metrology", con, if_exists="replace", index=False)
    pumps_raw.to_sql("pump_telemetry", con, if_exists="replace", index=False)

    total = len(runs_raw) + len(metro_raw) + len(pumps_raw)
    print(f"  wafer_runs: {len(runs_raw):,} rows | metrology: {len(metro_raw):,} | "
          f"pump_telemetry: {len(pumps_raw):,} | total {total:,} records")
    for tool, t in trips.items():
        print(f"  {tool} pump tripped at day {t:.2f}")
    return {"rows_wafer_runs": len(runs_raw), "rows_metrology": len(metro_raw),
            "rows_pump_telemetry": len(pumps_raw), "rows_total": total,
            "pump_trips_day": {k: round(v, 2) for k, v in trips.items()}}


if __name__ == "__main__":
    main()
