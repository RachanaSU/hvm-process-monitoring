"""Extract, clean and shape the raw SQL logs for analysis.

Replaces the manual spreadsheet step: one SQL query per table, then a fixed set
of cleaning rules, with every row that is removed or repaired counted in a
data-quality report.
"""
import numpy as np
import pandas as pd

from .config import N_DAYS
from .db import get_connection
from .simulate import START

# Physically plausible ranges. Anything outside is a sensor fault, not process data.
VALID_RANGES = {
    "rf_power_W": (1000, 2000),
    "pressure_mTorr": (10, 60),
    "chuck_temp_C": (20, 100),
    "gas_flow_sccm": (50, 200),
    "etch_rate_nm_min": (100, 400),
    "nonuniformity_pct": (0, 20),
    "pump_vibration_mm_s": (0, 30),
    "pump_current_A": (0, 40),
}

QUERY_METROLOGY = """
SELECT m.wafer_id, m.lot_id, m.tool_id, m.measured_at,
       m.etch_rate_nm_min, m.nonuniformity_pct,
       r.timestamp AS processed_at, r.rf_power_W, r.pressure_mTorr,
       r.chuck_temp_C, r.gas_flow_sccm
FROM metrology AS m
LEFT JOIN (SELECT DISTINCT * FROM wafer_runs) AS r USING (wafer_id)
"""


def _apply_ranges(df, report, table):
    for col, (lo, hi) in VALID_RANGES.items():
        if col in df.columns:
            bad = df[col].notna() & ~df[col].between(lo, hi)
            report.append((table, f"{col}: out-of-range / sentinel values set to missing", int(bad.sum())))
            df.loc[bad, col] = np.nan
    return df


def load_metrology(con, report):
    df = pd.read_sql(QUERY_METROLOGY, con, parse_dates=["measured_at", "processed_at"])
    n0 = len(df)
    df = df.drop_duplicates(subset="wafer_id", keep="first")
    report.append(("metrology", "duplicate wafer records removed", n0 - len(df)))
    df = _apply_ranges(df, report, "metrology")
    missing_ctq = df["etch_rate_nm_min"].isna()
    report.append(("metrology", "wafers dropped: no valid etch-rate reading", int(missing_ctq.sum())))
    df = df[~missing_ctq].copy()
    df["day"] = (df["processed_at"].fillna(df["measured_at"]) - START).dt.total_seconds() / 86400
    return df.sort_values(["tool_id", "processed_at"]).reset_index(drop=True)


def load_wafer_runs(con, report):
    """Chamber sensor log for every wafer (about 48,000 rows)."""
    df = pd.read_sql("SELECT * FROM wafer_runs", con, parse_dates=["timestamp"])
    n0 = len(df)
    df = df.drop_duplicates(subset="wafer_id", keep="first")
    report.append(("wafer_runs", "duplicate wafer records removed", n0 - len(df)))
    df = _apply_ranges(df, report, "wafer_runs")
    df["day"] = (df["timestamp"] - START).dt.total_seconds() / 86400
    return df.sort_values(["tool_id", "timestamp"]).reset_index(drop=True)


def daily_sensor_means(runs):
    """Daily mean of each chamber sensor: used to trace an SPC signal back to its cause."""
    runs = runs.assign(day_int=runs["day"].astype(int))
    return (runs.groupby(["tool_id", "day_int"])[["rf_power_W", "pressure_mTorr", "chuck_temp_C", "gas_flow_sccm"]]
            .mean().reset_index().rename(columns={"day_int": "day"}))


def load_pumps(con, report):
    df = pd.read_sql("SELECT * FROM pump_telemetry", con, parse_dates=["timestamp"])
    n0 = len(df)
    df = df.drop_duplicates()
    report.append(("pump_telemetry", "duplicate rows removed", n0 - len(df)))
    df = _apply_ranges(df, report, "pump_telemetry")
    df = df.sort_values(["tool_id", "timestamp"]).reset_index(drop=True)
    # Short gaps in a slow signal: interpolate within each pump (at most 3 hours)
    n_missing = int(df["pump_vibration_mm_s"].isna().sum())
    df["pump_vibration_mm_s"] = df.groupby("tool_id")["pump_vibration_mm_s"].transform(
        lambda s: s.interpolate(limit=3, limit_area="inside"))
    report.append(("pump_telemetry", "vibration gaps filled by interpolation", n_missing - int(df["pump_vibration_mm_s"].isna().sum())))
    df = df.dropna(subset=["pump_vibration_mm_s"])
    df["day"] = (df["timestamp"] - START).dt.total_seconds() / 86400
    return df


def lot_subgroups(metro, n=5):
    """One row per lot: subgroup mean and range of etch rate (complete subgroups only)."""
    g = metro.groupby(["tool_id", "lot_id"])
    sub = g.agg(day=("day", "min"), n=("etch_rate_nm_min", "size"),
                xbar=("etch_rate_nm_min", "mean"),
                r=("etch_rate_nm_min", lambda s: s.max() - s.min()),
                nu_mean=("nonuniformity_pct", "mean"),
                pressure=("pressure_mTorr", "mean"),
                chuck_temp=("chuck_temp_C", "mean")).reset_index()
    return sub[sub["n"] == n].sort_values(["tool_id", "day"]).reset_index(drop=True)


def run():
    con = get_connection()
    report = []
    runs = load_wafer_runs(con, report)
    metro = load_metrology(con, report)
    pumps = load_pumps(con, report)
    subgroups = lot_subgroups(metro)
    sensors = daily_sensor_means(runs)
    dq = pd.DataFrame(report, columns=["table", "rule", "rows_affected"])
    dq = dq[dq["rows_affected"] > 0].reset_index(drop=True)
    return metro, pumps, subgroups, dq, sensors, len(runs)
