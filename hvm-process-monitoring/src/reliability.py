"""Vacuum pump degradation model and predictive-maintenance forecast.

Model (same as matlab/pump_degradation.m):

    dD/dt = k * u(t) * (1 + D / Dc)      damage D grows faster as the pump wears
    V(t)  = V0 + D(t)                    measured vibration velocity (mm/s)

u(t) is the pump duty (load) cycle. k, Dc and V0 are fitted to the vibration
history by least squares, the ODE is integrated forward, and the forecast
crossing of the trip limit gives the predicted failure time. Because the ODE
is linear in (1 + D/Dc) it also has a closed form, used for the bootstrap:

    D(t) = Dc * (exp(k/Dc * U(t)) - 1),   U(t) = integral of u from 0 to t
"""
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .config import (C_CRITICAL, C_INK2, C_MUTED, C_WARNING, N_DAYS,
                     TOOL_COLOURS, TOOLS, VIB_ALERT, VIB_TRIP)
from .plotting import plt, save
from .simulate import pump_load

FORECAST_DAY = 40          # "today" for the headline forecast
HORIZON_DAYS = 90         # planning horizon: beyond this, no action is scheduled


def load_integral(t):
    """U(t) = integral of pump_load from 0 to t (pump_load = 0.85 + 0.15 sin 2*pi*t)."""
    t = np.asarray(t, dtype=float)
    return 0.85 * t + 0.15 * (1 - np.cos(2 * np.pi * t)) / (2 * np.pi)


def vib_closed_form(params, t):
    k, dc, v0 = params
    return v0 + dc * (np.exp(k / dc * load_integral(t)) - 1)


def vib_ode(params, t):
    k, dc, v0 = params
    sol = solve_ivp(lambda tt, d: k * pump_load(tt) * (1 + d / dc), (0, float(t[-1])), [0.0],
                    t_eval=t, rtol=1e-8, atol=1e-10)
    return v0 + sol.y[0]


def fit(t, v):
    res = least_squares(lambda p: vib_closed_form(p, t) - v, x0=[0.01, 5.0, float(np.median(v[:48]))],
                        bounds=([1e-5, 0.05, 0.0], [1.0, 500.0, 10.0]))
    return res.x, res.fun


def effective_trip(rmse):
    """The protective interlock trips on the first single reading over 7.1 mm/s, so the
    pump is at risk once the *mean* vibration is within two noise SDs of the limit."""
    return VIB_TRIP - 2 * rmse


def crossing_day(params, threshold, t_from):
    t = np.arange(t_from, t_from + HORIZON_DAYS, 1 / 24)
    v = vib_closed_form(params, t)
    idx = np.argmax(v >= threshold)
    return float(t[idx]) if v[idx] >= threshold else None


def first_life(pumps, tool):
    """Telemetry for the original pump only (stop at the first trip)."""
    p = pumps[pumps["tool_id"] == tool].sort_values("day")
    trip = p[p["pump_vibration_mm_s"] >= VIB_TRIP]
    actual = float(trip["day"].iloc[0]) if len(trip) else None
    if actual is not None:
        p = p[p["day"] < actual]
    return p, actual


def bootstrap_trip(params, t, v, resid, threshold, n=300, seed=7):
    """Residual bootstrap: refit on resampled noise to get a range for the failure day."""
    rng = np.random.default_rng(seed)
    fitted = vib_closed_form(params, t)
    days = []
    for _ in range(n):
        v_star = fitted + rng.choice(resid, size=len(resid), replace=True)
        p_star, _ = fit(t, v_star)
        d = crossing_day(p_star, threshold, t[-1])
        days.append(d if d is not None else np.inf)
    return np.array(days)


def analyse(pumps):
    rows, detail = [], {}
    for tool in TOOLS:
        p, actual = first_life(pumps, tool)
        train = p[p["day"] <= FORECAST_DAY]
        t, v = train["day"].values, train["pump_vibration_mm_s"].values
        params, resid = fit(t, v)
        # Check the closed form against a numerical ODE solve (ode45-style)
        ode_gap = float(np.max(np.abs(vib_ode(params, t) - vib_closed_form(params, t))))
        rmse = float(np.sqrt(np.mean(resid ** 2)))
        trip_level = effective_trip(rmse)
        pred_trip = crossing_day(params, trip_level, t[-1])
        pred_alert = crossing_day(params, VIB_ALERT, 0)
        row = dict(tool_id=tool, k=params[0], Dc=params[1], V0=params[2],
                   rmse=rmse, trip_level_used=trip_level, ode_check_max_gap=ode_gap,
                   current_vib=float(np.mean(v[-24:])),
                   predicted_trip_day=pred_trip, actual_trip_day=actual,
                   predicted_alert_day=pred_alert)
        if pred_trip is not None:
            boot = bootstrap_trip(params, t, v, resid, trip_level)
            row["trip_p05"], row["trip_p50"], row["trip_p95"] = np.percentile(boot, [5, 50, 95])
            row["recommended_pm_day"] = float(np.floor(row["trip_p05"] - 1))
            # Naive comparison: straight-line trend over the last 7 days
            last = train[train["day"] > FORECAST_DAY - 7]
            slope, icpt = np.polyfit(last["day"], last["pump_vibration_mm_s"], 1)
            row["linear_trend_trip_day"] = float((trip_level - icpt) / slope) if slope > 0 else None
            detail[tool] = dict(params=params, train=train, life=p, boot=boot)
        rows.append(row)
    summary = pd.DataFrame(rows)
    rolling = rolling_forecast(pumps, "ETCH04")
    return summary, detail, rolling


def rolling_forecast(pumps, tool, start=20):
    """Refit every day: how does the predicted failure day evolve as data arrives?"""
    p, actual = first_life(pumps, tool)
    out = []
    for day in range(start, int(np.floor(actual)) + 1):
        train = p[p["day"] <= day]
        params, resid = fit(train["day"].values, train["pump_vibration_mm_s"].values)
        level = effective_trip(np.sqrt(np.mean(resid ** 2)))
        pred = crossing_day(params, level, day)
        last = train[train["day"] > day - 7]
        slope, icpt = np.polyfit(last["day"], last["pump_vibration_mm_s"], 1)
        lin = (level - icpt) / slope if slope > 0 else np.nan
        out.append(dict(forecast_made_day=day, predicted_trip_day=pred, linear_trend_trip_day=lin,
                        actual_trip_day=actual))
    return pd.DataFrame(out)


def figures(pumps, summary, detail, rolling):
    paths = []
    tool = "ETCH04"
    d = detail[tool]
    row = summary.set_index("tool_id").loc[tool]
    full = pumps[pumps["tool_id"] == tool]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.plot(full["day"], full["pump_vibration_mm_s"], color="#b7d3f6", lw=0.8, label="measured vibration (hourly)")
    ax.plot(d["train"]["day"], d["train"]["pump_vibration_mm_s"], color=TOOL_COLOURS[tool], lw=0.8,
            label=f"data used for fit (up to day {FORECAST_DAY})")
    t_fc = np.arange(0, row["trip_p95"] + 2, 1 / 24)
    ax.plot(t_fc, vib_closed_form(d["params"], t_fc), color="#0b0b0b", lw=1.6, ls="--", label="fitted wear model + forecast")
    for y, c, name in ((VIB_ALERT, C_WARNING, "alert 4.5"), (VIB_TRIP, C_CRITICAL, "trip 7.1")):
        ax.axhline(y, color=c, lw=1, ls=":")
        ax.text(0.5, y + 0.12, f"{name} mm/s", color=c, fontsize=8)
    ax.axvspan(row["trip_p05"], row["trip_p95"], color="#f0efec", zorder=0)
    ax.axvline(row["actual_trip_day"], color=C_CRITICAL, lw=1.2)
    ax.text(row["actual_trip_day"] + 0.3, 1.2, f"actual trip\nday {row['actual_trip_day']:.1f}", color=C_CRITICAL, fontsize=8)
    ax.axvline(row["recommended_pm_day"], color=C_INK2, lw=1.2)
    ax.text(row["recommended_pm_day"] - 0.3, 1.2, f"planned PM\nday {row['recommended_pm_day']:.0f}",
            color=C_INK2, fontsize=8, ha="right")
    ax.axvline(FORECAST_DAY, color=C_MUTED, lw=0.8, ls="--")
    ax.text(FORECAST_DAY - 0.3, 7.9, f"forecast made\non day {FORECAST_DAY}", color=C_MUTED, fontsize=8, ha="right")
    ax.set_xlim(0, N_DAYS)
    ax.set_ylim(1, 8.6)
    ax.set_xlabel("day")
    ax.set_ylabel("pump vibration (mm/s)")
    ax.set_title(f"{tool} dry pump: wear model forecast vs actual failure "
                 f"(90% band: days {row['trip_p05']:.1f}-{row['trip_p95']:.1f})")
    ax.legend(loc="upper left")
    fig.tight_layout()
    paths.append(save(fig, "pump_forecast_ETCH04.png"))

    fig, ax = plt.subplots(figsize=(8, 3.8))
    ax.plot(rolling["forecast_made_day"], rolling["predicted_trip_day"], color="#0b0b0b", marker="o", ms=3,
            label="wear model (ODE) forecast")
    ax.plot(rolling["forecast_made_day"], rolling["linear_trend_trip_day"], color=TOOL_COLOURS["ETCH02"], marker="s", ms=3,
            label="7-day linear trend forecast")
    ax.axhline(rolling["actual_trip_day"].iloc[0], color=C_CRITICAL, lw=1.2)
    ax.text(rolling["forecast_made_day"].max(), rolling["actual_trip_day"].iloc[0] - 3.5,
            f"actual failure, day {rolling['actual_trip_day'].iloc[0]:.1f}", color=C_CRITICAL, fontsize=8, ha="right")
    ax.set_ylim(40, 85)
    ax.set_xlabel("day the forecast was made")
    ax.set_ylabel("predicted failure day")
    ax.set_title("Forecast stability: wear model vs straight-line trend, ETCH04 pump")
    ax.legend(loc="upper right")
    fig.tight_layout()
    paths.append(save(fig, "pump_forecast_convergence.png"))
    return paths
