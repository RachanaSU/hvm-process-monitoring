"""Statistical process control on etch rate.

* X-bar and R charts per chamber, subgroup = 5 measured wafers per lot.
* Phase I: control limits from the first BASELINE_DAYS of data.
* Phase II: every later lot is checked against four Western Electric rules.
* Process capability (Cp, Cpk, Ppk) and in-spec yield, baseline vs last 10 days.
* One-way ANOVA + Tukey HSD: do the chambers still match each other?
"""
import numpy as np
import pandas as pd
from scipy import stats

from .config import (BASELINE_DAYS, C_CRITICAL, C_INK2, C_MUTED, ER_LSL,
                     ER_TARGET, ER_USL, N_DAYS, NU_USL, TOOL_COLOURS, TOOLS)
from .plotting import plt, save

# Control chart constants for subgroup size n = 5
A2, D3, D4, D2 = 0.577, 0.0, 2.114, 2.326
RECENT_DAYS = 10


def limits(baseline):
    """Three-way chart limits (X-bar via moving range of lot means, plus R chart).

    Etch rate varies lot to lot as well as wafer to wafer. Classic X-bar limits
    (X-double-bar +/- A2*R-bar) only see the within-lot spread, so they are too
    tight and a healthy chamber raises a false alarm every couple of days. The
    X-bar limits here use the moving range of consecutive lot means instead,
    which includes the lot-to-lot component. The R chart still watches the
    within-lot spread.
    """
    xbb, rbar = baseline["xbar"].mean(), baseline["r"].mean()
    mr_bar = baseline["xbar"].diff().abs().mean()
    sigma_xbar = mr_bar / 1.128
    return dict(center=xbb, ucl=xbb + 3 * sigma_xbar, lcl=xbb - 3 * sigma_xbar,
                sigma_xbar=sigma_xbar, rbar=rbar, r_ucl=D4 * rbar, r_lcl=D3 * rbar,
                sigma_within=rbar / D2,
                classic_ucl=xbb + A2 * rbar, classic_lcl=xbb - A2 * rbar)


def classic_false_alarms(sg, lim):
    """Signals a healthy chamber would raise with textbook A2*R-bar limits (for comparison)."""
    classic = dict(lim, ucl=lim["classic_ucl"], lcl=lim["classic_lcl"], sigma_xbar=A2 * lim["rbar"] / 3)
    return sum(1 for f in western_electric(sg["xbar"], classic) if f)


def western_electric(x, lim):
    """Return a list of rule labels (or '') for each point in x."""
    z = (np.asarray(x) - lim["center"]) / lim["sigma_xbar"]
    flags = [""] * len(z)
    for i in range(len(z)):
        hits = []
        if abs(z[i]) > 3:
            hits.append("R1")                                    # 1 point beyond 3 sigma
        # Rules 2 and 3 flag the point that completes the pattern, so it must be on that side too
        w = z[max(0, i - 2): i + 1]
        if len(w) == 3 and ((z[i] > 2 and (w > 2).sum() >= 2) or (z[i] < -2 and (w < -2).sum() >= 2)):
            hits.append("R2")                                    # 2 of 3 beyond 2 sigma, same side
        w = z[max(0, i - 4): i + 1]
        if len(w) == 5 and ((z[i] > 1 and (w > 1).sum() >= 4) or (z[i] < -1 and (w < -1).sum() >= 4)):
            hits.append("R3")                                    # 4 of 5 beyond 1 sigma, same side
        w = z[max(0, i - 7): i + 1]
        if len(w) == 8 and (np.all(w > 0) or np.all(w < 0)):
            hits.append("R4")                                    # 8 in a row on one side
        flags[i] = ",".join(hits)
    return flags


def capability(values, sigma_within):
    mu, s_overall = values.mean(), values.std(ddof=1)
    return dict(mean=mu,
                cp=(ER_USL - ER_LSL) / (6 * sigma_within),
                cpk=min(ER_USL - mu, mu - ER_LSL) / (3 * sigma_within),
                ppk=min(ER_USL - mu, mu - ER_LSL) / (3 * s_overall))


def analyse(metro, subgroups):
    charts, cap_rows, signals = {}, [], []
    for tool in TOOLS:
        sg = subgroups[subgroups["tool_id"] == tool].copy()
        lim = limits(sg[sg["day"] < BASELINE_DAYS])
        sg["rule"] = western_electric(sg["xbar"], lim)
        sg["r_ooc"] = sg["r"] > lim["r_ucl"]
        charts[tool] = (sg, lim)

        mon = sg[(sg["day"] >= BASELINE_DAYS) & (sg["rule"] != "")]
        for _, row in mon.iterrows():
            signals.append((tool, row["lot_id"], round(row["day"], 2), round(row["xbar"], 2), row["rule"]))

        w = metro[metro["tool_id"] == tool]
        base_w = w[w["day"] < BASELINE_DAYS]
        rec_w = w[w["day"] >= N_DAYS - RECENT_DAYS]
        rec_sg = sg[sg["day"] >= N_DAYS - RECENT_DAYS]
        for period, wafers, sigma in (("baseline (days 0-15)", base_w, lim["sigma_within"]),
                                      (f"last {RECENT_DAYS} days", rec_w, rec_sg["r"].mean() / D2)):
            c = capability(wafers["etch_rate_nm_min"], sigma)
            in_spec = (wafers["etch_rate_nm_min"].between(ER_LSL, ER_USL)
                       & (wafers["nonuniformity_pct"].fillna(0) <= NU_USL))
            cap_rows.append(dict(tool_id=tool, period=period, n_wafers=len(wafers),
                                 mean_er=round(c["mean"], 2), cp=round(c["cp"], 2),
                                 cpk=round(c["cpk"], 2), ppk=round(c["ppk"], 2),
                                 yield_pct=round(100 * in_spec.mean(), 2)))

    signals = pd.DataFrame(signals, columns=["tool_id", "lot_id", "day", "xbar", "rules"])
    cap = pd.DataFrame(cap_rows)
    return charts, signals, cap


def first_detection(signals, event_day):
    """Days from the injected event to the first SPC signal on that chamber."""
    out = {}
    for tool, day0 in event_day.items():
        s = signals[(signals["tool_id"] == tool) & (signals["day"] >= day0)]
        before = signals[(signals["tool_id"] == tool) & (signals["day"] < day0)]
        out[tool] = dict(event_day=day0,
                         first_signal_day=float(s["day"].min()) if len(s) else None,
                         detection_lag_days=round(float(s["day"].min()) - day0, 2) if len(s) else None,
                         first_rules=s.iloc[0]["rules"] if len(s) else None,
                         signals_before_event=int(len(before)))
    return out


def anova(metro, start_day, end_day):
    w = metro[(metro["day"] >= start_day) & (metro["day"] < end_day)]
    groups = [w.loc[w["tool_id"] == t, "etch_rate_nm_min"].values for t in TOOLS]
    f, p = stats.f_oneway(*groups)
    res = stats.tukey_hsd(*groups)
    pairs = []
    for i in range(len(TOOLS)):
        for j in range(i + 1, len(TOOLS)):
            pairs.append(dict(pair=f"{TOOLS[i]} vs {TOOLS[j]}",
                              mean_diff=round(groups[i].mean() - groups[j].mean(), 2),
                              p_value=float(res.pvalue[i, j]),
                              significant=bool(res.pvalue[i, j] < 0.05)))
    return dict(window=f"days {start_day}-{end_day}", F=round(float(f), 2), p_value=float(p),
                means={t: round(float(g.mean()), 2) for t, g in zip(TOOLS, groups)}), pd.DataFrame(pairs)


# ----------------------------------------------------------------- figures
def _xbar_panel(ax, sg, lim, tool, show_spec=True):
    ax.plot(sg["day"], sg["xbar"], color=TOOL_COLOURS[tool], lw=1.0, marker="o", ms=2.2)
    for y, ls in ((lim["center"], "-"), (lim["ucl"], "--"), (lim["lcl"], "--")):
        ax.axhline(y, color=C_INK2, lw=0.9, ls=ls)
    ooc = sg[(sg["rule"] != "") & (sg["day"] >= BASELINE_DAYS)]
    ax.scatter(ooc["day"], ooc["xbar"], s=26, facecolor="none", edgecolor=C_CRITICAL, lw=1.2, zorder=3,
               label="Western Electric signal")
    ax.axvspan(0, BASELINE_DAYS, color="#f0efec", zorder=0)
    ax.text(0.5, 1.0, "Phase I baseline", transform=ax.get_xaxis_transform(), fontsize=7.5,
            color=C_MUTED, va="top")
    ax.text(N_DAYS + 0.5, lim["ucl"], "UCL", fontsize=7.5, color=C_INK2, va="center")
    ax.text(N_DAYS + 0.5, lim["lcl"], "LCL", fontsize=7.5, color=C_INK2, va="center")
    ax.set_xlim(0, N_DAYS)


def figures(charts, metro, events):
    paths = []
    fig, axes = plt.subplots(4, 1, figsize=(9, 9), sharex=True, sharey=True)
    for ax, tool in zip(axes, TOOLS):
        sg, lim = charts[tool]
        _xbar_panel(ax, sg, lim, tool)
        ax.set_title(f"{tool}  ·  {events.get(tool, 'no injected fault')}", fontsize=9.5)
        ax.set_ylabel("lot mean ER\n(nm/min)")
    axes[0].legend(loc="upper right")
    axes[-1].set_xlabel("day")
    fig.suptitle("X-bar charts of etch rate, all chambers (limits from Phase I baseline)",
                 x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    paths.append(save(fig, "spc_overview.png"))

    for tool in ("ETCH02", "ETCH03"):
        sg, lim = charts[tool]
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 5.2), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
        _xbar_panel(a1, sg, lim, tool)
        a1.set_ylabel("X-bar (nm/min)")
        a1.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), borderaxespad=0.2)
        a1.set_title(f"{tool}: X-bar and R chart, etch rate  ({events[tool]})", pad=22)
        a2.plot(sg["day"], sg["r"], color=TOOL_COLOURS[tool], lw=1.0, marker="o", ms=2.2)
        a2.axhline(lim["rbar"], color=C_INK2, lw=0.9)
        a2.axhline(lim["r_ucl"], color=C_INK2, lw=0.9, ls="--")
        a2.text(N_DAYS + 0.5, lim["r_ucl"], "UCL", fontsize=7.5, color=C_INK2, va="center")
        a2.set_ylabel("R (nm/min)")
        a2.set_xlabel("day")
        fig.tight_layout()
        paths.append(save(fig, f"xbar_r_{tool}.png"))
    return paths


def root_cause_figure(sensors):
    """Trace the two SPC signals to chamber sensors: ETCH02 pressure, ETCH03 chuck temperature."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    for ax, col, tool, unit, event in ((axes[0], "pressure_mTorr", "ETCH02", "mTorr", 32),
                                       (axes[1], "chuck_temp_C", "ETCH03", "C", 20)):
        for t in TOOLS:
            s = sensors[sensors["tool_id"] == t]
            hl = t == tool
            ax.plot(s["day"], s[col], color=TOOL_COLOURS[t] if hl else "#c3c2b7",
                    lw=1.8 if hl else 1.0, zorder=3 if hl else 2, label=t if hl else None)
        ax.axvline(event, color=C_MUTED, lw=0.8, ls="--")
        ax.text(event + 0.5, ax.get_ylim()[1], f"day {event}", color=C_MUTED, fontsize=8, va="top")
        ax.set_xlabel("day")
        ax.set_ylabel(f"daily mean ({unit})")
        ax.set_title(f"{tool} {col.split('_')[0] if col.startswith('pressure') else 'chuck temperature'} "
                     f"vs other chambers (grey)", fontsize=10)
    fig.tight_layout()
    return save(fig, "root_cause_sensors.png")


def capability_figure(cap):
    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    base = cap[cap["period"].str.startswith("baseline")].set_index("tool_id")
    rec = cap[~cap["period"].str.startswith("baseline")].set_index("tool_id")
    y = np.arange(len(TOOLS))[::-1]
    for yi, tool in zip(y, TOOLS):
        ax.plot([base.loc[tool, "cpk"], rec.loc[tool, "cpk"]], [yi, yi], color="#c3c2b7", lw=2, zorder=1)
        ax.scatter(base.loc[tool, "cpk"], yi, s=60, color="#86b6ef", zorder=2)
        ax.scatter(rec.loc[tool, "cpk"], yi, s=60, color=TOOL_COLOURS[tool], zorder=3)
        ax.text(rec.loc[tool, "cpk"], yi + 0.22, f"{rec.loc[tool, 'cpk']:.2f}", ha="center", fontsize=8, color=C_INK2)
    ax.axvline(1.33, color=C_CRITICAL, lw=1, ls="--")
    ax.text(1.34, y.max() + 0.45, "Cpk 1.33 target", color=C_CRITICAL, fontsize=8)
    ax.set_yticks(y, TOOLS)
    ax.set_xlabel("Cpk (etch rate, spec 235-265 nm/min)")
    ax.set_ylim(-0.6, y.max() + 0.8)
    ax.scatter([], [], s=50, color="#86b6ef", label="baseline (days 0-15)")
    ax.scatter([], [], s=50, color=C_INK2, label=f"last {RECENT_DAYS} days")
    ax.legend(loc="upper left")
    ax.set_title("Process capability: baseline vs last 10 days")
    fig.tight_layout()
    return save(fig, "capability_cpk.png")


def anova_figure(metro, start_day, end_day, result):
    w = metro[(metro["day"] >= start_day) & (metro["day"] < end_day)]
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    data = [w.loc[w["tool_id"] == t, "etch_rate_nm_min"] for t in TOOLS]
    bp = ax.boxplot(data, tick_labels=TOOLS, widths=0.45, patch_artist=True, showfliers=False,
                    medianprops=dict(color="#0b0b0b", lw=1.4))
    for patch, tool in zip(bp["boxes"], TOOLS):
        patch.set(facecolor=TOOL_COLOURS[tool], alpha=0.35, edgecolor=TOOL_COLOURS[tool])
    for y in (ER_LSL, ER_USL):
        ax.axhline(y, color=C_CRITICAL, lw=1, ls="--")
    ax.axhline(ER_TARGET, color=C_MUTED, lw=0.8)
    ax.text(4.55, ER_USL, "USL", color=C_CRITICAL, fontsize=8, va="center")
    ax.text(4.55, ER_LSL, "LSL", color=C_CRITICAL, fontsize=8, va="center")
    ax.set_ylabel("etch rate (nm/min)")
    p_txt = "p < 0.001" if result["p_value"] < 0.001 else f"p = {result['p_value']:.3f}"
    ax.set_title(f"Chamber matching, days {start_day}-{end_day}: one-way ANOVA "
                 f"F = {result['F']:.1f}, {p_txt}")
    fig.tight_layout()
    return save(fig, "anova_chamber_matching.png")
