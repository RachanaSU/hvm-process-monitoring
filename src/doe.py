"""Design of Experiments: tune the etch recipe for etch rate on target and lower non-uniformity.

* 2^3 full factorial on RF power, pressure and chuck temperature, 3 replicates,
  plus 4 centre points (28 test wafers), run on the simulated chamber.
* Effects estimated by least squares on coded factors (-1 / +1), with
  t-tests using the replicate error. Centre points test for curvature.
* The fitted models are searched for the setting that keeps etch rate within
  +/-1 nm/min of target while minimising non-uniformity.
* Confirmation run: 40 lots at the old and new recipe, compared on Cpk and yield.
"""
import itertools

import numpy as np
import pandas as pd
from scipy import stats

from . import process_model as pm
from .config import (C_INK2, C_MUTED, ER_LSL, ER_TARGET, ER_USL, NU_USL, RECIPE,
                     SEED, C_BLUE, C_ORANGE)
from .plotting import plt, save

FACTORS = ["rf_power_W", "pressure_mTorr", "chuck_temp_C"]
LEVELS = {"rf_power_W": (1400, 1600), "pressure_mTorr": (28, 32), "chuck_temp_C": (55, 65)}
LABELS = {"rf_power_W": "RF power", "pressure_mTorr": "Pressure", "chuck_temp_C": "Chuck temp"}
REPLICATES = 3
CENTRE_POINTS = 4


def to_physical(x):
    return {f: RECIPE[f] + x[i] * (LEVELS[f][1] - LEVELS[f][0]) / 2 for i, f in enumerate(FACTORS)}


def run_wafer(rng, settings):
    lot = rng.normal(0, pm.ER_LOT_SD)
    er = pm.mean_etch_rate(**settings) + lot + rng.normal(0, pm.ER_WAFER_SD)
    nu = pm.mean_nonuniformity(**settings) + rng.normal(0, pm.NU_SD)
    return er, nu


def design():
    runs = [x for x in itertools.product([-1, 1], repeat=3) for _ in range(REPLICATES)]
    runs += [(0, 0, 0)] * CENTRE_POINTS
    return np.array(runs, dtype=float)


def model_matrix(X):
    a, b, c = X[:, 0], X[:, 1], X[:, 2]
    cols = [np.ones(len(X)), a, b, c, a * b, a * c, b * c, a * b * c]
    names = ["Intercept", "RF", "P", "T", "RF x P", "RF x T", "P x T", "RF x P x T"]
    return np.column_stack(cols), names


def fit(X, y):
    M, names = model_matrix(X)
    beta, *_ = np.linalg.lstsq(M, y, rcond=None)
    resid = y - M @ beta
    dof = len(y) - M.shape[1]
    mse = resid @ resid / dof
    se = np.sqrt(np.diag(mse * np.linalg.inv(M.T @ M)))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    table = pd.DataFrame({"term": names, "coef": beta, "effect": 2 * beta, "t": t, "p_value": p})
    table.loc[0, "effect"] = np.nan
    return beta, table


def curvature_test(X, y):
    """Compare the centre-point mean with the factorial mean (pure-error t-test)."""
    centre = np.all(X == 0, axis=1)
    fac, ctr = y[~centre], y[centre]
    diff = ctr.mean() - fac.mean()
    # pooled pure error from replicate groups
    keys = ["".join(f"{v:+.0f}" for v in r) for r in X[~centre]]
    groups = pd.Series(y[~centre]).groupby(keys)
    ss = sum(((g - g.mean()) ** 2).sum() for _, g in groups) + ((ctr - ctr.mean()) ** 2).sum()
    dof = sum(len(g) - 1 for _, g in groups) + len(ctr) - 1
    s = np.sqrt(ss / dof)
    t = diff / (s * np.sqrt(1 / len(fac) + 1 / len(ctr)))
    return float(diff), float(2 * stats.t.sf(abs(t), dof))


def predict(beta, x):
    M, _ = model_matrix(np.atleast_2d(x))
    return M @ beta


def optimise(b_er, b_nu):
    grid = np.linspace(-1, 1, 41)
    pts = np.array(list(itertools.product(grid, grid, grid)))
    er, nu = predict(b_er, pts), predict(b_nu, pts)
    ok = np.abs(er - ER_TARGET) <= 1.0
    i = np.argmin(np.where(ok, nu, np.inf))
    return pts[i], float(er[i]), float(nu[i])


def confirmation(rng, settings, n_lots=40, wafers=5):
    er, nu = [], []
    for _ in range(n_lots):
        lot = rng.normal(0, pm.ER_LOT_SD)
        for _ in range(wafers):
            er.append(pm.mean_etch_rate(**settings) + lot + rng.normal(0, pm.ER_WAFER_SD))
            nu.append(pm.mean_nonuniformity(**settings) + rng.normal(0, pm.NU_SD))
    er, nu = np.array(er), np.array(nu)
    s_er, s_nu = er.std(ddof=1), nu.std(ddof=1)
    return dict(er_mean=er.mean(), nu_mean=nu.mean(),
                ppk_er=min(ER_USL - er.mean(), er.mean() - ER_LSL) / (3 * s_er),
                ppk_nu=(NU_USL - nu.mean()) / (3 * s_nu),
                yield_pct=100 * np.mean((er >= ER_LSL) & (er <= ER_USL) & (nu <= NU_USL)),
                n_wafers=len(er))


def run():
    rng = np.random.default_rng(SEED + 1)
    X = design()
    results = [run_wafer(rng, to_physical(x)) for x in X]
    er = np.array([r[0] for r in results])
    nu = np.array([r[1] for r in results])
    runs = pd.DataFrame(X, columns=["x_rf", "x_p", "x_t"])
    for f in FACTORS:
        runs[f] = [to_physical(x)[f] for x in X]
    runs["etch_rate_nm_min"], runs["nonuniformity_pct"] = er.round(2), nu.round(3)

    fac = ~np.all(X == 0, axis=1)
    b_er, t_er = fit(X[fac], er[fac])
    b_nu, t_nu = fit(X[fac], nu[fac])
    curv = {"etch_rate": curvature_test(X, er), "nonuniformity": curvature_test(X, nu)}

    x_opt, er_pred, nu_pred = optimise(b_er, b_nu)
    new_recipe = {f: round(float(v), 1) for f, v in to_physical(x_opt).items()}
    rng_c = np.random.default_rng(SEED + 2)
    before = confirmation(rng_c, dict(RECIPE))
    first = confirmation(rng_c, new_recipe)

    # The first confirmation run landed off target, so re-centre etch rate with RF
    # power (the strongest etch-rate knob), using the fitted local slope.
    a, b, c = x_opt
    slope = b_er[1] + b_er[4] * b + b_er[5] * c + b_er[7] * b * c      # d(ER)/d(x_rf)
    x_final = np.array([np.clip(a - (first["er_mean"] - ER_TARGET) / slope, -1, 1), b, c])
    final_recipe = {f: round(float(v), 1) for f, v in to_physical(x_final).items()}
    after = confirmation(rng_c, final_recipe)
    return dict(runs=runs, effects_er=t_er, effects_nu=t_nu, curvature=curv,
                x_opt=x_opt, doe_recipe=new_recipe, pred_er=er_pred, pred_nu=nu_pred,
                x_final=x_final, new_recipe=final_recipe,
                before=before, first_confirmation=first, after=after, b_nu=b_nu)


def figures(res):
    paths = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=False)
    for ax, key, title in ((axes[0], "effects_er", "Etch rate"), (axes[1], "effects_nu", "Non-uniformity")):
        t = res[key].iloc[1:].copy()
        t["abs_t"] = t["t"].abs()
        t = t.sort_values("abs_t")
        crit = stats.t.ppf(0.975, 24 - 8)
        colours = [C_BLUE if v > crit else "#c3c2b7" for v in t["abs_t"]]
        ax.barh(t["term"], t["abs_t"], color=colours, height=0.6)
        ax.axvline(crit, color=C_INK2, lw=1, ls="--")
        ax.text(crit, len(t) - 0.45, f" p = 0.05 (|t| = {crit:.2f})", color=C_INK2, fontsize=7.5)
        ax.set_title(f"{title}: standardised effects (|t|)")
        ax.set_xlabel("|t|")
        ax.grid(axis="y", visible=False)
    fig.tight_layout()
    paths.append(save(fig, "doe_pareto_effects.png"))

    # RF x Temp interaction on non-uniformity (pressure held at its optimum)
    fig, ax = plt.subplots(figsize=(6, 3.6))
    xs = np.array([-1, 1])
    for xt, colour, name in ((-1, C_BLUE, "chuck temp 55 C"), (1, C_ORANGE, "chuck temp 65 C")):
        y = [predict(res["b_nu"], [x, res["x_opt"][1], xt])[0] for x in xs]
        ax.plot([1400, 1600], y, color=colour, marker="o", ms=6, label=name)
        ax.text(1605, y[1], name, color=C_INK2, fontsize=8, va="center")
    ax.set_xlabel("RF power (W)")
    ax.set_ylabel("predicted non-uniformity (%)")
    ax.set_xlim(1380, 1680)
    ax.set_title("RF power x chuck temperature interaction")
    fig.tight_layout()
    paths.append(save(fig, "doe_interaction_rf_temp.png"))
    return paths
