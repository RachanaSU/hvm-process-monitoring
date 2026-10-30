"""'True' behaviour of the simulated etch process.

The simulator uses these functions to generate data. The analysis modules never
import them: SPC, ANOVA, reliability and DOE all work only from measured data,
just as an engineer in a fab would.
"""
import numpy as np

from .config import RECIPE


def coded(rf_power_W, pressure_mTorr, chuck_temp_C):
    """Convert physical settings to coded units (-1 .. +1 across the DOE range)."""
    x_rf = (np.asarray(rf_power_W) - RECIPE["rf_power_W"]) / 100.0
    x_p = (np.asarray(pressure_mTorr) - RECIPE["pressure_mTorr"]) / 2.0
    x_t = (np.asarray(chuck_temp_C) - RECIPE["chuck_temp_C"]) / 5.0
    return x_rf, x_p, x_t


def mean_etch_rate(rf_power_W, pressure_mTorr, chuck_temp_C):
    """Expected etch rate (nm/min) for given chamber conditions."""
    x_rf, x_p, x_t = coded(rf_power_W, pressure_mTorr, chuck_temp_C)
    return 250.0 + 9.0 * x_rf + 4.0 * x_p - 3.0 * x_t + 1.5 * x_rf * x_p


def mean_nonuniformity(rf_power_W, pressure_mTorr, chuck_temp_C):
    """Expected within-wafer non-uniformity (%, lower is better)."""
    x_rf, x_p, x_t = coded(rf_power_W, pressure_mTorr, chuck_temp_C)
    return 3.4 + 0.35 * x_rf - 0.45 * x_p + 0.30 * x_t + 0.60 * x_rf * x_t


ER_WAFER_SD = 3.0   # wafer-to-wafer noise (nm/min)
ER_LOT_SD = 1.5     # lot-to-lot noise (nm/min)
NU_SD = 0.30        # wafer-to-wafer noise on non-uniformity (%)
