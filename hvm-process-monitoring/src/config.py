"""Shared settings for the HVM process-monitoring pipeline.

All data in this project is simulated. The numbers below describe a generic
plasma-etch module so the analysis has realistic units and spec limits.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "outputs"
FIG_DIR = OUT_DIR / "figures"
DB_PATH = DATA_DIR / "fab_telemetry.db"

SEED = 2026

# Fab layout and schedule
TOOLS = ["ETCH01", "ETCH02", "ETCH03", "ETCH04"]
N_DAYS = 60
LOTS_PER_DAY = 8          # lots per chamber per day
WAFERS_PER_LOT = 25
SAMPLED_PER_LOT = 5       # wafers measured on metrology -> SPC subgroup size n = 5

# Nominal recipe (the "process of record")
RECIPE = {"rf_power_W": 1500.0, "pressure_mTorr": 30.0, "chuck_temp_C": 60.0}

# Critical-to-quality characteristic: etch rate (nm/min)
ER_TARGET = 250.0
ER_LSL = 235.0
ER_USL = 265.0
# Within-wafer non-uniformity (%), one-sided spec
NU_USL = 4.0

# SPC
BASELINE_DAYS = 15        # Phase I: limits are set from the first 15 days

# Pump vibration velocity thresholds (mm/s), in the spirit of ISO 10816 zones
VIB_ALERT = 4.5
VIB_TRIP = 7.1

# Colours (reference data-viz palette)
C_BLUE = "#2a78d6"
C_ORANGE = "#eb6834"
C_AQUA = "#1baf7a"
C_VIOLET = "#6250d6"
C_CRITICAL = "#d03b3b"
C_WARNING = "#c98500"
C_INK = "#0b0b0b"
C_INK2 = "#52514e"
C_MUTED = "#898781"
C_GRID = "#e1e0d9"
C_SURFACE = "#fcfcfb"
TOOL_COLOURS = dict(zip(TOOLS, [C_BLUE, C_ORANGE, C_AQUA, C_VIOLET]))
