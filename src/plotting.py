"""Common chart styling: thin marks, recessive grid, text in ink colours."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config import C_GRID, C_INK, C_INK2, C_MUTED, C_SURFACE, FIG_DIR

plt.rcParams.update({
    "figure.facecolor": C_SURFACE,
    "axes.facecolor": C_SURFACE,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": C_INK2,
    "axes.titlecolor": C_INK,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": C_GRID,
    "grid.linewidth": 0.6,
    "xtick.color": C_MUTED,
    "ytick.color": C_MUTED,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "lines.linewidth": 1.6,
    "font.family": "DejaVu Sans",
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
})


def save(fig, name):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / name
    fig.savefig(path)
    plt.close(fig)
    return path
