# Supplementary figure -- variance partitioning of community composition
# (Site x Season x Year, db-RDA / varpart in R). The three predictors are drawn as
# overlapping circles and the adjusted-R2 fractions written into each region.
# Values below are taken from the R output (supp_table_elev_vs_veget.xlsx).

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

panels = {
    "Bacteria": dict(site=0.627, season=0.001, year=0.000, total=0.619, resid=0.381),
    "Fungi":    dict(site=0.504, season=0.001, year=0.001, total=0.498, resid=0.502),
}
colour = {"Site": "#67B9A3", "Season": "#D98A8A", "Year": "#7C93C5"}
pos = {"Site": (-0.5, 0.35), "Season": (0.5, 0.35), "Year": (0.0, -0.5)}

fig, axes = plt.subplots(1, 2, figsize=(11, 6))
for ax, (name, d) in zip(axes, panels.items()):
    for k, (x, y) in pos.items():
        ax.add_patch(Circle((x, y), 1.0, facecolor=colour[k], edgecolor="#555555", lw=1.1, alpha=0.5))
    ax.text(-1.05, 0.35, f"[a]\n{d['site']:.3f} ***", ha="center", va="center", fontsize=9)
    ax.text(1.05, 0.35, f"[b]\n{d['season']:.3f} ns", ha="center", va="center", fontsize=9)
    ax.text(0.0, -0.95, f"[c]\n{d['year']:.3f} ns", ha="center", va="center", fontsize=9)
    for xy in [(0, 0.72), (-0.42, -0.15), (0.42, -0.15), (0, 0.12)]:
        ax.text(*xy, "0.000", ha="center", va="center", fontsize=8.5)
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.0, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(name, fontweight="bold")

fig.tight_layout()
fig.savefig("FigS_varpart.png", dpi=300, bbox_inches="tight")
