# Supplementary figure -- Spearman correlations between functional-guild relative
# abundances and soil enzyme activities, points coloured by vegetation. Reads the
# per-sample guild table (guilds + enzymes in one file).

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

import config as cfg

d = pd.read_excel(cfg.GUILD_TABLE)
d = d[~d["site"].astype(str).str.startswith("ctrl")]

pairs = [("F_ectomycorrhizal", "oxidase"), ("F_soil_saprotroph", "HYD"),
         ("B_Recalcitrant_C_decomp", "HYD"), ("F_litter_saprotroph", "BG"),
         ("B_Oligotrophic_C_processing", "PHOS"), ("F_ectomycorrhizal", "HYD")]

fig, axes = plt.subplots(2, 3, figsize=(11, 7))
for ax, (guild, enz) in zip(axes.ravel(), pairs):
    sub = d.dropna(subset=[guild, enz])
    for veg, c in cfg.veg_col.items():
        s = sub[sub["vegetation"] == veg]
        ax.scatter(s[guild], s[enz], c=c, s=25, edgecolor="white", linewidth=0.3)
    rho, p = spearmanr(sub[guild], sub[enz])
    ax.set_xlabel(f"{guild} (%)", fontsize=8)
    ax.set_ylabel(enz, fontsize=8)
    ax.set_title(f"rho = {rho:.2f}, p = {p:.3f}", fontsize=8)

fig.tight_layout()
fig.savefig("FigS_enzyme_correlations.png", dpi=300, bbox_inches="tight")
