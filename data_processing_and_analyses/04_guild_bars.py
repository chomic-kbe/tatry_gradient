# Guild composition of the dwarf-pine community, split by sharing category.
# Sharing categories come from the SH/ASV overlap (>= 2 samples per zone);
# the guild make-up of each category is computed from the per-SH lifestyle table,
# so the bars match the per-sample guild definitions used elsewhere.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import config as cfg
import phyloseq_io as ps

MIN_SAMPLES = 2

lifestyle_to_guild = {
    "ectomycorrhizal": "ECM", "soil_saprotroph": "Soil saprotroph",
    "litter_saprotroph": "Litter saprotroph", "wood_saprotroph": "Wood saprotroph",
    "unspecified_saprotroph": "Unspec. saprotroph", "dung_saprotroph": "Unspec. saprotroph",
    "nectar/tap_saprotroph": "Unspec. saprotroph", "animal_parasite": "Animal parasite",
    "lichenized": "Lichenised", "root_endophyte": "Root endophyte", "foliar_endophyte": "Root endophyte",
}
guild_order = list(cfg.FUN_GUILD_COLS.values())
guild_col = {"ECM": "#D44A32", "Soil saprotroph": "#E8845A", "Litter saprotroph": "#F5B88A",
             "Wood saprotroph": "#E89CC5", "Unspec. saprotroph": "#D466A8", "Animal parasite": "#5B7FBF",
             "Lichenised": "#2CA089", "Root endophyte": "#7A9E3A", "Unknown/other": "#BDC3C7"}

counts, taxa, samples, meta, tax = ps.load(cfg.PHY_FUN)
veg = meta["vegetation"].values
dp = np.isin(veg, ["dwarf_pine"])
occurs = lambda z: set(taxa[(counts[:, np.isin(veg, z)] > 0).sum(1) >= MIN_SAMPLES])
M, D, F = occurs(["meadow"]), occurs(["dwarf_pine"]), occurs(["forest_A", "forest_O"])
categories = {"DP∩Meadow only": (D & M) - F, "DP∩Forest only": (D & F) - M,
              "DP∩All three": D & M & F, "DP unique": D - M - F}

life = pd.read_excel(cfg.LIFESTYLE, sheet_name="working")
sh_life = dict(zip(life["Unnamed: 0"].astype(str), life["style"].astype(str)))
guild_of = lambda sh: lifestyle_to_guild.get(sh_life.get(sh, ""), "Unknown/other")

dp_total = counts[:, dp].sum(0)
comp = {c: dict.fromkeys(guild_order, 0.0) for c in categories}
for cat, ids in categories.items():
    sel = np.isin(taxa, list(ids))
    mean_ra = (counts[sel][:, dp] / dp_total).mean(1) * 100
    for sh, ra in zip(taxa[sel], mean_ra):
        comp[cat][guild_of(sh)] += ra

cats = list(categories)
fig, ax = plt.subplots(figsize=(5.5, 5))
bottom = np.zeros(len(cats))
for guild in guild_order:
    vals = np.array([comp[c][guild] for c in cats])
    ax.bar(range(len(cats)), vals, bottom=bottom, width=0.62, label=guild,
           color=guild_col[guild], edgecolor="white", linewidth=0.3)
    bottom += vals
ax.set_xticks(range(len(cats)))
ax.set_xticklabels(cats, fontsize=8)
ax.set_ylabel("Mean relative abundance in dwarf pine (%)")
ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig("Fig_dwarfpine_guilds.png", dpi=300, bbox_inches="tight")
