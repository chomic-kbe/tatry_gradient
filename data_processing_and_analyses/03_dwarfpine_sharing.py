# Dwarf pine as a meadow-forest transition: which taxa does it share, and with whom?
# A taxon counts as 'present' in a zone if it occurs in >= 2 samples of that zone.
# Panel A/B: Venn of meadow / dwarf pine / forest (bacteria ASVs, fungi SHs).
# Panel C/D: guild make-up of the dwarf-pine taxa, split by sharing category.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

import config as cfg
import phyloseq_io as ps

MIN_SAMPLES = 2


def zone_sets(path):
    counts, taxa, samples, meta, tax = ps.load(path)
    veg = meta["vegetation"].values
    def occurs(zones):
        cols = np.isin(veg, zones)
        return set(taxa[(counts[:, cols] > 0).sum(1) >= MIN_SAMPLES])
    meadow = occurs(["meadow"])
    dpine = occurs(["dwarf_pine"])
    forest = occurs(["forest_A", "forest_O"])
    return meadow, dpine, forest


def venn_tuple(m, d, f):
    return (len(m - d - f), len(d - m - f), len((m & d) - f), len(f - m - d),
            len((m & f) - d), len((d & f) - m), len(m & d & f))


fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
for ax, path, title in [(axes[0], cfg.PHY_BAC, "Bacteria (ASVs)"),
                        (axes[1], cfg.PHY_FUN, "Fungi (SHs)")]:
    subsets = venn_tuple(*zone_sets(path))
    v = venn3(subsets, set_labels=("Meadow", "Dwarf pine", "Forest"), ax=ax)
    for region, colour in zip(["100", "010", "001"],
                              [cfg.veg_col["meadow"], cfg.veg_col["dwarf_pine"], cfg.veg_col["forest_A"]]):
        patch = v.get_patch_by_id(region)
        if patch:
            patch.set_color(colour)
            patch.set_alpha(0.55)
    venn3_circles(subsets, ax=ax, linewidth=1.0, color="#555555")
    ax.set_title(title, fontweight="bold", loc="left")

fig.tight_layout()
fig.savefig("Fig_dwarfpine_venn.png", dpi=300, bbox_inches="tight")

# report the counts to the console so they can be checked against the figure
for path, name in [(cfg.PHY_BAC, "bacteria"), (cfg.PHY_FUN, "fungi")]:
    m, d, f = zone_sets(path)
    print(name, "venn (M-only,D-only,M&D,F-only,M&F,D&F,all):", venn_tuple(m, d, f))
