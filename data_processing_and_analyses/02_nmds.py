# NMDS ordination of community composition (Bray-Curtis).
# Bacteria: site scores exported from vegan::metaMDS in R.
# Fungi: recomputed here from the SH table (same Bray-Curtis, k = 2).


import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
import rdata

import config as cfg
import phyloseq_io as ps

# --- bacteria (stored scores) ---
bac = rdata.read_rds(cfg.NMDS_BAC)
bac.columns = [str(c) for c in bac.columns]
bac_xy = {"x": bac["MDS1"].values, "y": -bac["MDS2"].values, "veg": bac["vegetation"].values}

# --- fungi (recompute) ---
counts, taxa, samples, meta, tax = ps.load(cfg.PHY_FUN)
rel = (counts / counts.sum(0)).T
bc = squareform(pdist(rel, metric="braycurtis"))
mds = MDS(n_components=2, dissimilarity="precomputed", metric_mds=False,
          n_init=10, max_iter=500, random_state=1, normalized_stress="auto")
pts = mds.fit_transform(bc)
print(f"fungal NMDS stress = {mds.stress_:.3f}")
fun_xy = {"x": -pts[:, 0], "y": pts[:, 1], "veg": meta["vegetation"].values}

markers = {"scree": "D", "meadow": "o", "dwarf_pine": "s", "forest_A": "^", "forest_O": "v"}

fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
for ax, data, title in [(axes[0], bac_xy, "Bacteria"), (axes[1], fun_xy, "Fungi")]:
    ax.axhline(0, color="#cccccc", lw=0.8)
    ax.axvline(0, color="#cccccc", lw=0.8)
    for veg in markers:
        m = data["veg"] == veg
        ax.scatter(np.asarray(data["x"])[m], np.asarray(data["y"])[m],
                   c=cfg.veg_col[veg], marker=markers[veg], s=55,
                   edgecolor="white", linewidth=0.5)
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("NMDS axis 1")
axes[0].set_ylabel("NMDS axis 2")

fig.tight_layout()
fig.savefig("Fig_NMDS.png", dpi=300, bbox_inches="tight")
