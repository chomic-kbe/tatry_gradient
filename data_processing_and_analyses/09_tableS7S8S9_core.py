# Tables S7-S9 -- core microbiome.
# Core = taxa present in >= 80% of all samples, defined at ASV (bacteria) and SH (fungi)
# level so that 'shared' means the same organisms, not the same broad lineages.
# S7: core bacterial ASVs; S8: core fungal SHs; S9: summed core abundance per site.

import numpy as np
import pandas as pd

import config as cfg
import phyloseq_io as ps

THRESHOLD = 0.80


def core(path):
    counts, taxa, samples, meta, tax = ps.load(path)
    rel = counts / counts.sum(0)
    is_core = (counts > 0).mean(1) >= THRESHOLD
    members = pd.DataFrame({
        "taxon": taxa[is_core],
        "mean_rel_abund_%": (rel[is_core].mean(1) * 100).round(3),
    }).sort_values("mean_rel_abund_%", ascending=False)
    sk = ps.site_of(meta)
    per_site = {s: round(rel[is_core][:, sk == s].sum(0).mean() * 100, 1) for s in cfg.sites}
    return members, per_site


bac_members, bac_site = core(cfg.PHY_BAC)
fun_members, fun_site = core(cfg.PHY_FUN)

bac_members.to_csv("TableS7_core_bacteria.csv", index=False)
fun_members.to_csv("TableS8_core_fungi.csv", index=False)
pd.DataFrame({"Site": cfg.sites,
              "Core bacteria (%)": [bac_site[s] for s in cfg.sites],
              "Core fungi (%)": [fun_site[s] for s in cfg.sites]}).to_csv("TableS9_core_per_site.csv", index=False)

print(f"core bacteria: {len(bac_members)} ASVs; core fungi: {len(fun_members)} SH")
print("wrote TableS7/S8/S9")
