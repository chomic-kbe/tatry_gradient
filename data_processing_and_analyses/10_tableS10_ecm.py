# Table S10 -- ectomycorrhizal fungi per site.
# SHs assigned ECM by FungalTraits (via the per-SH lifestyle table), refined at SH level.
# For each SH: genus, species, UNITE SH code, prevalence (% of samples), and mean
# relative abundance per site. Only SHs reaching >= 0.05% at some site are listed.

import numpy as np
import pandas as pd

import config as cfg
import phyloseq_io as ps

counts, taxa, samples, meta, tax = ps.load(cfg.PHY_FUN)
genus, species, sh_code = tax[:, 5], tax[:, 6], tax[:, 7]
sk = ps.site_of(meta)
rel = counts / counts.sum(0)

life = pd.read_excel(cfg.LIFESTYLE, sheet_name="working")
ecm_ids = set(life.loc[life["style"] == "ectomycorrhizal", "SH"].astype(str)) | \
          set(life.loc[life["style"] == "ectomycorrhizal", "Unnamed: 0"].astype(str))

rows = []
for i in np.where(np.isin(sh_code, list(ecm_ids)))[0]:
    per_site = {s: rel[i][sk == s].mean() * 100 for s in cfg.sites}
    if max(per_site.values()) < 0.05:
        continue
    rows.append({"Genus": genus[i], "Species": species[i], "SH": sh_code[i],
                 "Prevalence_%": round((counts[i] > 0).mean() * 100),
                 **{s: round(per_site[s], 2) for s in cfg.sites}})

ecm = pd.DataFrame(rows).sort_values(["Genus", "Species"])
ecm.to_csv("TableS10_ECM.csv", index=False)
print(f"wrote TableS10_ECM.csv ( {len(ecm)} ECM SHs )")
