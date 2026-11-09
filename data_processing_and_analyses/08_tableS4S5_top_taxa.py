# Tables S4 / S5 -- mean relative abundance (%) of the top fungal genera (S4) and
# top bacterial orders (S5) per site. Taxa listed alphabetically; genus names are
# italicised in the Word/manuscript version.

import pandas as pd

import config as cfg
import phyloseq_io as ps

GENUS_RANK = 5
ORDER_RANK = 3


def top_table(path, rank, n):
    counts, taxa, samples, meta, tax = ps.load(path)
    sk = ps.site_of(meta)
    rel = counts / counts.sum(0)
    names = tax[:, rank]
    per_taxon = pd.DataFrame(rel, index=names).groupby(level=0).sum()
    per_site = per_taxon.T.groupby(sk).mean().T[cfg.sites] * 100
    keep = per_site.max(axis=1).sort_values(ascending=False).head(n).index
    return per_site.loc[sorted(keep)].round(2)


top_table(cfg.PHY_FUN, GENUS_RANK, 30).to_csv("TableS4_top_fungal_genera.csv")
top_table(cfg.PHY_BAC, ORDER_RANK, 29).to_csv("TableS5_top_bacterial_orders.csv")
print("wrote TableS4 (30 genera) and TableS5 (29 orders)")
