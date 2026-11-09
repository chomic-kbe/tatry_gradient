# Tables S2 / S3 -- mean relative abundance (%) of fungal guilds (S2) and bacterial
# functional groups (S3) per site, from the per-sample guild table. Kruskal-Wallis
# across sites with Dunn post-hoc (Bonferroni); letters set so 'a' is the highest.

import pandas as pd

import config as cfg
from stats_letters import kw_dunn

g = pd.read_excel(cfg.GUILD_TABLE)
g = g[~g["site"].astype(str).str.startswith("ctrl")].copy()
g["sitekey"] = g["vegetation"] + "_" + g["elevetaion"].astype(int).astype(str)


def guild_table(col_map, outfile):
    rows = []
    for col, name in col_map.items():
        means, ses, letters, p = kw_dunn(g, col, "sitekey", cfg.sites)
        row = {"Guild": name}
        for s in cfg.sites:
            row[s] = "" if pd.isna(means[s]) else f"{means[s]:.2f} {letters[s]}".strip()
        row["KW p"] = "<0.001" if p < 0.001 else f"{p:.3f}"
        rows.append(row)
    pd.DataFrame(rows).to_csv(outfile, index=False)
    print("wrote", outfile)


guild_table(cfg.FUN_GUILD_COLS, "TableS2_fungal_guilds.csv")
guild_table(cfg.BAC_GROUP_COLS, "TableS3_bacterial_groups.csv")
