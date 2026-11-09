# stats_letters.py -- Kruskal-Wallis + Dunn post-hoc, with compact letter display.
# Letters are ordered so that 'a' marks the site with the highest value
# (as requested during review).

import numpy as np
import string
from scipy.stats import kruskal
import scikit_posthocs as sp


def kw_dunn(df, value, site_col, order):
    """Return per-site mean, SE, compact letters, and the KW p-value."""
    groups = {s: df.loc[df[site_col] == s, value].dropna().values for s in order}
    present = [s for s in order if len(groups[s])]
    H, p = kruskal(*[groups[s] for s in present])

    means = {s: groups[s].mean() if len(groups[s]) else np.nan for s in order}
    ses = {s: groups[s].std(ddof=1) / np.sqrt(len(groups[s])) if len(groups[s]) > 1 else np.nan
           for s in order}

    letters = {s: "" for s in order}
    if p < 0.05:
        long = df[df[site_col].isin(present)][[site_col, value]].dropna()
        pmat = sp.posthoc_dunn(long, val_col=value, group_col=site_col, p_adjust="bonferroni")
        letters = _cld(present, pmat, means)
    return means, ses, letters, p


def _cld(groups, pmat, means):
    """Compact letter display: groups that do NOT differ (p > 0.05) share a letter.
    Letters are then renamed so 'a' sits on the highest mean."""
    n = len(groups)
    same = np.zeros((n, n), bool)
    for i, a in enumerate(groups):
        for j, b in enumerate(groups):
            same[i, j] = (i == j) or (pmat.loc[a, b] > 0.05)

    # greedy assignment of letters to maximal 'not different' cliques
    assigned = {g: set() for g in groups}
    col = 0
    remaining = list(range(n))
    order_by_mean = sorted(range(n), key=lambda k: -means[groups[k]])
    while remaining:
        clique = []
        for k in order_by_mean:
            if k in remaining and all(same[k, m] for m in clique):
                clique.append(k)
        for k in clique:
            assigned[groups[k]].add(col)
        remaining = [k for k in remaining if k not in clique or len(assigned[groups[k]]) == 0]
        remaining = [k for k in range(n) if not assigned[groups[k]]]
        col += 1
        if col > n:
            break

    # rename columns so highest mean -> 'a'
    colmax = {c: max(means[g] for g in groups if c in assigned[g]) for c in range(col)}
    rank = {c: string.ascii_lowercase[i] for i, c in enumerate(sorted(colmax, key=lambda c: -colmax[c]))}
    out = {g: "" for g in groups}
    for g in groups:
        out[g] = "".join(sorted(rank[c] for c in assigned[g]))
    return out
