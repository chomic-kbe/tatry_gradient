# phyloseq_io.py -- read phyloseq .rds objects in Python and attach the site key
# We keep the raw counts and compute relative abundance per sample (not rarefied),
# which matches how the published tables were made.

import numpy as np
import pandas as pd
import rdata

# catchment code -> elevation label (VW is the uppermost, RA the lowest forest)
_elev = {"VW": "2160", "NW": "2130", "FU": "1780", "JA": "1470", "RA": "1330"}


def load(path):
    """Read a phyloseq .rds and return counts (taxa x samples), taxa ids,
    sample ids, the sample_data as a DataFrame, and the taxonomy matrix."""
    ps = rdata.read_rds(path)
    otu = ps.otu_table
    d0 = [str(x) for x in otu.coords["dim_0"].values]
    d1 = [str(x) for x in otu.coords["dim_1"].values]

    # samples are labelled s001, s002, ...; orientation of otu_table varies
    sample_like = lambda v: v[:1] == "s" and v[1:].isdigit()
    if all(sample_like(x) for x in d0[:3]):
        samples, taxa, counts = d0, d1, otu.values.T
    else:
        taxa, samples, counts = d0, d1, otu.values

    sd = ps.sam_data
    meta = pd.DataFrame({k: list(sd[k]) for k in sd},
                        index=[str(s) for s in sd["sample_id"]]).reindex(samples)

    tax = None
    tt = np.array(ps.tax_table)
    if tt.size and tt.size % len(taxa) == 0:
        tax = tt.reshape((len(taxa), tt.size // len(taxa)), order="F")

    return counts.astype(float), np.asarray(taxa), np.asarray(samples), meta, tax


def site_of(meta):
    """vegetation-elevation site key for each sample, e.g. 'forest_A_1330'."""
    keys = []
    for s in meta.index:
        catch, veg = str(meta.loc[s, "site"]).split("_", 1)
        keys.append(f"{veg}_{_elev[catch]}")
    return np.asarray(keys)
