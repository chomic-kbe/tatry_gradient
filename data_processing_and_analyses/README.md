# Analysis code

Code to reproduce the figures and tables for the manuscript on soil bacterial and
fungal communities along an alpine-to-forest vegetation gradient (High Tatra Mts, Slovakia).

Sequence processing (DADA2) and the multivariate statistics (PERMANOVA, db-RDA,
variance partitioning; vegan in R) were done in R. The phyloseq objects and the
PERMANOVA/db-RDA results were exported from R and are read here; the scripts below
produce the figures and the per-site tables.

## Requirements

Python 3.10+ and:

    pip install numpy pandas scipy scikit-learn matplotlib matplotlib-venn \
                rdata geopandas python-docx openpyxl scikit-posthocs

## Data

Place the following in the same folder as the scripts (paths are set in `config.py`):

    tatry_gradient_v4_phy_ASV.rds          bacterial phyloseq (ASV)
    tatry_gradient_its_phy_sh.rds          fungal phyloseq (SH)
    bacteria_nmds_data.rds                 bacterial NMDS site scores (vegan::metaMDS)
    map_tatry_gradient_251131.xlsx         soil chemistry, enzymes, qPCR per sample
    tatry_gradient_BLAST_10_best_hits_with_lifestyle.xlsx   fungal lifestyles (sheet 'working')
    supp_table_elev_vs_veget.xlsx          PERMANOVA / db-RDA results (from R)
    map_tatry_gradient_with_guilds.xlsx    per-sample guild table (fungal F_, bacterial B_, enzymes)
    ne_50m_admin_0_countries.geojson       Natural Earth 50 m (Fig. 1 inset only)

## Run

    python run_all.py        # reproduces all figures and tables

## Files

    run_all.py          run every script below in order
    config.py           paths, site order, colours
    phyloseq_io.py      read .rds phyloseq objects; per-sample relative abundance
    stats_letters.py    Kruskal-Wallis + Dunn post-hoc, compact letters (a = highest)

    01_map.py                     Fig. 1   study-area map / inset
    02_nmds.py                    Fig.     NMDS ordination (bacteria, fungi)
    03_dwarfpine_sharing.py       Fig.     meadow / dwarf pine / forest sharing (Venn)
    04_guild_bars.py              Fig.     dwarf-pine guild composition by sharing group
    05_table2_microbial.py        Table 2  biomass / enzyme activities per site
    06_tableS1_permanova.py       Table S1 multivariate analyses
    07_tableS2S3_guilds.py        Tables S2/S3  functional guilds per site
    08_tableS4S5_top_taxa.py      Tables S4/S5  top genera / orders per site
    09_tableS7S8S9_core.py        Tables S7-S9  core microbiome
    10_tableS10_ecm.py            Table S10  ectomycorrhizal fungi per site
    11_figS_varpart.py            Fig. S   variance partitioning
    12_figS_enzyme_corr.py        Fig. S   guild-enzyme correlations

## Notes

- Relative abundance is computed per sample from raw counts (not rarefied).
- Core microbiome: present in >= 80% of samples, at ASV (bacteria) / SH (fungi) level.
- Dwarf-pine sharing: a taxon is 'present' in a zone if found in >= 2 of its samples.
- Post-hoc letters are ordered so that 'a' marks the highest value.
- Scripts 07, 08 also write the bacterial counterpart (S3, S5) with the same function;
  the bacterial order-to-group map (Methods) is applied the same way as the fungal lifestyles.
