# config.py -- file paths and site ordering used across the analysis
# Tatry gradient: soil bacterial & fungal communities, scree -> forest

# phyloseq objects (exported from R as .rds, read here with the rdata package)
PHY_BAC = "tatry_gradient_v4_phy_ASV.rds"   # bacteria, ASV level
PHY_FUN = "tatry_gradient_its_phy_sh.rds"   # fungi, SH level
NMDS_BAC = "bacteria_nmds_data.rds"         # bacterial NMDS site scores (from vegan::metaMDS)

# tables
MAP_XLSX   = "map_tatry_gradient_251131.xlsx"                      # soil chemistry + enzymes + qPCR
LIFESTYLE  = "tatry_gradient_BLAST_10_best_hits_with_lifestyle.xlsx"  # fungal lifestyles (sheet 'working')
PERMANOVA  = "supp_table_elev_vs_veget.xlsx"                       # multivariate results from R

# sites along the vegetation sequence (used to order columns everywhere)
sites = ["scree_2160", "meadow_2160", "meadow_2130", "meadow_1780", "dwarf_pine_1780",
         "forest_A_1470", "forest_O_1470", "forest_A_1330", "forest_O_1330"]

site_labels = ["Scree\n2160", "Meadow\n2160", "Meadow\n2130", "Meadow\n1780", "Dw.pine\n1780",
               "Forest A\n1470", "Forest O\n1470", "Forest A\n1330", "Forest O\n1330"]

# plotting colours (sampled from Fig. 1)
veg_col = {"scree": "#67B9A3", "meadow": "#FF875A", "dwarf_pine": "#7C93C5",
           "forest_A": "#B1DB66", "forest_O": "#5B8C2A"}

# per-sample guild table: master chemistry/enzymes + fungal (F_) and bacterial (B_)
# functional-group relative abundances, one row per sample. This is the table the
# guild figures and guild tables are built from, so results match the manuscript exactly.
GUILD_TABLE = "map_tatry_gradient_with_guilds.xlsx"

# fungal guild columns (as % per sample) and display names
FUN_GUILD_COLS = {
    "F_ectomycorrhizal": "ECM", "F_soil_saprotroph": "Soil saprotroph",
    "F_litter_saprotroph": "Litter saprotroph", "F_wood_saprotroph": "Wood saprotroph",
    "F_unspecified_saprotroph": "Unspec. saprotroph", "F_animal_parasite": "Animal parasite",
    "F_lichenized": "Lichenised", "F_root_endophyte": "Root endophyte",
    "F_unknown_other": "Unknown/other",
}
# bacterial functional-group columns and display names
BAC_GROUP_COLS = {
    "B_Recalcitrant_C_decomp": "Recalcitrant C", "B_Labile_C_heterotrophy": "Labile C",
    "B_Acidophilic_C_cycling": "Acidophilic C", "B_C1_metabolism": "C1 metabolism",
    "B_Oligotrophic_C_processing": "Oligotrophic C", "B_Predation_cell_lysis": "Predation",
    "B_Chemolithotrophy": "Chemolithotrophy", "B_PVC_complex": "PVC",
}
