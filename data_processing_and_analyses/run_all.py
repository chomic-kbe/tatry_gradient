# run_all.py -- reproduce every figure and table in one go.
# Expects the data files listed in README.md in this folder.

import runpy, sys

scripts = [
    "01_map.py", "02_nmds.py", "03_dwarfpine_sharing.py", "04_guild_bars.py",
    "05_table2_microbial.py", "06_tableS1_permanova.py", "07_tableS2S3_guilds.py",
    "08_tableS4S5_top_taxa.py", "09_tableS7S8S9_core.py", "10_tableS10_ecm.py",
    "11_figS_varpart.py", "12_figS_enzyme_corr.py",
]

for s in scripts:
    print(f"\n=== {s} ===")
    try:
        runpy.run_path(s, run_name="__main__")
    except FileNotFoundError as e:
        print(f"  skipped ({e})"); 
    except Exception as e:
        print(f"  ERROR: {e}", file=sys.stderr)
