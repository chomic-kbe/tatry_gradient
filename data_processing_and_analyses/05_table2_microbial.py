# Table 2 -- microbial biomass, diversity and enzyme activities per site.
# Mean +/- SE, Kruskal-Wallis across sites with Dunn post-hoc (Bonferroni);
# letters are set so 'a' marks the highest value. Rows grouped biomass / diversity / activity.
# (Chao1 richness is read from the alpha-diversity output; here we show the chemistry/enzyme block.)

import pandas as pd
from docx import Document
from docx.shared import Pt

import config as cfg
from stats_letters import kw_dunn

df = pd.read_excel(cfg.MAP_XLSX)
df = df[~df["site"].astype(str).str.startswith("ctrl")].copy()
df["sitekey"] = df["vegetation"] + "_" + df["elevetaion"].astype(int).astype(str)

blocks = {
    "BIOMASS": ["MBC", "MBN", "MBP", "Bacteria _genes", "Fungi _genes"],
    "ACTIVITY": ["Respiration", "HYD", "BG", "CBH", "PHOS", "LAP", "CHI",
                 "oxidase", "PhOx", "PerOx"],
}

doc = Document()
doc.add_paragraph("Table 2. Soil microbial biomass and enzyme activities along the sites "
                  "(mean ± SE; different letters: Kruskal–Wallis + Dunn, Bonferroni, p < 0.05).")
tbl = doc.add_table(rows=1, cols=len(cfg.sites) + 2)
tbl.style = "Table Grid"
hdr = tbl.rows[0].cells
hdr[0].text = "Parameter"
for i, lab in enumerate(cfg.site_labels):
    hdr[i + 1].text = lab.replace("\n", " ")
hdr[-1].text = "KW p"

for block, params in blocks.items():
    r = tbl.add_row().cells
    r[0].text = block
    for p in params:
        if p not in df.columns:
            continue
        means, ses, letters, kwp = kw_dunn(df, p, "sitekey", cfg.sites)
        # gene counts (qPCR) are reported scaled, as in the manuscript
        scale, suffix = (1e9, " (×10⁹)") if p == "Bacteria _genes" else \
                        (1e7, " (×10⁷)") if p == "Fungi _genes" else (1.0, "")
        row = tbl.add_row().cells
        row[0].text = p.replace("_genes", "").strip() + suffix
        for i, s in enumerate(cfg.sites):
            m, se, L = means[s], ses[s], letters[s]
            if pd.isna(m):
                row[i + 1].text = ""
            else:
                m, se = m / scale, se / scale
                fmt = (lambda x: f"{x:,.0f}") if abs(m) >= 100 else (lambda x: f"{x:.1f}")
                row[i + 1].text = f"{fmt(m)}±{fmt(se)} {L}".strip()
        row[-1].text = "<0.001" if kwp < 0.001 else f"{kwp:.3f}"

doc.save("Table2_microbial.docx")
print("wrote Table2_microbial.docx")
