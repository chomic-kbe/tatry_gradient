# Table S1 -- multivariate analyses (sequential & marginal PERMANOVA, partial db-RDA,
# variance partitioning, within-elevation tests). These were run in R (vegan::adonis2,
# dbrda, varpart); here we just format the exported results into the supplementary table.

import pandas as pd
from docx import Document

res = pd.read_excel(cfg_path := __import__("config").PERMANOVA, sheet_name="Sheet1", header=None)

doc = Document()
doc.add_paragraph("Table S1. Multivariate analyses of bacterial and fungal community composition "
                  "(Bray–Curtis distances, 999 permutations).")
tbl = doc.add_table(rows=1, cols=6)
tbl.style = "Table Grid"
for i, h in enumerate(["Analysis", "Effect / comparison", "Bacteria R²", "p", "Fungi R²", "p"]):
    tbl.rows[0].cells[i].text = h

for _, row in res.iterrows():
    cells = [str(x) for x in row.tolist() if str(x) != "nan"]
    if len(cells) >= 6:                       # skip caption / blank rows
        out = tbl.add_row().cells
        for i, v in enumerate(cells[:6]):
            out[i].text = v

doc.save("TableS1_PERMANOVA.docx")
print("wrote TableS1_PERMANOVA.docx")
