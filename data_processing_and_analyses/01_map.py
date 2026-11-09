# Figure 1 -- study area map with a Slovakia inset.
import geopandas as gpd
import matplotlib.pyplot as plt

world = gpd.read_file("ne_50m_admin_0_countries.geojson")
slovakia = world[world["ADMIN"] == "Slovakia"]
neighbours = world[world["ADMIN"].isin(["Czechia", "Poland", "Ukraine", "Hungary", "Austria"])]

# Furkotska valley, High Tatras
site = (20.0, 49.17)

fig, ax = plt.subplots(figsize=(3.4, 2.0))
neighbours.plot(ax=ax, facecolor="#f5f5f5", edgecolor="#b0b0b0", linewidth=0.5)
slovakia.plot(ax=ax, facecolor="#c9c9c9", edgecolor="#333333", linewidth=1.0)
ax.scatter(*site, marker="*", s=250, c="#d6272b", edgecolor="k", linewidth=0.5, zorder=5)
ax.set_xlim(16.2, 23.0)
ax.set_ylim(47.5, 50.2)
ax.set_aspect("equal")
ax.axis("off")

fig.tight_layout()
fig.savefig("Fig1_inset.png", dpi=300, bbox_inches="tight")
