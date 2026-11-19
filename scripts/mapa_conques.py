#!/usr/bin/env python3
"""Mapa de conques de Catalunya (HydroBASINS, nivells 7-10) per decidir com partir el MDT i el monitoratge.

Baixa HydroBASINS Europa (data.hydrosheds.org), es queda amb les conques que toquen Catalunya,
i escriu per a cada nivell:
  data/catalunya/conques/conques_nivNN.geojson  (simplificat, amb km² i id HydroBASINS)
  data/catalunya/conques/conques_nivNN.png       (mapa)
  data/catalunya/conques/resum.json               (nombre de conques i mides)"""
import io, json, zipfile
from pathlib import Path
import requests, geopandas as gpd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SORT = Path("data/catalunya/conques"); SORT.mkdir(parents=True, exist_ok=True)
NE = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_admin_1_states_provinces.geojson"
HB = "https://data.hydrosheds.org/file/hydrobasins/standard/hybas_eu_lev{:02d}_v1c.zip"
POBLES = {"Platja d'Aro": (3.067, 41.817), "Granollers": (2.288, 41.608), "Cabrera de Mar": (2.396, 41.527)}

prov = gpd.read_file(NE)
cat = prov[prov["iso_3166_2"].isin(["ES-B", "ES-GI", "ES-L", "ES-T"])].dissolve().to_crs(25831)
cat_geom = cat.geometry.iloc[0]
resum = {}
for niv in (7, 8, 9, 10):
    r = requests.get(HB.format(niv), timeout=600); r.raise_for_status()
    d = SORT / f"_hb{niv}"; zipfile.ZipFile(io.BytesIO(r.content)).extractall(d)
    shp = next(d.rglob("*.shp"))
    g = gpd.read_file(shp, bbox=(0.0, 40.4, 3.4, 42.9)).to_crs(25831)
    dins = g.geometry.intersection(cat_geom).area
    g = g[dins > 0.2 * g.geometry.area].copy()          # almenys un 20% dins de Catalunya
    g["km2"] = (g.geometry.area / 1e6).round(1)
    g = g[["HYBAS_ID", "NEXT_DOWN", "MAIN_BAS", "km2", "geometry"]]
    resum[niv] = {"conques": len(g), "km2_mediana": float(g.km2.median()), "km2_min": float(g.km2.min()), "km2_max": float(g.km2.max())}
    g2 = g.copy(); g2["geometry"] = g2.geometry.simplify(50)
    g2.to_crs(4326).to_file(SORT / f"conques_niv{niv:02d}.geojson", driver="GeoJSON", COORDINATE_PRECISION=5)
    fig, ax = plt.subplots(figsize=(10, 10), dpi=110)
    g.plot(ax=ax, column="HYBAS_ID", cmap="tab20", categorical=True, alpha=0.55, edgecolor="#333", linewidth=0.4)
    cat.boundary.plot(ax=ax, color="black", linewidth=1.2)
    for nom, (lon, lat) in POBLES.items():
        p = gpd.GeoSeries(gpd.points_from_xy([lon], [lat]), crs=4326).to_crs(25831).iloc[0]
        ax.plot(p.x, p.y, "o", color="red", ms=6, mec="white"); ax.annotate(nom, (p.x, p.y), xytext=(5, 5), textcoords="offset points", fontsize=9, weight="bold")
    ax.set_title(f"Conques HydroBASINS nivell {niv}: {len(g)} conques, mediana {resum[niv]['km2_mediana']:.0f} km²")
    ax.set_axis_off(); fig.tight_layout(); fig.savefig(SORT / f"conques_niv{niv:02d}.png"); plt.close(fig)
    print(niv, resum[niv])
(SORT / "resum.json").write_text(json.dumps(resum, indent=1))
for d in SORT.glob("_hb*"):
    import shutil; shutil.rmtree(d)
