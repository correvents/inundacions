#!/usr/bin/env python3
"""Prepara les conques per a conques.html: afegeix a cada conca la conca d'aigües avall (camp "avall")
i arrodoneix coordenades. Llegeix data/catalunya/conques/*.geojson i escriu data/catalunya/conques/web/*.json.

La conca d'avall és, d'entre les que toquen el punt de sortida, la que drena més superfície (km2_amunt)."""
import json
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.strtree import STRtree

D = Path("data/catalunya/conques"); W = D / "web"; W.mkdir(parents=True, exist_ok=True)

def rnd(c):
    return [round(c[0], 4), round(c[1], 4)] if isinstance(c[0], (int, float)) else [rnd(x) for x in c]

for n in ("25", "100"):
    fc = json.loads((D / f"conques_{n}km2.geojson").read_text())
    F = fc["features"]; geoms = [shape(f["geometry"]) for f in F]; tree = STRtree(geoms)
    for i, f in enumerate(F):
        p = f["properties"]; pt = Point(p["sortida_lon"], p["sortida_lat"]); best = None
        for j in tree.query(pt.buffer(0.0018)):
            q = F[j]["properties"]
            if j != i and geoms[j].distance(pt) < 0.0018 and q["km2_amunt"] > p["km2_amunt"] and \
               (best is None or q["km2_amunt"] > F[best]["properties"]["km2_amunt"]):
                best = j
        p["avall"] = F[best]["properties"]["id"] if best is not None else None
        f["geometry"]["coordinates"] = rnd(f["geometry"]["coordinates"])
    (W / f"conques_{n}.json").write_text(json.dumps(fc, separators=(",", ":")))
r = json.loads((D / "rius.geojson").read_text())
for f in r["features"]:
    f["geometry"]["coordinates"] = rnd(f["geometry"]["coordinates"])
(W / "rius.json").write_text(json.dumps(r, separators=(",", ":")))
