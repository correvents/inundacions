#!/usr/bin/env python3
"""Creua les zones inundables de l'ACA (T10, T50, T100, T500) amb els elements d'OpenStreetMap i els punts de risc del SNCZI.
Per a cada element es guarda el període de retorn més petit que l'afecta (T10 < T50 < T100 < T500).

Escriu data/<conca>/afectats.geojson (per al visor) i data/<conca>/afectats_resum.json.
Ús:  python3 scripts/creua_zones.py --config config/ridaura.json"""
import argparse, json
from pathlib import Path
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from shapely.strtree import STRtree

TS = [10, 50, 100, 500]


def categoria(p):
    if p.get("amenity") in ("school", "kindergarten", "college", "university"): return "educació"
    if p.get("amenity") in ("hospital", "clinic", "nursing_home", "social_facility"): return "sanitari / residència"
    if p.get("tourism") in ("camp_site", "caravan_site"): return "càmping"
    if p.get("tourism") in ("hotel", "hostel", "apartment"): return "allotjament"
    if p.get("amenity") == "parking": return "pàrquing soterrat"
    if p.get("amenity") in ("fire_station", "police", "townhall"): return "servei d'emergència / públic"
    if p.get("amenity") or p.get("leisure"): return "concurrència pública"
    if p.get("ford") == "yes": return "gual"
    if p.get("highway"): return "pont" if p.get("bridge") == "yes" else "carrer / carretera"
    if p.get("building"): return "edifici"
    return "altres"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"]
    codis = set(cfg.get("codis_zones_aca") or [])     # zones de l'ACA d'aquesta conca (el bbox també en baixa de veïnes)
    zones = {T: unary_union([shape(f["geometry"]) for f in json.loads((d / "zones" / f"aca_T{T}.geojson").read_text())["features"]
                             if not codis or f["properties"].get("CODI") in codis]).buffer(0) for T in TS}
    # una zona de període més gran inclou les de període més petit (les capes no sempre ho garanteixen)
    acum = None
    for T in TS:
        acum = zones[T] if acum is None else acum.union(zones[T])
        zones[T] = acum

    elements = []
    for f in json.loads((d / "zones" / "osm.geojson").read_text())["features"]:
        p = f["properties"]
        elements.append((shape(f["geometry"]), {"font": "OSM", "osm": p["osm"], "nom": p.get("name"), "categoria": categoria(p)}))
    for f in json.loads((d / "zones" / "snczi_punts_risc_T500.geojson").read_text())["features"]:
        p = f["properties"]
        elements.append((shape(f["geometry"]), {"font": "SNCZI", "nom": p.get("desc_elto"), "categoria": f"{p.get('tipo_elto')} · {p.get('subtipo_el')}",
                                                 "gravetat": p.get("clasif_afe")}))
    arbre = STRtree([g for g, _ in elements])
    afectats = {}
    for T in reversed(TS):                       # de T500 a T10: el que queda és el T mínim
        for i in arbre.query(zones[T], predicate="intersects"):
            afectats[int(i)] = T
    feats, resum = [], {f"T{T}": {} for T in TS}
    carrers = {f"T{T}": set() for T in TS}
    for i, T in afectats.items():
        g, p = elements[i]
        if p["categoria"] == "carrer / carretera":
            g = g.intersection(zones[500])       # només el tram dins la zona inundable
            if p["nom"]:
                carrers[f"T{T}"].add(p["nom"])
        feats.append({"type": "Feature", "properties": {**p, "T": T}, "geometry": mapping(g.simplify(0.00002) if g.geom_type != "Point" else g)})
        resum[f"T{T}"][p["categoria"]] = resum[f"T{T}"].get(p["categoria"], 0) + 1
    out = {"type": "FeatureCollection", "features": feats}
    (d / "afectats.geojson").write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")))
    noms = {f"T{T}": sorted({f["properties"]["nom"] for f in feats if f["properties"]["T"] == T and f["properties"]["nom"]
                             and f["properties"]["categoria"] not in ("carrer / carretera", "edifici")}) for T in TS}
    zz = {f"T{T}": mapping(zones[T].simplify(0.00003)) for T in TS}
    (d / "zones_acumulades.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"T": T}, "geometry": zz[f"T{T}"]} for T in reversed(TS)]}, separators=(",", ":")))
    (d / "afectats_resum.json").write_text(json.dumps({"per_categoria": resum, "llocs": noms,
                                                         "carrers": {k: sorted(v) for k, v in carrers.items()}}, ensure_ascii=False, indent=1))
    for T in TS:
        print(f"T{T}:", resum[f"T{T}"], "| llocs:", noms[f"T{T}"][:15], "| carrers:", sorted(carrers[f"T{T}"])[:25])


if __name__ == "__main__":
    main()
