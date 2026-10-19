#!/usr/bin/env python3
"""Baixa (un cop) les capes estàtiques per saber què s'inunda a cada conca:
  - zones inundables de l'ACA (T10, T50, T100, T500) — WFS aplicacions.aca.gencat.cat
  - punts i població en risc del SNCZI (MITECO) — WFS gis.miteco.gob.es
  - inundació costanera del SNCZI (T100, T500)
  - elements sensibles d'OpenStreetMap (escoles, càmpings, hospitals, residències, pàrquings soterrats, edificis...)
Escriu GeoJSON (EPSG:4326) a data/<conca>/zones/.

Ús:  python3 scripts/prepara_zones.py --config config/ridaura.json"""
import argparse, json, sys
from pathlib import Path
import requests
from pyproj import Transformer

UA = {"User-Agent": "inundacions/0.1 (github.com/correvents/inundacions)"}
ACA = "https://aplicacions.aca.gencat.cat/geoserver/wfs"
MITECO = "https://gis.miteco.gob.es/geoserver/wfs"
CAPES = {
    "aca_T10": (ACA, "vishid:Zones_inundables_T10"),
    "aca_T50": (ACA, "vishid:Zones_inundables_T50"),
    "aca_T100": (ACA, "vishid:Zones_inundables_T100"),
    "aca_T500": (ACA, "vishid:Zones_inundables_T500"),
    "snczi_punts_risc_T10": (MITECO, "agua:ZIF_RiesgoInundacion_Q10"),
    "snczi_punts_risc_T100": (MITECO, "agua:Zif_riesgoinundacion_q100"),
    "snczi_punts_risc_T500": (MITECO, "agua:Zif_riesgoinundacion_q500"),
    "snczi_poblacio_T10": (MITECO, "agua:Zif_poblacion_afect_q10"),
    "snczi_poblacio_T100": (MITECO, "agua:Zif_riesgo_pob_q100"),
    "snczi_poblacio_T500": (MITECO, "agua:Zif_riesgo_pob_q500"),
    "snczi_costa_T100": (MITECO, "costas:zim_laminas_q100"),
    "snczi_costa_T500": (MITECO, "costas:zim_laminas_q500"),
}
OVERPASS = ["https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter"]


def wfs(url, capa, bbox):
    t = Transformer.from_crs(4326, 25831, always_xy=True)
    x0, y0 = t.transform(bbox[0], bbox[1]); x1, y1 = t.transform(bbox[2], bbox[3])
    p = {"service": "WFS", "version": "1.1.0", "request": "GetFeature", "typeName": capa, "outputFormat": "application/json",
         "srsName": "EPSG:4326", "bbox": f"{x0:.0f},{y0:.0f},{x1:.0f},{y1:.0f},EPSG:25831"}
    r = requests.get(url, params=p, headers=UA, timeout=180)
    if r.status_code != 200 or not r.text.lstrip().startswith("{"):
        print(f"  {capa}: HTTP {r.status_code} {r.text[:200]!r}"); return None
    return r.json()


def osm(bbox):
    s, w, n, e = bbox[1], bbox[0], bbox[3], bbox[2]
    q = f"""[out:json][timeout:120];
(
  nwr["amenity"~"^(school|kindergarten|college|university|hospital|clinic|nursing_home|social_facility|fire_station|police|townhall|community_centre|place_of_worship|fuel)$"]({s},{w},{n},{e});
  nwr["tourism"~"^(camp_site|caravan_site|hotel|hostel|apartment)$"]({s},{w},{n},{e});
  nwr["amenity"="parking"]["parking"~"underground|multi-storey"]({s},{w},{n},{e});
  nwr["leisure"~"^(sports_centre|stadium|water_park)$"]({s},{w},{n},{e});
  way["building"]({s},{w},{n},{e});
  way["highway"~"^(motorway|trunk|primary|secondary|tertiary|residential|unclassified|living_street|pedestrian)$"]({s},{w},{n},{e});
  way["bridge"="yes"]["highway"]({s},{w},{n},{e});
  node["ford"="yes"]({s},{w},{n},{e});
  way["ford"="yes"]({s},{w},{n},{e});
);
out geom tags;"""
    for u in OVERPASS:
        try:
            r = requests.post(u, data={"data": q}, headers=UA, timeout=300)
            if r.ok:
                return r.json()
            print(f"  overpass {u}: HTTP {r.status_code}")
        except requests.RequestException as e:
            print(f"  overpass {u}: {type(e).__name__}")
    return None


def osm_a_geojson(j):
    feats = []
    for el in j.get("elements", []):
        tags = el.get("tags", {})
        if el["type"] == "node":
            g = {"type": "Point", "coordinates": [el["lon"], el["lat"]]}
        elif "geometry" in el:
            c = [[p["lon"], p["lat"]] for p in el["geometry"]]
            tancat = len(c) > 3 and c[0] == c[-1] and ("building" in tags or "amenity" in tags or "tourism" in tags or "leisure" in tags)
            g = {"type": "Polygon", "coordinates": [c]} if tancat else {"type": "LineString", "coordinates": c}
        elif "bounds" in el:
            b = el["bounds"]
            g = {"type": "Point", "coordinates": [(b["minlon"] + b["maxlon"]) / 2, (b["minlat"] + b["maxlat"]) / 2]}
        else:
            continue
        feats.append({"type": "Feature", "properties": {"osm": f"{el['type']}/{el['id']}", **tags}, "geometry": g})
    return {"type": "FeatureCollection", "features": feats}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    bbox = cfg["bbox_zones"]
    dest = Path("data") / cfg["conca"] / "zones"; dest.mkdir(parents=True, exist_ok=True)
    errors = 0
    for nom, (url, capa) in CAPES.items():
        j = wfs(url, capa, bbox)
        if j is None:
            errors += 1; continue
        (dest / f"{nom}.geojson").write_text(json.dumps(j, ensure_ascii=False, separators=(",", ":")))
        print(f"{nom}: {len(j.get('features', []))} elements")
    j = osm(bbox)
    if j is None:
        errors += 1
    else:
        g = osm_a_geojson(j)
        (dest / "osm.geojson").write_text(json.dumps(g, ensure_ascii=False, separators=(",", ":")))
        print(f"osm: {len(g['features'])} elements")
    return 1 if errors == len(CAPES) + 1 else 0


if __name__ == "__main__":
    sys.exit(main())
