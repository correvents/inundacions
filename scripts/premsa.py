#!/usr/bin/env python3
"""Episodis d'inundació a Castell-Platja d'Aro recollits de la premsa, amb els carrers citats situats a l'OpenStreetMap.

Escriu data/<conca>/episodis/index.json (llista d'episodis per a pluja.html) i premsa.geojson (carrers i llocs citats).
Els carrers es busquen pel nom dins la plana de Platja d'Aro i Castell d'Aro; el que no es troba queda només com a text.

Ús:  python3 scripts/premsa.py --config config/ridaura.json"""
import argparse, json, sys
from pathlib import Path
from shapely.geometry import shape, mapping
from shapely.ops import unary_union

CAIXA = (3.02, 41.79, 3.085, 41.84)

EPISODIS = [
    {"id": "2026-10-04", "titol": "3-4 d'octubre de 2026: el Ridaura es desborda", "radar": "2026-10-03_radar.json",
     "episodi": "2026-10-04.json", "ia_mm": 150,
     "resum": "173 mm a Platja d'Aro, 198 mm a l'aforador. El Ridaura passa de 30 a 300 m³/s entre les 9:30 i les 11:30 i es desborda "
              "al tram final (Aro Gran): soterranis i baixos inundats, esquerda al pont. Centre negat des de les 9-10 h.",
     "carrers": ["Carrer de Juli Garreta", "Carrer de la Diagonal", "Avinguda de s'Agaró", "Avinguda del Cavall Bernat",
                 "Carrer Mestral", "Carrer de la Tramuntana", "Passeig del Ridaura", "Avinguda de París"],
     "llocs": ["Aro Gran (soterranis i baixos)", "Pont sobre el Ridaura (esquerda a l'asfalt)", "Parc dels Estanys",
               "Riera de la Coma desbordada a la part alta (zona rústica)"],
     "fonts": ["https://news.rpa.cat/2026/10/04/balanc-local-de-lemergencia-173-litres-per-metre-quadrat-nombroses-afectacions-i-una-recuperacio-progressiva-de-la-normalitat/",
               "https://www.gerio.cat/noticia/1757796/reobren-el-pont-de-platja-d-aro-afectat-per-una-esquerda-arran-del-desbordament-del-ridaura",
               "https://es.ara.cat/sociedad/habiamos-visto-litoral-recuperarse-estragos-trombas-agua_1_5869549.html"]},
    {"id": "2026-05-18", "titol": "18 de maig de 2026: xàfec al centre", "radar": "2026-05-18_radar.json", "ia_mm": 60,
     "resum": "Pluja intensa a migdia. L'aigua inunda l'avinguda (s'Agaró / Platja d'Aro) fins als genolls en alguns punts; comerços tancats.",
     "carrers": ["Avinguda de s'Agaró"], "llocs": [],
     "fonts": ["https://www.tvcostabrava.com/redaccio/medi-ambient/el-temps/la-pluja-intensa-provoca-les-primeres-incidencies-al-baix-emporda"]},
    {"id": "2026-01-20", "titol": "20 de gener de 2026: temporal", "radar": "2026-01-19_radar.json", "ia_mm": 40,
     "resum": "Temporal de llevant. Es tallen els carrers propers al Ridaura, el barri de la Coma i el Rec del Molí (Castell d'Aro) i la "
              "riera nord a la zona d'Aquadiver; passeig marítim tancat.",
     "carrers": ["Passeig del Ridaura", "Carrer de la Coma", "Passeig Marítim"], "llocs": ["Aquadiver (riera nord)", "Rec del Molí"],
     "fonts": ["https://www.tvcostabrava.com/redaccio/societat/seguretat-/castell-daro-platja-daro-i-sagaro-suspenen-activitats-i-tanquen-carrers-pel-temporal"]},
    {"id": "2024-09-23", "titol": "23 de setembre de 2024: xàfec de 25 mm en 30 min", "radar": "2024-09-23_radar.json", "ia_mm": 60,
     "resum": "25,3 mm en mitja hora a Castell d'Aro (18,6 en 10 min), més de 60 mm en menys d'una hora. Uns 60 comerços del centre afectats.",
     "carrers": ["Avinguda de s'Agaró"],
     "llocs": ["Rotonda de l'antic Magic Park", "Sortida 314 de la C-31 (giratori dels Xais)"],
     "fonts": ["https://www.radiocapital.cat/un-aiguat-inunda-carrers-de-platja-daro-i-fa-tallar-una-sortida-de-la-c-31-2/",
               "https://www.gerio.cat/noticia/un-fort-xafec-inunda-el-centre-de-platja-daro-afectant-baixos-i-vehicles",
               "https://www.tvcostabrava.com/laia-campos/societat/meteo/platja-daro-fa-balanc-de-les-inundacions-dahir-la-tarda"]},
    {"id": "2020-01-21", "titol": "19-23 de gener de 2020: temporal Gloria", "radar": "2020-01-19_radar.json", "ia_mm": 30,
     "resum": "Temporal marítim i de pluja de diversos dies. Les crescudes del Ridaura erosionen el talús del marge dret a Santa Cristina "
              "(camí del Golf), reparat el 2021.",
     "carrers": [], "llocs": ["Talús del Ridaura al Golf Costa Brava (Santa Cristina)"],
     "fonts": ["https://www.gerio.cat/noticia/1349055/santa-cristina-daro-inicia-les-obres-per-reparar-lescullera-del-riu-ridaura-danyat-pel-glo"]},
    {"id": "2019-10-22", "titol": "22 d'octubre de 2019: llevantada", "radar": "2019-10-22_radar.json", "ia_mm": 60,
     "resum": "A les 16:30 la Policia Local tanca tots els carrers del centre. Carrer Doctor Fleming i el centre de s'Agaró, inundats.",
     "carrers": ["Avinguda del Doctor Fleming"], "llocs": ["Centre de s'Agaró"],
     "fonts": ["https://www.gerio.cat/noticia/854586/platja-daro-tanca-els-carrers-del-centre-inundats-per-la-llevantada"]},
    {"id": "2018-11-18", "titol": "18 de novembre de 2018: temporal i crescuda del Ridaura", "radar": "2018-11-17_radar.json", "ia_mm": 60,
     "resum": "Més de 37 mm de matinada. Inundacions al nucli urbà, a l'entorn del Magic Park i de la comissaria de la Policia Local; "
              "forta pujada del Ridaura i accessos a la llera tancats; barri de la Coma (Castell d'Aro).",
     "carrers": ["Carrer de la Coma"], "llocs": ["Policia Local", "Entorn del Magic Park", "Accessos al Ridaura"],
     "fonts": ["https://news.rpa.cat/2018/11/18/el-temporal-de-pluja-provoca-inundacions-talls-de-carreteres-i-una-pujada-forta-del-cabal-del-ridaura-a-castell-platja-daro-i-sagaro/"]},
    {"id": "2025-11-06", "titol": "6 de novembre de 2025: xàfec", "radar": "2025-11-06_radar.json", "ia_mm": 60,
     "resum": "Més de 34 mm a Platja d'Aro (30 mm en 1 h a Castell d'Aro). La Policia Local talla carrers per acumulació d'aigua; sense incidents greus.",
     "carrers": ["Carrer de Juli Garreta", "Carrer de Luxemburg"], "llocs": ["Avinguda central de s'Agaró"],
     "fonts": ["https://www.radiocapital.cat/el-temporal-descarrega-amb-forca-al-baix-emporda-i-deixa-acumulacions-daigua-pero-sense-incidencies-greus/"]},
    {"id": "2025-05-05", "titol": "5 de maig de 2025: tempesta sense danys", "radar": "2025-05-05_radar.json", "ia_mm": 60, "inundacio": False,
     "resum": "23 mm en mitja hora a Castell d'Aro (31 mm en 1 h). Protocol preventiu activat, però la premsa no cita inundacions a Platja d'Aro.",
     "carrers": [], "llocs": [], "fonts": ["https://www.gerio.cat/noticia/1712351/una-forta-tempesta-inunda-els-carrers-de-palamos", "https://news.rpa.cat/2025/05/05/lajuntament-activa-el-protocol-de-prevencio-dinundacions/"]},
    {"id": "2024-02-26", "titol": "26 de febrer de 2024: tempesta amb calamarsa", "radar": "2024-02-26_radar.json", "ia_mm": 60,
     "resum": "Tempesta amb pluja i calamarsa al vespre (21 mm en 30 min a Castell d'Aro). Carrers inundats que dificulten el trànsit.",
     "carrers": [], "llocs": [], "fonts": ["https://cronicaglobal.elespanol.com/vida/20240227/un-tormenta-sacude-girona-granizo-algunas-inundaciones/835916439_0.html"]},
    {"id": "2023-09-15", "titol": "15-16 de setembre de 2023: xàfec sense danys", "radar": "2023-09-15_radar.json", "ia_mm": 60, "inundacio": False,
     "resum": "54 mm en menys d'una hora a Castell d'Aro (40 mm en 30 min) i 55 mm a Platja d'Aro, però l'episodi acaba \"sense gaires incidències\".",
     "carrers": [], "llocs": [], "fonts": ["https://govern.cat/salapremsa/notes-premsa/539362/lepisodi-pluja-deixa-30-50-mm-al-litoral-mes-100-mm-punts-del-baix-ebre", "https://www.elpuntavui.cat/societat/article/5-societat/2334172-l-episodi-de-pluges-acaba-sense-gaires-incidencies-i-amb-precipitacions-intenses-a-platja-d-aro-55l-i-palafrugell-48l.html"]},
    {"id": "2019-09-21", "titol": "21 de setembre de 2019: aiguat mortal", "radar": "2019-09-20_radar.json", "ia_mm": 60,
     "resum": "93 mm (51 mm en 1 h). Carrers, baixos, locals i aparcaments soterranis inundats fins a 2-3 pams; embornals tapats per la pinassa. "
              "Un home mor ofegat en uns baixos del carrer Ciutat de Girona.",
     "carrers": ["Avinguda Ciutat de Girona"], "llocs": ["Palau d'Esports i Congressos", "Entorn del Magic Park", "Zona comercial"],
     "fonts": ["http://www.elpuntavui.cat/societat/article/5-societat/1668851-un-home-mor-ofegat-en-els-baixos-d-un-local-inundat-per-les-fortes-tempestes-a-platja-d-aro.html", "http://www.elpuntavui.cat/punt-divers/article/4-divers/1668954-inundacions-a-l-emporda-amb-una-victima-mortal-i-el-vallvalles-oriental.html"]},
    {"id": "2018-10-10", "titol": "10 d'octubre de 2018: xàfec", "radar": "2018-10-09_radar.json", "ia_mm": 60,
     "resum": "84 mm a Castell d'Aro, uns 50 mm en menys d'una hora al matí. Inundacions i afectacions a la via pública; cotxes i aigua a s'Agaró.",
     "carrers": [], "llocs": [], "fonts": ["https://www.radiocapital.cat/una-gota-freda-responsable-daiguats-i-inundacions-al-baix-emporda/"]},
    {"id": "2016-09-14", "titol": "14 de setembre de 2016: xàfec sense danys", "radar": "2016-09-13_radar.json", "ia_mm": 60, "inundacio": False,
     "resum": "21 mm en 30 min a Castell-Platja d'Aro (33 mm en 1 h); la premsa no hi esmenta danys.", "carrers": [], "llocs": [], "fonts": ["https://bibgirona.cat/assets/documents/000/343/649/15_set_2016_-_Un_aiguat_escombra_el_territori___%C3%92._Pinilla___E._Agull%C3%B3___ACN___girona___Punt_divers___El_Punt_Avui___Not%C3%ADcia.pdf"]},
]
SENSE_RADAR = [
    {"id": "2005-10-13", "titol": "11-14 d'octubre de 2005: el Ridaura es desborda",
     "resum": "259 mm a Santa Cristina d'Aro entre l'11 i el 13 (130 mm en 24 h). Últim desbordament del Ridaura abans del 2026, "
              "danys molt greus a Castell-Platja d'Aro. Sense radar OPERA (l'arxiu comença el 2011).",
     "fonts": ["https://agora.ub.edu/wp-content/uploads/2023/04/20051011InformeTecnic_final.pdf"]},
    {"id": "2004-03-30", "titol": "30 de març de 2004: crescuda del Ridaura",
     "resum": "La crescuda malmet un pilar del pont entre la urbanització del Golf i el barri Vilà (Santa Cristina d'Aro).",
     "fonts": ["https://www.3cat.cat/3catinfo/les-fortes-pluges-enfonsen-parcialment-un-pont-de-santa-cristina-daro/noticia/62182/"]},
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"]
    osm = json.loads((d / "zones" / "osm.geojson").read_text())["features"]
    w, s, e, n = CAIXA
    per_nom = {}
    for f in osm:
        nom = f["properties"].get("name")
        if not nom:
            continue
        g = shape(f["geometry"]); c = g.centroid
        if w < c.x < e and s < c.y < n:
            per_nom.setdefault(nom, []).append(g)
    feats, falten = [], []
    for ep in EPISODIS:
        for nom in ep["carrers"] + [l.split(" (")[0] for l in ep["llocs"]]:
            if nom in per_nom:
                feats.append({"type": "Feature", "properties": {"episodi": ep["id"], "nom": nom},
                              "geometry": mapping(unary_union(per_nom[nom]).simplify(0.00002))})
            else:
                falten.append(f"{ep['id']}: {nom}")
    (d / "episodis" / "premsa.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                              ensure_ascii=False, separators=(",", ":")))
    def te_radar(ep):                     # el fitxer existeix i no és buit (si el format de l'arxiu no s'ha pogut llegir, tot és 0)
        f = d / "episodis" / ep["radar"]
        return f.exists() and any(v > 0 for fila in json.loads(f.read_text())["mm"] for v in fila)
    disponibles = [ep for ep in EPISODIS if te_radar(ep)]
    (d / "episodis" / "index.json").write_text(json.dumps({"episodis": disponibles, "sense_radar": SENSE_RADAR,
                                                           "pendents_de_radar": [ep["id"] for ep in EPISODIS if ep not in disponibles]},
                                                          ensure_ascii=False, indent=1))
    print(f"{len(feats)} carrers i llocs situats; sense situar: {falten}")
    print(f"episodis amb radar: {[ep['id'] for ep in disponibles]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
