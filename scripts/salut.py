"""Control de salut de la cadena en directe (Ridaura): avisa quan una font deixa d'arribar, quan torna, i dels buits.

Comprova a Supabase l'última dada de cada font (aforador ACA, radar, radar corregit, previsió) i els buits del radar
de les últimes 24 h. Avisa només quan l'estat canvia (falla / es recupera) i de cada buit una sola vegada; l'estat es
desa en un fitxer JSON (--estat) que el workflow conserva entre execucions amb la memòria cau d'Actions.
Si el mateix Supabase no respon, també avisa.
"""
import argparse, json, sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
import requests

sys.path.insert(0, str(Path(__file__).parent))
import sb
from avisa import envia

PAS_RADAR = 5   # minuts
BUIT_MIN = 30   # un forat de radar d'aquesta durada o més es notifica

# nom, taula, filtres, columna de temps, retard màxim acceptable (min)
FONTS = [
    ("Aforador ACA (cabal)", "observacions", {"var": "eq.cabal_m3s"}, "temps", 60),
    ("Radar OPERA", "pluja_conca", {"font": "eq.radar", "subconca": "eq.aforador"}, "temps", 60),
    ("Radar corregit", "pluja_conca", {"font": "eq.radar_corregit", "subconca": "eq.aforador"}, "temps", 120),
    ("Previsió cada 10 min", "previsions_cabal", {"punt": "eq.aforador"}, "run", 60),
]


def hora(t):
    """UTC → hora local de Catalunya (text curt)."""
    try:
        from zoneinfo import ZoneInfo
        return t.astimezone(ZoneInfo("Europe/Madrid")).strftime("%d/%m %H:%M")
    except Exception:
        return t.strftime("%d/%m %H:%M UTC")


def darrera(taula, filtres, col, conca):
    f = sb.llegeix(taula, {"select": col, "conca": f"eq.{conca}", **filtres, "order": f"{col}.desc", "limit": 1})
    return datetime.fromisoformat(f[0][col]) if f else None


def buits_radar(conca, ara):
    ini = ara - timedelta(hours=24)
    f = sb.llegeix("pluja_conca", {"select": "temps", "conca": f"eq.{conca}", "font": "eq.radar", "subconca": "eq.aforador",
                                   "temps": f"gte.{ini.isoformat()}", "order": "temps.asc", "limit": 2000})
    t = [datetime.fromisoformat(r["temps"]) for r in f]
    return [(a, b) for a, b in zip(t, t[1:]) if (b - a) >= timedelta(minutes=BUIT_MIN + PAS_RADAR)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--conca", default="ridaura")
    ap.add_argument("--estat", default="salut_estat.json")
    ap.add_argument("--ara", help="hora UTC (YYYY-MM-DDTHH:MM), per provar")
    a = ap.parse_args()
    ara = datetime.fromisoformat(a.ara).replace(tzinfo=timezone.utc) if a.ara else datetime.now(timezone.utc)
    pe = Path(a.estat)
    estat = json.loads(pe.read_text()) if pe.exists() else {}
    fallen, buits_avisats = estat.get("fallen", {}), set(estat.get("buits", []))
    nous_fallen, recuperats, linies = {}, [], []

    # sb.llegeix fa sys.exit si Supabase falla: ho capturem com a una font més
    try:
        for nom, taula, filtres, col, max_min in FONTS:
            t = darrera(taula, filtres, col, a.conca)
            retard = (ara - t).total_seconds() / 60 if t else None
            ok = retard is not None and retard <= max_min
            linies.append(f"{nom}: {'OK' if ok else 'FALLA'} · última {hora(t) if t else 'cap'}")
            if not ok:
                nous_fallen[nom] = fallen.get(nom) or (t.isoformat() if t else "")
            elif nom in fallen:
                recuperats.append(f"{nom} (aturat des de {hora(datetime.fromisoformat(fallen[nom])) if fallen[nom] else '?'})")
        buits = buits_radar(a.conca, ara)
    except (SystemExit, requests.RequestException) as e:
        nous_fallen = {"Supabase": fallen.get("Supabase") or ara.isoformat()}
        linies.append(f"Supabase no respon: {e}")
        buits = []
        for nom in fallen:       # no sabem res de la resta: es queden com estaven
            nous_fallen.setdefault(nom, fallen[nom])
    else:
        if "Supabase" in fallen:
            recuperats.append("Supabase")

    nous_buits = [(x, y) for x, y in buits if x.isoformat() not in buits_avisats]
    print("\n".join(linies))

    noves = [n for n in nous_fallen if n not in fallen]
    if noves:
        envia("S'ha aturat: " + ", ".join(noves) + ".\n\n" + "\n".join(linies), "⚠️ Inundacions: la cadena en directe falla")
    if recuperats:
        envia("Torna a funcionar: " + ", ".join(recuperats) + ".", "✅ Inundacions: cadena recuperada")
    if nous_buits:
        envia("Buits al radar (24 h): " + "; ".join(f"{hora(x)} → {hora(y)} ({(y - x).total_seconds() / 60 - PAS_RADAR:.0f} min)"
                                                   for x, y in nous_buits) + ".", "⚠️ Inundacions: forats al radar")

    pe.write_text(json.dumps({"fallen": nous_fallen,
                              "buits": sorted({x.isoformat() for x, _ in buits} | {b for b in buits_avisats
                                              if datetime.fromisoformat(b) > ara - timedelta(hours=48)}),
                              "comprovat": ara.isoformat()}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
