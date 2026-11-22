#!/usr/bin/env python3
"""Passos de torrentets sota vies i ponts dins la simulació 2D (opció --passos de simula_2d.py).

L'MDT de l'IGN és el terreny nu: un terraplè de carretera o de tren hi queda com un dic massís i l'obra de drenatge no
s'hi veu. Per això, sense aquesta opció, la simulació 2D tracta tots els passos com a TAPATS (és el comportament per
defecte i no canvia). Amb l'opció, cada pas detectat per scripts/passos_cat.py (data/catalunya/passos/cNNN.geojson,
precalculat per a tot Catalunya) es converteix en un conducte que travessa el terraplè:

  - Entrada: la cel·la més baixa del llac que es faria si l'obra es tapés (contorn del llac de passos_cat.py).
  - Sortida: la cel·la més propera a l'altra banda del terraplè, fora del llac, més baixa que l'entrada, dins
    l'amplada de la base del terraplè: calçada de 10 m (carril 3,5 m × 2 + vorals, norma 3.1-IC) + talussos 3H:2V a
    banda i banda (1,5 × fondària del llac cada un).
  - Cabal: el d'un orifici amb la diferència de nivell d'aigua entre les dues bandes (Δh):
      · si es coneix la secció real de l'obra (propietat seccio_m2 o diametre_m del pas): Q = Cd·A·√(2gΔh), Cd = 0,6;
      · si no, amb el cabal de disseny de la norma 5.2-IC que passos_cat.py ja calcula (escenaris.obert.q_m3s):
        l'obra porta aquest cabal quan el llac és ple (Δh = fondària del llac) i Q = Qd·√(Δh/fondària), que és la
        mateixa llei d'orifici sense haver de saber la mida;
      · sense cap de les dues (fora de les conques pilot, fins que hi hagi la mida real), el pas queda tapat i es diu.
  - Taponament: Q × (1 − bloqueig). Modes:
      · oberts: bloqueig 0 a tots els passos;
      · arrossegalls: bloqueig de disseny de l'ARR 2019 (scripts/arrossegalls.py, data/catalunya/arrossegalls),
        segons el potencial d'arrossegalls de la conca del pas i l'obertura (petita / mitjana, o la real si se sap).
        També als ponts d'OSM: les cares de la cel·la del pont perden la fracció de pas bloquejada.
    (tapats: cap conducte; és el que fa la simulació sense l'opció.)

L'aigua que passa pel conducte surt d'una cel·la i entra a l'altra (el volum es conserva) i mai més del que cal per
igualar el nivell a banda i banda en un pas de temps."""
import json, math
from pathlib import Path
import numpy as np
from numba import njit

G = 9.81
CD = 0.6                       # coeficient de descàrrega d'un orifici (entrada de conducte)
CALCADA = 10.0                 # m: 2 carrils de 3,5 m + vorals (3.1-IC)
TALUS = 1.5                    # talús 3H:2V
L10 = 1.5                      # ARR 2019: llargada típica dels arrossegalls (arrossegalls.py)
PASSOS = Path("data/catalunya/passos")
ARROSS = Path("data/catalunya/arrossegalls")


def _classe_obertura(w):
    return "petita" if w < L10 else ("mitjana" if w <= 3 * L10 else "gran")


def prepara(mode, obertura, tr, nf, nc, crs, z, mar, psx, psy, dir_passos=PASSOS, dir_arr=ARROSS):
    """Conductes dels passos dins el domini. Modifica psx/psy als ponts (mode arrossegalls).
    Torna (taula de conductes per a transfereix(), llista d'informació per pas, resum)."""
    from pyproj import Transformer
    from rasterio import features
    from scipy import ndimage
    from shapely.geometry import shape
    from shapely.ops import transform as stransf
    to_c = Transformer.from_crs(4326, crs, always_xy=True)
    dx = float(tr.a)
    x0, y1 = tr.c, tr.f; x1, y0 = x0 + nc * dx, y1 - nf * dx
    ll = Transformer.from_crs(crs, 4326, always_xy=True)
    lw, ls = ll.transform(x0, y0); le, ln = ll.transform(x1, y1)

    def dins(lon, lat):
        return lw - 0.01 <= lon <= le + 0.01 and ls - 0.01 <= lat <= ln + 0.01

    passos, llacs = [], {}
    for f in sorted(Path(dir_passos).glob("c[0-9]*.geojson")):
        for ft in json.loads(f.read_text())["features"]:
            g, p = ft["geometry"], ft["properties"]
            if g["type"] == "Point":
                if p.get("tipus") == "terraplè" and dins(*g["coordinates"]):
                    passos.append((g["coordinates"], p))
            elif "llac_de" in p:
                llacs.setdefault(p["llac_de"], []).append(g)
    arr = {}
    if mode == "arrossegalls":
        fitxers = sorted(Path(dir_arr).glob("c[0-9]*.geojson"))
        if not fitxers:
            raise SystemExit(f"--passos arrossegalls: no hi ha {dir_arr}/cNNN.geojson (scripts/arrossegalls.py)")
        for f in fitxers:
            for ft in json.loads(f.read_text())["features"]:
                p = ft["properties"]
                if p.get("tipus") in ("terraplè", "pont") and dins(*ft["geometry"]["coordinates"]):
                    arr[p["id"]] = (ft["geometry"]["coordinates"], p)

    def cella(lon, lat):
        x, y = to_c.transform(lon, lat)
        return int((y1 - y) // dx), int((x - x0) // dx)

    taula, info = [], []
    cnt = {"passos": len(passos), "conductes": 0, "sense_capacitat": 0, "ja_passa_a_la_reixa": 0, "fora": 0}
    for (lon, lat), p in passos:
        r, c = cella(lon, lat)
        if not (1 <= r < nf - 1 and 1 <= c < nc - 1) or mar[r, c]:
            cnt["fora"] += 1; continue
        b, ob = 0.0, None
        A = float(p.get("seccio_m2") or 0) or (math.pi * float(p["diametre_m"]) ** 2 / 4 if p.get("diametre_m") else 0.0)
        qd = ((p.get("escenaris") or {}).get("obert") or {}).get("q_m3s")
        if mode == "arrossegalls":
            ob = _classe_obertura(float(p["diametre_m"])) if p.get("diametre_m") else obertura
            bl = (arr.get(p["id"], (None, {}))[1].get("bloqueig") or {})
            b = float(bl.get(ob, bl.get("petita", 0.0)) or 0.0)
        reg = {"id": p["id"], "via": p.get("via"), "nom": p.get("nom"), "conca_km2": p.get("conca_km2"), "bloqueig": b, "obertura": ob}
        if not A and not qd:
            cnt["sense_capacitat"] += 1; reg["estat"] = "tapat (sense mida ni cabal de disseny)"; info.append(reg); continue
        # entrada: cel·la més baixa del llac (a la reixa de la simulació)
        m = np.zeros((nf, nc), bool)
        if p["id"] in llacs:
            geoms = [stransf(to_c.transform, shape(g)) for g in llacs[p["id"]]]
            m = features.rasterize([(g, 1) for g in geoms], out_shape=(nf, nc), transform=tr, all_touched=True).astype(bool)
        m[r, c] = True
        m &= ~mar
        rr, cc = np.nonzero(m)
        k = int(np.argmin(z[rr, cc])); ru, cu = int(rr[k]), int(cc[k])
        H = max(float(p.get("prof_m") or 0), 0.1)
        K = max(int(math.ceil((CALCADA + 2 * TALUS * H) / dx)), 1)
        ra, rb, ca, cb = max(ru - K - 1, 1), min(ru + K + 2, nf - 1), max(cu - K - 1, 1), min(cu + K + 2, nc - 1)
        sub_m = m[ra:rb, ca:cb]
        dist = ndimage.distance_transform_edt(~sub_m) * dx                 # distància al llac
        cand = (~sub_m) & (~mar[ra:rb, ca:cb]) & (z[ra:rb, ca:cb] < z[ru, cu] - 1e-3) & (dist <= K * dx + 1e-6)
        if not cand.any():
            cnt["ja_passa_a_la_reixa"] += 1
            reg["estat"] = "sense desnivell a la reixa (el terraplè ja no fa de dic a aquesta mida de cel·la)"; info.append(reg); continue
        clau = np.where(cand, dist * 1e4 + z[ra:rb, ca:cb], np.inf)
        kd = np.unravel_index(np.argmin(clau), clau.shape); rd, cd = ra + int(kd[0]), ca + int(kd[1])
        taula.append((ru, cu, rd, cd, A, float(qd or 0), H, b))
        reg.update({"estat": "conducte", "entrada": [ru, cu], "sortida": [rd, cd], "seccio_m2": A or None, "q_disseny_m3s": qd,
                    "fondaria_m": H})
        info.append(reg)
        cnt["conductes"] += 1
    # ponts: el bloqueig treu pas a les cares de la cel·la del pont
    ponts = 0
    for pid, ((lon, lat), p) in arr.items():
        if p.get("tipus") != "pont":
            continue
        b = float(next(iter((p.get("bloqueig") or {"x": 0.0}).values())) or 0.0)
        r, c = cella(lon, lat)
        if b <= 0 or not (1 <= r < nf - 1 and 1 <= c < nc - 1):
            continue
        f = 1 - b
        psx[r, c] *= f; psx[r, c - 1] *= f; psy[r, c] *= f; psy[r - 1, c] *= f
        ponts += 1
        info.append({"id": pid, "tipus": "pont", "nom": p.get("nom"), "llum_m": p.get("llum_m"), "bloqueig": b, "cella": [r, c]})
    cnt["ponts_amb_bloqueig"] = ponts
    t = np.array(taula, dtype=np.float64).reshape(-1, 8)
    print(f"Passos ({mode}{'' if mode != 'arrossegalls' else ', obertura ' + obertura}): {cnt}")
    return t, info, cnt


@njit(cache=True)
def transfereix(t, z, h, phi, dx, dt, qmax, vol):
    """Un pas de temps dels conductes: orifici entre l'entrada i la sortida, en qualsevol sentit."""
    a = dx * dx
    for k in range(t.shape[0]):
        ru, cu, rd, cd = int(t[k, 0]), int(t[k, 1]), int(t[k, 2]), int(t[k, 3])
        A, qd, H, b = t[k, 4], t[k, 5], t[k, 6], t[k, 7]
        if b >= 1.0:
            continue
        dh = (z[ru, cu] + h[ru, cu]) - (z[rd, cd] + h[rd, cd])
        if dh > 0:
            rs, cs, rt, ct = ru, cu, rd, cd
        else:
            rs, cs, rt, ct = rd, cd, ru, cu
            dh = -dh
        if dh < 1e-4 or h[rs, cs] <= 1e-4:
            continue
        q = CD * A * np.sqrt(2 * G * dh) if A > 0 else qd * np.sqrt(dh / H)
        q *= 1.0 - b
        As, At = phi[rs, cs] * a, phi[rt, ct] * a
        v = min(q * dt, h[rs, cs] * As, 0.5 * dh / (1.0 / As + 1.0 / At))   # mai més del que iguala els nivells
        h[rs, cs] -= v / As
        h[rt, ct] += v / At
        if v / dt > qmax[k]:
            qmax[k] = v / dt
        vol[k] += v if rs == ru and cs == cu else -v
