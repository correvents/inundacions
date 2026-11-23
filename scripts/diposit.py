#!/usr/bin/env python3
"""Edificis com a dipòsit: quanta aigua hi pot entrar (planta baixa i soterranis) i per on, a partir del Cadastre.

Opció de la simulació 2D (simula_2d.py --edificis-diposit). Per defecte els edificis només són obstacles (porositat) i
l'aigua no hi entra; amb l'opció, cada edifici és un dipòsit connectat al carrer per les seves obertures:
  - Capacitat: soterranis (àrea de cada part × plantes sota rasant × ALCADA_UTIL) i planta baixa (àrea de les parts amb
    plantes sobre rasant × ALCADA_UTIL), totes dues per la fracció útil F_UTIL (la resta són murs, pilars, envans i trasters
    tancats). L'aigua omple primer el soterrani (rampes i escales) i després puja a la planta baixa.
  - Obertures: OBERTURA_PORTAL m per edifici (portal o porta de botiga) i OBERTURA_RAMPA m més si té soterrani (rampa de
    garatge). Els edificis grans (indústria, comerç) tenen una porta més per cada PERIMETRE_PER_PORTA m de façana.
  - Entrada: abocador de llavi gruixut, Q = CD · b · √(2g) · H^1,5, amb H el calat del carrer per sobre del terra de
    l'edifici més el LLINDAR del portal; si a dins ja hi ha aigua a la planta baixa, factor d'abocador negat (Villemonte 1947).
    L'aigua de la planta baixa torna a sortir quan el carrer baixa; la del soterrani s'hi queda (cal bombar-la).
Tot es precalcula per edifici i per a tot Catalunya (Release `diposits-cat`); a la simulació només s'hi afegeix la cota
del terra (mediana del MDT sota la planta) i la cel·la de carrer on obre.

Ús:  python3 scripts/diposit.py conca data/ridaura/cadastre/edificis.json.gz          → data/ridaura/cadastre/diposits.json.gz
     python3 scripts/diposit.py cat /tmp/cadastre_cat/edificis_*.json.gz --sortida /tmp/diposits_cat"""
import argparse, gzip, json, sys
from pathlib import Path
import numpy as np

ALCADA_UTIL = 2.5            # m d'alçada lliure per planta (CTE DB-SUA: ≥ 2,2 m garatges, ≥ 2,5 m habitatge)
F_UTIL = 0.8                 # superfície útil / construïda (habitual en habitatge i garatge, 0,75-0,85)
OBERTURA_PORTAL = 1.0        # m (porta de portal o de botiga)
OBERTURA_RAMPA = 3.0         # m (porta de garatge, CTE DB-SI: rampa ≥ 3 m)
PERIMETRE_PER_PORTA = 40.0   # m de façana per porta addicional en edificis grans
LLINDAR = 0.10               # m entre la vorera i el terra de la planta baixa (graó del portal)
CD = 0.40                    # coeficient d'abocador de llavi gruixut (Q = CD b √(2g) H^1,5)
CAMPS = ["ref", "x", "y", "area_pb_m2", "vol_soterrani_m3", "obertura_m", "plantes_sota"]


def per_edifici(parts):
    """[ref, x, y, àrea útil planta baixa, volum útil soterrani, obertura, plantes sota] per edifici (ref cadastral)."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    agr = {}
    for ref, sobre, sota, an in parts:
        if len(an[0]) < 4:
            continue
        p = Polygon(an[0], an[1:]).buffer(0)
        if p.is_empty:
            continue
        r = agr.setdefault(ref, {"pb": [], "sot": 0.0, "nsot": 0, "tot": []})
        r["tot"].append(p)
        if sobre >= 1:
            r["pb"].append(p)
        if sota > 0:
            r["sot"] += p.area * sota
            r["nsot"] = max(r["nsot"], sota)
    out = []
    for ref, r in agr.items():
        tot = unary_union(r["tot"])
        pb = unary_union(r["pb"]).area if r["pb"] else 0.0
        if pb <= 0 and r["sot"] <= 0:
            continue
        c = tot.representative_point()
        per = tot.length
        ob = OBERTURA_PORTAL * max(1, int(per // PERIMETRE_PER_PORTA)) + (OBERTURA_RAMPA if r["sot"] > 0 else 0)
        out.append([ref, round(c.x, 1), round(c.y, 1), round(pb * F_UTIL, 1), round(r["sot"] * ALCADA_UTIL * F_UTIL, 1),
                    round(ob, 1), r["nsot"]])
    return out


def desa(fitxer, files, extra=None):
    with gzip.open(fitxer, "wt", compresslevel=9) as g:
        json.dump({"crs": "EPSG:25831", "camps": CAMPS, **(extra or {}), "hipotesis": hipotesis(), "edificis": files},
                  g, separators=(",", ":"), ensure_ascii=False)


def hipotesis():
    return {"alcada_util_m": ALCADA_UTIL, "fraccio_util": F_UTIL, "obertura_portal_m": OBERTURA_PORTAL,
            "obertura_rampa_m": OBERTURA_RAMPA, "perimetre_per_porta_m": PERIMETRE_PER_PORTA, "llindar_m": LLINDAR, "cd": CD,
            "font": "Cadastre (INSPIRE Buildings), plantes sobre i sota rasant"}


def llegeix(fitxer):
    d = json.loads(gzip.open(fitxer, "rt").read())
    return d["edificis"]


class Diposits:
    """Estat dels dipòsits durant la simulació. Cada edifici obre a una cel·la de carrer `cel` (índex pla de la reixa)."""

    def __init__(self, files, tr, nf, nc, phi, desti, z_mdt, tr_mdt, llindar=LLINDAR, cd=CD):
        from scipy.ndimage import distance_transform_edt
        a = np.array([f[1:6] for f in files], dtype=float)
        self.ref = [f[0] for f in files]
        x, y, self.apb, self.vsot, self.ob = a.T
        dx = tr.a
        c = ((x - tr.c) // dx).astype(int); r = ((tr.f - y) // dx).astype(int)
        dins = (c >= 0) & (c < nc) & (r >= 0) & (r < nf)
        self.ref = [k for k, d in zip(self.ref, dins) if d]
        x, y, c, r = x[dins], y[dins], c[dins], r[dins]
        self.apb, self.vsot, self.ob = self.apb[dins], self.vsot[dins], self.ob[dins]
        cel = r * nc + c
        if desti is not None:
            cel = desti[cel]                                 # cel·la tancada per edificis → el carrer més proper
        self.cel = cel
        self.uc, self.inv = np.unique(cel, return_inverse=True)   # cel·les de carrer amb algun edifici
        # cota del terra: MDT original (5 m) al punt representatiu de la planta, mediana 3×3
        mx = ~tr_mdt
        cc, rr = mx * (x, y)
        cc, rr = cc.astype(int), rr.astype(int)
        H, W = z_mdt.shape
        mostres = []
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                mostres.append(z_mdt[np.clip(rr + di, 0, H - 1), np.clip(cc + dj, 0, W - 1)])
        with np.errstate(all="ignore"):
            self.terra = np.nanmedian(np.array(mostres, dtype=float), axis=0)
        ok = np.isfinite(self.terra)
        self.terra[~ok] = 1e9                                # sense MDT (mar): mai hi entra aigua
        self.llindar, self.cd = llindar, cd
        self.cap = self.vsot + self.apb * ALCADA_UTIL
        self.V = np.zeros(len(self.cel))                    # m³ dins de cada edifici
        self.vmax = np.zeros(len(self.cel))
        self.t_entra = np.full(len(self.cel), np.inf)       # minut en què hi comença a entrar aigua
        self.entrat = 0.0; self.sortit = 0.0                # m³ acumulats
        self.n = nf * nc
        print(f"Dipòsits: {len(self.cel)} edificis, capacitat {self.cap.sum() / 1e6:.3f} hm³ "
              f"(soterranis {self.vsot.sum() / 1e6:.3f}), obertures {self.ob.sum() / 1000:.1f} km")

    def pas(self, z, h, phi, dx, dt, t_min):
        """Intercanvi carrer ↔ edifici en un pas de temps. Modifica h (calat al carrer) al lloc."""
        zf, hf, pf = z.ravel(), h.ravel(), phi.ravel()
        cota = self.terra + self.llindar
        Ho = zf[self.cel] + hf[self.cel] - cota              # càrrega de fora sobre el llindar
        hpb = np.maximum(self.V - self.vsot, 0) / np.maximum(self.apb, 1.0)
        Hi = np.where(self.V > self.vsot, hpb - self.llindar, -1.0)   # nivell de dins sobre el llindar (soterrani ple)
        k = self.cd * self.ob * np.sqrt(2 * 9.81)
        Q = np.zeros(len(self.cel))
        ent = Ho > 0.005
        Hmx = np.maximum(Ho, Hi)
        sub = np.clip(np.minimum(np.maximum(Ho, 0), np.maximum(Hi, 0)) / np.maximum(Hmx, 1e-6), 0, 1)
        f = (1 - sub ** 1.5) ** 0.385                        # abocador negat (Villemonte)
        dentro = ent & (Ho > Hi)
        Q[dentro] = k[dentro] * Ho[dentro] ** 1.5 * f[dentro]
        surt = (Hi > np.maximum(Ho, 0) + 0.005)
        Q[surt] = -k[surt] * Hi[surt] ** 1.5 * f[surt]
        dV = Q * dt
        dV = np.minimum(dV, self.cap - self.V)              # no passa de la capacitat
        dV = np.maximum(dV, -np.maximum(self.V - self.vsot, 0))   # només surt la de la planta baixa
        if not dV.any():
            return
        # l'aigua que entra no pot ser més de la meitat de la que hi ha al carrer en aquest pas
        disp = hf[self.uc] * pf[self.uc] * dx * dx * 0.5
        dem = np.bincount(self.inv, weights=np.maximum(dV, 0), minlength=len(self.uc))
        fac = np.where(dem > disp, disp / np.maximum(dem, 1e-12), 1.0)
        dV = np.where(dV > 0, dV * fac[self.inv], dV)
        nou = (dV > 0) & ~np.isfinite(self.t_entra)
        self.t_entra[nou] = t_min
        self.V += dV
        np.maximum(self.vmax, self.V, out=self.vmax)
        self.entrat += float(dV[dV > 0].sum()); self.sortit += float(-dV[dV < 0].sum())
        hf[self.uc] = np.maximum(hf[self.uc] - np.bincount(self.inv, weights=dV, minlength=len(self.uc)) / (pf[self.uc] * dx * dx), 0)

    def resultat(self):
        """{ref: [cm màx a la planta baixa, % soterrani ple, minut d'entrada, m³ màx]} dels edificis on ha entrat aigua."""
        out = {}
        for i in np.nonzero(self.vmax > 0.5)[0]:
            pb = max(self.vmax[i] - self.vsot[i], 0) / max(self.apb[i], 1.0)
            sot = 100 * min(self.vmax[i], self.vsot[i]) / self.vsot[i] if self.vsot[i] > 0 else None
            out[self.ref[i]] = [int(round(pb * 100)), None if sot is None else int(round(sot)),
                                int(self.t_entra[i]) if np.isfinite(self.t_entra[i]) else None, int(round(self.vmax[i]))]
        return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    a1 = sub.add_parser("conca"); a1.add_argument("edificis"); a1.add_argument("--sortida")
    a2 = sub.add_parser("cat"); a2.add_argument("fitxers", nargs="+"); a2.add_argument("--sortida", default="/tmp/diposits_cat")
    a = ap.parse_args()
    if a.ordre == "conca":
        files = per_edifici(json.loads(gzip.open(a.edificis, "rt").read())["parts"])
        sort = Path(a.sortida or Path(a.edificis).with_name("diposits.json.gz"))
        desa(sort, files)
        print(f"{len(files)} edificis, planta baixa {sum(f[3] for f in files) / 1e6:.2f} km² útils, "
              f"soterranis {sum(f[4] for f in files) / 1e6:.3f} hm³ → {sort}")
        return 0
    dest = Path(a.sortida); dest.mkdir(parents=True, exist_ok=True)
    for f in a.fitxers:
        d = json.loads(gzip.open(f, "rt").read())
        codi = Path(f).name.split("_")[1].split(".")[0]
        files = per_edifici(d["parts"])
        desa(dest / f"diposits_{codi}.json.gz", files, {"municipi": d.get("municipi", "")})
        print(f"{codi} {d.get('municipi', '')}: {len(files)} edificis, soterranis {sum(x[4] for x in files) / 1e6:.3f} hm³")
    return 0


if __name__ == "__main__":
    sys.exit(main())
