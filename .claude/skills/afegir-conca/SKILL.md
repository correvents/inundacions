---
name: afegir-conca
description: Afegir una conca o població nova al sistema (configuració, pluja, relleu, terreny precalculat, episodis i calibratge). Fes-la servir quan es demani ampliar el model a un poble o riu nou.
---

# Afegir una conca nova

## 1. Triar i configurar

- Millor si té **aforador** (ACA), episodis recents i premsa. Sense aforador, la validació és només per carrers.
- Busca a quina conca del tall de 100 km² (`data/catalunya/conques`, `mapa_conques.py`) cau: d'aquí surten els fitxers
  precalculats (MDT, MET, cobertes, sòls, Cadastre, passos, clavegueram) a les Releases. **No els recalculis**: baixa'ls.
- Crea `config/<conca>.json` copiant `config/granollers.json`: `conca`, `nom`, `centre`, `caixa`, `caixa_radar`,
  `aforador` (codi ACA), `punts`, cabals Q10-Q500 del tram si n'hi ha (ACA/SNCZI), `estacio_api` (buit = automàtica).
  Cada valor amb la seva font a `nota`.

## 2. Pluja, relleu i episodis

- Llista d'episodis amb data i font (premsa, Bombers, estacions XEMA): `[{"id","titol","resum","des","fins"}]` en UTC.
- Llança `conca_nova.yml` (`conca`, `episodis`): estacions XEMA properes, radar i correcció per episodi, índex,
  relleu i graella del visor. Resultat a `data/<conca>/`.
- Si hi ha aforador, comprova que surt a l'arxiu mensual (`aca_arxiu.yml` desa tots els sensors; l'ACA només en guarda
  ~100 dies) i a `observacions` de Supabase.

## 3. Calibrar (per conca, mai global)

- Absorció: `python3 scripts/humitat.py calibra --config config/<conca>.json` → `config/<conca>_humitat.json`.
- Calibra amb **un** episodi i valida amb **un altre** que no s'hagi fet servir. Si només n'hi ha un, digues-ho.
- Paràmetres que toquen els límits del rang (p. ex. C = 1) són senyal d'error en dades, no de calibratge bo.

## 4. Simulació 2D i web

- `sim2d.yml` amb la conca nova i els seus episodis (vegeu la skill `reconstruir-episodi`). Amb MET, pas 2-5 m.
- Afegeix la conca a la web només amb fitxers estàtics de `data/<conca>/` (la web no calcula). Les dades grans van a
  Releases (Pages té un límit d'1 GB).
- Actualitza `docs/replicar.md` si canvia algun pas de la cadena.
