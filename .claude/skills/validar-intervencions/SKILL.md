---
name: validar-intervencions
description: Validar la simulació o les regles de risc d'una conca amb dades independents (premsa, actuacions de Bombers en dades obertes, carrers citats). Fes-la servir quan es demani comprovar si el model encerta on i quan es va inundar.
---

# Validar amb intervencions i premsa

## Fonts (independents del model)

- **Premsa**: `scripts/premsa.py` → `data/<conca>/episodis/premsa.geojson` i `model/premsa_episodis.json` (carrers citats
  situats a OSM).
- **Bombers, dades obertes** (dataset g2ay-3vnj): `scripts/bombers.py --config ...` → `model/bombers_aigua_dia.csv`
  (actuacions per aigua per dia i municipi). És una etiqueta de "dia amb inundació", no un punt exacte.
- **Dades no obertes**: no entren mai al repo, a Pages ni a Releases. Es processen fora del repo i només se'n publiquen
  resultats agregats (comptes per dia, per municipi o per tipus de lloc), sense adreces, text lliure ni noms.

## Passos

1. **Per carrers**: `python3 scripts/valida_sim2d.py --config ...` → `model/sim2d_carrers.csv`. Per a cada carrer de la
   premsa, el rang percentil del seu calat (100 = el més inundat). Bo si els citats queden al percentil ≥ 80.
2. **Per dies**: creua els dies amb actuacions per aigua amb la pluja (màx. 30 min, 1 h, episodi) i amb les regles de la
   ciutat (`config/<conca>_llindars_*.json`). Taula de contingència: encerts, falses alarmes, inundacions no avisades.
3. **Per tipus de lloc** (si hi ha coordenades): passos inferiors (`punts_baixos`), soterranis i plantes baixes
   (`soterranis.json` de la sim) vs. llocs amb actuacions. Distància al punt simulat més proper amb > 15 cm.
4. **Hora**: compara l'hora d'arribada simulada (`minut_15cm.u2`, via `resum.py`) amb l'hora de les primeres actuacions.

## Regles

- No ajustis el model per encertar la validació: si falla, documenta on i per què (relleu, clavegueram, pluja).
- Separa sempre episodis de calibratge i de validació.
- Lliura una taula curta (MD/CSV) per episodi: encerts, errors, percentils, i la causa probable de cada error.
