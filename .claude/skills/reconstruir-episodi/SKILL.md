---
name: reconstruir-episodi
description: Reconstruir un episodi de pluja passat en una conca (pluja corregida, absorció, cabal, simulació 2D) i comparar-lo amb aforadors, premsa i intervencions. Fes-la servir quan es demani "reconstruir", "simular" o "validar" un episodi concret.
---

# Reconstruir un episodi

Regla: **res ajustable a mà**. Tots els valors surten de la configuració de la conca o del càlcul; si un resultat no quadra,
s'explica per què, no es retoca un paràmetre per fer-lo quadrar.

## Passos

1. **Situa l'episodi**: conca (`config/<conca>.json`), inici i final en UTC. Mira si ja existeix a
   `data/<conca>/episodis/` (`*_radar.json`, `*_estacio.json`) o a `data/<conca>/sim2d/<id>/`. No ho tornis a calcular si hi és.
2. **Pluja**:
   - Mes corregit a tot Catalunya: `data/pluja_cat/index.json` diu si el mes ja és corregit (`correccio_cat.py --mes AAAA-MM`).
   - Radar a la conca: `python3 scripts/radar_graella.py --config config/<conca>.json --des <ini> --fins <fi>`.
   - Correcció amb estacions: `python3 scripts/correccio_estacio.py --config ... --radar <AAAA-MM-DD>_radar.json --id <id>`.
   - Comprova el biaix (radar vs. pluviòmetres) i anota'l; si és > ±30 %, digues-ho abans de seguir.
3. **Absorció inicial**: `python3 scripts/humitat.py ia --config config/<conca>.json --t0 <ini>` (dipòsit d'humitat, pluja
   dels 60 dies anteriors). No posis una Ia fixa.
4. **Cabal**: si hi ha aforador, compara el pic i l'hora del pic (`observacions` a Supabase o arxiu ACA amb
   `extreu_aforador.py`). Reprodueix la previsió com si fos en directe:
   `python3 scripts/previsio.py --config ... --ara <hora> --sense-desar` (antelació de l'avís de T10).
5. **Simulació 2D**: llança `sim2d.yml` (workflow_dispatch) amb `episodis=[{"id","radar","episodi","ia"}]`. Opcions noves
   d'una en una i sempre comparades amb la sim per defecte. Desa cada prova a la Release `sim2d-proves-*`.
6. **Resum**: `python3 scripts/resum.py data/<conca>/sim2d/<id>` (o llegeix `resum.md`). No obris els `.u2` ni els PNG.
7. **Validació**: segueix la skill `validar-intervencions` (premsa, Bombers, carrers).

## Què lliurar

Un informe curt (MD) amb: pluja total i màxima en 30 min/1 h, Ia, pic simulat vs. observat (m³/s i hora), període de
retorn, km² inundats > 15 cm, soterranis en risc, encerts i errors amb premsa/intervencions, i què falla i per què.
Taules, no sortides de consola.
