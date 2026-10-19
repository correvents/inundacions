# Inundacions: previsió de crescudes per conca

Sistema de previsió de crescudes a les rieres i rius de Catalunya. Combina la pluja observada (radar), la pluja prevista (models) i el nivell actual dels aforadors per estimar on serà l'aigua a cada moment.

Inspirat en el visor ONOFRE del GRAF (pauguarquec/echotops-viewer), que mostra les dades però no fa modelització hidrològica.

## Conca pilot: Ridaura (Platja d'Aro)

| | A l'aforador (Santa Cristina d'Aro, ACA 171812-001) | A la desembocadura (Platja d'Aro) |
|---|---|---|
| Superfície | 35,3 km² | 70,7 km² |
| Temps de concentració (Témez) | ~4,2 h | ~6,9 h |

L'episodi del 3-4/10/2026 va fer passar el riu de 30 a 300 m³/s en 2 h, de manera que la resposta real en episodis intensos és més ràpida que la que dona Témez.

## Pla

1. ✅ Conca i xarxa de drenatge (MDE Copernicus 30 m; pendent passar a l'MDT de 5 m de l'ICGC)
2. Dades: radar OPERA, nowcast, models AROME/WRF/ICON, aforador i pluviòmetre de l'ACA
3. Capa 1: perill per pluja acumulada a la conca (1/3/6/24 h) comparada amb llindars
4. Capa 2: model de cabal (`scripts/model_cabal.py`): llindar de pluja + dipòsits no lineals, calibrat amb el 3-4/10/2026
   (amb el pluviòmetre NSE 0,87, pic 206 de 304 m³/s; amb el radar NSE 0,63, perquè OPERA infravalora els nuclis convectius).
   Pendent: fer-lo córrer cada 10 min, nowcast/AROME, propagació a la desembocadura i correcció amb l'aforador.
5. ✅ Visor web (GitHub Pages), amb el període de retorn equivalent al cabal (cabals oficials de l'ACA per tram)
6. Validació amb episodis històrics

## Funcionament

- **GitHub Actions** (`.github/workflows/recull.yml`) recull les dades cada 10 minuts i fa els càlculs.
- **Supabase** (projecte minairo, esquema `inundacions`) guarda els resultats: `observacions`, `pluja_conca`, `previsions_cabal`, `avisos`. Lectura pública; escriptura només amb la clau secreta.
- **GitHub Pages** serveix el visor, que llegeix de Supabase.

Secret del repo (Settings → Secrets and variables → Actions): `SUPABASE_SECRET_KEY`.

## Estructura

- `index.html`: visor
- `config/ridaura.json`: aforador, sensors, llindars i cabals de càlcul de l'ACA (Q10–Q500) de la conca
- `config/ridaura_model_*.json`, `data/ridaura/calibratge_*.png`, `data/ridaura/episodi_*.csv`: calibratge del model
- `scripts/radar_conca.py`: pluja del radar OPERA sobre la conca (temps real o, amb `--des/--fins`, de l'arxiu des de 2011)
- `scripts/aca_aforador.py`: recull el nivell, el cabal i la pluja de l'aforador de l'ACA
- `data/ridaura/conca.geojson`: conca a la desembocadura, conca a l'aforador i xarxa de drenatge
- `scripts/retalla_mde.py`, `scripts/delimita_conca.py`, `scripts/dibuixa_conca.py`: delimitació de la conca
- `docs/`: recerca de models i rieres litorals de risc
