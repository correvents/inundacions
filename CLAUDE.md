# CLAUDE.md — Inundacions

Sistema propi, no oficial i en desenvolupament per saber **quins punts dels pobles s'inunden amb pluges fortes i quan**.
Les zones inundables de l'ACA (T10-T500) ja existeixen; aquí es busca quina pluja les activa i es comprova amb episodis reals.
Guia llarga: `docs/replicar.md`. Treballa i escriu en **català**, de manera concisa.

## Regles (no negociables)

1. **Res ajustable a mà en reconstruir episodis.** Cada paràmetre surt d'una dada o d'una norma (5.2-IC, TR-55, Chow,
   ARR 2019) i diu d'on surt. Si cal un valor, es calibra amb un episodi i es valida amb un altre.
2. **Paràmetres calibrats per conca** (`config/<conca>.json`, `config/<conca>_humitat.json`), mai un valor únic per a tot Catalunya.
3. **Precalcular tot el que no depèn de la pluja** (terreny, cobertes, sòls, edificis, passos, clavegueram) per a tot
   Catalunya i desar-ho a Releases. En directe només s'hi aplica la pluja.
4. **Càlculs al servidor** (GitHub Actions). La web només llegeix fitxers estàtics i Supabase; mai calcula.
5. **Cada millora és una opció desactivada per defecte.** La simulació per defecte ha de quedar idèntica bit a bit fins
   que l'opció demostri amb episodis reals que millora. Desa els resultats de cada prova per comparar-los.
6. **Ordre de treball**: primer afinar el terreny, després la pluja i els episodis.
7. **Repo públic**: cap dada personal, cap credencial, cap dada d'intervencions que no sigui oberta (ni text lliure,
   adreces o noms). Només resultats agregats. No anomenis fonts internes o no publicades als docs, commits ni web.

## La cadena

```
radar OPERA + pluviòmetres XEMA/ACA → pluja corregida (biaix horari + IDW) → pluja per conca (5 min)
→ absorció inicial segons la pluja dels 60 dies anteriors (dipòsit d'humitat) + pèrdua contínua
→ cabal (model agregat, validat amb aforadors ACA) → període de retorn (Q10-Q500 del tram)
→ zones ACA + simulació 2D sobre MDT/MET amb edificis, soterranis, passos i clavegueram
→ punts afectats: hora d'arribada, calat, velocitat, temps d'assecat
```

Tres vies: desbordament de riu, xàfec amb clavegueram saturat, rieres litorals amb temporal de mar.

## On és cada cosa

| Què | On |
|---|---|
| Scripts | `scripts/` (un per pas; el docstring diu l'ús). Accés a Supabase: `scripts/sb.py` |
| Configuració per conca | `config/ridaura.json` (pilot, Platja d'Aro), `granollers.json` (Congost), `cabrera.json` |
| Resultats lleugers (web) | `data/<conca>/` (episodis, sim2d, model, zones) |
| Workflows | `.github/workflows/` (`vigila.yml` = directe cada 10 min; la resta, `workflow_dispatch`) |
| MDT 5 m IGN per conca | Release `mdt5-cat` (metres enters: compte amb clots falsos) |
| MET ICGC | Releases `met21-cat` (1 m, 2021-23) i `met2-cat` (2 m, 2008-11), només zones urbanes |
| Cobertes ICGC 2024 / sòls IGME | Releases `cobertes-cat` / `sols-cat` |
| Edificis i soterranis Cadastre | Releases `cadastre-cat`, `diposits-cat`, `planta-baixa-cat` |
| Radar OPERA retallat mensual | Releases `opera-cat-AAAA` (uint8, dBZ = v/2 − 32) |
| Proves de simulació 2D | Release `sim2d-proves-*` (un `.tar.gz` per prova) |
| Sèries vives | Supabase, esquema `inundacions`: `observacions` (aforadors), `pluja_conca` (camp `font`: radar, `radar_corregit`...), `avisos`, `punts_baixos` |
| Aforadors | API Sentilo de l'ACA (~100 dies d'històric; arxiu mensual amb `aca_arxiu.py`) |

Secrets: `SUPABASE_SECRET_KEY` (escriptura). La web fa servir només la clau publicable (lectura amb RLS).

## Ordres habituals

```bash
# Previsió (directe o reproduint una hora passada sense escriure a Supabase)
python3 scripts/previsio.py --config config/ridaura.json [--ara 2026-10-04T07:00 --sense-desar]

# Absorció inicial (dipòsit d'humitat): calibrar Smax per conca
python3 scripts/humitat.py calibra --config config/ridaura.json

# Simulació 2D d'un episodi (normalment via el workflow sim2d.yml)
python3 scripts/simula_2d.py --config config/ridaura.json --mdt <mdt.tif> --radar <radar.json> --id 2026-10-04 --pas 20
#   opcions (totes desactivades per defecte): --edificis --edificis-diposit --passos oberts|arrossegalls
#   --escenari-clavegueram funciona|antic|collapsat --clavegueram-desti riu --grup-sol --planta-baixa --velocitat

# Validació per carrers amb la premsa
python3 scripts/valida_sim2d.py --config config/ridaura.json
```

## Estalviar context

- Llegeix el docstring del script (`head -30`) abans d'obrir-lo sencer; molts fan 500+ línies.
- No bolquis sortides llargues: `| tail -20`, o llegeix el `resum.md` que escriu cada simulació 2D
  (`python3 scripts/resum.py data/<conca>/sim2d/<id>`).
- Mai llegeixis fitxers binaris (`*.u2`, `*.tif`, `*.h5`) ni GeoJSON grans sencers: resumeix-los amb Python.
- Logs de GitHub Actions: busca-hi només l'error o l'última línia de cada pas.

## Skills del projecte

`.claude/skills/`: `reconstruir-episodi`, `afegir-conca`, `validar-intervencions`.
