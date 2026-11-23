# Com replicar el sistema en una altra comunitat

Guia de com s'ha construït aquest sistema (Catalunya, 2026) i què caldria canviar per fer-lo en una altra comunitat autònoma. És un sistema propi i en desenvolupament, no oficial.

**Objectiu**: saber quins punts concrets dels pobles s'inunden amb pluges fortes i quan, a partir de la pluja real. Les zones inundables oficials (T10-T500) ja existeixen; el que falta és saber quina pluja les activa, i comprovar-ho amb episodis passats.

**Principis que s'han mantingut i que convé copiar**

1. **Res ajustable a mà en reconstruir episodis.** Tot paràmetre surt d'una dada pública o d'una norma (5.2-IC, TR-55, Chow, ARR 2019...), i cada valor diu d'on surt.
2. **Precalcular tot el que no depèn de la pluja** (terreny, cobertes, passos, edificis, clavegueram) per a tota la comunitat, i desar-ho. En directe només s'hi aplica la pluja.
3. **Cada millora del model és una opció desactivada per defecte.** La simulació per defecte no canvia bit a bit fins que una opció demostra amb episodis reals que millora.
4. **Validar amb dades independents**: aforadors per al cabal; premsa i intervencions de bombers per als carrers.
5. **Calibrar per conca**, no un valor per a tota la comunitat.

---

## 1. La cadena

```
radar OPERA ─┐
             ├─> pluja corregida ─> pluja per conca ─> infiltració ─> cabal ─> període de retorn
estacions ───┘     (biaix + IDW)       (cada 5 min)    (60 dies abans)  (aforador)   (Q10-Q500)
                                                                                      │
              terreny precalculat (MDT/MET, cobertes, edificis, passos, clavegueram)   │
                                                                                      v
                                  simulació 2D (pluja sobre la malla) + zones inundables oficials
                                                                                      │
                                                                                      v
                               punts afectats: hora d'arribada, calat, velocitat, temps d'assecat
```

Tres vies d'inundació: **desbordament de riu**, **xàfec amb clavegueram saturat** i **rieres litorals amb temporal de mar**.

### 1.1 Pluja

| Pas | Com s'ha fet | Script |
|---|---|---|
| Radar | Compostos OPERA (EUMETNET Open Radar Data, CC BY 4.0), DBZH cada 5 min a ~1 km des de 2024 (15 min i 2 km abans). Lectura per HTTP Range: Catalunya cau en un sol bloc del HDF5 (~200 kB de 1,9 MB). | `radar_conca.py`, `radar_graella.py`, `pluja_catalunya.py` |
| Z-R | Marshall-Palmer per defecte; relació ajustada per tipus de pluja (estratiforme/convectiva) amb estacions locals. | `calibra_zr.py` |
| Correcció | Per hora: factor de biaix F = Σ pluviòmetres / Σ radar (límits 0,25-4, només amb ≥3 pluviòmetres i ≥1 mm), després residus interpolats per IDW (radi 30 km) i sostre 1,5 × màxim del pluviòmetre. Octubre 2026: biaix -26% → -1%, error diari 16 → 11 mm. | `correccio_cat.py`, `radar_corregit.py`, `pluviometres_cat.py` |
| Arxiu | Retall mensual de la comunitat en HDF5 uint8 (dBZ = v/2 − 32), 2012-2026 en 3,3 GB. | Releases `opera-cat-AAAA` |
| Previsió | AROME France HD (1,5 km, 15 min) via Open-Meteo. | — |

### 1.2 Infiltració i cabal

- **Absorció inicial** segons la pluja dels **60 dies anteriors** (índex de pluja antecedent amb l'estació més propera), després una pèrdua contínua (~12 mm/h a la Ridaura; 65 mm i 10 mm/h al Congost).
- **Model agregat** a pas de 5 min (`model_cabal.py`): llindar Ia + fracció d'escolament + dipòsits no lineals. Ridaura 4/10/2026: amb pluviòmetre NSE 0,87; amb radar sense corregir 0,63.
- **Validació amb aforadors** (API Sentilo de l'ACA, sense clau, només ~100 dies → cal arxivar-ho cada mes: `aca_arxiu.py`).
- **Període de retorn**: el cabal es compara amb els Q10/Q50/Q100/Q500 oficials del tram (`config/<conca>.json`).
- **Previsió en directe** (`previsio.py`, cada 10 min): busca l'inici de l'episodi (última finestra de 48 h amb <1 mm) i propaga a la desembocadura.

### 1.3 Terreny precalculat (tota la comunitat)

| Capa | Font a Catalunya | Mètode | Script / Release |
|---|---|---|---|
| Conques | MDT propi | Arbre de confluències; talls de 100 km² (per partir l'MDT, 365 conques) i 25 km² (per vigilar, 711) | `mapa_conques.py` |
| MDT 5 m | IGN WCS `Elevacion25830_5` (LiDAR PNOA) | Per conca + 500 m, int32 cm | `mdt_cat.py` / `mdt5-cat` |
| MET 1-2 m | ICGC (COG a datacloud) | Només zones urbanes, lectura per finestres | `met2_cat.py` / `met21-cat`, `met2-cat` |
| Drenatge i zones baixes | MDT | pysheds, D8, rieretes ≥ llindar | `relleu.py` |
| Cobertes del sòl | ICGC 2024 (41 classes) via WMS ràster | Color → classe amb GetFeatureInfo; n de Manning (Chow 1959) i CN (TR-55) per classe | `cobertes.py` / `cobertes-cat` |
| Grup hidrològic de sòl | IGME Permeabilitats 1:200.000 | Permeabilitat → A-D | `sols.py` / `sols-cat` |
| Edificis i soterranis | Cadastre INSPIRE Buildings | Plantes sobre/sota rasant; obstacles i dipòsits | `cadastre_cat.py`, `diposit.py` / `cadastre-cat`, `diposits-cat` |
| Ús de planta baixa | Cadastre INSPIRE + OVC Consulta_DNPRC | Habitatge, local, garatge... | `planta_baixa.py` / `planta-baixa-cat` |
| Passos sota vies | MDT + vies OSM (Geofabrik) | Depressió que toca una via = pas probable; escenaris obert/mig tapat/tapat | `passos_cat.py` |
| Secció dels passos | ACA PSE (WFS + fitxes), catàleg de ponts, OSM | Orifici Cd 0,6; si no hi ha mida, cabal de disseny 5.2-IC | `obres_cat.py`, `passos_2d.py` |
| Cabal de disseny | Màxims diaris XEMA | Gumbel + 5.2-IC (P0 per cobertes) | `xema_maxims.py` |
| Arrossegalls | Cobertes + incendis (Agricultura, EFFIS) | ARR 2019 llibre 6 cap. 6 | `arrossegalls.py` |
| Clavegueram | IDF de les estacions + cobertes | Mètode racional C·I(T10, tc) per subconca urbana; escenaris funciona/antic/col·lapsat | `clavegueram_cat.py` |
| Estrenyiments de llera | MET 2 m + cobertes | Seccions cada 20 m; capacitat Manning vs Q2-Q500 | `estrenyiments.py` |
| Zones inundables | ACA (T10-T500), SNCZI | Creuades amb elements OSM | `prepara_zones.py`, `creua_zones.py` |

### 1.4 Simulació 2D

`simula_2d.py`: pluja sobre la malla (rain-on-grid), esquema inercial local de Bates et al. 2010 (el de LISFLOOD-FP), Manning per cel·la, pas 5-20 m (2-5 m amb MET). Opcions: `--edificis`, `--edificis-diposit`, `--passos`, `--escenari-clavegueram`, `--grup-sol`, `--planta-baixa`, `--velocitat` (perill AIDR H1-H6). Corre a GitHub Actions (`sim2d.yml`), ~20 min a 20 m per a la Ridaura, ~3 h a 10 m.

Sortides: calat màxim, hora d'arribada, temps d'assecat, velocitat, soterranis plens. Validació per carrer amb la premsa (`valida_sim2d.py`).

### 1.5 Regles per a la ciutat

A més del model físic, una regla simple per poble ajustada amb el catàleg d'episodis (estacions + premsa). Platja d'Aro: risc si ≥15 mm en 30 min, o ≥70 mm d'episodi, o llevant amb onatge ≥2,8 m i ≥40 mm. Encerta 23 de 24 episodis. L'estat del mar surt d'Open-Meteo Marine (`mar.py`).

---

## 2. Fonts de dades i equivalents fora de Catalunya

| Necessitat | A Catalunya | Fora de Catalunya |
|---|---|---|
| Radar | OPERA (europeu) | **Igual**: OPERA cobreix tot l'Estat. AEMET té radars propis (OpenData, amb clau). |
| Pluviòmetres oficials | XEMA (Meteocat, dades obertes Socrata, des de 2009) + ACA | **AEMET OpenData** (clau gratuïta) + serveis autonòmics: Euskalmet, MeteoGalicia, AVAMET (xarxa ciutadana, CV), Meteoclimatic, SAIH de cada confederació |
| Aforadors en temps real | ACA (API Sentilo, sense clau) | **SAIH** de cada confederació (Ebre, Xúquer, Segura, Guadalquivir, Tajo, Duero, Miño-Sil, Cantàbric, Guadiana) o agència autonòmica (URA, Augas de Galicia, Agencia Andaluza del Agua). Formats i accés molt diferents: és el pas que més costa. |
| Cabals de càlcul Q10-Q500 | ACA (estudis per tram) | **SNCZI** (MITECO): zones inundables i cabals dels ARPSI de tot l'Estat |
| Zones inundables | ACA T10-T500 | **SNCZI** T10/T100/T500 (fluvial i costanera) |
| MDT 5 m | IGN | **Igual** (IGN cobreix tot l'Estat). MDT02 del CNIG on n'hi hagi (no cobreix Catalunya). |
| MDT 1-2 m | ICGC | CNIG MDT02/MDT01, o l'institut cartogràfic autonòmic (IGCV a València, IECA a Andalusia, GeoEuskadi...) |
| Cobertes del sòl | ICGC 2024 | **SIOSE** (IGN) o CORINE / ESA WorldCover (10 m, tot el món) |
| Sòls | IGME | **Igual** (IGME estatal) |
| Edificis | Cadastre INSPIRE | **Igual**, excepte País Basc i Navarra (cadastres forals propis) |
| Vies | OSM (Geofabrik) | **Igual** |
| Ponts i obres de drenatge | ACA PSE, catàleg de ponts Generalitat | Inventario de Obstáculos Transversales (MITECO, només conques intercomunitàries), inventaris de carreteres autonòmics i de diputacions |
| Intervencions d'emergència | Bombers Generalitat (dades obertes) | 112 / serveis autonòmics; sovint cal demanar-ho |
| Premsa | Hemeroteca local | Igual |
| Estat del mar | Open-Meteo Marine, Puertos del Estado | **Igual** |
| Previsió | AROME HD (Open-Meteo) | **Igual** a la península; HARMONIE-AROME d'AEMET |

**Normes que valen igual a tot l'Estat**: 5.2-IC (drenatge de carreteres), mètode racional, TR-55, Manning/Chow.

---

## 3. Infraestructura

Tot gratuït o gairebé:

| Peça | Ús |
|---|---|
| **GitHub** (repo públic) | Codi i capes web lleugeres |
| **GitHub Actions** | Tots els càlculs. Un workflow `vigila.yml` fa un bucle de ~5 h 40 min cada 10 min (aforador → radar → correcció → previsió) i es rellança sol; cron horari de reserva. La resta, a mà (`workflow_dispatch`) per conca o per mes. |
| **GitHub Releases** | Fitxers grans: MDT, MET, cobertes, Cadastre, arxiu OPERA (un asset per conca o mes). Límit: 1.000 assets per Release i 2 GB per fitxer. |
| **GitHub Pages** | Web estàtica (`index.html`, `previsio.html`, `conques.html`...). La web només baixa resultats; no calcula. |
| **Supabase** (Postgres + PostgREST) | Sèries vives: observacions, pluja per conca, previsions, avisos. Lectura pública amb RLS; escriptura amb la clau secreta. Bucket privat per a l'arxiu de l'ACA. |
| Secrets | `SUPABASE_SECRET_KEY` (i `AEMET_API_KEY` quan n'hi hagi). |

Configuració per conca: `config/<conca>.json` (caixa, aforador, punts, cabals Q10-Q500, estació per a la pluja antecedent). Conca nova: workflow `conca_nova.yml`.

---

## 4. Passos en ordre per a una comunitat nova

1. **Repo i infra**: fork, Supabase amb l'esquema `inundacions`, secrets, Pages activat.
2. **Conques**: MDT de 25 m de l'IGN → `mapa_conques.py` (talls de 100 i 25 km²).
3. **MDT 5 m per conca** a una Release (`mdt_cat.py`, canviant la caixa i l'EPSG si cal: Canàries 32628, la resta 25830).
4. **Pluja**: arxiu OPERA retallat a la comunitat; llista de pluviòmetres (AEMET + autonòmics); correcció de radar mensual. Comprovar el biaix amb un episodi conegut.
5. **Aforadors**: trobar l'API del SAIH, arxivar-ne les dades cada mes des del primer dia (es perden).
6. **Conca pilot**: triar-ne una amb aforador, episodis recents i premsa. Calibrar absorció i pèrdua amb un episodi, validar amb un altre.
7. **Terreny precalculat**: cobertes, sòls, Cadastre, passos, clavegueram (en aquest ordre; cada capa depèn de les anteriors).
8. **Simulació 2D** de la conca pilot amb un episodi real; comparar amb aforador i premsa. Activar opcions d'una en una.
9. **Regles de ciutat** i catàleg d'episodis amb premsa i intervencions d'emergència.
10. **Directe**: `vigila.yml` per a la conca pilot; després ampliar.

**Ordre que ha funcionat**: primer afinar el terreny, després la pluja i els episodis.

---

## 5. Lliçons apreses

**Dades**
- **El WCS de l'IGN torna l'MDT en metres enters.** Terraplens i clotades de menys d'1 m no s'hi veuen i hi ha aigua que queda atrapada en clots falsos. Per a zones urbanes cal un MDT en decimals (MET 1-2 m).
- **El terreny nu no té obres de drenatge**: cada terraplè fa de dic. Cal detectar els passos i obrir-los; i a OSM gairebé cap té mida (183 a tot Catalunya).
- **Les APIs d'aforadors guarden poc** (l'ACA ~100 dies). Arxivar des del primer dia.
- **El radar s'equivoca en tots dos sentits** (gener 2026: 29 mm de radar per 162 a l'estació). Sense correcció amb pluviòmetres no serveix per a cabal.
- **OPERA pot canviar d'adreça sense avís** (11/11/2026 el directe va passar a un servidor de reserva). Cal un avís quan la pluja deixa d'entrar.
- **SoilGrids dona sòl C a gairebé tot arreu**: no serveix. La litologia de l'IGME sí, però és permeabilitat de roca.
- **Mapes de sòls i descàrregues oficials** sovint només són accessibles des de certes xarxes: alguns servidors tallen GitHub Actions per IP o per User-Agent (no s'ha de falsejar), d'altres no responen des de fora d'Espanya. Cal tenir un ordinador propi per a aquests casos i desar-ne una instantània al repo.
- **El Cadastre limita molt les consultes per parcel·la** (~0,1/s des de GitHub): fer-ho de nit i de forma incremental, prioritzant soterranis.

**Model**
- **No comptar dues vegades les pèrdues**: absorció i pèrdua contínua s'han de calibrar juntes.
- **El clavegueram no fa desaparèixer l'aigua**: cal portar-la a la riera o al mar.
- **Manning de llera real ~0,08-0,09** (vegetació), no 0,035: amb el valor de manual la Ridaura no desbordava ni a T500.
- **La forma de mesurar compta**: un transsecte infra-mesurava el cabal un 20-25% respecte d'una caixa.
- **L'absorció inicial depèn de l'estat del sòl**: 150-160 mm a la Ridaura després de l'estiu, 65 mm al Congost. Un llindar fix calibrat amb un sol episodi avisa tard quan el sòl és moll: cal lligar-lo a la pluja dels dies anteriors i calibrar-lo amb diversos episodis.
- **Validar la previsió com si fos en directe**: `previsio.py --ara <hora passada> --sense-desar` i mesurar amb quina antelació avisa. A la Ridaura el 4/10/2026 l'escenari «plou 1 h més» avisava de T10 ~1 h 25 min abans; el «no plou més», gairebé en el moment.
- **Les opcions noves s'han de poder comparar**: desar els resultats de cada prova (aquí, un `.tar.gz` per prova a una Release).

**Infraestructura**
- GitHub Actions basta per a tot, amb 6 h per tasca; les tasques llargues es parteixen per conca o per mes.
- Un workflow nou en una branca no es pot llançar a mà fins que el fitxer és a `main`: afegir-hi un `push` a la branca per provar-lo.
- Releases: 1.000 assets màxim; separar en Releases per tema.
- La web llegeix fitxers estàtics i Supabase; mai calcula.
- GitHub Pages publica com a màxim 1 GB: les dades que la web no llegeix (rasters, clavegueram...) han d'anar a Releases.

---

## 6. Referències

- Bates, Horritt i Fewtrell (2010), *J. Hydrol.* 387: esquema inercial local.
- Norma 5.2-IC Drenaje superficial (2016).
- USDA-NRCS TR-55 i National Engineering Handbook part 630.
- Chow (1959), *Open-channel hydraulics*.
- Australian Rainfall and Runoff 2019, llibre 6 cap. 6 (bloqueig per arrossegalls).
- ASCE/WEF (1992), coeficients del mètode racional.
- Beneyto et al. (2026), *J. Hydrol. Reg. Stud.* 67:103793 (DANA València, edificis i plantes baixes).
- Documents d'aquest repo: `docs/models-previsio-inundacions.md`, `docs/rieres-litorals-risc.md`, `docs/estrenyiments.md`.
