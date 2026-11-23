# Com s'ha fet i com replicar-ho

Guia per muntar el mateix sistema en una altra conca de Catalunya o en una altra comunitat. Explica la cadena, d'on surt
cada dada, en quin ordre s'executa i què hem après pel camí. Estat: novembre de 2026.

## 1. Què fa el sistema

Vol respondre **quins punts d'un poble s'inunden amb una pluja forta, i quan**. La cadena és:

```
radar OPERA (5 min, 1 km) ──► corregit amb pluviòmetres (XEMA + ACA)
        │
        ▼
pluja per conca ──► absorció del terreny (llindar Ia + infiltració) ──► cabal (model agregat, validat amb aforadors ACA)
        │                                                                     │
        ▼                                                                     ▼
simulació 2D sobre el relleu (MDT/MET, edificis, passos,          període de retorn (Q10…Q500 oficials de l'ACA)
clavegueram) ──► calat, hora d'arribada i temps d'assecat                     │
        │                                                                     ▼
        └──────────────► punts afectats (carrers, soterranis, plantes baixes, passos inferiors) ◄── zones inundables ACA T10-T500
```

Tres vies d'inundació: desbordament de riu, xàfec amb clavegueram saturat, i riera + temporal de mar.

Dos principis de disseny:
- **Res ajustable a mà** quan es reconstrueix un episodi: tot surt de dades reals i es valida amb aforadors, premsa i
  intervencions de Bombers.
- **Tot el que es pot precalcular, es precalcula** per a tot el territori (relleu, cobertes, sòls, passos, edificis,
  clavegueram) i es desa. En directe només s'hi aplica la pluja.

## 2. Infraestructura (cost zero)

| Peça | Ús | Notes |
|---|---|---|
| GitHub (repo públic) | codi, configuració, capes lleugeres per a la web | repo públic = minuts d'Actions il·limitats |
| GitHub Actions | tots els càlculs; `vigila.yml` cada 10 min | una execució de ~5 h 40 min amb un bucle de 10 min que es torna a llançar sola; un `cron` horari de reserva |
| GitHub Releases | dades grosses (MDT, MET, cobertes, cadastre, radar mensual) | fins a 2 GB per fitxer; no compten per a Pages |
| GitHub Pages | visor web estàtic | límit 1 GB publicat: no hi posis rasters ni dades que la web no llegeix |
| Supabase (pla gratuït) | sèries en directe: `observacions`, `pluja_conca`, `previsions_cabal`, `avisos` | RLS amb només lectura per a `anon`; escriptura amb la clau secreta (secret `SUPABASE_SECRET_KEY`) |

## 3. Dades i equivalents fora de Catalunya

| Dada | Catalunya | Resta d'Espanya / altres |
|---|---|---|
| Radar | OPERA (EUMETNET, CC BY 4.0), compost DBZH 5 min ~1 km. Directe: bucket de reserva d'ECMWF des de l'11/11/2026; arxiu a S3 | El mateix OPERA cobreix tot Europa. AEMET té radar propi (OpenData, cal clau) |
| Pluviòmetres | XEMA (Meteocat, dades obertes Socrata, des de 2009) i pluviòmetres ACA (API Sentilo) | AEMET OpenData; SAIH de cada confederació (Ebre, Xúquer, Segura, Guadalquivir…) |
| Cabal en directe | Aforadors ACA (API Sentilo pública, ~100 dies d'històric) | SAIH de la confederació; ROEA (anuari d'aforaments, MITECO) per a històric |
| Cabals de període de retorn | ACA, tramificació de cabals (Q10…Q500 per tram) | SNCZI (MITECO): mapes de caudales máximos (CAUMAX) i estudis ARPSI |
| Zones inundables | ACA (T10, T50, T100, T500) | SNCZI: zones inundables T10/T100/T500 de tot Espanya |
| Relleu | MDT 5 m de l'IGN (WCS) i MET 1-2 m de l'ICGC | MDT 5 m i 2 m de l'IGN/CNIG (PNOA LiDAR) a tot Espanya |
| Cobertes del sòl | Mapa de cobertes de l'ICGC 2024 (WMS) | SIOSE (IGN) o CORINE / ESA WorldCover |
| Sòls | Litologia IGME (permeabilitat) → grups hidrològics A-D | Igual (IGME és estatal). SoilGrids dona gairebé tot grup C i no serveix |
| Edificis | Cadastre INSPIRE (edificis, plantes sobre/sota rasant, ús) | Igual a tot Espanya excepte País Basc i Navarra (cadastres forals) |
| Vies i passos | OpenStreetMap (Geofabrik) + inventari d'obres de l'ACA | OSM + inventaris de carreteres de cada administració |
| Mar | Puertos del Estado (onatge, nivell) | Igual |
| Validació | Intervencions de Bombers (dades obertes), premsa local | Serveis d'emergència de cada comunitat, Consorcio de Compensación de Seguros |

## 4. Ordre de construcció

Els workflows són a `.github/workflows/` i els scripts a `scripts/`. Cada script té l'ús explicat a la capçalera.

**A. Base per a tot el territori (una sola vegada, després només si canvien les fonts)**
1. Conques de treball: `mapa_conques.yml` (tall de 100 km² per partir l'MDT; tall de 25 km² per vigilar).
2. Relleu: `mdt_cat.yml` (MDT 5 m per conca → Release `mdt5-cat`), `met2_cat.py` (MET 2 m i 1 m per a la 2D),
   `relleu.yml` (rieretes i zones baixes).
3. Cobertes i sòls: `cobertes_llegenda.yml` → `cobertes_cat.yml`; `sols_cat.yml` (grups hidrològics).
4. Pluja de disseny: `xema_maxims.yml` (màxims anuals XEMA), `clavegueram_cat.py idf`.
5. Passos i obres: `passos_cat.yml` (torrentets sota terraplens, amb cabal de disseny 5.2-IC), `obres_cat.yml` (secció i
   capacitat de cada obra), `arrossegalls_cat.yml`.
6. Ciutat: `cadastre_cat.yml`, `planta_baixa.yml`, `diposits_cat.yml`, `clavegueram_cat.yml`.

**B. Conca nova**
1. Fitxer `config/<conca>.json` (copiar-ne un): `conca`, `nom`, **`centre` [lat, lon]**, aforador i sensors, `geometria`,
   `cabals_T` (Q del tram de l'aforador i de la desembocadura), `codis_zones_aca`, `municipis_cadastre`.
2. `conca_nova.yml` (radar, estacions, relleu i graella) i `prepara.yml` (zones inundables creuades amb OSM i SNCZI).
3. Històric: `pluja_xema.yml`, `radar_episodis.yml` / `radar_arxiu.yml`, `aca_arxiu.yml` (l'ACA només guarda ~100 dies:
   cal arxivar-ho cada mes).
4. Calibratge: `model.yml` (`model_cabal.py`, un episodi amb aforador). Comprovar-lo amb altres episodis (vegeu §6).
5. Simulació 2D d'episodis: `sim2d.yml` (`simula_2d.py`), validada amb l'hidrograma de l'aforador i amb premsa/Bombers.

**C. Directe**
- Afegir la conca al bucle de `vigila.yml`: `aca_aforador.py`, `radar_conca.py`, `radar_corregit.py`, `previsio.py`.
- La web (`previsio.html`) llegeix Supabase amb la clau publicable.

## 5. Mètodes principals (resum)

- **Radar → pluja**: Z-R ajustada amb estacions en directe (`calibra_zr.py`), Marshall-Palmer a l'arxiu; tall a 125 mm/h.
  Correcció per biaix mitjà horari amb els pluviòmetres a 40 km (≥ 3, factor 0,25-4) i, a l'històric mensual, també
  residus interpolats (`correccio_cat.py`, validació deixant cada estació fora: biaix -26 % → -1 % a l'octubre de 2026).
- **Cabal en directe**: llindar d'absorció Ia + coeficient C + dos dipòsits no lineals Q = (S/K)^(1/m), pas de 5 min;
  s'assimila el cabal mesurat ara (la diferència s'esvaeix en 2 h). Dos escenaris: no plou més / plou 1 h més igual.
  Desembocadura: mateix període de retorn amb 45 min de retard. Rieres sense aforador: mètode racional.
- **Període de retorn**: interpolació en log T entre els Q oficials del tram.
- **Simulació 2D**: esquema inercial local (Bates et al. 2010, com LISFLOOD-FP), Manning per coberta, absorció per cel·la
  repartida amb el número de corba, edificis com a porositat o dipòsit, clavegueram per subconca amb capacitat C·I(T10),
  passos sota vies com a conductes.
- **Cabal de disseny de petites conques**: 5.2-IC (Témez, Kt, Ka, I1/Id) amb Gumbel regional dels màxims XEMA.

## 6. Lliçons apreses (el que ens hauria estalviat temps)

1. **L'MDT del WCS de l'IGN ve en metres enters**: crea graons i falsos clots. Fer servir el MET 2 m/1 m on importi.
2. **Els terraplens sense passos fan de dic** a l'MDT: 1,4-1,6 hm³ quedaven en clots a la Ridaura. Cal detectar els passos.
3. **No comptar les pèrdues dues vegades** (llindar d'un calibratge + infiltració d'un altre).
4. **Mesurar el cabal simulat amb una caixa al voltant de l'aforador**, no amb un transsecte (infra-mesura un 20-25 %).
5. **El radar infravalora els nuclis convectius** (-26 %): sempre corregir-lo amb pluviòmetres.
6. **Les fonts canvien d'adreça sense avís** (OPERA va passar a un bucket de reserva l'11/11/2026 i es van perdre 12 h).
   Cal vigilar la frescor de les dades, no només que el workflow acabi bé.
7. **L'absorció inicial depèn de l'estat del sòl**: 150-160 mm a la Ridaura després de l'estiu, 65 mm al Congost.
   Un sol episodi no n'hi ha prou per calibrar-la.
8. **SoilGrids dona gairebé tot grup C** a Catalunya: la litologia de l'IGME discrimina molt millor.
9. **Les dades grosses a Releases**, no al repo: Pages té un límit d'1 GB.
10. **Opcions noves desactivades per defecte**: cada millora del model és una opció que es prova amb episodis abans
    de canviar la simulació per defecte (que ha de donar el mateix resultat bit a bit).

## 7. Validació mínima abans de fer-ho servir

- Cabal: NSE i error de pic contra l'aforador en **episodis que no s'han fet servir per calibrar**.
- Previsió: reproduir-la amb `previsio.py --ara <hora passada> --sense-desar` en diversos episodis i mesurar amb quina
  antelació avisa. A la Ridaura el 4/10/2026 (episodi de calibratge) l'escenari «plou 1 h més» avisava de T10 ~1 h 25 min
  abans que l'aforador hi arribés; l'escenari «no plou més», pràcticament en el moment.
- Punts afectats: comparar amb premsa i intervencions de Bombers de cada episodi (encerts, falses alarmes, omissions).
