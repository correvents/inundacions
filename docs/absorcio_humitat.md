# Absorció inicial segons la pluja dels 60 dies anteriors (23/11/2026)

## Què canvia
Abans, la previsió de la Ridaura feia servir sempre Ia = 159 mm (calibrat el 4/10, després d'un estiu sec).
Ara el terreny és un **dipòsit d'humitat del sòl** (`scripts/humitat.py`):

- la pluja l'omple i l'evapotranspiració l'asseca (ET0 FAO-56 d'ERA5, en climatologia de cada dia de l'any, precalculada),
  cada cop menys com més sec és;
- **Ia = el dèficit del dipòsit en començar l'episodi**, calculat amb la pluja dels 60 dies anteriors
  (pluviòmetre ACA en directe; XEMA per als episodis antics). L'estat de fa 60 dies és el típic d'aquell dia de l'any
  (dipòsit fet córrer amb la sèrie XEMA des del 2009);
- després d'omplir Ia, el terreny encara perd una quantitat contínua **phi**.

Res no s'ajusta a mà. Els paràmetres surten del cabal observat:

| Conca | Smax | phi | Com |
|---|---|---|---|
| Ridaura (Santa Cristina) | 195 mm | 7 mm/h | Hidrograma del 4/10/2026 (NSE) i pic del 17-21/1/2026 (28,6 m³/s), junts |
| Congost (La Garriga) | 130 mm | — | Ia = 65 mm el 3/10/2026 (el calibratge que ja hi havia) |

Ara mateix (23/11) l'absorció és d'uns **24 mm a la Ridaura** (abans 159) i uns 9 mm al Congost.

## Validació (cabal observat, ACA)

**Ridaura**

| Prova | Ia fixa 159 mm | Dipòsit |
|---|---|---|
| 4/10/2026, calibratge (NSE) | 0,87 (pic 206 de 304) | 0,86 (pic 195), Ia 153 mm |
| 17-21/1/2026, independent: pic obs. 28,6 m³/s | 0 m³/s (no l'hauria vist) | 27 m³/s, Ia 22 mm |
| Sèrie contínua 11/8-23/11/2026 (NSE) | 0,69 | **0,83** |
| 17/10 (20 mm) i 12/11 (28 mm): sense crescuda | cap crescuda | cap crescuda (Ia 20 i 36 mm) |
| 89 episodis XEMA 2009-2026 per sobre de QMCO (37 m³/s), fora del 4/10 | 2 (nov. 2018, Glòria) | 3 (nov. 2018 a 95 m³/s, Glòria 40, feb. 2010 39) |

La premsa diu que el novembre del 2018 hi va haver crescuda del Ridaura i la riera de Malvet es va desbordar, i que
amb el Glòria el Ridaura va fer danys al marge dret a Santa Cristina. Per tant, aquests casos no són falses alarmes.
Del febrer del 2010 no se'n troba res.

**Congost** (model ràpid ajustat al 3-5/10)

| Prova | Ia fixa 65 mm | Dipòsit |
|---|---|---|
| Sèrie contínua jul-nov 2026 (NSE) | 0,80 | 0,81 |
| 17/10 (53 mm), obs. 1,6 m³/s | 0,1 | 2,2 |
| 16/9 (31 mm), obs. 2,8 m³/s | 0 | 0 (cap dels dos) |

## Límits
- Només hi ha dos episodis amb cabal a la Ridaura i un al Congost. L'arxiu de l'ACA comença a l'agost del 2026.
  Cada episodi nou que es mesuri reforça la calibració.
- La pluja de gener i dels episodis antics és la de les estacions XEMA (no de la conca), sense cap factor.
- En directe, si falten dades del pluviòmetre de l'ACA, aquelles hores compten com a seques, i per tant l'absorció surt més alta.
- Els episodis de la web (`premsa.py`, `index_episodis.py`) i la simulació 2D encara fan servir l'índex antic.
  Es poden passar al dipòsit en un pas següent.
- Cabrera no té aforador i, per tant, no té Smax.
