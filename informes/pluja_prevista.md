# Pluja prevista a la previsió de la Ridaura: validació (24/11/2026)

Opció `previsio.py --pluja-prevista tot` (PR #37). Fonts:
- **Nowcast del radar OPERA (0-2 h):** extrapolació de l'últim compost amb el moviment estimat entre dos compostos separats 15 min.
- **AROME France HD** (Météo-France, via Open-Meteo, fins a 6 h). A cada hora passada es fa servir l'execució que hi hauria hagut disponible (≥ 3 h abans).

Validació amb `scripts/valida_pluja_prevista.py`, executada a Actions. Període 11/08 - 11/11/2026, una previsió cada hora. Pluja mesurada: pluviòmetre de l'ACA a l'aforador.

## Resum

| | Primer avís de T10 el 4/10 | Antelació (l'aforador hi arriba a les 08:35 UTC) |
|---|---|---|
| Ara (sense pluja prevista) | 07:10 | 85 min |
| Només nowcast | 06:30 | 125 min |
| Només AROME | 06:30 | 125 min |
| **Nowcast + AROME** | **04:20** | **255 min (4 h 15 min)** |
| Amb la pluja que va caure (sostre) | 02:30 | 365 min |

Per a l'avís de crescuda ordinària, l'antelació passa de 120 min (06:30) a 380-390 min (02:00-02:10).

- **Falses alarmes:** hi ha 8 hores amb una falsa alarma, repartides en 2 dies (6/10 i 8/10). Totes cauen després del 4/10, quan el model donava per esgotada l'absorció inicial (vegeu més avall).
- **La pluja d'AROME** segueix la mesurada:
  - biaix de +5 % a +10 % i correlació de 0,72 (1 h), 0,82 (3 h) i 0,83 (6 h);
  - de les hores amb ≥ 40 mm en 6 h, n'encerta 5 de 6, amb 4 falses alarmes.
- **El nowcast** és més feble:
  - sobreestima un 69 % a 1 h i té una correlació de 0,45;
  - al 4/10, el moviment estimat era 0 km/h, de manera que va fer de persistència: el valor que hi aporta és el de dir "continua plovent així".

## Detall

### A. Pluja prevista vs mesurada (281 hores amb pluja a prop; AROME disponible en 175)

| Font | Termini | Llindar | Encerts | Falses alarmes | No previstos |
|---|---|---|---|---|---|
| nowcast | 1 h | ≥ 10 mm | 5 | 11 | 3 |
| arome | 1 h | ≥ 10 mm | 4 | 4 | 1 |
| arome | 3 h | ≥ 20 mm | 5 | 6 | 1 |
| arome | 3 h | ≥ 40 mm | 3 | 2 | 1 |
| arome | 6 h | ≥ 40 mm | 5 | 4 | 1 |
| arome | 6 h | ≥ 70 mm | 3 | 0 | 2 |
| tot | 6 h | ≥ 40 mm | 5 | 2 | 1 |

Totals:

| Font | Termini | Prevista | Mesurada | Correlació |
|---|---|---|---|---|
| arome | 1 h | 203 mm | 187 mm | 0,72 |
| arome | 3 h | 596 mm | 577 mm | 0,82 |
| arome | 6 h | 1.151 mm | 1.046 mm | 0,83 |
| nowcast | 1 h | 498 mm | 295 mm | 0,45 |

### C. Falses alarmes de l'avís

S'han avaluat 110 hores, les que tenien ≥ 5 mm previstos o mesurats en 6 h. La pluja prevista fa pujar el nivell en 12:
- **4 són bones:** el 4/10 a les 02:00 i les 03:00, i el 5/10 a les 09:00 i les 11:00. En aquestes, la previsió feta amb la pluja real dona el mateix nivell.
- **8 són falses:** el 6/10 a les 16, 18, 20, 21 i 22 h, i el 8/10 a les 11, 17 i 18 h.
  - El 6/10, AROME preveia ~46 mm en 3 h i en van caure 5.
  - El nivell arriba a T10-T50 perquè el model encara compta l'episode del 4/10. La regla "inici d'episodi després de 48 h seques" no es reinicia, i l'absorció de 159 mm continua gastada.
  - Amb l'absorció calculada a partir dels 60 dies anteriors (l'altra feina en curs), aquests nivells baixarien. Aquest efecte no és de la pluja prevista, sinó del model de cabal: el 5/10, la previsió amb la pluja real ja dona T10 i l'aforador va marcar 1 m³/s.

## Limitacions

- Només hi ha **un episodi amb crescuda** (4/10), i és el mateix que es va fer servir per calibrar el model. L'antelació s'ha de tornar a mirar amb el pròxim episodi.
- **Límit d'Open-Meteo:**
  - La versió gratuïta limita les peticions. Executant trossos en paral·lel es van perdre execucions d'AROME, i per això en surten 175 hores de 281.
  - La taula de l'antelació és d'una execució sola (sense límit). En paral·lel, nowcast + AROME donava T10 a les 06:30.
  - Ara el codi torna a provar-ho quan Open-Meteo diu que hi ha massa peticions, i el workflow va d'un en un.
- **Nowcast:**
  - Fa servir Marshall-Palmer, sense correcció amb pluviòmetres.
  - El moviment per correlació de fase sovint surt 0. Una millora possible és el flux òptic de pysteps.
- **AROME:**
  - Com a retard de publicació s'han suposat 3 h, una xifra prudent.
  - En directe es fa servir l'última execució publicada.
