# Estrenyiments de llera que desborden primer (opció)

Generat per `scripts/estrenyiments.py informe` amb les dades de `data/catalunya/estrenyiments/`. El mètode és a la capçalera de l'script: seccions cada 20 m amb el MET 2 m de l'ICGC (LiDAR 2021-2023), nivell de desbordament quan l'aigua arriba a conreus, carrers o edificis (cobertes ICGC), Manning amb la n de les cobertes, cabals T2-T500 de la 5.2-IC amb la XEMA o els oficials de l'ACA on n'hi ha, i capacitat real dels ponts i passos coneguts. Res ajustat a mà.

## Ridaura (Santa Cristina d'Aro, Castell d'Aro) (c125)

4408 seccions; 451 amb la llera desbordant abans de T10. Pluja: XEMA de la conca pilot ridaura (51 anys-estació), P0 = 50 mm (episodi més moll del model). 44 estrenyiments, 10 dins o a tocar de nucli.

| T desborda | Què | Riu (ACA) | Capacitat m³/s | Q10 | Q100 | Conca km² | Lloc |
|---|---|---|---|---|---|---|---|
| ≤2 | llera, 180 m, surt pel — |  | 0.5 | 5.0 | 12.0 | 3.12 | [41.80638, 3.01863](https://www.openstreetmap.org/?mlat=41.80638&mlon=3.01863#map=18/41.80638/3.01863) |
| <10 | llera, 60 m, surt pel — | Torrent de Canyet | 0.9 | 13.4 | 33.4 | 1.83 | [41.81577, 3.0184](https://www.openstreetmap.org/?mlat=41.81577&mlon=3.0184#map=18/41.81577/3.0184) |
| <10 | llera, 80 m, surt pel dret | Torrent de la Coma | 15.2 | 29.0 | 63.4 | 2.05 | [41.81569, 3.03168](https://www.openstreetmap.org/?mlat=41.81569&mlon=3.03168#map=18/41.81569/3.03168) |
| 10 | pont ACA: Obra de drenatge (Φ≥2m) GI-8521 sobre Torrent de C | Torrent de Canyet | 14.0 | 14.0 | 34.9 | 1.94 | [41.81284, 3.02101](https://www.openstreetmap.org/?mlat=41.81284&mlon=3.02101#map=18/41.81284/3.02101) |
| 10 | pas: GI-662 | Torrent de la Coma | 29.4 | 29.4 | 64.1 | 2.09 | [41.81464, 3.0321](https://www.openstreetmap.org/?mlat=41.81464&mlon=3.0321#map=18/41.81464/3.0321) |
| 45 | pont ACA: Obra de drenatge (Φ≥2m) C-31 sobre Curs innominat | Curs innominat | 11.8 | 6.9 | 15.7 | 0.53 | [41.8163, 3.0431](https://www.openstreetmap.org/?mlat=41.8163&mlon=3.0431#map=18/41.8163/3.0431) |
| 53 | pont ACA: Obra de drenatge (Φ≥2m) C-250 sobre Riera de Bell- | Riera de Bell-lloc | 45.8 | 18.0 | 65.5 | 3.52 | [41.82364, 2.97628](https://www.openstreetmap.org/?mlat=41.82364&mlon=2.97628#map=18/41.82364/2.97628) |
| 53 | pas: motorway_link | Riera de Bell-lloc | 45.8 | 18.0 | 65.5 | 3.52 | [41.82345, 2.97621](https://www.openstreetmap.org/?mlat=41.82345&mlon=2.97621#map=18/41.82345/2.97621) |
| 96 | pont ACA: Pont de trànsit rodat GIV-6611 sobre el Ridaura | el Ridaura | 193.8 | 51.6 | 199.0 | 36.14 | [41.81805, 2.97955](https://www.openstreetmap.org/?mlat=41.81805&mlon=2.97955#map=18/41.81805/2.97955) |
| 134 | pont ACA: Pont de trànsit rodat Carrer de la Teulera sobre e | el Ridaura | 322.9 | 73.1 | 292.6 | 54.24 | [41.80993, 3.0014](https://www.openstreetmap.org/?mlat=41.80993&mlon=3.0014#map=18/41.80993/3.0014) |

Ponts de l'inventari ACA que no passen el T500 (7; capacitat de la fitxa oficial):

| T desborda | Pont | Capacitat m³/s | Q10 / Q100 / Q500 |
|---|---|---|---|
| 10 | Obra de drenatge (Φ≥2m) GI-8521 sobre Torrent de Canyet | 14.0 | 14.0 / 34.9 / 53.0 |
| 10 | Cobertura de llera o soterrament GIV-6621 sobre Torrent de la Coma | 29.4 | 29.4 / 64.2 / 92.8 |
| 45 | Obra de drenatge (Φ≥2m) C-31 sobre Curs innominat | 11.78 | 6.9 / 15.67 / 23.68 |
| 53 | Obra de drenatge (Φ≥2m) C-250 sobre Riera de Bell-lloc | 45.76 | 18.0 / 65.5 / 112.8 |
| 96 | Pont de trànsit rodat GIV-6611 sobre el Ridaura | 193.76 | 51.6 / 199.0 / 341.7 |
| 134 | Pont de trànsit rodat Carrer de la Teulera sobre el Ridaura | 322.9 | 73.11 / 292.6 / 500.39 |
| 312 | Pont de trànsit rodat C-65 sobre el Ridaura | 361.13 | 63.4 / 251.6 / 419.2 |

## Congost (Granollers) (c188)

14311 seccions; 779 amb la llera desbordant abans de T10. Pluja: XEMA a 20 km (107 anys-estació), P0 = 65 mm (episodi més moll del model). 142 estrenyiments, 63 dins o a tocar de nucli.

| T desborda | Què | Riu (ACA) | Capacitat m³/s | Q10 | Q100 | Conca km² | Lloc |
|---|---|---|---|---|---|---|---|
| <10 | llera, 60 m, surt pel esquerre | El Congost | 25.9 | 59.1 | 133.6 | 18.34 | [41.79571, 2.23308](https://www.openstreetmap.org/?mlat=41.79571&mlon=2.23308#map=18/41.79571/2.23308) |
| <10 | llera, 60 m, surt pel esquerre | Riera d'Avencó | 27.6 | 70.6 | 191.7 | 35.5 | [41.76511, 2.26087](https://www.openstreetmap.org/?mlat=41.76511&mlon=2.26087#map=18/41.76511/2.26087) |
| <10 | pont ACA: Pont de trànsit rodat Pont del Poncic sobre El Con | El Congost | 53.4 | 134.5 | 338.0 | 63.08 | [41.76797, 2.24899](https://www.openstreetmap.org/?mlat=41.76797&mlon=2.24899#map=18/41.76797/2.24899) |
| <10 | llera, 140 m, surt pel esquerre | El Congost | 54.2 | 134.2 | 337.3 | 63.53 | [41.76635, 2.251](https://www.openstreetmap.org/?mlat=41.76635&mlon=2.251#map=18/41.76635/2.251) |
| <10 | llera, 60 m, surt pel dret | El Congost | 66.6 | 351.6 | 926.7 | 214.72 | [41.57626, 2.27465](https://www.openstreetmap.org/?mlat=41.57626&mlon=2.27465#map=18/41.57626/2.27465) |
| <10 | llera, 60 m, surt pel — | El Congost | 76.8 | 214.5 | 571.1 | 118.66 | [41.73085, 2.2679](https://www.openstreetmap.org/?mlat=41.73085&mlon=2.2679#map=18/41.73085/2.2679) |
| <10 | llera, 180 m, surt pel esquerre | El Congost | 77.0 | 249.0 | 676.4 | 150.31 | [41.68426, 2.28235](https://www.openstreetmap.org/?mlat=41.68426&mlon=2.28235#map=18/41.68426/2.28235) |
| <10 | llera, 80 m, surt pel esquerre | El Congost | 107.4 | 213.8 | 569.4 | 116.8 | [41.73746, 2.26648](https://www.openstreetmap.org/?mlat=41.73746&mlon=2.26648#map=18/41.73746/2.26648) |
| <10 | llera, 120 m, surt pel esquerre | Torrent de Can Bosquets | 139.1 | 257.5 | 695.6 | 155.16 | [41.66904, 2.2803](https://www.openstreetmap.org/?mlat=41.66904&mlon=2.2803#map=18/41.66904/2.2803) |
| 4 | llera, 100 m, surt pel — |  | 0.1 | 0.5 | 1.3 | 0.3 | [41.60653, 2.28321](https://www.openstreetmap.org/?mlat=41.60653&mlon=2.28321#map=18/41.60653/2.28321) |
| <10 | pont ACA: Pont de trànsit rodat C-17 sobre El Congost | El Congost | 156.7 | 243.8 | 658.4 | 132.51 | [41.71884, 2.27267](https://www.openstreetmap.org/?mlat=41.71884&mlon=2.27267#map=18/41.71884/2.27267) |
| <10 | llera, 140 m, surt pel esquerre | El Congost | 170.9 | 249.1 | 676.7 | 149.0 | [41.68996, 2.28488](https://www.openstreetmap.org/?mlat=41.68996&mlon=2.28488#map=18/41.68996/2.28488) |
| 5 | llera, 80 m, surt pel — |  | 0.4 | 0.8 | 2.2 | 0.42 | [41.61299, 2.28033](https://www.openstreetmap.org/?mlat=41.61299&mlon=2.28033#map=18/41.61299/2.28033) |
| 5 | llera, 140 m, surt pel dret |  | 2.6 | 8.3 | 25.1 | 11.38 | [41.80291, 2.22907](https://www.openstreetmap.org/?mlat=41.80291&mlon=2.22907#map=18/41.80291/2.22907) |
| 5 | llera, 80 m, surt pel esquerre |  | 3.6 | 8.3 | 25.3 | 11.37 | [41.80424, 2.22803](https://www.openstreetmap.org/?mlat=41.80424&mlon=2.22803#map=18/41.80424/2.22803) |

Ponts de l'inventari ACA que no passen el T500 (46; capacitat de la fitxa oficial):

| T desborda | Pont | Capacitat m³/s | Q10 / Q100 / Q500 |
|---|---|---|---|
| <10 | Pont de trànsit rodat C. del Gorg Negre sobre El Congost | 47.63 | 59.1 / 133.6 / 197.0 |
| <10 | Pont de trànsit rodat Ctra. De Ribes sobre El Congost | 101.91 | 134.5 / 338.0 / 519.8 |
| <10 | Pont de trànsit rodat Pont del Poncic sobre El Congost | 53.35 | 134.5 / 338.0 / 519.8 |
| <10 | Pont de trànsit rodat C-17 sobre El Congost | 156.67 | 243.8 / 658.4 / 1041.0 |
| <10 | Pont de trànsit rodat Camí Fluvial del Congost sobre El Congost | 166.16 | 243.8 / 658.4 / 1041.0 |
| <10 | Pont de trànsit rodat C. Moranta sobre El Congost | 209.45 | 249.0 / 676.4 / 1068.3 |
| <10 | Pont de trànsit rodat C. de Can Noguera sobre El Congost | 209.38 | 249.0 / 676.4 / 1068.3 |
| <10 | Pont de trànsit rodat C. de Josep Maria Ruera sobre El Congost | 216.07 | 249.0 / 676.4 / 1068.3 |
| <10 | Pont de trànsit rodat Camí de Can Bassa sobre El Congost | 336.0 | 351.6 / 926.7 / 1437.4 |
| 10 | Pont de trànsit rodat C-17 sobre El Congost | 214.37 | 214.0 / 569.8 / 895.2 |
| 11 | Pont de trànsit rodat N-152z sobre El Congost | 252.74 | 243.8 / 658.4 / 1041.0 |
| 12 | Pont de trànsit rodat C-17 sobre Torrent de Can Borrell | 7.01 | 6.68 / 14.09 / 20.94 |
| 13 | Pont de trànsit rodat C-17 sobre El Congost | 276.06 | 243.8 / 658.4 / 1041.0 |
| 16 | Pont de trànsit rodat C. de la Indústra sobre El Congost | 161.42 | 134.5 / 338.0 / 519.8 |
| 16 | Pont de trànsit rodat Camí de Querol sobre El Congost | 303.62 | 249.0 / 676.4 / 1068.3 |
| 18 | Pont de trànsit rodat C. del Camí Antic de Vic sobre El Congost | 275.93 | 214.0 / 569.8 / 895.2 |
| 18 | Pont de trànsit rodat C-352 sobre Torrent de Fangues | 27.61 | 22.13 / 52.57 / 82.66 |
| 19 | Pont de trànsit rodat Pont de Centelles sobre El Congost | 74.7 | 59.1 / 133.6 / 197.0 |
| 20 | Pont de trànsit rodat C. del Camí Antic de Vic sobre El Congost | 285.1 | 214.0 / 569.8 / 895.2 |
| 21 | Passarel·la per a vianants C. de la Rompuda sobre El Congost | 336.11 | 243.8 / 658.4 / 1041.0 |

## Cabrera de Mar (c199)

360 seccions; 10 amb la llera desbordant abans de T10. Pluja: XEMA a 20 km (88 anys-estació), P0 = 50 mm (episodi més moll del model). 2 estrenyiments, 1 dins o a tocar de nucli.

| T desborda | Què | Riu (ACA) | Capacitat m³/s | Q10 | Q100 | Conca km² | Lloc |
|---|---|---|---|---|---|---|---|
| 500 | pas: R1;RG1 | Torrent de Ca n'Amat | 18.7 | 9.2 | 15.0 | 1.08 | [41.50737, 2.40299](https://www.openstreetmap.org/?mlat=41.50737&mlon=2.40299#map=18/41.50737/2.40299) |

## Platja d'Aro (nucli) (c124)

136 seccions; 40 amb la llera desbordant abans de T10. Pluja: XEMA de la conca pilot ridaura (51 anys-estació), P0 = 50 mm (episodi més moll del model). 1 estrenyiments, 1 dins o a tocar de nucli.

| T desborda | Què | Riu (ACA) | Capacitat m³/s | Q10 | Q100 | Conca km² | Lloc |
|---|---|---|---|---|---|---|---|
| <10 | llera, 80 m, surt pel dret | el Ridaura | 39.7 | 87.7 | 275.5 | 0.86 | [41.81046, 3.05392](https://www.openstreetmap.org/?mlat=41.81046&mlon=3.05392#map=18/41.81046/3.05392) |

