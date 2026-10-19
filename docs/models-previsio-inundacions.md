# Models de previsió d'inundacions: què es fa al món

*Recerca inicial, 19/10/2026. Objectiu: triar l'enfocament per a un visor nou (tipus ONOFRE) que previngui crescudes per conca a Catalunya.*

## 1. Sistemes operatius de referència

| Sistema | On | Model | Què fa | Interès per a nosaltres |
|---|---|---|---|---|
| **EFAS** (Copernicus) | Europa | LISFLOOD (distribuït, ~1,5 km) + indicadors de crescuda sobtada (ERIC) | Previsió de cabal a 10-15 dies amb conjunts ECMWF/DWD; avís de pluja extrema en conques petites | Cobreix Catalunya. Es pot fer servir de referència, però va lent per a rieres |
| **GloFAS** | Món | LISFLOOD + ECMWF | El mateix a escala global, rius grans | Poc útil per a conques petites |
| **Google Flood Hub** | ~100 països | LSTM (IA) entrenada amb milers d'aforadors | Previsió a 7 dies, també en conques sense aforador | Demostra que la IA funciona en conques no aforades. Tancat, no integrable |
| **NOAA National Water Model** | EUA | WRF-Hydro (física completa, 2,7 M trams de riu) | Cabal a cada tram, de 18 h a 30 dies | Molt potent però massa pesat (supercomputador) |
| **Vigicrues Flash** (França) | França | Mètode AIGA: model GR distribuït + pluja radar | Avís de crescuda sobtada a ~35.000 punts sense aforador | **El més semblant al que volem**, i en clima mediterrani |
| **ERICHA / CRAHI-UPC** | Catalunya i Europa | Pluja radar (i nowcast) acumulada a la conca, comparada amb llindars de període de retorn | Avís de perill de crescuda sobtada a tota la xarxa de drenatge | Fet aquí, provat en episodis catalans |
| **FloodPROOFS** (CIMA, Itàlia) | Itàlia | Continuum (distribuït) + Delft-FEWS | Previsió operativa de Protecció Civil italiana | Codi obert, bona referència d'arquitectura |
| **Delft-FEWS** (Deltares) | Molts països | Plataforma, no model | Orquestra dades, models i assimilació | Massa gran per començar |

## 2. Famílies de models

1. **Índex de perill per pluja sobre la conca** (ERICHA, ERIC, AIGA simplificat). Es calcula la pluja mitjana acumulada aigües amunt de cada punt de riu (1, 3, 6, 24 h) i es compara amb llindars de període de retorn. No dona cabal, però funciona a cada riera sense aforador i és ràpid. Ideal per a crescudes sobtades mediterrànies.
2. **Models conceptuals agregats per subconca** (GR4J/GR5H, HBV, SCS-CN + hidrograma unitari, a HEC-HMS). Converteixen pluja en cabal amb pocs paràmetres, calibrats amb aforadors. La propagació riu avall es fa amb Muskingum i el temps de viatge entre subconques. Es poden corregir amb el nivell actual dels aforadors (assimilació).
3. **Models distribuïts amb física** (LISFLOOD, WRF-Hydro, Continuum, TOPKAPI). Són més realistes, però necessiten molta dada (sòls, usos), molta computació i molta calibració.
4. **IA / LSTM** (Google, NeuralHydrology, de codi obert). Aprenen pluja → cabal directament de l'històric. Funcionen molt bé quan hi ha anys de dades d'aforadors. L'ACA només en guarda 3 mesos, així que caldria arxiu propi o dades històriques.
5. **Hidràulics / mapes d'inundació** (HEC-RAS 2D, LISFLOOD-FP). Calculen l'extensió de l'aigua desbordada. Es fan servir offline per fer mapes; en temps real normalment s'usen mapes precalculats per cabal (a Catalunya, les zones inundables de l'ACA / INUNCAT).

## 3. Recomanació per al projecte

Fer-ho en tres capes, de més simple a més complexa:

1. **Capa 1, perill de crescuda sobtada (com ERICHA/AIGA).** Xarxa de drenatge i conques a partir del MDT de l'ICGC (pysheds o WhiteboxTools). A cada tram es calcula la pluja acumulada aigües amunt amb el radar (OPERA, que ONOFRE ja baixa), el nowcast i els models (AROME/WRF/ICON). Es compara amb llindars de període de retorn. Cobreix totes les rieres i es pot fer amb el que ja té ONOFRE.
2. **Capa 2, previsió de cabal a les conques amb aforador** (Llobregat, Ter, Besòs, Tordera, Francolí...). Es fa amb un model GR4J/GR5H (paquet airGR, de codi obert) per subconca, amb propagació Muskingum, i es corregeix amb el nivell actual dels aforadors de l'ACA i del SAIH. Dona hidrogrames previstos i hora estimada d'arribada del pic a cada estació, comparats amb els llindars de cabal que ja existeixen.
3. **Capa 3, més endavant.** Afegir un LSTM entrenat amb l'històric (quan n'hi hagi prou) i creuar el cabal previst amb les zones inundables de l'ACA per mostrar quines zones es mullarien.

Conca pilot suggerida: **Besòs o Tordera** (resposta ràpida, prou aforadors, impacte urbà), o el Llobregat si es vol una conca gran.

## Fonts
- [EFAS: indicadors de crescuda sobtada](https://european-flood.emergency.copernicus.eu/en/node/475)
- [Google Research: model d'IA de previsió d'inundacions global](https://research.google/blog/a-flood-forecasting-ai-model-trained-and-evaluated-globally/)
- [Google Research: IA per ampliar l'accés a previsions](https://research.google/blog/using-ai-to-expand-global-access-to-reliable-flood-forecasts/)
- [NOAA National Water Model](https://water.noaa.gov/about/nwm)
- [Vigicrues Flash (UNDRR/PreventionWeb)](https://www.preventionweb.net/publications/view/65750)
- [Mètode AIGA (INRAE)](https://hal.inrae.fr/hal-02608804)
- [CRAHI-UPC: ERICHA](https://www.crahi.upc.edu/en/ericha-3)
- [Comparació de dos sistemes d'alerta de crescudes sobtades (arXiv)](https://arxiv.org/pdf/1904.07213)
- [Plataformes operatives de previsió (WMO APFM)](https://www.floodmanagement.info/TT-E2/documents/info/Operational-platforms-for-flood-forecasting_TG2.docx)
