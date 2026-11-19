// Mapes de fons comuns per a les pestanyes amb Leaflet.
// fonsMapa(mapa, 'Gris') afegeix el selector de fons (amb opció «Sense mapa») i una capa de noms
// de llocs que es pinta per sobre de les dades (radar, conques), i que també es pot treure.
function fonsMapa(mapa, defecte = 'Topogràfic', atrib = '') {
  const ICGC = c => `https://geoserveis.icgc.cat/servei/catalunya/mapa-base/wmts/${c}/MON3857NW/{z}/{x}/{y}.png`;
  const A = (t) => t + (atrib ? ' · ' + atrib : '');
  const fons = {
    'Topogràfic': L.tileLayer(ICGC('topografic'), { maxZoom: 19, attribution: A('© ICGC') }),
    'Gris': L.tileLayer('https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}{r}.png', { maxZoom: 19, subdomains: 'abcd', attribution: A('© OpenStreetMap · © CARTO') }),
    'Satèl·lit': L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', { maxZoom: 19, attribution: A('© Esri, Maxar, Earthstar Geographics') }),
    'Sense mapa': L.layerGroup(),
  };
  mapa.createPane('etiquetes');
  const pane = mapa.getPane('etiquetes');
  pane.style.zIndex = 650; pane.style.pointerEvents = 'none';
  const noms = L.tileLayer('https://{s}.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}{r}.png',
    { maxZoom: 19, subdomains: 'abcd', pane: 'etiquetes', attribution: '© OpenStreetMap · © CARTO' });
  (fons[defecte] || fons['Topogràfic']).addTo(mapa);
  noms.addTo(mapa);
  const ctl = L.control.layers(fons, { 'Noms de pobles i carrers (a sobre)': noms }, // desplegat en ordinador; en mòbil, plegat (ocupa molt)
    { collapsed: window.matchMedia('(max-width: 760px)').matches }).addTo(mapa);
  return { fons, noms, ctl, senseMapa: () => mapa.hasLayer(fons['Sense mapa']) };
}
