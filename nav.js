// Barra de navegació comuna de totes les pestanyes.
// Per afegir una pestanya nova: posa-la a PAGINES i afegeix <script src="nav.js"></script> al final del <body>.
// Les pàgines amb `opcional: true` s'amaguen si encara no existeixen (p. ex. abans de fusionar-ne la PR).
(function () {
  const PAGINES = [
    { grup: 'Catalunya', items: [
      { href: 'catalunya.html', text: 'Pluja' },
      { href: 'conques.html', text: 'Conques' },
    ] },
    { grup: 'Casos d\'estudi', items: [
      { href: 'evolucio.html', text: 'Episodis (4 mapes)' },
      { href: 'episodi.html', text: 'Reproducció 4/10' },
      { href: 'pluja.html', text: 'Pluja sobre el relleu' },
    ] },
    { grup: 'En directe', items: [
      { href: 'previsio.html', text: 'Previsió Ridaura' },
    ] },
  ];

  const actual = (location.pathname.split('/').pop() || 'index.html');
  const css = `
    #nav-inund { position: fixed; top: 0; left: 0; right: 0; height: 40px; z-index: 2000; display: flex; align-items: center; gap: 4px;
      padding: 0 10px; background: #14212e; color: #e6edf3; font: 13px/1 system-ui, sans-serif; overflow-x: auto; white-space: nowrap;
      box-shadow: 0 1px 4px rgba(0,0,0,.25); scrollbar-width: none; }
    #nav-inund::-webkit-scrollbar { display: none; }
    #nav-inund a { color: #c9d6e2; text-decoration: none; padding: 7px 9px; border-radius: 6px; }
    #nav-inund a:hover { background: #22364a; color: #fff; }
    #nav-inund a.actiu { background: #2b7bba; color: #fff; font-weight: 600; }
    #nav-inund .marca { font-weight: 700; color: #fff; margin-right: 6px; padding-left: 4px; }
    #nav-inund .marca span { color: #6bb0e8; }
    #nav-inund .grup { color: #7f93a6; font-size: 11px; text-transform: uppercase; letter-spacing: .05em; margin: 0 2px 0 10px;
      padding-left: 10px; border-left: 1px solid #2c3f52; }
    #nav-inund .amaga { display: none; }
    html body { padding-top: 40px !important; box-sizing: border-box; }
    @media (max-width: 760px) { #nav-inund .grup { display: none; } }`;
  const st = document.createElement('style'); st.textContent = css; document.head.appendChild(st);

  const nav = document.createElement('nav');
  nav.id = 'nav-inund';
  nav.setAttribute('aria-label', 'Pestanyes');
  let html = `<a class="marca${actual === 'index.html' ? ' actiu' : ''}" href="./">≋ <span>Inundacions</span></a>`;
  for (const g of PAGINES) {
    html += `<span class="grup">${g.grup}</span>`;
    for (const p of g.items) {
      const cls = [p.href === actual ? 'actiu' : '', p.opcional ? 'amaga' : ''].filter(Boolean).join(' ');
      html += `<a href="${p.href}"${cls ? ` class="${cls}"` : ''}${p.opcional ? ' data-opcional' : ''}>${p.text}</a>`;
    }
  }
  nav.innerHTML = html;
  document.body.prepend(nav);

  // mostra les pestanyes opcionals que ja existeixen; amaga els grups buits
  const opc = [...nav.querySelectorAll('a[data-opcional]')];
  const grups = () => nav.querySelectorAll('.grup').forEach(g => {
    let n = g.nextElementSibling, alguna = false;
    while (n && !n.classList.contains('grup')) { if (n.tagName === 'A' && !n.classList.contains('amaga')) alguna = true; n = n.nextElementSibling; }
    g.classList.toggle('amaga', !alguna);
  });
  grups();
  Promise.all(opc.map(a => (a.getAttribute('href') === actual ? Promise.resolve(true)
      : fetch(a.getAttribute('href'), { method: 'HEAD' }).then(r => r.ok).catch(() => false))
    .then(ok => { if (ok) a.classList.remove('amaga'); }))).then(grups);

  // els mapes Leaflet s'han de recol·locar perquè l'alçada útil ha canviat
  window.addEventListener('load', () => setTimeout(() => window.dispatchEvent(new Event('resize')), 50));
})();
