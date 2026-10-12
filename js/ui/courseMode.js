// Moduli del corso di Chimica Generale. Solo formule e dati di riferimento:
// masse IUPAC del sito, N_A del SI, scattering coulombiano di Rutherford.
// Nessuna traiettoria o reazione inventata.

import { ELEMENTS, element } from '../physics/elements.js';

const NA = 6.02214076e23; // mol⁻¹, definizione SI dal 2019

// Isotopi usati nel corso (abbondanze IUPAC approssimate alle cifre delle slide).
// Masse: numero di massa A, come nel calcolo della slide sul cloro.
const ISOTOPES = {
  Cl: [
    { A: 35, Z: 17, n: 18, frac: 0.7576 },
    { A: 37, Z: 17, n: 20, frac: 0.2424 },
  ],
  Ne: [
    { A: 20, Z: 10, n: 10, frac: 0.9048 },
    { A: 21, Z: 10, n: 11, frac: 0.0027 },
    { A: 22, Z: 10, n: 12, frac: 0.0925 },
  ],
  Hg: [
    { A: 196, Z: 80, n: 116, frac: 0.0015 },
    { A: 198, Z: 80, n: 118, frac: 0.1004 },
    { A: 199, Z: 80, n: 119, frac: 0.1694 },
    { A: 200, Z: 80, n: 120, frac: 0.2313 },
    { A: 201, Z: 80, n: 121, frac: 0.1322 },
    { A: 202, Z: 80, n: 122, frac: 0.2986 },
    { A: 204, Z: 80, n: 124, frac: 0.0687 },
  ],
};

// Reazioni bilanciate del corso. Coefficienti fissi, non un bilanciatore automatico.
const REACTIONS = {
  combustione: {
    name: 'Combustione del metano',
    eq: 'CH₄ + 2 O₂ → CO₂ + 2 H₂O',
    reactants: [
      { id: 'CH4', label: 'CH₄', M: 16.0426, coeff: 1 },
      { id: 'O2', label: 'O₂', M: 31.9988, coeff: 2 },
    ],
    products: [
      { id: 'CO2', label: 'CO₂', M: 44.0095, coeff: 1 },
      { id: 'H2O', label: 'H₂O', M: 18.0153, coeff: 2 },
    ],
  },
  termite: {
    name: 'Termite (Fe₂O₃ + Al)',
    eq: 'Fe₂O₃ + 2 Al → 2 Fe + Al₂O₃',
    reactants: [
      { id: 'Fe2O3', label: 'Fe₂O₃', M: 159.688, coeff: 1 },
      { id: 'Al', label: 'Al', M: 26.9815, coeff: 2 },
    ],
    products: [
      { id: 'Fe', label: 'Fe', M: 55.845, coeff: 2 },
      { id: 'Al2O3', label: 'Al₂O₃', M: 101.961, coeff: 1 },
    ],
  },
  so2: {
    name: 'Ossidazione di SO₂',
    eq: '2 SO₂ + O₂ + 2 H₂O → 2 H₂SO₄',
    reactants: [
      { id: 'SO2', label: 'SO₂', M: 64.066, coeff: 2 },
      { id: 'O2', label: 'O₂', M: 31.9988, coeff: 1 },
      { id: 'H2O', label: 'H₂O', M: 18.0153, coeff: 2 },
    ],
    products: [
      { id: 'H2SO4', label: 'H₂SO₄', M: 98.079, coeff: 2 },
    ],
  },
};

const OX_RULES = [
  'Elemento non combinato: n.o. = 0 (es. H₂, O₂, Fe).',
  'La somma dei n.o. in un composto neutro è 0; in uno ione è la carica.',
  'Metalli del gruppo 1: +1. Gruppo 2: +2. Fluoro: −1.',
  'Idrogeno: +1 salvo negli idruri metallici (−1).',
  'Ossigeno: −2, salvo perossidi (−1) e superossidi (−1/2).',
  'Alogeni in binari: −1. Gruppo 16 in binari: −2. Gruppo 15 in binari: −3.',
];

function fmt(n, d = 4) {
  if (!Number.isFinite(n)) return '—';
  return n.toLocaleString('it-IT', { maximumFractionDigits: d });
}

function atomicMass(list) {
  return list.reduce((s, iso) => s + iso.A * iso.frac, 0);
}

function rutherfordAngle(bFm, EMeV, Z1, Z2) {
  // b = (Z1 Z2 e² / (8 π ε0 E)) cot(θ/2)
  // costante: e²/(4πε0) = 1.43996 MeV·fm
  const k = 1.439964548;
  const half = Math.atan((Z1 * Z2 * k) / (2 * EMeV * bFm));
  return 2 * half * (180 / Math.PI);
}

export function mountCourse(root) {
  root.innerHTML = `
    <div class="course">
      <nav class="course-tabs" role="tablist">
        <button type="button" data-tab="iso" class="on">Isotopi e massa</button>
        <button type="button" data-tab="mole">Mole</button>
        <button type="button" data-tab="stoich">Stechiometria</button>
        <button type="button" data-tab="ox">Numeri di ossidazione</button>
        <button type="button" data-tab="ruth">Rutherford</button>
      </nav>
      <div id="course-body"></div>
    </div>`;
  const body = root.querySelector('#course-body');
  const show = (tab) => {
    root.querySelectorAll('.course-tabs button').forEach(b => b.classList.toggle('on', b.dataset.tab === tab));
    if (tab === 'iso') renderIso(body);
    if (tab === 'mole') renderMole(body);
    if (tab === 'stoich') renderStoich(body);
    if (tab === 'ox') renderOx(body);
    if (tab === 'ruth') renderRuth(body);
  };
  root.querySelectorAll('.course-tabs button').forEach(b => b.addEventListener('click', () => show(b.dataset.tab)));
  show('iso');
}

function renderIso(body) {
  const opts = Object.keys(ISOTOPES).map(s => `<option value="${s}">${s}</option>`).join('');
  body.innerHTML = `
    <h2>Massa atomica media</h2>
    <p class="course-note">M = Σ Aᵢ · xᵢ. Le abbondanze sono quelle delle slide del corso (Ne, Cl) o i valori standard (Hg). Il numero di massa A è l’approssimazione usata a lezione; la massa IUPAC del sito è più precisa.</p>
    <label>Elemento <select id="iso-el">${opts}</select></label>
    <div id="iso-out"></div>`;
  const draw = () => {
    const key = body.querySelector('#iso-el').value;
    const list = ISOTOPES[key];
    const M = atomicMass(list);
    const el = ELEMENTS.find(e => e.symbol === key);
    const rows = list.map(iso => `<tr><td>${key}-${iso.A}</td><td>${iso.Z}</td><td>${iso.n}</td><td>${iso.A}</td><td>${fmt(iso.frac * 100, 2)} %</td></tr>`).join('');
    body.querySelector('#iso-out').innerHTML = `
      <table class="course-table">
        <thead><tr><th>Nuclide</th><th>Z</th><th>N</th><th>A</th><th>Abbondanza</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
      <p>M calcolata = <b>${fmt(M, 3)} u</b>${el ? ` · massa IUPAC nel sito: ${el.mass} u` : ''}</p>
      <p class="course-note">N = A − Z. Isotopi: stesso Z, N diverso. Isobari: stesso A. Isotoni: stesso N.</p>`;
  };
  body.querySelector('#iso-el').addEventListener('change', draw);
  draw();
}

function renderMole(body) {
  const opts = ELEMENTS.filter(e => !e.mass.startsWith('[')).map(e => `<option value="${e.Z}">${e.symbol} · ${e.mass} u</option>`).join('');
  body.innerHTML = `
    <h2>Mole</h2>
    <p class="course-note">1 mol = numero di atomi in 12 g di ¹²C. N<sub>A</sub> = 6,02214076×10²³ mol⁻¹ (esatto, SI 2019). n = m / M, N = n · N<sub>A</sub>.</p>
    <label>Elemento <select id="mole-el">${opts}</select></label>
    <label>Massa m (g) <input id="mole-m" type="number" step="any" value="12"></label>
    <div id="mole-out"></div>`;
  const draw = () => {
    const el = element(Number(body.querySelector('#mole-el').value));
    const m = Number(body.querySelector('#mole-m').value);
    const M = Number(el.mass);
    const n = m / M;
    const N = n * NA;
    body.querySelector('#mole-out').innerHTML = `
      <p>M(${el.symbol}) = ${el.mass} g/mol</p>
      <p>n = m/M = ${fmt(m, 4)} / ${el.mass} = <b>${fmt(n, 6)} mol</b></p>
      <p>N = n · N<sub>A</sub> = <b>${N.toExponential(4).replace('.', ',')}</b> atomi</p>`;
  };
  body.querySelector('#mole-el').addEventListener('change', draw);
  body.querySelector('#mole-m').addEventListener('input', draw);
  draw();
}

function renderStoich(body) {
  const opts = Object.entries(REACTIONS).map(([k, r]) => `<option value="${k}">${r.name}</option>`).join('');
  body.innerHTML = `
    <h2>Stechiometria e reagente limitante</h2>
    <p class="course-note">Reazioni già bilanciate. Il fattore stechiometrico è il rapporto dei coefficienti. Il reagente limitante è quello con il minor numero di moli di reazione.</p>
    <label>Reazione <select id="st-rx">${opts}</select></label>
    <div id="st-inputs"></div>
    <div id="st-out"></div>`;
  const draw = () => {
    const rx = REACTIONS[body.querySelector('#st-rx').value];
    const box = body.querySelector('#st-inputs');
    if (!box.dataset.rx || box.dataset.rx !== body.querySelector('#st-rx').value) {
      box.dataset.rx = body.querySelector('#st-rx').value;
      box.innerHTML = rx.reactants.map(r => `<label>${r.label} (g) <input data-id="${r.id}" type="number" step="any" value="10"></label>`).join('');
      box.querySelectorAll('input').forEach(i => i.addEventListener('input', draw));
    }
    const masses = {};
    box.querySelectorAll('input').forEach(i => { masses[i.dataset.id] = Number(i.value) || 0; });
    const extents = rx.reactants.map(r => ({ ...r, n: masses[r.id] / r.M, xi: (masses[r.id] / r.M) / r.coeff }));
    const lim = extents.reduce((a, b) => (b.xi < a.xi ? b : a));
    const products = rx.products.map(p => ({ ...p, g: lim.xi * p.coeff * p.M }));
    body.querySelector('#st-out').innerHTML = `
      <p>${rx.eq}</p>
      <table class="course-table">
        <thead><tr><th>Reagente</th><th>moli</th><th>ξ = n/coeff</th></tr></thead>
        <tbody>${extents.map(e => `<tr><td>${e.label}</td><td>${fmt(e.n, 4)}</td><td>${fmt(e.xi, 4)}</td></tr>`).join('')}</tbody>
      </table>
      <p>Reagente limitante: <b>${lim.label}</b> (ξ = ${fmt(lim.xi, 4)} mol)</p>
      <p>Resa teorica: ${products.map(p => `${p.label} ${fmt(p.g, 3)} g`).join(' · ')}</p>
      <p class="course-note">Resa percentuale = (resa reale / resa teorica) × 100. Qui si calcola solo la teorica.</p>`;
  };
  body.querySelector('#st-rx').addEventListener('change', () => { body.querySelector('#st-inputs').dataset.rx = ''; draw(); });
  draw();
}

function renderOx(body) {
  body.innerHTML = `
    <h2>Numeri di ossidazione</h2>
    <p class="course-note">Regole del corso, applicate a esempi standard. Non è un parser di formule: i casi sono quelli della lezione.</p>
    <ul class="course-rules">${OX_RULES.map(r => `<li>${r}</li>`).join('')}</ul>
    <label>Esempio
      <select id="ox-ex">
        <option value="h2o">H₂O</option>
        <option value="h2o2">H₂O₂ (perossido)</option>
        <option value="licl">LiCl</option>
        <option value="fe2o3">Fe₂O₃</option>
        <option value="so4">SO₄²⁻</option>
        <option value="nh4">NH₄⁺</option>
        <option value="pdh2">PdH₂ (idruro)</option>
      </select>
    </label>
    <div id="ox-out"></div>`;
  const cases = {
    h2o: 'H: +1, O: −2. Somma: 2·(+1) + (−2) = 0.',
    h2o2: 'Perossido: O ha n.o. −1. H: +1. Somma: 2·(+1) + 2·(−1) = 0.',
    licl: 'Gruppo 1: Li +1. Cl −1. Somma 0.',
    fe2o3: 'O: −2 × 3 = −6. Due Fe devono dare +6, quindi Fe +3 ciascuno.',
    so4: 'O: −2 × 4 = −8. Carica dello ione −2, quindi S = −2 − (−8) = +6.',
    nh4: 'H: +1 × 4 = +4. Carica +1, quindi N = +1 − 4 = −3.',
    pdh2: 'Idruro metallico: H −1. Due H danno −2, quindi Pd +2.',
  };
  const draw = () => { body.querySelector('#ox-out').innerHTML = `<p>${cases[body.querySelector('#ox-ex').value]}</p>`; };
  body.querySelector('#ox-ex').addEventListener('change', draw);
  draw();
}

function renderRuth(body) {
  body.innerHTML = `
    <h2>Scattering di Rutherford</h2>
    <p class="course-note">Deflessione coulombiana di una particella α su un nucleo. Formula: cot(θ/2) = 2Eb / (Z₁Z₂ k), con k = e²/(4πε₀) = 1,440 MeV·fm. Non è una traiettoria animata: solo l’angolo dalla formula.</p>
    <label>Nucleo bersaglio
      <select id="ru-z">
        <option value="79">Au (Z = 79)</option>
        <option value="47">Ag (Z = 47)</option>
        <option value="13">Al (Z = 13)</option>
      </select>
    </label>
    <label>Energia α (MeV) <input id="ru-e" type="number" step="any" value="5"></label>
    <label>Parametro d’impatto b (fm) <input id="ru-b" type="number" step="any" value="10"></label>
    <div id="ru-out"></div>`;
  const draw = () => {
    const Z2 = Number(body.querySelector('#ru-z').value);
    const E = Number(body.querySelector('#ru-e').value);
    const b = Number(body.querySelector('#ru-b').value);
    const theta = rutherfordAngle(b, E, 2, Z2);
    body.querySelector('#ru-out').innerHTML = `
      <p>Z₁ = 2 (α), Z₂ = ${Z2}, E = ${fmt(E, 2)} MeV, b = ${fmt(b, 2)} fm</p>
      <p>θ = <b>${fmt(theta, 2)}°</b></p>
      <p class="course-note">Nel modello a panettone di Thomson la carica è diffusa: deflessioni grandi sono rare. Nel modello nucleare di Rutherford un urto ravvicinato (b piccolo) dà θ grande, fino a retrodiffusione. L’esperimento della lamina d’oro (Geiger–Marsden, 1909) ha misurato quelle deflessioni.</p>`;
  };
  body.querySelectorAll('input, select').forEach(el => el.addEventListener('input', draw));
  draw();
}
