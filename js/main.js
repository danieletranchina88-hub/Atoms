// Atlante Orbitale: collega calcolo quantistico, scena 3D, grafici e interfaccia.

import { ELEMENTS, element, CATEGORY_LABELS } from './physics/elements.js';
import {
  configurationParts, isMadelungException, hundBoxes, L_LETTERS, subshellLabel,
} from './physics/configuration.js';
import { loadAtom, getOrbital } from './physics/atomStore.js';
import { term, samplePoints, sampleGrid, sampleAtom, evaluate, RadialFunction } from './physics/wavefunction.js';
import { orbitalName, mOrder, angularNodesDescription } from './physics/harmonics.js';
import { Viewer, formatPm } from './render/viewer.js';
import { buildPeriodicTable, blockOf, ramp } from './ui/periodicTable.js';
import { ATOM_SUMMARY } from './physics/atomSummary.js';
import { PAULING, ALLEN, ELECTRON_AFFINITY, OXIDATION_STATES } from './chem/elementData.js';
import { drawLineChart, drawLevels, drawMODiagram, drawSlice, superscript } from './ui/charts.js';
import { initSandboxMode, activateSandbox, deactivateSandbox, sandboxRedraw } from './ui/sandboxMode.js';

const HARTREE_EV = 27.211386245988;
const BOHR_PM = 52.917721090;

const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// ---------------------------------------------------------------------------
// Stato
// ---------------------------------------------------------------------------

const state = {
  Z: 6,
  atom: null,
  mode: 'atom',
  relativistic: true,       // equazione di Koelling–Harmon (ScRLDA del NIST); false: Schrödinger (LDA)
  atomView: 'density',      // 'density' | 'orbitals' — la densità totale è |ψ|², non un colore di sottolivello
  render: 'both',            // 'cloud' | 'surface' | 'both'
  enclosed: 0.9,
  points: 60000,
  quality: 72,
  orbital: { n: 2, l: 1, m: 0 },
  hidden: new Set(),
  bond: { kind: 'mo', mo: 'sigma-p', hybrid: 'sp3', which: 'all', distancePm: null },
  radialLog: false,
  slicePlane: 'auto',
  token: 0,
};

// ---------------------------------------------------------------------------
// Formattazione
// ---------------------------------------------------------------------------

const nf = (v, d = 2) => Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d });
const nfp = (v, p = 4) => {
  const a = Math.abs(v);
  // separatore delle migliaia con spazio sottile (stile SI), per non confonderlo con la virgola decimale
  if (a >= 1000) return String(Math.round(v)).replace(/\B(?=(\d{3})+(?!\d))/g, '\u202F');
  return Number(v.toPrecision(p)).toLocaleString('it-IT', { maximumFractionDigits: 6 });
};
const pm = (bohr) => `${nfp(bohr * BOHR_PM, 3)} pm`;
const eV = (ha) => `${nfp(ha * HARTREE_EV, 4)} eV`;

function orbitalHTML(n, l, m) {
  const { base, sub } = orbitalName(n, l, m);
  if (!sub) return base;
  if (sub.startsWith('m=')) return `${base}<sub>(${sub})</sub>`;
  return `${base}<sub>${sub}</sub>`;
}

function configHTML(Z, config) {
  const { core, subshells } = configurationParts(Z, config);
  const parts = subshells.map(s => `${subshellLabel(s.n, s.l)}<sup>${s.occ}</sup>`);
  return `${core ? `<span class="core">[${core}]</span> ` : ''}${parts.join(' ')}`;
}

// Colori per sottolivello: tinta dal guscio n (come uno spettro), luminosità da l.
const SHELL_HUES = [0, 8, 32, 50, 140, 190, 222, 272];
function subshellColor(n, l) {
  const hue = SHELL_HUES[n] ?? 0;
  const light = cssVar('--scene-mode') === 'light';
  const L = light ? 38 + l * 6 : 58 + l * 7;
  const S = light ? 75 : 85;
  return `hsl(${hue} ${S}% ${L}%)`;
}

function colorToRGB(color) {
  const c = document.createElement('canvas').getContext('2d');
  c.fillStyle = color;
  const hex = c.fillStyle; // normalizzato in #rrggbb
  const v = parseInt(hex.slice(1), 16);
  return [((v >> 16) & 255) / 255, ((v >> 8) & 255) / 255, (v & 255) / 255];
}

const HYBRID_COLORS = ['#ff6b57', '#3fa7ff', '#62d48f', '#f0b429', '#c77dff', '#4dd6c9'];
const SET_ENCLOSED = 0.55;

// ---------------------------------------------------------------------------
// Inizializzazione
// ---------------------------------------------------------------------------

const viewer = new Viewer($('viewport'));
const table = buildPeriodicTable($('periodic-table'), (Z) => selectElement(Z));

document.querySelectorAll('.modes button').forEach(b => {
  b.addEventListener('click', () => setMode(b.dataset.mode));
});

$('btn-reset').addEventListener('click', () => viewer.resetView());
$('btn-rotate').addEventListener('click', (e) => {
  const on = e.currentTarget.getAttribute('aria-pressed') !== 'true';
  e.currentTarget.setAttribute('aria-pressed', String(on));
  viewer.setAutoRotate(on);
});
$('btn-clip').addEventListener('click', (e) => {
  const on = e.currentTarget.getAttribute('aria-pressed') !== 'true';
  e.currentTarget.setAttribute('aria-pressed', String(on));
  viewer.setClipping(on);
});
const PT_PROPS = {
  chi: { label: 'χ Pauling', unit: '', get: (Z) => PAULING[Z] ?? null },
  ie: { label: 'energia di ionizzazione', unit: 'eV', get: (Z) => element(Z).ionizationEV },
  iecalc: { label: 'energia di ionizzazione calcolata (ΔSCF)', unit: 'eV', get: (Z) => ATOM_SUMMARY[Z]?.ie ?? null },
  ea: { label: 'affinità elettronica', unit: 'eV', get: (Z) => ELECTRON_AFFINITY[Z] ?? null },
  radius: { label: 'raggio dell\'orbitale più esterno', unit: 'pm', get: (Z) => ATOM_SUMMARY[Z]?.rmax ?? null },
};
$('pt-color').addEventListener('change', (e) => {
  const key = e.target.value;
  const prop = PT_PROPS[key];
  const box = $('pt-key');
  if (!prop) {
    table.heatmap(null);
    box.innerHTML = ['s', 'p', 'd', 'f'].map(b => `<span><i class="sw block-${b}"></i>blocco ${b}</span>`).join('');
    return;
  }
  const values = new Map(ELEMENTS.map(el => [el.Z, prop.get(el.Z)]));
  const range = table.heatmap(values);
  const stops = [0, 0.25, 0.5, 0.75, 1].map(t => `rgb(${ramp(t).join(',')}) ${t * 100}%`).join(', ');
  const fmtv = (v) => v.toLocaleString('it-IT', { maximumFractionDigits: 2 });
  box.innerHTML = `<span>${fmtv(range.min)} ${prop.unit}</span><i class="heat-bar" style="background: linear-gradient(90deg, ${stops})"></i><span>${fmtv(range.max)} ${prop.unit}</span><span class="pt-color-label">${key === 'radius' ? 'raggio di massima probabilità, calcolato' : ''}${key === 'iecalc' ? 'calcolo DFT-LDA relativistico scalare' : ''} · caselle grigie: dato non disponibile</span>`;
});

$('radial-log').addEventListener('change', (e) => { state.radialLog = e.target.checked; drawCharts(); });
$('slice-plane').addEventListener('change', (e) => { state.slicePlane = e.target.value; drawCharts(); });

// Tema: scuro di default; il pulsante alterna chiaro e scuro.
function currentTheme() {
  const t = document.documentElement.dataset.theme;
  if (t) return t;
  return window.matchMedia('(prefers-color-scheme: light)').matches ? 'light' : 'dark';
}
$('theme-toggle').addEventListener('click', () => {
  document.documentElement.dataset.theme = currentTheme() === 'dark' ? 'light' : 'dark';
  try { localStorage.setItem('atlante-theme', document.documentElement.dataset.theme); } catch { /* storage non disponibile */ }
  onThemeChange();
});
try {
  const saved = localStorage.getItem('atlante-theme');
  if (saved) document.documentElement.dataset.theme = saved;
} catch { /* storage non disponibile */ }
window.matchMedia('(prefers-color-scheme: light)').addEventListener('change', onThemeChange);
function onThemeChange() {
  viewer.applyTheme();
  if (state.mode === 'sandbox') { viewer.setAxesVisible(false); sandboxRedraw(); return; }
  if (state.atom) { render3D(); drawCharts(); }
}

let resizeTimer = null;
new ResizeObserver(() => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (state.mode === 'sandbox') sandboxRedraw();
    else if (state.atom) drawCharts();
  }, 120);
}).observe(document.querySelector('.charts'));

window.addEventListener('hashchange', () => {
  const h = location.hash.slice(1);
  if (h === 'sandbox' && state.mode !== 'sandbox') { setMode('sandbox'); return; }
  const Z = zFromHash();
  if (Z && Z !== state.Z) selectElement(Z);
});

function zFromHash() {
  const h = decodeURIComponent(location.hash.slice(1));
  if (!h) return null;
  const el = ELEMENTS.find(e => e.symbol.toLowerCase() === h.toLowerCase());
  return el ? el.Z : null;
}

// ---------------------------------------------------------------------------
// Selezione dell'elemento
// ---------------------------------------------------------------------------

async function selectElement(Z) {
  state.Z = Z;
  const token = ++state.token;
  table.select(Z);
  const el = element(Z);
  if (!isChemMode() && location.hash.slice(1) !== el.symbol) history.replaceState(null, '', `#${el.symbol}`);
  $('busy').hidden = false;
  $('busy-text').textContent = `Calcolo autoconsistente ${state.relativistic ? 'relativistico ' : ''}di ${el.name} (${Z} elettroni)…`;
  try {
    const atom = await loadAtom(Z, state.relativistic);
    if (token !== state.token) return;
    state.atom = atom;
    state.hidden.clear();
    // orbitale predefinito: il sottolivello occupato di energia più alta
    const homo = atom.orbitals.reduce((a, b) => (b.e > a.e ? b : a));
    state.orbital = { n: homo.n, l: homo.l, m: homo.l <= 2 ? 0 : -homo.l };
    state.bond.distancePm = null;
    if (!isChemMode()) renderAll();
  } catch (err) {
    $('busy-text').textContent = `Errore nel calcolo: ${err.message}`;
    return;
  }
  $('busy').hidden = true;
}

const CHEM_MODES = new Set(['sandbox']);
const isChemMode = () => CHEM_MODES.has(state.mode);

function setMode(mode) {
  const prev = state.mode;
  state.mode = mode;
  const models = {
    atom: 'Densità |ψ|² dal DFT-LDA relativistico scalare (Koelling–Harmon, validato sui dati NIST ScRLDA) · senza spin–orbita · nucleo non in scala',
    orbital: 'Densità di probabilità |ψ|² · orbitali del modello a campo centrale; nessuna traiettoria elettronica',
    sandbox: 'Modello scelto per la scatola: Lennard–Jones per i gas nobili · GFN2-xTB quantistico per tutti gli elementi fino al radon (fino a 120 atomi), anche in acqua implicita · MINDO/3 e Hartree–Fock a scelta · campo classico qualitativo solo su scelta esplicita · verifica la deriva energetica',
  };
  $('model-status').textContent = models[mode] ?? '';
  document.querySelectorAll('.modes button').forEach(b => b.setAttribute('aria-selected', String(b.dataset.mode === mode)));
  if (prev === 'sandbox' && mode !== 'sandbox') deactivateSandbox();
  $('analysis').hidden = mode !== 'sandbox';
  try { history.replaceState(null, '', mode === 'sandbox' ? '#sandbox' : `#${element(state.Z).symbol}`); } catch { /* ignora */ }
  viewer.renderer.domElement.style.visibility = '';
  if (mode === 'sandbox') { activateSandbox(); return; }
  document.body.dataset.mode = 'atom';
  viewer.setAxesVisible(true);
  $('busy').hidden = true;
  if (state.atom) renderAll();
}

/** Riporta le scelte dei legami a opzioni valide per l'elemento corrente. */
function sanitizeBond() { /* bond mode removed */ }

function renderAll() {
  renderElementCard();
  renderControls();
  renderPanelBody();
  render3D();
  drawCharts();
}
