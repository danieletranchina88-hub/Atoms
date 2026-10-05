// Atlante Orbitale: collega calcolo quantistico, scena 3D, grafici e interfaccia.

import { ELEMENTS, element, CATEGORY_LABELS } from './physics/elements.js';
import {
  configurationParts, isMadelungException, hundBoxes, L_LETTERS, subshellLabel,
} from './physics/configuration.js';
import { loadAtom, getOrbital } from './physics/atomStore.js';
import { term, samplePoints, sampleGrid, sampleAtom, evaluate, RadialFunction } from './physics/wavefunction.js';
import { orbitalName, mOrder, angularNodesDescription } from './physics/harmonics.js';
import {
  BOND_LENGTH_PM, MO_TYPES, moTerms, overlapIntegral, valenceMODiagram,
  HYBRIDS, hybridTerms, hybridEnergy,
} from './physics/bonds.js';
import { Viewer, formatPm } from './render/viewer.js';
import { buildPeriodicTable, blockOf, ramp } from './ui/periodicTable.js';
import { ATOM_SUMMARY } from './physics/atomSummary.js';
import { PAULING, ALLEN, ELECTRON_AFFINITY, OXIDATION_STATES } from './chem/elementData.js';
import { drawLineChart, drawLevels, drawMODiagram, drawSlice, superscript } from './ui/charts.js';
import { initMoleculeMode, activateMolecule, deactivateMolecule, moleculeThemeChanged, moleculeResize, openMoleculeSmiles } from './ui/moleculeMode.js';
import { initReactionMode, activateReaction, deactivateReaction, reactionRedraw } from './ui/reactionMode.js';
import { activateLab, deactivateLab, labRedraw } from './ui/labMode.js';
import { initSandboxMode, activateSandbox, deactivateSandbox, sandboxRedraw } from './ui/sandboxMode.js';
import { activatePhase, deactivatePhase, phaseRedraw } from './ui/phaseMode.js';
import { activateBeaker, deactivateBeaker, beakerRedraw } from './ui/beakerMode.js';
import { activateKinetics, deactivateKinetics, kineticsRedraw } from './ui/kineticsMode.js';
import { activateFormation, deactivateFormation, formationRedraw, initFormationMode } from './ui/bondFormationMode.js';

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
  if (state.mode === 'molecule') { viewer.setAxesVisible(false); moleculeThemeChanged(); return; }
  if (state.mode === 'reaction') { viewer.setAxesVisible(false); reactionRedraw(); return; }
  if (state.mode === 'lab') { labRedraw(); return; }
  if (state.mode === 'sandbox') { viewer.setAxesVisible(false); sandboxRedraw(); return; }
  if (state.mode === 'phase') { phaseRedraw(); return; }
  if (state.mode === 'beaker') { beakerRedraw(); return; }
  if (state.mode === 'kinetics') { kineticsRedraw(); return; }
  if (state.mode === 'formation') { formationRedraw(); return; }
  if (state.atom) { render3D(); drawCharts(); }
}

let resizeTimer = null;
new ResizeObserver(() => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (state.mode === 'molecule') moleculeResize();
    else if (state.mode === 'reaction') reactionRedraw();
    else if (state.mode === 'lab') labRedraw();
    else if (state.mode === 'sandbox') sandboxRedraw();
    else if (state.mode === 'phase') phaseRedraw();
    else if (state.mode === 'beaker') beakerRedraw();
    else if (state.mode === 'kinetics') kineticsRedraw();
    else if (state.mode === 'formation') formationRedraw();
    else if (state.atom) drawCharts();
  }, 120);
}).observe(document.querySelector('.charts'));

window.addEventListener('hashchange', () => {
  const h = location.hash.slice(1);
  if (h === 'molecole' && state.mode !== 'molecule') { setMode('molecule'); return; }
  if (h === 'reazioni' && state.mode !== 'reaction') { setMode('reaction'); return; }
  if (h === 'laboratorio' && state.mode !== 'lab') { setMode('lab'); return; }
  if (h === 'sandbox' && state.mode !== 'sandbox') { setMode('sandbox'); return; }
  if (h === 'becher' && state.mode !== 'beaker') { setMode('beaker'); return; }
  if (h === 'cinetica' && state.mode !== 'kinetics') { setMode('kinetics'); return; }
  if (h === 'legame' && state.mode !== 'formation') { setMode('formation'); return; }
  if (h === 'fasi' && state.mode !== 'phase') { setMode('phase'); return; }
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

const CHEM_MODES = new Set(['molecule', 'reaction', 'lab', 'sandbox', 'phase', 'beaker', 'kinetics', 'formation']);
const isChemMode = () => CHEM_MODES.has(state.mode);

function setMode(mode) {
  const prev = state.mode;
  state.mode = mode;
  const models = {
    atom: 'Densità |ψ|² dal DFT-LDA relativistico scalare (Koelling–Harmon, validato sui dati NIST ScRLDA) · senza spin–orbita · nucleo non in scala',
    orbital: 'Densità di probabilità |ψ|² · orbitali del modello a campo centrale; nessuna traiettoria elettronica',
    bond: 'LCAO: σ e π come combinazione degli orbitali calcolati · il bastoncino non è il legame',
    molecule: 'Hartree–Fock / basi gaussiane · molecole isolate; correlazione e solvente limitano l’accuratezza',
    reaction: 'Gas ideali · HF/MP2 + rotore rigido e oscillatore armonico · equilibrio distinto dalla cinetica',
    lab: 'Modelli termodinamici e cinetici · attività, unità e condizioni esplicite',
    sandbox: 'Modello scelto per la scatola: Lennard–Jones per i gas nobili · GFN2-xTB quantistico per tutti gli elementi fino al radon (fino a 120 atomi), anche in acqua implicita · MINDO/3 e Hartree–Fock a scelta · campo classico qualitativo solo su scelta esplicita · verifica la deriva energetica',
    phase: 'Dinamica Lennard–Jones · gas nobile modello · NVE / NVT · transizioni emergenti',
    beaker: 'Equilibrio in acqua (database MINTEQ v4) · van \'t Hoff per la temperatura · attività Debye–Hückel esteso/Davies, affidabilità limitata ad alta forza ionica · nessuna cinetica',
    kinetics: 'Legge di azione di massa · Arrhenius · integratore implicito per sistemi rigidi · isoterma, a volume costante',
    formation: 'Equazione di Schrödinger risolta a ogni distanza: Hartree–Fock RHF/UHF, FCI per H₂, GFN2-xTB · nuclei fermi (Born–Oppenheimer), non relativistico · confronto con le costanti spettroscopiche (Huber–Herzberg, NIST)',
  };
  $('model-status').textContent = models[mode] ?? '';
  document.querySelectorAll('.modes button').forEach(b => b.setAttribute('aria-selected', String(b.dataset.mode === mode)));
  if (prev === 'molecule' && mode !== 'molecule') deactivateMolecule();
  if (prev === 'reaction' && mode !== 'reaction') deactivateReaction();
  if (prev === 'lab' && mode !== 'lab') deactivateLab();
  if (prev === 'sandbox' && mode !== 'sandbox') deactivateSandbox();
  if (prev === 'phase' && mode !== 'phase') deactivatePhase();
  if (prev === 'beaker' && mode !== 'beaker') deactivateBeaker();
  if (prev === 'kinetics' && mode !== 'kinetics') deactivateKinetics();
  if (prev === 'formation' && mode !== 'formation') deactivateFormation();
  $('analysis').hidden = !CHEM_MODES.has(mode);
  const hashes = { molecule: 'molecole', reaction: 'reazioni', lab: 'laboratorio', sandbox: 'sandbox', beaker: 'becher', phase: 'fasi', kinetics: 'cinetica', formation: 'legame' };
  try { history.replaceState(null, '', CHEM_MODES.has(mode) ? `#${hashes[mode]}` : `#${element(state.Z).symbol}`); } catch { /* ignora */ }
  if (mode !== 'lab' && mode !== 'beaker' && mode !== 'kinetics' && mode !== 'formation') viewer.renderer.domElement.style.visibility = '';
  if (mode === 'molecule') { activateMolecule(); return; }
  if (mode === 'reaction') { activateReaction(); return; }
  if (mode === 'lab') { viewer.clear(); activateLab(); return; }
  if (mode === 'sandbox') { activateSandbox(); return; }
  if (mode === 'phase') { viewer.clear(); activatePhase(); return; }
  if (mode === 'beaker') { viewer.clear(); activateBeaker(); return; }
  if (mode === 'kinetics') { viewer.clear(); activateKinetics(); return; }
  if (mode === 'formation') { viewer.clear(); activateFormation(); return; }
  document.body.dataset.mode = 'atom';
  viewer.setAxesVisible(true);
  $('busy').hidden = true;
  if (state.atom) renderAll();
}

/** Riporta le scelte dei legami a opzioni valide per l'elemento corrente. */
function sanitizeBond() {
  const { bond, atom } = state;
  const n = valenceShell(atom).n;
  if (atom.Z <= 2) {
    bond.kind = 'mo';
    if (!MO_TYPES.find(t => t.id === bond.mo && t.ao === 's')) bond.mo = 'sigma-s';
  }
  if (HYBRIDS[bond.hybrid].needsD && n < 3) { bond.hybrid = 'sp3'; bond.which = 'all'; }
}

function renderAll() {
  sanitizeBond();
  renderElementCard();
  renderControls();
  renderPanelBody();
  render3D();
  drawCharts();
}

// ---------------------------------------------------------------------------
// Scheda dell'elemento
// ---------------------------------------------------------------------------

function renderElementCard() {
  const { atom } = state;
  const el = element(atom.Z);
  const block = blockOf(el);
  const t = atom.term;
  const termHTML = `<sup>${t.multiplicity}</sup>${t.letter}${t.odd ? '<sup>o</sup>' : ''}<sub>${t.Jlabel}</sub>`;
  const ieCalc = atom.ionization.eV;
  const ieExp = el.ionizationEV;
  let ieDelta = '';
  if (ieExp) {
    const d = ieCalc - ieExp;
    ieDelta = `<span class="${Math.abs(d) < 0.6 ? 'delta-good' : 'delta-warn'}">${d >= 0 ? '+' : '−'}${nf(Math.abs(d))}</span>`;
  }
  const from = atom.ionization.from;
  const valence = atom.orbitals.reduce((a, b) => (b.e > a.e ? b : a));
  const exc = isMadelungException(atom.Z);
  $('element-card').innerHTML = `
    <div class="el-head">
      <div class="el-symbol block-${block}"><span class="z">${el.Z}</span><span class="sym">${el.symbol}</span></div>
      <div>
        <p class="el-name">${el.name}</p>
        <p class="el-cat">${CATEGORY_LABELS[el.category]}</p>
        <p class="el-cat">Periodo ${el.period}${el.group ? ` · gruppo ${el.group}` : ''} · blocco ${block}</p>
      </div>
    </div>
    <p class="config">${configHTML(atom.Z, atom.config)}</p>
    <p class="config-note">${exc ? '<b>Eccezione alla regola di Madelung</b> (configurazione sperimentale).' : 'Configurazione dello stato fondamentale.'}</p>
    <dl class="facts">
      <dt>Massa atomica</dt><dd>${el.mass} u</dd>
      <dt>χ di Pauling</dt><dd>${PAULING[el.Z] !== null && PAULING[el.Z] !== undefined ? nf(PAULING[el.Z]) : '—'}</dd>
      ${ALLEN[el.Z] ? `<dt>χ di Allen</dt><dd>${nf(ALLEN[el.Z], 3)}</dd>` : ''}
      ${ELECTRON_AFFINITY[el.Z] !== undefined ? `<dt>Affinità elettronica</dt><dd>${ELECTRON_AFFINITY[el.Z] > 0 ? `${nf(ELECTRON_AFFINITY[el.Z], 3)} eV` : 'anione instabile'}</dd>` : ''}
      <dt>Stati di ossidazione</dt><dd>${(OXIDATION_STATES[el.Z] ?? []).map(v => (v > 0 ? `+${v}` : v < 0 ? `−${-v}` : '0')).join(', ') || '—'}</dd>
      <dt>Termine di Hund</dt><dd class="term">${termHTML}</dd>
      <dt>Elettroni spaiati</dt><dd>${t.unpaired}</dd>
      <div class="sep"></div>
      <dt>Ionizzazione calc.</dt><dd>${nf(ieCalc)} eV</dd>
      <dt>Ionizzazione sper.</dt><dd>${ieExp ? `${nf(ieExp)} eV` : 'non misurata'}</dd>
      ${ieExp ? `<dt>Scarto</dt><dd>${ieDelta} eV</dd>` : ''}
      <dt>Elettrone rimosso</dt><dd>${subshellLabel(from.n, from.l)}</dd>
      <div class="sep"></div>
      <dt>Orbitale esterno</dt><dd>${valence.label}, ${eV(valence.e)}</dd>
      <dt>Raggio (90% e⁻)</dt><dd>${pm(atom.r90)}</dd>
      <dt>Energia totale</dt><dd>${nfp(atom.energy.total, 7)} Ha</dd>
      <dt>Teoria</dt><dd>${atom.relativistic ? 'DFT-LDA relativistica scalare' : 'DFT-LDA non relativistica'}</dd>
      <dt>Convergenza SCF</dt><dd>${atom.converged ? 'raggiunta' : 'NON raggiunta: risultato non validato'}</dd>
      <dt>Iterazioni SCF</dt><dd>${atom.iterations} · ${atom.elapsed} ms</dd>
    </dl>
    <div id="rel-effect" class="rel-effect"></div>`;
  relativisticEffect(atom);
}

/**
 * Effetto relativistico: confronto con il calcolo dell'altra teoria (calcolato in secondo piano e in cache).
 * Negli atomi pesanti gli elettroni s vicini al nucleo hanno velocità ≈ Zα·c: la massa relativistica li contrae
 * e abbassa di energia gli orbitali s e p₁/₂ esterni, mentre d e f, più schermati, si espandono.
 */
async function relativisticEffect(atom) {
  const box = $('rel-effect');
  if (!box) return;
  const token = state.token;
  box.innerHTML = '<p class="desc-muted">Calcolo il confronto relativistico…</p>';
  const other = await loadAtom(atom.Z, !atom.relativistic).catch(() => null);
  if (token !== state.token || !other || !$('rel-effect')) return;
  const [sr, nr] = atom.relativistic ? [atom, other] : [other, atom];
  const pick = (a, n, l) => a.orbitals.find(o => o.n === n && o.l === l);
  const outer = sr.orbitals.filter(o => o.l <= 2).sort((a, b) => b.e - a.e).slice(0, 2);
  const rows = outer.map(o => {
    const q = pick(nr, o.n, o.l);
    if (!q) return '';
    const dr = (o.rAvg / q.rAvg - 1) * 100;
    return `<tr><td>${o.label}</td><td class="num">${nfp(q.e * HARTREE_EV, 4)} → ${nfp(o.e * HARTREE_EV, 4)}</td><td class="num">${dr >= 0 ? '+' : '−'}${nf(Math.abs(dr), 1)}%</td></tr>`;
  }).join('');
  const dIE = sr.ionization.eV - nr.ionization.eV;
  $('rel-effect').innerHTML = `
    <h3 class="side-h">Effetto relativistico</h3>
    <table class="orb-table"><thead><tr><th>Orbitale</th><th class="num">ε senza → con (eV)</th><th class="num">⟨r⟩</th></tr></thead><tbody>${rows}</tbody></table>
    <p class="desc-muted">Ionizzazione: ${nf(nr.ionization.eV)} eV senza relatività, ${nf(sr.ionization.eV)} eV con (${dIE >= 0 ? '+' : '−'}${nf(Math.abs(dIE))} eV). Energia totale: ${nfp((sr.energy.total - nr.energy.total) * HARTREE_EV, 4)} eV.
    ${atom.Z >= 55 ? 'Negli atomi pesanti gli s interni si muovono a una frazione importante della velocità della luce: la contrazione degli s esterni spiega il colore dell\'oro e il mercurio liquido.' : 'Per gli atomi leggeri l\'effetto è piccolo.'} Senza spin–orbita (p₁/₂ e p₃/₂ mediati).</p>`;
}

// ---------------------------------------------------------------------------
// Controlli (dipendono dalla modalità)
// ---------------------------------------------------------------------------

function segHTML(id, options, current, cls = '') {
  return `<div class="seg ${cls}" id="${id}">${options.map(o =>
    `<button type="button" data-v="${o.v}" aria-pressed="${String(o.v) === String(current)}" ${o.disabled ? 'disabled' : ''} ${o.title ? `title="${o.title}"` : ''}>${o.label}</button>`).join('')}</div>`;
}

function bindSeg(id, fn) {
  const root = $(id);
  if (!root) return;
  root.querySelectorAll('button').forEach(b => b.addEventListener('click', () => fn(b.dataset.v)));
}

function renderControls() {
  const c = $('controls');
  const { mode, atom } = state;
  const renderOpts = [
    { v: 'cloud', label: 'Nuvola' },
    { v: 'surface', label: 'Superficie' },
    { v: 'both', label: 'Entrambe' },
  ];
  const pointsCtl = `
    <div class="ctl"><span class="lbl">Punti Monte Carlo</span>
      <div class="range-row"><input type="range" id="points" min="5000" max="200000" step="5000" value="${state.points}"><output id="points-out">${state.points.toLocaleString('it-IT')}</output></div>
    </div>`;
  const surfCtl = `
    <div class="ctl"><span class="lbl">Probabilità racchiusa dalla superficie</span>
      <div class="range-row"><input type="range" id="enclosed" min="0.5" max="0.98" step="0.01" value="${state.enclosed}"><output id="enclosed-out">${Math.round(state.enclosed * 100)}%</output></div>
    </div>
    <div class="ctl"><span class="lbl">Risoluzione della superficie</span>
      ${segHTML('quality', [{ v: 48, label: 'Bassa' }, { v: 72, label: 'Media' }, { v: 100, label: 'Alta' }], state.quality)}
    </div>`;

  const theoryCtl = `<div class="ctl"><span class="lbl">Teoria</span>
      ${segHTML('theory', [{ v: 'sr', label: 'Relativistica', title: 'Equazione di Koelling–Harmon + scambio di MacDonald–Vosko (NIST ScRLDA)' }, { v: 'nr', label: 'Non relativistica', title: 'Equazione di Schrödinger, LDA (NIST LDA)' }], state.relativistic ? 'sr' : 'nr')}</div>`;
  if (mode === 'atom') {
    c.innerHTML = `
      ${theoryCtl}
      <div class="ctl"><span class="lbl">Rappresentazione dell'atomo</span>
        ${segHTML('atom-view', [{ v: 'orbitals', label: 'Orbitali occupati' }, { v: 'density', label: 'Densità totale' }], state.atomView)}
      </div>
      ${pointsCtl}
      <p class="desc-muted">${state.atomView === 'orbitals'
        ? 'Gli orbitali sono funzioni di probabilità, non traiettorie. Visualizzazione nella base reale (p<sub>x</sub>, d<sub>xy</sub>…), riempito secondo la regola di Hund. Colori per sottolivello.'
        : 'Densità elettronica totale a simmetria sferica: si vedono i gusci K, L, M… Un sottolivello pieno è sempre sferico (teorema di Unsöld).'}</p>`;
    bindSeg('atom-view', v => { state.atomView = v; renderControls(); render3D(); drawCharts(); });
  } else if (mode === 'orbital') {
    const { n, l, m } = state.orbital;
    const nOpts = [1, 2, 3, 4, 5, 6, 7].map(v => ({ v, label: v }));
    const lOpts = [0, 1, 2, 3, 4, 5, 6].map(v => ({ v, label: L_LETTERS[v], disabled: v >= n }));
    const mOpts = mOrder(l).map(v => ({ v, label: orbitalHTML('', l, v) }));
    c.innerHTML = `
      ${theoryCtl}
      <div class="ctl"><span class="lbl">Numero quantico principale n</span>${segHTML('q-n', nOpts, n, 'qn')}</div>
      <div class="ctl"><span class="lbl">Numero quantico secondario l</span>${segHTML('q-l', lOpts, l, 'qn')}</div>
      <div class="ctl"><span class="lbl">Orbitale reale (m<sub>l</sub>)</span>${segHTML('q-m', mOpts, m)}</div>
      <div class="ctl"><span class="lbl">Rappresentazione</span>${segHTML('render', renderOpts, state.render)}</div>
      ${state.render !== 'cloud' ? surfCtl : ''}
      ${state.render !== 'surface' ? pointsCtl : ''}
      <div class="ctl"><label class="lbl" for="orb-probe">Sonda radiale r (a₀)</label>
        <input id="orb-probe" type="range" min="0" max="${getOrbital(atom, n, l).radial.radiusEnclosing(0.999)}" step="any" value="${getOrbital(atom, n, l).rAvg}">
        <output id="orb-probe-out"></output></div>`;
    const radial = getOrbital(atom, n, l).radial;
    const probe = () => {
      const r = +$('orb-probe').value;
      $('orb-probe-out').textContent = `${nf(r, 3)} a₀ · P(r′ ≤ r) = ${nf(radial.probabilityWithin(r) * 100, 2)}% · nodi: ${n-l-1} radiali, ${l} angolari`;
    };
    $('orb-probe').addEventListener('input', probe); probe();
    bindSeg('q-n', v => {
      const nn = +v;
      state.orbital.n = nn;
      if (state.orbital.l >= nn) { state.orbital.l = nn - 1; state.orbital.m = 0; }
      refreshOrbital();
    });
    bindSeg('q-l', v => { state.orbital.l = +v; state.orbital.m = 0; refreshOrbital(); });
    bindSeg('q-m', v => { state.orbital.m = +v; refreshOrbital(); });
    bindSeg('render', v => { state.render = v; renderControls(); render3D(); });
  } else {
    const { bond } = state;
    const el = element(atom.Z);
    const val = valenceShell(atom);
    const moOpts = MO_TYPES.map(t => ({ v: t.id, label: `${t.label}<sub>${val.n}${t.sub}</sub>`, disabled: t.ao !== 's' && atom.Z <= 2 }));
    const hybOpts = Object.entries(HYBRIDS).map(([k, h]) => ({ v: k, label: h.label, disabled: (h.needsD && val.n < 3) || atom.Z <= 2 }));
    const def = defaultBondPm(atom);
    const dist = bond.distancePm ?? def;
    const set = HYBRIDS[bond.hybrid];
    const whichOpts = [{ v: 'all', label: 'Tutti' }, ...set.set.map((_, i) => ({ v: i, label: `h${i + 1}` }))];
    c.innerHTML = `
      <div class="ctl"><span class="lbl">Tipo</span>
        ${segHTML('bond-kind', [{ v: 'mo', label: `Orbitali molecolari ${el.symbol}₂` }, { v: 'hybrid', label: 'Orbitali ibridi' }], bond.kind)}
      </div>
      ${bond.kind === 'mo' ? `
        <div class="ctl"><span class="lbl">Orbitale molecolare</span>${segHTML('mo-type', moOpts, bond.mo)}</div>
        <div class="ctl"><span class="lbl">Distanza tra i nuclei</span>
          <div class="range-row"><input type="range" id="bond-d" min="${Math.round(def * 0.4)}" max="${Math.round(def * 2.5)}" step="1" value="${Math.round(dist)}"><output id="bond-d-out">${Math.round(dist)} pm</output></div>
        </div>` : `
        <div class="ctl"><span class="lbl">Ibridazione</span>${segHTML('hyb-type', hybOpts, bond.hybrid)}</div>
        <div class="ctl"><span class="lbl">Ibrido mostrato</span>${segHTML('hyb-which', whichOpts, bond.which)}</div>`}
      <div class="ctl"><span class="lbl">Rappresentazione</span>${segHTML('render', renderOpts, state.render)}</div>
      ${state.render !== 'cloud' ? surfCtl : ''}
      ${state.render !== 'surface' ? pointsCtl : ''}`;
    bindSeg('bond-kind', v => { bond.kind = v; renderControls(); renderPanelBody(); render3D(); drawCharts(); });
    bindSeg('mo-type', v => { bond.mo = v; renderControls(); renderPanelBody(); render3D(); drawCharts(); });
    bindSeg('hyb-type', v => { bond.hybrid = v; bond.which = 'all'; renderControls(); renderPanelBody(); render3D(); drawCharts(); });
    bindSeg('hyb-which', v => { bond.which = v === 'all' ? 'all' : +v; renderControls(); render3D(); drawCharts(); });
    bindSeg('render', v => { state.render = v; renderControls(); render3D(); });
    const d = $('bond-d');
    if (d) {
      d.addEventListener('input', () => { $('bond-d-out').textContent = `${d.value} pm`; });
      d.addEventListener('change', () => { bond.distancePm = +d.value; renderPanelBody(); render3D(); drawCharts(); });
    }
  }

  bindSeg('theory', v => {
    state.relativistic = v === 'sr';
    const keep = { ...state.orbital };
    selectElement(state.Z).then(() => { if (state.atom?.orbitals.some(o => o.n === keep.n && o.l === keep.l) || state.mode === 'orbital') { state.orbital = keep; if (!isChemMode()) renderAll(); } });
  });
  const pts = $('points');
  if (pts) {
    pts.addEventListener('input', () => { $('points-out').textContent = (+pts.value).toLocaleString('it-IT'); });
    pts.addEventListener('change', () => { state.points = +pts.value; render3D(); });
  }
  const enc = $('enclosed');
  if (enc) {
    enc.addEventListener('input', () => { $('enclosed-out').textContent = `${Math.round(enc.value * 100)}%`; });
    enc.addEventListener('change', () => { state.enclosed = +enc.value; render3D(); });
  }
  bindSeg('quality', v => { state.quality = +v; renderControls(); render3D(); });
}

function refreshOrbital() {
  renderControls();
  renderPanelBody();
  render3D();
  drawCharts();
}

// ---------------------------------------------------------------------------
// Corpo del pannello: tabella degli orbitali / dettagli
// ---------------------------------------------------------------------------

function boxesHTML(l, occ) {
  return `<span class="boxes">${hundBoxes(l, occ).map(b =>
    `<span class="box">${b.up ? '↑' : ''}${b.down ? '↓' : ''}</span>`).join('')}</span>`;
}

function renderPanelBody() {
  const body = $('panel-body');
  const { atom, mode } = state;
  if (mode === 'atom' || mode === 'orbital') {
    const rows = atom.orbitals.map(o => {
      const active = mode === 'orbital' && o.n === state.orbital.n && o.l === state.orbital.l;
      const checked = !state.hidden.has(o.label);
      return `<tr class="clickable ${active ? 'active' : ''}" data-n="${o.n}" data-l="${o.l}">
        ${mode === 'atom' ? `<td><input type="checkbox" class="vis-toggle" data-label="${o.label}" ${checked ? 'checked' : ''} aria-label="Mostra ${o.label}"></td>` : ''}
        <td><span class="orb-label"><i style="background:${subshellColor(o.n, o.l)}"></i>${o.label}</span></td>
        <td>${o.l <= 3 ? boxesHTML(o.l, o.occ) : o.occ}</td>
        <td class="num">${nfp(o.e * HARTREE_EV, 4)}</td>
        <td class="num">${nf(o.zeffSlater)}</td>
        <td class="num">${nfp(o.rMaxProb * BOHR_PM, 3)}</td>
      </tr>`;
    }).join('');
    let details = '';
    if (mode === 'orbital') details = orbitalDetailsHTML();
    body.innerHTML = `
      ${details}
      <div>
        <h3>Orbitali occupati</h3>
        <div class="table-scroll">
        <table class="orb-table">
          <thead><tr>${mode === 'atom' ? '<th></th>' : ''}<th>Orbitale</th><th>Elettroni</th><th class="num">ε (eV)</th><th class="num" title="Carica nucleare efficace, regole di Slater">Z<sub>eff</sub></th><th class="num" title="Raggio di massima probabilità radiale">r<sub>max</sub> (pm)</th></tr></thead>
          <tbody>${rows}</tbody>
        </table>
        </div>
        <p class="desc-muted" style="margin-top:8px">ε: energia dell'orbitale calcolata${atom.relativistic ? ' con l\'equazione relativistica scalare (p e d: media sullo spin–orbita)' : ''}. Z<sub>eff</sub>: regole di Slater. Clicca una riga per vedere l'orbitale.</p>
      </div>`;
    body.querySelectorAll('tr.clickable').forEach(tr => tr.addEventListener('click', (ev) => {
      if (ev.target.classList.contains('vis-toggle')) return;
      const n = +tr.dataset.n;
      const l = +tr.dataset.l;
      state.orbital = { n, l, m: l <= 2 ? 0 : -l };
      setMode('orbital');
    }));
    body.querySelectorAll('.vis-toggle').forEach(cb => cb.addEventListener('change', () => {
      if (cb.checked) state.hidden.delete(cb.dataset.label); else state.hidden.add(cb.dataset.label);
      render3D();
      drawCharts();
    }));
  } else {
    body.innerHTML = bondDetailsHTML();
  }
}

function orbitalDetailsHTML() {
  const { n, l, m } = state.orbital;
  const orb = getOrbital(state.atom, n, l);
  const radialNodes = orb.nodes.map(r => nfp(r * BOHR_PM, 3)).join(' · ');
  const status = orb.occupied
    ? `<span class="pill occ">occupato · ${orb.occ} e⁻ nel sottolivello</span>`
    : `<span class="pill virt">non occupato (stato eccitato)</span>`;
  return `
    <div>
      <h3>Orbitale ${orbitalHTML(n, l, m)}</h3>
      <p style="margin:0 0 10px">${status}</p>
      <dl class="info-list">
        <dt>Numeri quantici</dt><dd>n = ${n}, l = ${l}, m<sub>l</sub> = ${m}</dd>
        <dt>Energia ε</dt><dd>${eV(orb.e)}</dd>
        <dt>Nodi radiali n−l−1</dt><dd>${n - l - 1}</dd>
        ${orb.nodes.length ? `<dt>Posizione dei nodi</dt><dd>${radialNodes} pm</dd>` : ''}
        <dt>Nodi angolari l</dt><dd>${l}</dd>
        <dt>Superfici nodali</dt><dd class="wrap">${angularNodesDescription(l, m)}</dd>
        <dt>Raggio più probabile</dt><dd>${pm(orb.rMaxProb)}</dd>
        <dt>Raggio medio ⟨r⟩</dt><dd>${pm(orb.rAvg)}</dd>
        <dt>Raggio al 90%</dt><dd>${pm(orb.r90)}</dd>
      </dl>
      <div class="phase-key" style="margin-top:10px">
        <span><i style="background:var(--phase-pos)"></i>ψ &gt; 0</span>
        <span><i style="background:var(--phase-neg)"></i>ψ &lt; 0</span>
      </div>
    </div>`;
}

/** Guscio di valenza: n più alto tra s/p occupati, oppure (n+1) dopo un d, (n+2) dopo un f (es. Pd: 5s vuoto). */
function valenceShell(atom) {
  let n = 1;
  for (const o of atom.orbitals) n = Math.max(n, o.l <= 1 ? o.n : o.n + o.l - 1);
  return { n };
}

function defaultBondPm(atom) {
  if (BOND_LENGTH_PM[atom.Z]) return BOND_LENGTH_PM[atom.Z];
  // Stima: due volte il raggio di massima probabilità dell'orbitale di valenza, ridotto del 25%
  const valence = atom.orbitals.reduce((a, b) => (b.e > a.e ? b : a));
  return Math.round(2 * valence.rMaxProb * BOHR_PM * 0.75);
}

/**
 * Per i legami si usa la convenzione chimica: il lobo radiale più esterno (quello di valenza) è positivo.
 * (Nella convenzione dei polinomi di Laguerre è positivo il lobo vicino al nucleo, e il 2s è negativo fuori.)
 */
function valenceSign(orb) {
  if (!orb.valenceRadial) {
    const R = orb.radial;
    const flip = R.at(orb.rMaxProb) < 0;
    orb.valenceRadial = flip ? new RadialFunction(state.atom.fullGrid, orb.u.map(v => -v), orb.l) : R;
  }
  return { ...orb, radial: orb.valenceRadial };
}

function bondRadials(atom) {
  const n = valenceShell(atom).n;
  const s = valenceSign(getOrbital(atom, n, 0));
  // per H e He non esiste 1p: si usa il 2p (non occupato)
  const p = valenceSign(getOrbital(atom, Math.max(n, 2), 1));
  const d = n >= 3 ? valenceSign(getOrbital(atom, n, 2)) : null;
  return { n, s, p, d };
}

function bondDetailsHTML() {
  const { atom, bond } = state;
  const el = element(atom.Z);
  const { n, s, p, d } = bondRadials(atom);
  if (bond.kind === 'mo') {
    const type = MO_TYPES.find(t => t.id === bond.mo);
    const dist = bond.distancePm ?? defaultBondPm(atom);
    const terms = moTerms(type, s.radial, p.radial, dist / BOHR_PM);
    const radial = type.ao === 's' ? s : p;
    const S = overlapIntegral(terms, radial.radial.radiusEnclosing(0.999));
    const effective = terms[0].coeff * terms[1].coeff * S;
    const valenceE = atom.config.filter(c => c.n === n && c.l <= 1).reduce((a, c) => a + c.occ, 0) * 2;
    const diag = valenceMODiagram(atom.Z, valenceE);
    const known = BOND_LENGTH_PM[atom.Z] !== undefined;
    return `
      <div>
        <h3>${type.label}<sub>${n}${type.sub}</sub> di ${el.symbol}₂</h3>
        <p class="desc">${type.text}</p>
        <dl class="info-list" style="margin-top:10px">
          <dt>Orbitali atomici combinati</dt><dd>${n}${type.ao === 's' ? 's' : type.ao === 'pz' ? 'p<sub>z</sub>' : 'p<sub>x</sub>'} + ${n}${type.ao === 's' ? 's' : type.ao === 'pz' ? 'p<sub>z</sub>' : 'p<sub>x</sub>'}</dd>
          <dt>Distanza tra i nuclei</dt><dd>${Math.round(dist)} pm${known && bond.distancePm === null ? ' (sperim.)' : ''}</dd>
          <dt>Sovrapposizione S = ⟨φ<sub>A</sub>|φ<sub>B</sub>⟩</dt><dd>${nf(Math.abs(S), 3)}</dd>
          <dt>Interferenza tra i nuclei</dt><dd>${effective >= 0 ? 'costruttiva' : 'distruttiva'}</dd>
          <dt>Normalizzazione</dt><dd>1/√(2 ${effective >= 0 ? '+' : '−'} 2S)</dd>
        </dl>
      </div>
      <div>
        <h3>Configurazione di valenza di ${el.symbol}₂</h3>
        <dl class="info-list">
          <dt>Elettroni di valenza</dt><dd>${valenceE}</dd>
          <dt>Ordine di legame</dt><dd>${nf(diag.bondOrder, 1)}</dd>
          <dt>Elettroni spaiati</dt><dd>${diag.unpaired} · ${diag.unpaired ? 'paramagnetica' : 'diamagnetica'}</dd>
        </dl>
        <p class="desc-muted" style="margin-top:8px">Ordine dei livelli per le molecole del 2° periodo: ${diag.sMixing ? 'fino all\'azoto π(2p) sta sotto σ(2p) per il mescolamento s–p.' : 'da O₂ in poi σ(2p) sta sotto π(2p).'} Per gli altri elementi lo schema è indicativo.</p>
      </div>`;
  }
  const h = HYBRIDS[bond.hybrid];
  const energies = { 0: s.e, 1: p.e, 2: d ? d.e : 0 };
  const eh = h.set.map(c => hybridEnergy(c, energies));
  const composition = (c) => {
    const w = [0, 0, 0];
    for (const [l, , k] of c) w[l] += k * k;
    return w;
  };
  const shown = bond.which === 'all' ? h.set : [h.set[bond.which]];
  const comps = [...new Set(shown.map(c => {
    const w = composition(c);
    return `${Math.round(w[0] * 100)} / ${Math.round(w[1] * 100)}${h.needsD ? ` / ${Math.round(w[2] * 100)}` : ''}%`;
  }))];
  return `
    <div>
      <h3>Ibridi ${h.label} di ${el.symbol}</h3>
      <p class="desc">Geometria ${h.geometry}, angolo tra gli ibridi ${h.angle}. Esempi: ${h.example}.</p>
      <dl class="info-list" style="margin-top:10px">
        <dt>Orbitali atomici usati</dt><dd>${n}s, ${n}p${h.needsD ? `, ${n}d` : ''}</dd>
        <dt>Carattere s / p${h.needsD ? ' / d' : ''}</dt><dd>${comps.join(' · ')}</dd>
        <dt>ε(${n}s)</dt><dd>${eV(s.e)}</dd>
        <dt>ε(${n}p)${p.occupied ? '' : ' (vuoto)'}</dt><dd>${eV(p.e)}</dd>
        ${h.needsD ? `<dt>ε(${n}d) (vuoto)</dt><dd>${eV(d.e)}</dd>` : ''}
        <dt>ε dell'ibrido = Σc²ε</dt><dd>${[...new Set(eh.map(e => eV(e)))].join(' · ')}</dd>
      </dl>
      ${h.needsD ? `<p class="desc" style="margin-top:8px">Nell'atomo libero il ${n}d calcolato ha raggio più probabile ${pm(d.rMaxProb)}, contro ${pm(p.rMaxProb)} del ${n}p: gli ibridi con carattere d risultano molto più estesi.</p>` : ''}
      ${h.needsD ? '<p class="desc-muted" style="margin-top:8px">Gli orbitali d vuoti degli elementi del blocco p sono molto più alti in energia e più diffusi di s e p: per questo i chimici oggi descrivono PCl₅ e SF₆ soprattutto con legami a più centri. Il modello sp³d/sp³d² resta utile per prevedere la geometria (VSEPR).</p>' : ''}
      <div class="phase-key" style="margin-top:10px">${h.set.map((_, i) => `<span><i style="background:${HYBRID_COLORS[i]}"></i>h${i + 1}</span>`).join('')}</div>
    </div>`;
}

// ---------------------------------------------------------------------------
// Scena 3D
// ---------------------------------------------------------------------------

function setTitle(html, sub) {
  $('viewport-title').innerHTML = `${html}${sub ? `<small>${sub}</small>` : ''}`;
}

function setLegend(items) {
  $('viewport-legend').innerHTML = items.map(([color, label]) => `<span><i style="background:${color}"></i>${label}</span>`).join('');
}

function render3D() {
  const { atom, mode } = state;
  if (!atom) return;
  viewer.clear();
  const el = element(atom.Z);
  const note = $('viewport-note');

  if (mode === 'atom') {
    const visible = atom.orbitals.filter(o => !state.hidden.has(o.label));
    const parts = [];
    visible.forEach((o, gi) => {
      if (state.atomView === 'density') {
        parts.push({ radial: o.radial, l: o.l, m: null, weight: o.occ, group: gi });
      } else {
        const ms = mOrder(o.l);
        hundBoxes(o.l, o.occ).forEach((b, i) => {
          const e = b.up + b.down;
          if (e) parts.push({ radial: o.radial, l: o.l, m: ms[i], weight: e, group: gi });
        });
      }
    });
    const extent = Math.max(...visible.map(o => o.radial.radiusEnclosing(0.985)), 0.05) * 1.05;
    viewer.frame(extent, true);
    if (parts.length) {
      const res = sampleAtom(parts, state.points, 11);
      const palette = state.atomView === 'density'
        ? [colorToRGB(cssVar('--scene-mode') === 'light' ? '#3d6ea8' : '#d7e6ff')]
        : visible.map(o => colorToRGB(subshellColor(o.n, o.l)));
      const colors = new Float32Array(res.count * 3);
      for (let i = 0; i < res.count; i++) colors.set(palette[state.atomView === 'density' ? 0 : res.groups[i]], 3 * i);
      viewer.addPoints(res.positions, colors, { size: state.atomView === 'density' ? 1.15 : 0.9, opacity: cssVar('--scene-mode') === 'light' ? 0.45 : 0.62 });
    }
    viewer.addNucleus();
    setTitle(`${el.name}`, state.atomView === 'orbitals' ? 'scomposizione in orbitali occupati, non l’aspetto' : 'densità di probabilità elettronica');
    setLegend(state.atomView === 'density'
      ? [[cssVar('--scene-mode') === 'light' ? '#3d6ea8' : '#d7e6ff', 'ρ = Σ nᵢ|ψᵢ|²']]
      : visible.map(o => [subshellColor(o.n, o.l), `${o.label}${superscript(o.occ)}`]));
    note.textContent = state.atomView === 'density'
      ? `Ogni punto è campionato con la regola di Born, probabilità ∝ |ψ|² degli orbitali Kohn–Sham occupati. Il nucleo è disegnato circa 100 000 volte più grande del reale: il raggio nucleare è dell’ordine dei fm, la nube degli Å. Tacche ogni ${formatPm(viewer.tickPm)}.`
      : `Scomposizione negli orbitali occupati, colori per sottolivello: non è l’aspetto dell’atomo. La densità totale è la somma. Nucleo non in scala. Tacche ogni ${formatPm(viewer.tickPm)}.`;
    return;
  }

  if (mode === 'orbital') {
    const { n, l, m } = state.orbital;
    const orb = getOrbital(atom, n, l);
    const terms = [term(orb.radial, l, m)];
    const extent = orb.radial.radiusEnclosing(0.995) * 1.02;
    viewer.frame(extent, true);
    drawTerms(terms, extent, { phase: true });
    viewer.addNucleus();
    setTitle(`${orbitalHTML(n, l, m)} · ${el.name}`, orb.occupied ? `ε = ${eV(orb.e)}` : `stato non occupato, ε = ${eV(orb.e)}`);
    setLegend([[cssVar('--phase-pos'), 'ψ > 0'], [cssVar('--phase-neg'), 'ψ < 0']]);
    note.textContent = `Isosuperficie di ψ e nuvola campionata con |ψ|². Il segno di ψ è una fase matematica, non un colore osservato. La superficie racchiude il ${Math.round(state.enclosed * 100)}% della probabilità. Sezione per i nodi. Tacche ogni ${formatPm(viewer.tickPm)}.`;
    return;
  }

  // Legami
  const { bond } = state;
  const { n, s, p, d } = bondRadials(atom);
  if (bond.kind === 'mo') {
    const type = MO_TYPES.find(t => t.id === bond.mo);
    const dist = (bond.distancePm ?? defaultBondPm(atom)) / BOHR_PM;
    const terms = moTerms(type, s.radial, p.radial, dist);
    const radial = type.ao === 's' ? s.radial : p.radial;
    const extent = dist / 2 + radial.radiusEnclosing(0.99);
    viewer.frame(extent, true);
    drawTerms(terms, extent, { phase: true });
    viewer.addNucleus([0, 0, -dist / 2], el.symbol);
    viewer.addNucleus([0, 0, dist / 2], el.symbol);
    setTitle(`${type.label}<sub>${n}${type.sub}</sub> · ${el.symbol}₂`, `asse di legame z, d = ${Math.round(dist * BOHR_PM)} pm`);
    setLegend([[cssVar('--phase-pos'), 'ψ > 0'], [cssVar('--phase-neg'), 'ψ < 0']]);
    note.textContent = `Legame come combinazione lineare degli orbitali atomici calcolati (LCAO). L’addensamento fra i nuclei è la densità di probabilità, non un bastoncino. Tacche ogni ${formatPm(viewer.tickPm)}.`;
    return;
  }
  const h = HYBRIDS[bond.hybrid];
  const radials = { 0: s.radial, 1: p.radial, 2: d ? d.radial : null };
  const indices = bond.which === 'all' ? h.set.map((_, i) => i) : [bond.which];
  const extent = Math.max(s.radial.radiusEnclosing(0.99), p.radial.radiusEnclosing(0.99), h.needsD ? d.radial.radiusEnclosing(0.95) : 0);
  viewer.frame(extent, true);
  // Con tutti gli ibridi insieme una superficie al 90% farebbe sovrapporre i lobi: si usa al massimo il 55%.
  const enclosed = indices.length > 1 ? Math.min(state.enclosed, SET_ENCLOSED) : state.enclosed;
  indices.forEach((i) => {
    const terms = hybridTerms(h.set[i], radials);
    drawTerms(terms, extent, {
      phase: indices.length === 1, color: HYBRID_COLORS[i], seed: 100 + i,
      points: Math.round(state.points / indices.length),
      quality: indices.length > 1 ? Math.min(state.quality, 64) : state.quality,
      enclosed,
    });
  });
  viewer.addNucleus([0, 0, 0], el.symbol);
  setTitle(`Ibridi ${h.label} · ${el.name}`, `geometria ${h.geometry}`);
  setLegend(indices.length === 1
    ? [[cssVar('--phase-pos'), `h${indices[0] + 1}, ψ > 0`], [cssVar('--phase-neg'), 'ψ < 0']]
    : indices.map(i => [HYBRID_COLORS[i], `h${i + 1}`]));
  note.textContent = `Ogni ibrido ha un lobo grande (ψ > 0) e un piccolo lobo di segno opposto${indices.length > 1 ? ' (più scuro)' : ''}.${indices.length > 1 && state.render !== 'cloud' && enclosed < state.enclosed ? ` Con tutti gli ibridi la superficie racchiude il ${Math.round(enclosed * 100)}% della probabilità, per non sovrapporre i lobi.` : ''} Tacche degli assi ogni ${formatPm(viewer.tickPm)}.`;
}

/** Disegna una funzione d'onda (somma di termini) come nuvola e/o superficie. */
function drawTerms(terms, extent, { phase, color, seed = 1, points = state.points, quality = state.quality, enclosed = state.enclosed }) {
  const pos = cssVar('--phase-pos');
  const neg = cssVar('--phase-neg');
  const light = cssVar('--scene-mode') === 'light';
  const posColor = phase ? pos : color;
  const negColor = phase ? neg : shade(color, light ? 0.45 : 0.32);
  if (state.render !== 'surface') {
    const pts = samplePoints(terms, points, seed);
    const cp = colorToRGB(posColor);
    const cn = colorToRGB(negColor);
    const colors = new Float32Array(pts.count * 3);
    for (let i = 0; i < pts.count; i++) colors.set(pts.signs[i] > 0 ? cp : cn, 3 * i);
    const both = state.render === 'both';
    viewer.addPoints(pts.positions, colors, { size: both ? 0.7 : 1, opacity: both ? 0.35 : (light ? 0.6 : 0.55) });
  }
  if (state.render !== 'cloud') {
    const g = sampleGrid(terms, quality, extent, enclosed);
    const opacity = state.render === 'both' ? 0.6 : 0.9;
    viewer.addSurface(g, g.iso, 1, posColor, { opacity });
    viewer.addSurface(g, g.iso, -1, negColor, { opacity });
  }
}

function shade(color, k) {
  const [r, g, b] = colorToRGB(color);
  const to = (v) => Math.round(v * k * 255).toString(16).padStart(2, '0');
  return `#${to(r)}${to(g)}${to(b)}`;
}

// ---------------------------------------------------------------------------
// Grafici
// ---------------------------------------------------------------------------

function drawCharts() {
  if (!state.atom) return;
  drawRadial();
  drawLevelChart();
  drawSliceChart();
}

function radialSeries(orb, weight, color, label, fill) {
  const { r, u } = { r: orb.radial.r, u: orb.u };
  const xs = [];
  const ys = [];
  const step = Math.max(1, Math.floor(r.length / 700));
  for (let i = 0; i < r.length; i += step) {
    xs.push(r[i] * BOHR_PM);
    ys.push(weight * u[i] * u[i] / BOHR_PM);
  }
  return { xs, ys, color, label, fill };
}

function drawRadial() {
  const { atom, mode } = state;
  const canvas = $('chart-radial');
  const accent = cssVar('--accent');
  const muted = cssVar('--muted');
  const orbitals = atom.orbitals.filter(o => !state.hidden.has(o.label));
  let log = state.radialLog;

  if (mode === 'atom') {
    const cutoff = Math.max(...orbitals.map(o => o.radial.radiusEnclosing(0.999)), 0.05) * BOHR_PM;
    const trim = (s) => {
      const keep = s.xs.map((x, i) => i).filter(i => s.xs[i] <= cutoff);
      return { ...s, xs: keep.map(i => s.xs[i]), ys: keep.map(i => s.ys[i]) };
    };
    // In scala logaritmica si usa r·D(r): elettroni per unità di ln r, così l'area sotto ogni picco resta proporzionale agli elettroni.
    const perLn = (sr) => (log ? { ...sr, ys: sr.ys.map((y, i) => y * sr.xs[i]) } : sr);
    const series = orbitals.map(o => perLn(trim(radialSeries(o, o.occ, subshellColor(o.n, o.l), o.label))));
    // densità radiale totale 4πr²ρ
    const total = trim(radialSeries(orbitals[0] ?? atom.orbitals[0], 0, accent, 'totale', `color-mix(in srgb, ${accent} 14%, transparent)`));
    total.ys = total.ys.map(() => 0);
    for (const s of series) s.ys.forEach((v, i) => { total.ys[i] += v; });
    total.width = 2.2;
    drawLineChart(canvas, { series: [total, ...series], log, xlabel: 'r (pm)', ylabel: log ? 'elettroni per unità di ln r' : 'elettroni per pm' });
    $('radial-title').textContent = 'Distribuzione radiale degli elettroni';
    $('radial-note').innerHTML = log
      ? 'Area gialla: r·4πr²ρ(r), elettroni per unità di ln r. Ogni picco è un guscio (K, L, M, …); l\'area sotto ciascun picco è proporzionale ai suoi elettroni.'
      : `Area gialla: numero totale di elettroni per unità di raggio, 4πr²ρ(r). I massimi corrispondono ai gusci. ${atom.Z > 10 ? 'Con la scala log i gusci interni si vedono meglio.' : ''}`;
    return;
  }
  if (mode === 'orbital') {
    const { n, l } = state.orbital;
    const orb = getOrbital(atom, n, l);
    const P = radialSeries(orb, 1, accent, 'P(r) = r²R²', `color-mix(in srgb, ${accent} 16%, transparent)`);
    // R(r) con segno, riscalata per stare nello stesso grafico
    const Rxs = [];
    const Rys = [];
    const step = Math.max(1, Math.floor(orb.radial.r.length / 700));
    for (let i = 0; i < orb.radial.r.length; i += step) {
      Rxs.push(orb.radial.r[i] * BOHR_PM);
      Rys.push(orb.radial.R[i]);
    }
    const pmax = Math.max(...P.ys);
    const cutoff = orb.radial.radiusEnclosing(0.9995) * BOHR_PM;
    const inRange = Rys.filter((_, i) => Rxs[i] > cutoff * 0.02 && Rxs[i] < cutoff);
    const rmaxAbs = Math.max(...inRange.map(Math.abs), 1e-12);
    const R = { xs: Rxs, ys: Rys.map(v => Math.max(-1.5, Math.min(1.5, v / rmaxAbs)) * pmax), color: muted, dash: [5, 4], label: 'R(r) (scala arbitraria)', width: 1.3 };
    const trimmed = (s) => {
      const keep = s.xs.map((x, i) => i).filter(i => s.xs[i] <= cutoff);
      return { ...s, xs: keep.map(i => s.xs[i]), ys: keep.map(i => s.ys[i]) };
    };
    const markers = [
      ...orb.nodes.map(r => ({ x: r * BOHR_PM, label: 'nodo', color: cssVar('--phase-neg') })),
      { x: orb.rMaxProb * BOHR_PM, label: 'r max', color: accent, dash: [], row: 1 },
    ];
    drawLineChart(canvas, { series: [trimmed(P), trimmed(R)], log, xlabel: 'r (pm)', ylabel: 'probabilità per pm', markers, signed: true });
    $('radial-title').textContent = `Funzione radiale ${n}${L_LETTERS[l]}`;
    $('radial-note').textContent = `${n - l - 1} nodi radiali: sono i punti in cui R(r) cambia segno e la probabilità si annulla. Raggio più probabile: ${pm(orb.rMaxProb)}.`;
    return;
  }
  const { n, s, p, d } = bondRadials(atom);
  const h = HYBRIDS[state.bond.hybrid];
  const series = [
    radialSeries(s, 1, cssVar('--block-s'), `${n}s`),
    radialSeries(p, 1, cssVar('--block-p'), `${n}p`),
  ];
  if (state.bond.kind === 'hybrid' && h.needsD && d) series.push(radialSeries(d, 1, cssVar('--block-d'), `${n}d`));
  const cutoff = Math.max(...series.map((_, i) => [s, p, d][i]?.radial.radiusEnclosing(0.999) ?? 0)) * BOHR_PM;
  const markers = state.bond.kind === 'mo' ? [{ x: (state.bond.distancePm ?? defaultBondPm(atom)) / 2, label: 'd/2', color: accent, dash: [] }] : [];
  drawLineChart(canvas, {
    series: series.map(sr => {
      const keep = sr.xs.map((x, i) => i).filter(i => sr.xs[i] <= cutoff);
      return { ...sr, xs: keep.map(i => sr.xs[i]), ys: keep.map(i => sr.ys[i]) };
    }),
    log, xlabel: 'r (pm)', ylabel: 'probabilità per pm', markers,
  });
  $('radial-title').textContent = 'Orbitali atomici di valenza';
  $('radial-note').textContent = state.bond.kind === 'mo'
    ? 'd/2 è metà della distanza di legame: gli orbitali si sovrappongono dove la loro probabilità radiale è ancora grande.'
    : 'Per gli ibridi si mescolano orbitali di estensione diversa: per questo i lobi hanno forme asimmetriche.';
}

function drawLevelChart() {
  const { atom, mode } = state;
  const canvas = $('chart-levels');
  if (mode === 'bond' && state.bond.kind === 'mo') {
    const { n, s, p } = bondRadials(atom);
    const valenceE = atom.config.filter(c => c.n === n && c.l <= 1).reduce((a, c) => a + c.occ, 0) * 2;
    const diag = valenceMODiagram(atom.Z, valenceE);
    drawMODiagram(canvas, { eS: s.e * HARTREE_EV, eP: p.e * HARTREE_EV, levels: diag.levels, selected: state.bond.mo, symbol: element(atom.Z).symbol });
    $('levels-title').textContent = `Diagramma degli orbitali molecolari`;
    $('levels-note').textContent = `Livelli atomici ${n}s e ${n}p calcolati; posizioni degli orbitali molecolari qualitative. Ordine di legame ${nf(diag.bondOrder, 1)}.`;
    return;
  }
  if (mode === 'bond') {
    const { n, s, p, d } = bondRadials(atom);
    const h = HYBRIDS[state.bond.hybrid];
    const energies = { 0: s.e, 1: p.e, 2: d ? d.e : 0 };
    const levels = [
      { label: `${n}s`, e: s.e * HARTREE_EV, occ: 0, color: cssVar('--block-s') },
      { label: `${n}p`, e: p.e * HARTREE_EV, occ: 0, color: cssVar('--block-p') },
    ];
    if (h.needsD) levels.push({ label: `${n}d`, e: d.e * HARTREE_EV, occ: 0, color: cssVar('--block-d') });
    const eh = [...new Set(h.set.map(c => hybridEnergy(c, energies).toFixed(8)))].map(Number);
    eh.forEach((e, i) => levels.push({ label: `${h.label}${eh.length > 1 ? ` (${i + 1})` : ''}`, e: e * HARTREE_EV, occ: 0, highlight: true }));
    drawLevels(canvas, levels);
    $('levels-title').textContent = 'Energia degli ibridi';
    $('levels-note').textContent = 'L\'energia di un ibrido è la media delle energie degli orbitali di partenza, pesata con i quadrati dei coefficienti.';
    return;
  }
  const levels = atom.orbitals.map(o => ({
    label: o.label, e: o.e * HARTREE_EV, occ: o.occ, cap: 2 * (2 * o.l + 1),
    color: subshellColor(o.n, o.l),
    highlight: mode === 'orbital' && o.n === state.orbital.n && o.l === state.orbital.l,
  }));
  if (mode === 'orbital') {
    const { n, l } = state.orbital;
    const orb = getOrbital(atom, n, l);
    if (!orb.occupied) levels.push({ label: orb.label, e: orb.e * HARTREE_EV, occ: 0, cap: 2 * (2 * l + 1), highlight: true });
  }
  drawLevels(canvas, levels);
  $('levels-title').textContent = 'Livelli energetici calcolati';
  $('levels-note').textContent = 'Energie degli orbitali (eV) in scala logaritmica. La barra colorata indica il riempimento del sottolivello.';
}

const PLANES = {
  xz: { axes: ['x', 'z'], map: (u, v) => [u, 0, v] },
  xy: { axes: ['x', 'y'], map: (u, v) => [u, v, 0] },
  yz: { axes: ['y', 'z'], map: (u, v) => [0, u, v] },
};

function bestPlane(fn, half) {
  let best = 'xz';
  let bestV = -1;
  for (const [k, P] of Object.entries(PLANES)) {
    let s = 0;
    for (let j = 0; j < 24; j++) {
      for (let i = 0; i < 24; i++) {
        const [x, y, z] = P.map(-half + (i + 0.5) * half / 12, -half + (j + 0.5) * half / 12);
        s += Math.abs(fn(x, y, z));
      }
    }
    if (s > bestV * 1.0001) { bestV = s; best = k; }
  }
  return best;
}

function drawSliceChart() {
  const { atom, mode } = state;
  const canvas = $('chart-slice');
  const N = 180;
  let fn;
  let signed = true;
  let half;
  let nuclei = [[0, 0]];
  let preferred = 'xz';
  let composite = false;

  if (mode === 'atom') {
    const visible = atom.orbitals.filter(o => !state.hidden.has(o.label));
    half = Math.max(...visible.map(o => o.radial.radiusEnclosing(0.98)), 0.05);
    signed = false;
    if (state.atomView === 'density') {
      fn = (x, y, z) => {
        const r = Math.sqrt(x * x + y * y + z * z);
        let rho = 0;
        for (const o of visible) { const R = o.radial.at(r); rho += o.occ * R * R; }
        return rho;
      };
    } else {
      const terms = [];
      for (const o of visible) {
        const ms = mOrder(o.l);
        hundBoxes(o.l, o.occ).forEach((b, i) => {
          const e = b.up + b.down;
          if (e) terms.push({ t: term(o.radial, o.l, ms[i]), e });
        });
      }
      fn = (x, y, z) => terms.reduce((s, { t, e }) => { const v = evaluate([t], x, y, z); return s + e * v * v; }, 0);
    }
  } else if (mode === 'orbital') {
    const { n, l, m } = state.orbital;
    const orb = getOrbital(atom, n, l);
    const terms = [term(orb.radial, l, m)];
    half = orb.radial.radiusEnclosing(0.99);
    fn = (x, y, z) => evaluate(terms, x, y, z);
  } else if (state.bond.kind === 'mo') {
    const { s, p } = bondRadials(atom);
    const type = MO_TYPES.find(t => t.id === state.bond.mo);
    const dist = (state.bond.distancePm ?? defaultBondPm(atom)) / BOHR_PM;
    const terms = moTerms(type, s.radial, p.radial, dist);
    half = dist / 2 + (type.ao === 's' ? s : p).radial.radiusEnclosing(0.98);
    fn = (x, y, z) => evaluate(terms, x, y, z);
    nuclei = [[0, -dist / 2 * BOHR_PM], [0, dist / 2 * BOHR_PM]];
    preferred = 'xz';
  } else {
    const { s, p, d } = bondRadials(atom);
    const h = HYBRIDS[state.bond.hybrid];
    const radials = { 0: s.radial, 1: p.radial, 2: d ? d.radial : null };
    const which = state.bond.which === 'all' ? h.set.map((_, i) => i) : [state.bond.which];
    const sets = which.map(i => hybridTerms(h.set[i], radials));
    half = Math.max(s.radial.radiusEnclosing(0.98), p.radial.radiusEnclosing(0.98));
    if (sets.length === 1) {
      fn = (x, y, z) => evaluate(sets[0], x, y, z);
    } else {
      // con più ibridi si mostra, punto per punto, l'ibrido con |ψ| maggiore
      fn = (x, y, z) => {
        let best = 0;
        for (const t of sets) { const v = evaluate(t, x, y, z); if (Math.abs(v) > Math.abs(best)) best = v; }
        return best;
      };
      composite = true;
    }
    preferred = state.bond.hybrid === 'sp2' ? 'xy' : 'xz';
  }

  let plane = state.slicePlane;
  if (plane === 'auto') plane = (mode === 'orbital' || (mode === 'bond' && state.bond.kind === 'hybrid')) ? bestPlane(fn, half) : preferred;
  const P = PLANES[plane];
  if (mode === 'bond' && state.bond.kind === 'mo') {
    // i nuclei stanno sull'asse z: visibili nei piani xz e yz
    nuclei = plane === 'xy' ? [[0, 0]] : nuclei;
  }
  const field = new Float32Array(N * N);
  for (let j = 0; j < N; j++) {
    const v = half - (j + 0.5) * 2 * half / N;
    for (let i = 0; i < N; i++) {
      const u = -half + (i + 0.5) * 2 * half / N;
      const [x, y, z] = P.map(u, v);
      field[j * N + i] = fn(x, y, z);
    }
  }
  drawSlice(canvas, { field, n: N, signed, extentPm: half * BOHR_PM, axes: P.axes, nuclei });
  $('slice-title').textContent = signed ? 'Sezione di ψ' : 'Sezione della densità';
  $('slice-note').textContent = composite
    ? 'In ogni punto è mostrato l\'ibrido con |ψ| maggiore: si vedono i lobi grandi (ψ > 0) e i piccoli lobi opposti.'
    : signed
    ? 'Colore: segno e ampiezza di ψ nel piano. Linee chiare: nodi (ψ = 0).'
    : 'Densità di probabilità nel piano, in scala logaritmica: si vedono i gusci.';
}

// ---------------------------------------------------------------------------
// Avvio
// ---------------------------------------------------------------------------

viewer.setAutoRotate(true);
initMoleculeMode(viewer);
initReactionMode(viewer);
initFormationMode();
initSandboxMode(viewer, {
  openMolecule: (smiles, name) => { setMode('molecule'); openMoleculeSmiles(smiles, name); },
});
document.body.dataset.mode = 'atom';
{
  const h = location.hash.slice(1);
  const startMode = { molecole: 'molecule', reazioni: 'reaction', laboratorio: 'lab', sandbox: 'sandbox', becher: 'beaker', fasi: 'phase', cinetica: 'kinetics', legame: 'formation' }[h] ?? null;
  if (startMode) state.mode = startMode;
  selectElement(zFromHash() ?? 6);
  if (startMode) { state.mode = 'atom'; setMode(startMode); }
}
