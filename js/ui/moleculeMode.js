// Modalità "Molecole": costruzione, calcolo quantistico e analisi chimica di una molecola.

import { MOLECULES, CATEGORIES, moleculeById } from '../chem/library.js';
import { LIBRARY_DATA } from '../chem/libraryData.js';
import * as THREE from 'three';
import { parseSmiles, hillFormula, SmilesError, writeSmiles, fillHydrogens } from '../chem/smiles.js';
import { analyzeStructure } from '../chem/structure.js';
import { embedMolecule } from '../chem/embed.js';
import { buildBasis } from '../chem/integrals.js';
import { thermochemistry, ISOTOPE_MASS } from '../chem/vibrations.js';
import { PAULING, ALLEN } from '../chem/elementData.js';
import { angle as angleOf } from '../chem/optimize.js';
import { request, cancelAll } from '../chem/moleculeClient.js';
import { marchingCubes } from '../render/marching.js';
import { buildMolecule, dipoleArrow, espColor, chargeColor } from '../render/moleculeView.js';
import { formatPm } from '../render/viewer.js';
import { drawLewis, drawMOLevels, drawIR } from './chemCharts.js';
import { drawLineChart } from './charts.js';
import { ELEMENTS } from '../physics/elements.js';

const BOHR_ANG = 0.52917721090;
const HARTREE_EV = 27.211386245988;
const HARTREE_KJ = 2625.4996394799;
// fattori di scala per le frequenze HF (Scott e Radom, J. Phys. Chem. 1996)
export const FREQ_SCALE = { 'STO-3G': 0.8929, '3-21G': 0.9085, '6-31G': 0.9054, '6-31G*': 0.8953, '6-31G**': 0.8992 };

const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const sym = (Z) => ELEMENTS[Z - 1].symbol;
const nf = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d }));
const sgn = (v, d = 2) => (v > 0 ? '+' : v < 0 ? '−' : '') + nf(Math.abs(v), d);

let viewer = null;
const M = {
  entry: null, smiles: '', name: '', charge: 0, multiplicity: 1,
  graph: null, analysis: null, atoms: null, basis: 'STO-3G', basisAuto: true,
  summary: null, freq: null, lib: null, scan: null,
  view: { overlay: 'model', mo: null, labels: 'formal', dipole: true, iso: 0.05, vib: null, chart2: 'ir' },
  status: '', busy: false, token: 0,
  cat: 'hydride',
  moRegions: [],
  sel: [],
  mol3d: null,
};

export function initMoleculeMode(v) {
  viewer = v;
  // selezione degli atomi con un clic (senza trascinamento) nel modello 3D
  const canvas = viewer.renderer.domElement;
  let down = null;
  canvas.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
  canvas.addEventListener('pointerup', (e) => {
    if (!active || !down || !M.mol3d) return;
    const moved = Math.hypot(e.clientX - down[0], e.clientY - down[1]);
    down = null;
    if (moved > 5) return;
    const rect = canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
    const ray = new THREE.Raycaster();
    ray.setFromCamera(ndc, viewer.camera);
    const hits = ray.intersectObjects(M.mol3d.atomMeshes, false);
    if (!hits.length) { M.sel = []; } else {
      const idx = M.mol3d.atomMeshes.indexOf(hits[0].object);
      if (M.sel.includes(idx)) M.sel = M.sel.filter(i => i !== idx);
      else M.sel = e.shiftKey || M.sel.length === 1 ? [...M.sel.slice(-1), idx] : [idx];
    }
    renderEditor();
    highlightSelection();
  });
  $('chart-radial').addEventListener('click', (ev) => {
    if (!isActive()) return;
    const rect = ev.currentTarget.getBoundingClientRect();
    const x = ev.clientX - rect.left, y = ev.clientY - rect.top;
    const hit = M.moRegions.find(r => x >= r.x0 && x <= r.x1 && y >= r.y0 && y <= r.y1);
    if (hit) selectMO(hit.index, hit.spin);
  });
}

let active = false;
function isActive() { return active; }

export function activateMolecule() {
  active = true;
  document.body.dataset.mode = 'molecule';
  if (!M.graph) loadEntry(moleculeById('H2O'));
  else renderAll();
}

export function deactivateMolecule() {
  active = false;
  stopVibration();
}

export function moleculeThemeChanged() {
  if (active && M.graph) { render3D(); drawCharts(); }
}

export function moleculeResize() {
  if (active && M.graph) drawCharts();
}

// ---------------------------------------------------------------------------
// Caricamento di una molecola
// ---------------------------------------------------------------------------

function autoBasis(atoms) {
  for (const b of ['6-31G*', '3-21G', 'STO-3G']) {
    try {
      if (buildBasis(atoms, b).nbf <= (b === '6-31G*' ? 72 : 90)) return b;
    } catch { /* elemento non disponibile in questa base */ }
  }
  return 'STO-3G';
}

function centerAtoms(atoms) {
  const m = atoms.map(a => ISOTOPE_MASS[a.Z] ?? 2 * a.Z);
  const M0 = m.reduce((s, x) => s + x, 0);
  const c = [0, 1, 2].map(k => atoms.reduce((s, a, i) => s + m[i] * a.xyz[k], 0) / M0);
  return atoms.map(a => ({ Z: a.Z, xyz: a.xyz.map((v, k) => v - c[k]) }));
}

function loadEntry(entry) {
  M.entry = entry;
  M.name = entry.name;
  M.cat = entry.cat;
  loadSmiles(entry.smiles, { multiplicity: entry.multiplicity, lib: LIBRARY_DATA[entry.id] ?? null });
}

function loadSmiles(smiles, { multiplicity, lib = null } = {}) {
  cancelAll();
  stopVibration();
  let graph;
  try {
    graph = parseSmiles(smiles);
  } catch (e) {
    M.status = e instanceof SmilesError ? `SMILES non valido: ${e.message}` : String(e.message);
    renderSide();
    return;
  }
  if (graph.atoms.some(a => a.Z > 36)) {
    M.status = 'Le basi gaussiane disponibili arrivano fino al kripton (Z = 36): la molecola contiene elementi più pesanti.';
    renderSide();
    return;
  }
  M.smiles = smiles;
  M.graph = graph;
  M.sel = [];
  M.resetView = true;
  M.charge = graph.atoms.reduce((s, a) => s + a.charge, 0);
  const nel = graph.atoms.reduce((s, a) => s + a.Z, 0) - M.charge;
  M.multiplicity = multiplicity ?? (nel % 2 === 0 ? 1 : 2);
  M.analysis = analyzeStructure({ atoms: graph.atoms, bonds: graph.bonds, charge: M.charge });
  M.lib = lib;
  if (lib && lib.geometry.length === graph.atoms.length) {
    M.atoms = lib.geometry.map(([Z, x, y, z]) => ({ Z, xyz: [x / BOHR_ANG, y / BOHR_ANG, z / BOHR_ANG] }));
  } else {
    M.atoms = embedMolecule(graph, M.analysis.atoms.map(a => a.lonePairs + a.radical));
  }
  M.atoms = centerAtoms(M.atoms);
  if (M.basisAuto) M.basis = lib?.basis ?? autoBasis(M.atoms);
  M.freq = lib && lib.freqs?.length && lib.basis === M.basis
    ? { modes: lib.freqs.map(f => ({ freq: f.freq, ir: f.ir, mu: f.mu, d: f.d })), source: 'libreria' }
    : null;
  M.summary = null;
  M.scan = null;
  M.view.mo = null;
  M.view.vib = null;
  if (M.view.overlay !== 'model' && M.view.overlay !== 'charges') M.view.overlay = 'model';
  M.status = lib ? 'Geometria ottimizzata HF (libreria precalcolata).' : 'Geometria iniziale dal modello VSEPR: puoi ottimizzarla con Hartree–Fock.';
  renderAll();
  runSCF();
}

async function runSCF() {
  const token = ++M.token;
  setBusy(`Calcolo Hartree–Fock ${M.multiplicity > 1 ? '(UHF) ' : ''}in base ${M.basis}…`);
  try {
    const summary = await request('scf', { atoms: M.atoms, opts: { basis: M.basis, charge: M.charge, multiplicity: M.multiplicity } });
    if (token !== M.token) return;
    M.summary = summary;
    // orbitale predefinito: HOMO
    M.view.mo = { index: summary.nalpha - 1, spin: summary.unrestricted ? 'α' : null };
    setBusy(null);
    renderAll();
  } catch (e) {
    if (token !== M.token || e.message === 'annullato') return;
    setBusy(null);
    M.status = `Calcolo non riuscito: ${e.message}`;
    renderSide();
  }
}

async function optimize() {
  const token = ++M.token;
  stopVibration();
  setBusy('Ottimizzazione della geometria…');
  try {
    const r = await request('optimize', { atoms: M.atoms, opts: { basis: M.basis, charge: M.charge, multiplicity: M.multiplicity } }, (p, info) => {
      if (token !== M.token || !info?.atoms) return;
      M.atoms = info.atoms;
      setBusy(`Ottimizzazione: passo ${info.step}, E = ${nf(info.energy, 6)} Ha, gradiente max ${info.gmax.toExponential(1)}`);
      render3D();
    });
    if (token !== M.token) return;
    M.atoms = centerAtoms(r.atoms);
    M.summary = r.summary;
    M.freq = null;
    M.view.mo = { index: r.summary.nalpha - 1, spin: r.summary.unrestricted ? 'α' : null };
    M.status = `Geometria ottimizzata in ${r.history.length - 1} passi (HF/${M.basis})${r.converged ? '' : ', convergenza non raggiunta'}.`;
    setBusy(null);
    // ricalcola con la geometria centrata per le griglie
    renderAll();
    runSCF();
  } catch (e) {
    if (token !== M.token || e.message === 'annullato') return;
    setBusy(null);
    M.status = `Ottimizzazione non riuscita: ${e.message}`;
    renderSide();
  }
}

async function computeFrequencies() {
  const token = ++M.token;
  stopVibration();
  setBusy('Calcolo della matrice hessiana (derivate seconde dell\'energia)…');
  try {
    const r = await request('frequencies', { atoms: M.atoms, opts: { basis: M.basis, charge: M.charge, multiplicity: M.multiplicity } }, (p) => {
      if (token === M.token && p !== null) setBusy(`Hessiana: ${Math.round(p * 100)}% (${M.atoms.length * 6} gradienti analitici)`);
    });
    if (token !== M.token) return;
    M.freq = { modes: r.modes, source: 'calcolo' };
    setBusy(null);
    renderAll();
    // riporta il worker sulla funzione d'onda della geometria corrente
    runSCF();
  } catch (e) {
    if (token !== M.token || e.message === 'annullato') return;
    setBusy(null);
    M.status = `Frequenze non riuscite: ${e.message}`;
    renderSide();
  }
}

async function scanBond() {
  if (M.atoms.length !== 2) return;
  const token = ++M.token;
  stopVibration();
  const [a, b] = M.atoms;
  const re = Math.hypot(a.xyz[0] - b.xyz[0], a.xyz[1] - b.xyz[1], a.xyz[2] - b.xyz[2]);
  const distances = [];
  for (let f = 0.6; f <= 3.01; f += f < 1.6 ? 0.05 : 0.2) distances.push(re * f);
  setBusy('Curva di energia potenziale…');
  try {
    const r = await request('scan', { Z1: a.Z, Z2: b.Z, distances, opts: { basis: M.basis, charge: M.charge, multiplicity: M.multiplicity } }, (p) => {
      if (token === M.token && p !== null) setBusy(`Curva di energia potenziale: ${Math.round(p * 100)}%`);
    });
    if (token !== M.token) return;
    M.scan = r.points;
    M.view.chart2 = 'scan';
    setBusy(null);
    renderAll();
    runSCF();
  } catch (e) {
    if (token !== M.token || e.message === 'annullato') return;
    setBusy(null);
    M.status = `Scansione non riuscita: ${e.message}`;
    renderSide();
  }
}

function setBusy(text) {
  M.busy = !!text;
  $('busy').hidden = !text;
  if (text) $('busy-text').textContent = text;
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

function renderAll() {
  if (!active) return;
  renderSide();
  renderControls();
  renderPanelBody();
  render3D();
  drawCharts();
  renderAnalysis();
}

function entryFormula(m) {
  if (m.formula) return m.formula;
  const g = parseSmiles(m.smiles);
  return hillFormula(g.atoms, g.atoms.reduce((s, a) => s + a.charge, 0));
}

function formulaOf() {
  if (M.entry?.formula && M.entry.smiles === M.smiles) return M.entry.formula;
  return hillFormula(M.graph.atoms, M.charge);
}

function renderSide() {
  if (!active) return;
  const S = M.summary;
  const exp = M.entry?.smiles === M.smiles ? M.entry.exp ?? {} : {};
  const catButtons = CATEGORIES.map(c => `<option value="${c.id}" ${c.id === M.cat ? 'selected' : ''}>${c.label}</option>`).join('');
  const list = MOLECULES.filter(m => m.cat === M.cat).map(m =>
    `<button type="button" class="mol-chip ${M.entry?.id === m.id && M.entry.smiles === M.smiles ? 'active' : ''}" data-id="${m.id}" title="${m.name}">${entryFormula(m)}</button>`).join('');
  const is2 = M.atoms?.length === 2;
  const mult = { 1: 'singoletto', 2: 'doppietto', 3: 'tripletto', 4: 'quartetto' }[M.multiplicity] ?? `${M.multiplicity}`;
  const muExp = exp.dipole;
  $('element-card').innerHTML = `
    <div class="mol-picker">
      <label class="lbl" for="mol-cat">Libreria</label>
      <select id="mol-cat">${catButtons}</select>
      <div class="mol-list">${list}</div>
      <label class="lbl" for="mol-smiles">Oppure scrivi uno SMILES</label>
      <form class="smiles-row" id="smiles-form">
        <input id="mol-smiles" type="text" spellcheck="false" autocomplete="off" value="${escapeHtml(M.smiles)}" placeholder="es. CC(=O)O">
        <button type="submit" class="btn">Costruisci</button>
      </form>
      <p class="hint">Esempi: <code>CCO</code> etanolo, <code>C=CC=C</code> butadiene, <code>c1ccccc1O</code> fenolo, <code>[NH4+]</code>, <code>O=C=O</code>.</p>
    </div>
    <div class="editor" id="mol-editor"></div>
    <div class="mol-head">
      <p class="mol-formula">${formulaOf()}</p>
      <p class="el-name">${escapeHtml(M.entry?.smiles === M.smiles ? M.name : 'Molecola personalizzata')}</p>
      <p class="el-cat">${M.analysis.valenceElectrons} elettroni di valenza · ${mult}${M.charge ? ` · carica ${sgn(M.charge, 0)}` : ''}</p>
    </div>
    ${M.entry?.note && M.entry.smiles === M.smiles ? `<p class="mol-note">${M.entry.note}</p>` : ''}
    <div class="calc-opts">
      <label class="lbl" for="mol-basis">Base gaussiana</label>
      <select id="mol-basis">
        ${['STO-3G', '3-21G', '6-31G', '6-31G*', '6-31G**'].map(b => `<option ${b === M.basis ? 'selected' : ''}>${b}</option>`).join('')}
      </select>
      <div class="btn-row">
        <button type="button" class="btn" id="btn-opt" ${M.busy ? 'disabled' : ''}>Ottimizza geometria</button>
        <button type="button" class="btn" id="btn-freq" ${M.busy ? 'disabled' : ''}>Calcola vibrazioni</button>
        ${is2 ? `<button type="button" class="btn" id="btn-scan" ${M.busy ? 'disabled' : ''}>Curva di dissociazione</button>` : ''}
      </div>
      <p class="hint">${escapeHtml(M.status)}</p>
    </div>
    ${S ? `
    <dl class="facts">
      <dt>Energia HF</dt><dd>${nf(S.energy, 6)} Ha</dd>
      ${S.mp2 !== null ? `<dt>Correlazione MP2</dt><dd>${nf(S.mp2, 6)} Ha</dd><dt>Energia MP2</dt><dd>${nf(S.energy + S.mp2, 6)} Ha</dd>` : ''}
      <dt>Gruppo puntuale</dt><dd>${S.pointGroup}</dd>
      <dt>Momento di dipolo</dt><dd>${nf(S.dipole.debye, 2)} D${muExp !== undefined ? ` <span class="exp">(sper. ${nf(muExp, 2)})</span>` : ''}</dd>
      <dt>HOMO / LUMO</dt><dd>${nf(S.homo * HARTREE_EV, 2)} / ${nf(S.lumo * HARTREE_EV, 2)} eV</dd>
      <dt>Gap HOMO–LUMO</dt><dd>${nf(S.gapEV, 2)} eV</dd>
      <dt>Ionizzazione (Koopmans)</dt><dd>${nf(S.koopmansIE, 2)} eV</dd>
      ${S.unrestricted ? `<dt>⟨S²⟩ (atteso ${nf((M.multiplicity - 1) / 2 * ((M.multiplicity - 1) / 2 + 1), 2)})</dt><dd>${nf(S.S2, 3)}</dd>` : ''}
      <dt>Funzioni di base</dt><dd>${S.nbf} (${S.basis})</dd>
      <dt>Iterazioni SCF</dt><dd>${S.iterations} · ${S.time} ms</dd>
    </dl>` : ''}`;
  $('mol-cat').addEventListener('change', (e) => { M.cat = e.target.value; renderSide(); });
  $('element-card').querySelectorAll('.mol-chip').forEach(b => b.addEventListener('click', () => loadEntry(moleculeById(b.dataset.id))));
  $('smiles-form').addEventListener('submit', (e) => {
    e.preventDefault();
    const s = $('mol-smiles').value.trim();
    const match = MOLECULES.find(m => m.smiles === s);
    if (match) loadEntry(match);
    else { M.entry = null; loadSmiles(s); }
  });
  $('mol-basis').addEventListener('change', (e) => {
    M.basis = e.target.value;
    M.basisAuto = false;
    M.freq = M.lib && M.lib.basis === M.basis && M.freq ? M.freq : null;
    M.summary = null;
    renderAll();
    runSCF();
  });
  renderEditor();
  $('btn-opt')?.addEventListener('click', optimize);
  $('btn-freq')?.addEventListener('click', computeFrequencies);
  $('btn-scan')?.addEventListener('click', scanBond);
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

function seg(id, options, current) {
  return `<div class="seg" id="${id}">${options.map(o => `<button type="button" data-v="${o.v}" aria-pressed="${String(o.v) === String(current)}" ${o.disabled ? 'disabled' : ''}>${o.label}</button>`).join('')}</div>`;
}

function renderControls() {
  const c = $('controls');
  const S = M.summary;
  c.innerHTML = `
    <div class="ctl"><span class="lbl">Visualizza</span>
      ${seg('mol-overlay', [
        { v: 'model', label: 'Modello' },
        { v: 'mo', label: 'Orbitale', disabled: !S },
        { v: 'density', label: 'Densità', disabled: !S },
        { v: 'esp', label: 'Potenziale el.', disabled: !S },
        { v: 'charges', label: 'Cariche', disabled: !S },
      ], M.view.overlay)}
    </div>
    <div class="ctl"><span class="lbl">Etichette sugli atomi</span>
      <select id="mol-labels">
        ${[['none', 'Nessuna'], ['symbol', 'Simbolo'], ['formal', 'Carica formale (Lewis)'], ['ox', 'Numero di ossidazione'], ['mulliken', 'Carica parziale di Mulliken'], ['lowdin', 'Carica parziale di Löwdin'], ['chi', 'Elettronegatività di Pauling'], ['hyb', 'Ibridazione']]
          .map(([v, l]) => `<option value="${v}" ${M.view.labels === v ? 'selected' : ''} ${(v === 'mulliken' || v === 'lowdin') && !S ? 'disabled' : ''}>${l}</option>`).join('')}
      </select>
    </div>
    ${M.view.overlay === 'mo' || M.view.overlay === 'density' ? `
    <div class="ctl"><span class="lbl">${M.view.overlay === 'mo' ? 'Isovalore |ψ| (a.u.)' : 'Isovalore ρ (e/bohr³)'}</span>
      <div class="range-row"><input type="range" id="mol-iso" min="${M.view.overlay === 'mo' ? 0.01 : 0.002}" max="${M.view.overlay === 'mo' ? 0.2 : 0.3}" step="0.002" value="${M.view.overlay === 'mo' ? M.view.iso : (M.view.isoRho ?? 0.05)}"><output id="mol-iso-out">${nf(M.view.overlay === 'mo' ? M.view.iso : (M.view.isoRho ?? 0.05), 3)}</output></div>
    </div>` : ''}
    <label class="mini-toggle"><input type="checkbox" id="mol-dipole" ${M.view.dipole ? 'checked' : ''} ${!S ? 'disabled' : ''}> freccia del momento di dipolo</label>`;
  bindSeg('mol-overlay', (v) => { M.view.overlay = v; renderControls(); render3D(); });
  $('mol-labels').addEventListener('change', (e) => { M.view.labels = e.target.value; render3D(); });
  $('mol-dipole').addEventListener('change', (e) => { M.view.dipole = e.target.checked; render3D(); });
  const iso = $('mol-iso');
  if (iso) {
    iso.addEventListener('input', () => { $('mol-iso-out').textContent = nf(+iso.value, 3); });
    iso.addEventListener('change', () => {
      if (M.view.overlay === 'mo') M.view.iso = +iso.value; else M.view.isoRho = +iso.value;
      render3D();
    });
  }
}

function bindSeg(id, fn) {
  $(id)?.querySelectorAll('button').forEach(b => b.addEventListener('click', () => fn(b.dataset.v)));
}

/** Nome dell'OM in base alla composizione (atomo e tipo di funzione dominanti). */
function moLabel(k, spin) {
  const S = M.summary;
  const comp = S.compositions[k]?.[spin === 'β' ? 'beta' : 'alpha'];
  if (!comp) return '';
  const byAtom = comp.byAtom.map((w, i) => ({ w, i })).sort((a, b) => b.w - a.w);
  const main = byAtom.filter(x => x.w > 0.15).slice(0, 3).map(x => `${sym(M.atoms[x.i].Z)}${x.i + 1}`);
  const L = comp.byL;
  const tot = L.reduce((s, x) => s + x, 0) || 1;
  const typ = ['s', 'p', 'd'].map((t, i) => (L[i] / tot > 0.2 ? `${t} ${Math.round(100 * L[i] / tot)}%` : null)).filter(Boolean).join(', ');
  const kind = comp.type ? `<b>${comp.type}</b> ` : '';
  const bond = comp.bonding ? `<span class="bchar ${comp.bonding === 'legante' ? 'b-pos' : comp.bonding === 'antilegante' ? 'b-neg' : ''}">${comp.bonding}</span> · ` : '';
  return `${kind}${bond}${main.join(' ')}${typ ? ` · ${typ}` : ''}`;
}

function renderPanelBody() {
  const body = $('panel-body');
  const S = M.summary;
  if (!S) { body.innerHTML = '<p class="desc-muted">Il calcolo quantistico è in corso…</p>'; return; }
  const n = S.nbf;
  const rows = [];
  const add = (k, spin) => {
    const eps = spin === 'β' ? S.epsB : S.epsA;
    const nocc = spin === 'β' ? S.nbeta : S.nalpha;
    const occ = S.unrestricted ? (k < nocc ? 1 : 0) : (k < S.nalpha ? 2 : 0);
    const tag = k === nocc - 1 ? 'HOMO' : k === nocc ? 'LUMO' : '';
    const sel = M.view.mo && M.view.mo.index === k && M.view.mo.spin === spin;
    rows.push(`<tr class="clickable ${sel ? 'active' : ''}" data-k="${k}" data-spin="${spin ?? ''}">
      <td class="num">${k + 1}${spin ? spin : ''}</td><td>${tag ? `<span class="pill">${tag}</span>` : ''}</td>
      <td class="num">${nf(eps[k] * HARTREE_EV, 2)}</td><td class="num">${occ}</td><td class="comp">${moLabel(k, spin)}</td></tr>`);
  };
  const nShow = Math.min(n, Math.max(S.nalpha, S.nbeta) + 6);
  for (let k = 0; k < nShow; k++) {
    add(k, S.unrestricted ? 'α' : null);
    if (S.unrestricted) add(k, 'β');
  }
  const scale = FREQ_SCALE[M.basis] ?? 0.9;
  const vib = M.freq ? M.freq.modes.map((m, i) => `
    <tr class="clickable ${M.view.vib === i ? 'active' : ''}" data-vib="${i}">
      <td class="num">${i + 1}</td><td class="num">${m.freq < 0 ? `${nf(-m.freq, 0)}i` : nf(m.freq * scale, 0)}</td><td class="num">${nf(m.freq, 0)}</td><td class="num">${nf(m.ir ?? 0, 1)}</td></tr>`).join('') : '';
  body.innerHTML = `
    <div>
      <h3>Orbitali molecolari</h3>
      <div class="table-scroll short">
      <table class="orb-table mo-table">
        <thead><tr><th class="num">#</th><th></th><th class="num">ε (eV)</th><th class="num">e⁻</th><th>carattere</th></tr></thead>
        <tbody>${rows.join('')}</tbody>
      </table></div>
      <p class="desc-muted" style="margin-top:6px">Clicca un orbitale per vederlo in 3D. σ/π: simmetria rispetto all'asse (molecole lineari) o al piano (molecole planari). Legante/antilegante: popolazione di sovrapposizione di Mulliken sui legami. Poi atomi e funzioni (s, p, d) con il peso maggiore.</p>
    </div>
    <div>
      <h3>Vibrazioni</h3>
      ${M.freq ? `
      <div class="table-scroll short">
      <table class="orb-table">
        <thead><tr><th class="num">modo</th><th class="num">ν̃ scalata</th><th class="num">ν̃ armonica</th><th class="num">IR (km/mol)</th></tr></thead>
        <tbody>${vib}</tbody>
      </table></div>
      <p class="desc-muted" style="margin-top:6px">Frequenze in cm⁻¹. La scalatura (×${nf(scale, 4)} per HF/${M.basis}) compensa l'anarmonicità e la mancanza di correlazione. Clicca un modo per animarlo.${M.view.vib !== null ? ' <button type="button" class="linkish" id="vib-stop">Ferma</button>' : ''}</p>` : '<p class="desc-muted">Non ancora calcolate: usa "Calcola vibrazioni".</p>'}
    </div>`;
  body.querySelectorAll('tr[data-k]').forEach(tr => tr.addEventListener('click', () => selectMO(+tr.dataset.k, tr.dataset.spin || null)));
  body.querySelectorAll('tr[data-vib]').forEach(tr => tr.addEventListener('click', () => {
    M.view.vib = +tr.dataset.vib;
    renderPanelBody();
    render3D();
    drawCharts();
  }));
  $('vib-stop')?.addEventListener('click', () => { M.view.vib = null; renderPanelBody(); render3D(); drawCharts(); });
}

function selectMO(index, spin) {
  M.view.mo = { index, spin: spin || null };
  if (M.view.overlay !== 'mo') M.view.overlay = 'mo';
  renderControls();
  renderPanelBody();
  render3D();
  drawCharts();
}

function stopVibration() {
  if (viewer) viewer.frameCallbacks.clear();
}

function extentOf(atoms) {
  let r = 0;
  for (const a of atoms) r = Math.max(r, Math.hypot(...a.xyz));
  return r + 3.2;
}

let gridToken = 0;

async function render3D() {
  if (!active || !M.atoms) return;
  const S = M.summary;
  const token = ++gridToken;
  viewer.clear();
  const ext = extentOf(M.atoms);
  viewer.frame(ext, !M.resetView);
  if (M.resetView) { viewer.resetView(); M.resetView = false; }
  viewer.setAxesVisible(false);
  const lw = M.analysis;
  const labels = M.atoms.map((a, i) => {
    const info = lw.atoms[i];
    switch (M.view.labels) {
      case 'symbol': return sym(a.Z);
      case 'formal': return info.formalCharge ? sgn(info.formalCharge, 0) : '';
      case 'ox': return (info.oxidation > 0 ? '+' : info.oxidation < 0 ? '−' : '') + Math.abs(info.oxidation);
      case 'mulliken': return S ? sgn(S.mulliken[i], 2) : '';
      case 'lowdin': return S ? sgn(S.lowdin[i], 2) : '';
      case 'chi': return PAULING[a.Z] !== null ? nf(PAULING[a.Z], 2) : '';
      case 'hyb': return info.hybridization !== '—' ? info.hybridization : '';
      default: return '';
    }
  });
  const bonds = M.graph.bonds.map((b, k) => ({ a: b.a, b: b.b, order: lw.bonds[k].averageOrder }));
  const atomColors = M.view.overlay === 'charges' && S ? S.mulliken.map(q => chargeColor(q)) : null;
  const mol = buildMolecule(M.atoms, bonds, { labels, atomColors });
  viewer.add(mol.group);
  M.mol3d = mol;
  highlightSelection();
  if (M.view.dipole && S) {
    const arrow = dipoleArrow(S.dipole.vector, S.dipole.debye, ext);
    if (arrow) viewer.add(arrow);
  }
  // animazione del modo normale selezionato
  if (M.view.vib !== null && M.freq) {
    const mode = M.freq.modes[M.view.vib];
    const x0 = M.atoms.map(a => a.xyz);
    const amp = 0.35; // bohr
    const period = 1.1;
    viewer.frameCallbacks.add((t) => {
      const s = amp * Math.sin(2 * Math.PI * t / period);
      mol.update(x0.map((p, i) => [p[0] + s * mode.d[3 * i], p[1] + s * mode.d[3 * i + 1], p[2] + s * mode.d[3 * i + 2]]));
    });
  }
  const title = `${formulaOf()}${M.entry?.smiles === M.smiles ? ` · ${escapeHtml(M.name)}` : ''}`;
  let sub = S ? `HF/${M.basis} · ${S.pointGroup}` : 'calcolo in corso';
  let legend = [];
  let note = '';
  $('viewport-title').innerHTML = `${title}<small>${sub}</small>`;

  if (S && (M.view.overlay === 'mo' || M.view.overlay === 'density' || M.view.overlay === 'esp')) {
    const N = M.view.overlay === 'esp' ? 56 : 60;
    const half = M.view.overlay === 'esp' ? ext + 1.5 : ext;
    try {
      if (M.view.overlay === 'mo') {
        const { index, spin } = M.view.mo;
        const g = await request('grid', { kind: 'mo', index, spin, res: N, half, center: [0, 0, 0] });
        if (token !== gridToken) return;
        viewer.addSurface(g, M.view.iso, 1, cssVar('--phase-pos'), { opacity: 0.78 });
        viewer.addSurface(g, M.view.iso, -1, cssVar('--phase-neg'), { opacity: 0.78 });
        const eps = (spin === 'β' ? S.epsB : S.epsA)[index];
        const nocc = spin === 'β' ? S.nbeta : S.nalpha;
        const tag = index === nocc - 1 ? 'HOMO' : index === nocc ? 'LUMO' : index < nocc ? 'occupato' : 'virtuale';
        sub = `OM ${index + 1}${spin ?? ''} (${tag}) · ε = ${nf(eps * HARTREE_EV, 2)} eV`;
        legend = [[cssVar('--phase-pos'), 'ψ > 0'], [cssVar('--phase-neg'), 'ψ < 0']];
        note = `Isosuperficie |ψ| = ${nf(M.view.iso, 3)} a.u. dell'orbitale molecolare calcolato (combinazione lineare delle ${S.nbf} funzioni di base).`;
      } else if (M.view.overlay === 'density') {
        const g = await request('grid', { kind: 'density', res: N, half, center: [0, 0, 0] });
        if (token !== gridToken) return;
        const iso = M.view.isoRho ?? 0.05;
        viewer.addSurface(g, iso, 1, cssVar('--accent'), { opacity: 0.55 });
        sub = `densità elettronica ρ = ${nf(iso, 3)} e/bohr³`;
        note = 'Superficie di densità elettronica costante. A valori alti (≈ 0,2) resta attorno ai nuclei; attorno a 0,05 si vedono i legami; a 0,002 la superficie di van der Waals.';
      } else {
        setBusy('Potenziale elettrostatico sulla superficie di van der Waals…');
        const g = await request('grid', { kind: 'density', res: N, half, center: [0, 0, 0] });
        if (token !== gridToken) return;
        const { positions, normals } = marchingCubes(g, 0.002, 1);
        // ESP esatto su un reticolo grossolano vicino alla superficie, poi interpolato ai vertici
        const cell = 0.6;
        const key = (i, j, k) => `${i},${j},${k}`;
        const lattice = new Map();
        for (let v = 0; v < positions.length; v += 3) {
          const i = Math.floor(positions[v] / cell), j = Math.floor(positions[v + 1] / cell), k = Math.floor(positions[v + 2] / cell);
          for (let a = 0; a <= 1; a++) for (let b = 0; b <= 1; b++) for (let c = 0; c <= 1; c++) lattice.set(key(i + a, j + b, k + c), [i + a, j + b, k + c]);
        }
        const pts = new Float64Array(lattice.size * 3);
        const keys = [...lattice.keys()];
        keys.forEach((kk, q) => { const [i, j, k] = lattice.get(kk); pts[3 * q] = i * cell; pts[3 * q + 1] = j * cell; pts[3 * q + 2] = k * cell; });
        const r = await request('esp', { points: pts }, (p) => { if (p !== null) setBusy(`Potenziale elettrostatico: ${Math.round(p * 100)}% (${keys.length} punti)`); });
        if (token !== gridToken) return;
        setBusy(null);
        const val = new Map(keys.map((kk, q) => [kk, r.values[q]]));
        const vv = new Float32Array(positions.length / 3);
        for (let v = 0; v < vv.length; v++) {
          const x = positions[3 * v] / cell, y = positions[3 * v + 1] / cell, z = positions[3 * v + 2] / cell;
          const i = Math.floor(x), j = Math.floor(y), k = Math.floor(z);
          const fx = x - i, fy = y - j, fz = z - k;
          let s = 0;
          for (let a = 0; a <= 1; a++) for (let b = 0; b <= 1; b++) for (let c = 0; c <= 1; c++) {
            s += (a ? fx : 1 - fx) * (b ? fy : 1 - fy) * (c ? fz : 1 - fz) * (val.get(key(i + a, j + b, k + c)) ?? 0);
          }
          vv[v] = s;
        }
        let vmax = 0;
        for (const x of vv) vmax = Math.max(vmax, Math.abs(x));
        vmax = M.charge === 0 ? Math.min(vmax, 0.08) : vmax;
        const colors = new Float32Array(positions.length);
        for (let v = 0; v < vv.length; v++) espColor(vv[v], vmax, colors, 3 * v);
        viewer.addColoredSurface(positions, normals, colors, { opacity: 0.88 });
        sub = `potenziale elettrostatico sulla superficie ρ = 0,002 e/bohr³`;
        legend = [['#e8322a', `negativo (${nf(-vmax * HARTREE_KJ, 0)} kJ/mol)`], ['#f3f3f3', '0'], ['#2d63e8', `positivo (+${nf(vmax * HARTREE_KJ, 0)} kJ/mol)`]];
        note = 'Energia di una carica di prova +1 sulla superficie della molecola: in rosso le zone ricche di elettroni (attaccate dagli elettrofili), in blu quelle povere.';
      }
    } catch (e) {
      if (e.message !== 'annullato') note = `Griglia non calcolata: ${e.message}`;
      setBusy(null);
    }
  } else if (M.view.overlay === 'charges' && S) {
    legend = [[chargeColor(-0.8), 'δ− (Mulliken)'], [chargeColor(0), '0'], [chargeColor(0.8), 'δ+']];
    note = 'Atomi colorati secondo la carica parziale di Mulliken calcolata dalla funzione d\'onda.';
  } else {
    note = 'Modello a sfere e bastoncini: i legami doppi e tripli sono disegnati con 2 e 3 bastoncini; il tratto trasparente indica un legame delocalizzato (risonanza).';
  }
  if (M.view.vib !== null && M.freq) {
    const m = M.freq.modes[M.view.vib];
    sub = `modo ${M.view.vib + 1}: ${nf(m.freq * (FREQ_SCALE[M.basis] ?? 0.9), 0)} cm⁻¹ (scalata)`;
  }
  if (token !== gridToken) return;
  $('viewport-title').innerHTML = `${title}<small>${sub}</small>`;
  $('viewport-legend').innerHTML = legend.map(([c, l]) => `<span><i style="background:${c}"></i>${l}</span>`).join('');
  if (M.view.dipole && S && S.dipole.debye > 0.02) note += ` Freccia: momento di dipolo (${nf(S.dipole.debye, 2)} D), verso l'estremità negativa.`;
  $('viewport-note').textContent = note;
}

// ---------------------------------------------------------------------------
// Grafici
// ---------------------------------------------------------------------------

function drawCharts() {
  if (!active || !M.analysis) return;
  // struttura di Lewis
  $('slice-title').textContent = 'Struttura di Lewis';
  const lw = M.analysis;
  drawLewis($('chart-slice'), {
    atoms: M.atoms, bonds: M.graph.bonds, orders: lw.orders, lone: lw.lone,
    formalCharges: lw.formalCharges, resonanceCount: lw.resonanceCount,
  });
  $('slice-note').textContent = `${lw.sigmaBonds} legami σ e ${lw.piBonds} legami π. Coppie di punti: coppie solitarie; cerchi: cariche formali.${lw.unpaired ? ' Punto giallo: elettrone spaiato (radicale).' : ''}`;

  const S = M.summary;
  if (S) {
    // diagramma degli OM (valenza): si omettono gli orbitali di core molto profondi
    const levels = [];
    const nShow = Math.min(S.nbf, Math.max(S.nalpha, S.nbeta) + 3);
    let core = 0;
    for (const a of M.atoms) core += a.Z > 18 ? 9 : a.Z > 10 ? 5 : a.Z > 2 ? 1 : 0;
    const start = Math.min(core, Math.max(0, S.nbeta - 1));
    const push = (eps, spin, nocc) => {
      for (let k = start; k < nShow; k++) {
        levels.push({ index: k, spin, e: eps[k] * HARTREE_EV, occ: S.unrestricted ? (k < nocc ? 1 : 0) : (k < S.nalpha ? 2 : 0), label: k === nocc - 1 ? 'HOMO' : k === nocc ? 'LUMO' : '' });
      }
    };
    push(S.epsA, S.unrestricted ? 'α' : null, S.nalpha);
    if (S.unrestricted) push(S.epsB, 'β', S.nbeta);
    M.moRegions = drawMOLevels($('chart-radial'), {
      levels, selected: M.view.mo, unrestricted: S.unrestricted,
      ncoreShown: start ? `${start} orbitali di core non mostrati` : '',
    });
    $('radial-title').textContent = 'Orbitali molecolari calcolati';
    $('radial-note').textContent = S.unrestricted
      ? 'Calcolo a spin non ristretto (UHF): gli orbitali α e β hanno energie diverse. Le frecce indicano gli elettroni.'
      : 'Energie degli orbitali canonici Hartree–Fock. Livelli degeneri affiancati. Clicca un livello per vederne la forma.';
  } else {
    const c = $('chart-radial');
    c.getContext('2d').clearRect(0, 0, c.width, c.height);
    $('radial-note').textContent = '';
  }

  // spettro IR o curva di dissociazione
  if (M.view.chart2 === 'scan' && M.scan) {
    drawScan();
  } else if (M.freq) {
    const exp = M.entry?.smiles === M.smiles ? M.entry.exp?.freq ?? [] : [];
    drawIR($('chart-levels'), { modes: M.freq.modes, scale: FREQ_SCALE[M.basis] ?? 0.9, exp, selected: M.view.vib });
    $('levels-title').textContent = 'Spettro infrarosso calcolato';
    $('levels-note').textContent = `Frequenze armoniche HF/${M.basis} scalate; intensità dalle derivate del momento di dipolo.${exp.length ? ' Triangoli: bande fondamentali sperimentali.' : ''} I modi con intensità nulla non si vedono in IR (ma possono essere attivi in Raman).`;
  } else {
    const c = $('chart-levels');
    const ctx = c.getContext('2d');
    ctx.clearRect(0, 0, c.width, c.height);
    $('levels-title').textContent = 'Spettro infrarosso';
    $('levels-note').textContent = 'Calcola le vibrazioni per ottenere lo spettro IR.';
  }
}

function drawScan() {
  const pts = M.scan;
  const rhf = { xs: pts.map(p => p.d * BOHR_ANG * 100), ys: [], color: cssVar('--phase-neg'), label: 'RHF' };
  const uhf = { xs: pts.map(p => p.d * BOHR_ANG * 100), ys: [], color: cssVar('--accent'), label: 'UHF' };
  // energie relative al minimo, in kJ/mol
  const all = pts.flatMap(p => [p.rhf, p.uhf].filter(v => v !== null));
  const emin = Math.min(...all);
  rhf.ys = pts.map(p => (p.rhf - emin) * HARTREE_KJ);
  uhf.ys = pts.map(p => (p.uhf === null ? NaN : (p.uhf - emin) * HARTREE_KJ));
  const cap = (s) => ({ ...s, ys: s.ys.map(v => Math.min(v, 1500)), keepZero: true });
  // minimo, costante di forza e dissociazione (curva UHF)
  const E = pts.map(p => Math.min(p.rhf, p.uhf ?? Infinity));
  let im = 0;
  E.forEach((e, i) => { if (e < E[im]) im = i; });
  let text = '';
  if (im > 0 && im < pts.length - 1) {
    const [x0, x1, x2] = [pts[im - 1].d, pts[im].d, pts[im + 1].d];
    const [y0, y1, y2] = [E[im - 1], E[im], E[im + 1]];
    // parabola per tre punti: k = E''(r_e)
    const k = 2 * ((y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)) / (x2 - x0);
    const re = x1 - ((y2 - y0) / (x2 - x0)) / k;
    const [a, b] = M.atoms;
    const mu = ISOTOPE_MASS[a.Z] * ISOTOPE_MASS[b.Z] / (ISOTOPE_MASS[a.Z] + ISOTOPE_MASS[b.Z]);
    const omega = Math.sqrt(k / mu) * 5140.48;
    const De = (E[E.length - 1] - E[im]) * HARTREE_KJ;
    text = `Minimo a r_e ≈ ${nf(re * BOHR_ANG, 3)} Å; costante di forza k = ${nf(k * 15.569, 1)} N/cm; ω = √(k/μ) ≈ ${nf(omega, 0)} cm⁻¹; energia di dissociazione D_e ≈ ${nf(De, 0)} kJ/mol (UHF). `;
  }
  drawLineChart($('chart-levels'), { series: [cap(rhf), cap(uhf)], log: false, xlabel: 'r (pm)', ylabel: 'E − E_min (kJ/mol)', legend: true });
  $('levels-title').textContent = 'Curva di energia potenziale';
  $('levels-note').textContent = `${text}A grande distanza RHF sbaglia: impone elettroni appaiati nello stesso orbitale anche quando il legame è rotto. UHF separa correttamente gli atomi.`;
}

// ---------------------------------------------------------------------------
// Tabelle di analisi
// ---------------------------------------------------------------------------

function renderAnalysis() {
  const el = $('analysis');
  const lw = M.analysis;
  const S = M.summary;
  const exp = M.entry?.smiles === M.smiles ? M.entry.exp ?? {} : {};
  const label = (i) => `${sym(M.atoms[i].Z)}${i + 1}`;
  const dist = (i, j) => Math.hypot(...[0, 1, 2].map(k => M.atoms[i].xyz[k] - M.atoms[j].xyz[k])) * BOHR_ANG;
  const atomRows = lw.atoms.map((a, i) => `
    <tr>
      <td><b>${label(i)}</b></td>
      <td class="num">${a.chi !== null ? nf(a.chi, 2) : '—'}</td>
      <td class="num">${ALLEN[a.Z] ? nf(ALLEN[a.Z], 2) : '—'}</td>
      <td class="num">${a.valence}</td>
      <td class="num">${a.lonePairs}${a.radical ? ' + 1•' : ''}</td>
      <td class="num">${a.formalCharge ? sgn(a.formalCharge, 0) : '0'}</td>
      <td class="num">${a.oxidation > 0 ? '+' : a.oxidation < 0 ? '−' : ''}${Math.abs(a.oxidation)}</td>
      <td class="num">${S ? sgn(S.mulliken[i], 3) : '…'}</td>
      <td class="num">${S ? sgn(S.lowdin[i], 3) : '…'}</td>
      ${S?.unrestricted ? `<td class="num">${nf(S.spin[i], 3)}</td>` : ''}
      <td class="num">${S ? nf(S.valence[i], 2) : '…'}</td>
      <td>${a.hybridization}</td>
      <td>${a.vsepr ? `${a.vsepr.notation} · ${a.vsepr.molecularGeometry}` : '—'}</td>
    </tr>`).join('');
  const bondRows = M.graph.bonds.map((b, k) => {
    const info = lw.bonds[k];
    const d = dist(b.a, b.b);
    const keyName = [sym(M.atoms[b.a].Z), sym(M.atoms[b.b].Z)].sort().join('-');
    const dexp = exp.r?.[keyName];
    const mayer = S ? S.bondOrder[b.a][b.b] : null;
    const ord = info.averageOrder;
    const ordTxt = Math.abs(ord - Math.round(ord)) < 0.01 ? ['', 'singolo', 'doppio', 'triplo'][Math.round(ord)] ?? ord : `${nf(ord, 2)} (risonanza)`;
    return `<tr>
      <td><b>${label(b.a)}–${label(b.b)}</b></td>
      <td>${ordTxt}</td>
      <td class="num">${info.sigma} σ${info.pi ? ` + ${info.pi} π` : ''}</td>
      <td class="num">${mayer !== null ? nf(mayer, 2) : '…'}</td>
      <td class="num">${nf(d, 3)}</td>
      <td class="num">${dexp ? nf(dexp, 3) : '—'}</td>
      <td class="num">${nf(info.predictedLength / 100, 2)}</td>
      <td class="num">${nf(info.deltaChi, 2)}</td>
      <td class="num">${nf(info.ionicCharacter * 100, 0)}%</td>
      <td>${info.type}</td>
    </tr>`;
  }).join('');
  // angoli di legame attorno agli atomi centrali
  const x = M.atoms.flatMap(a => a.xyz);
  const angleRows = [];
  lw.atoms.forEach((a, c) => {
    if (a.neighbors < 2) return;
    const nb = M.graph.bonds.filter(b => b.a === c || b.b === c).map(b => (b.a === c ? b.b : b.a));
    const seen = new Set();
    for (let i = 0; i < nb.length; i++) for (let j = i + 1; j < nb.length; j++) {
      const ang = angleOf(x, nb[i], c, nb[j]) * 180 / Math.PI;
      const name = `${sym(M.atoms[nb[i]].Z)}-${sym(M.atoms[c].Z)}-${sym(M.atoms[nb[j]].Z)}`;
      const k = `${name}:${Math.round(ang)}`;
      if (seen.has(k)) continue;
      seen.add(k);
      const rev = name.split('-').reverse().join('-');
      const aexp = exp.angle?.[name] ?? exp.angle?.[rev];
      angleRows.push(`<tr><td><b>${label(nb[i])}–${label(c)}–${label(nb[j])}</b></td><td class="num">${nf(ang, 1)}°</td><td class="num">${aexp !== undefined ? `${nf(aexp, 1)}°` : '—'}</td><td>${a.vsepr ? `${a.vsepr.idealAngle} (${a.vsepr.notation})` : '—'}</td></tr>`);
    }
  });
  // termodinamica
  let thermoHTML = '<p class="desc-muted">Calcola le vibrazioni per le grandezze termodinamiche.</p>';
  if (M.freq && S) {
    const T = M.view.T ?? 298.15;
    const scale = FREQ_SCALE[M.basis] ?? 0.9;
    const modes = M.freq.modes.map(m => ({ freq: m.freq * scale }));
    const th = thermochemistry(M.atoms, modes, S.energy, { T, multiplicity: M.multiplicity });
    const Emp2 = S.mp2 !== null ? S.energy + S.mp2 : null;
    thermoHTML = `
      <div class="range-row" style="max-width:420px"><label for="mol-T" class="lbl" style="min-width:110px">Temperatura</label><input type="range" id="mol-T" min="50" max="2000" step="5" value="${T}"><output>${nf(T, 0)} K</output></div>
      <dl class="facts thermo">
        <dt>Energia di punto zero (ZPE)</dt><dd>${nf(th.zpe, 2)} kJ/mol</dd>
        <dt>Entalpia H − E<sub>el</sub></dt><dd>${nf(th.Hcorr * HARTREE_KJ, 2)} kJ/mol</dd>
        <dt>Entropia S°</dt><dd>${nf(th.S, 2)} J/(mol·K)</dd>
        <dt>&nbsp;&nbsp;traslazionale / rotazionale</dt><dd>${nf(th.parts.St, 1)} / ${nf(th.parts.Sr, 1)}</dd>
        <dt>&nbsp;&nbsp;vibrazionale / elettronica</dt><dd>${nf(th.parts.Sv, 1)} / ${nf(th.parts.Se, 1)}</dd>
        <dt>Capacità termica C<sub>v</sub></dt><dd>${nf(th.Cv, 2)} J/(mol·K)</dd>
        <dt>Energia libera G − E<sub>el</sub></dt><dd>${nf(th.Gcorr * HARTREE_KJ, 2)} kJ/mol</dd>
        <dt>H totale (HF${Emp2 !== null ? ' / MP2' : ''})</dt><dd>${nf(S.energy + th.Hcorr, 6)}${Emp2 !== null ? ` / ${nf(Emp2 + th.Hcorr, 6)}` : ''} Ha</dd>
        <dt>Numero di simmetria σ</dt><dd>${th.sigma} (${th.pointGroup})</dd>
        <dt>Costanti rotazionali</dt><dd>${th.rotationalConstants.filter(Boolean).map(b => nf(b, 3)).join(' · ')} cm⁻¹</dd>
      </dl>
      <p class="desc-muted">Gas ideale a 1 atm, rotore rigido, oscillatore armonico con frequenze scalate (×${nf(scale, 4)}).</p>`;
  }
  el.innerHTML = `
    <h2>Analisi chimica</h2>
    <div class="analysis-grid">
      <article class="wide">
        <h3>Atomi: elettronegatività, cariche, stati di ossidazione, VSEPR</h3>
        <div class="table-wrap"><table class="data-table">
          <thead><tr><th>Atomo</th><th class="num" title="Scala di Pauling">χ Pauling</th><th class="num" title="Scala di Allen (energia media degli elettroni di valenza)">χ Allen</th><th class="num">e⁻ di valenza</th><th class="num">coppie solitarie</th><th class="num">carica formale</th><th class="num">n. ossidazione</th><th class="num">q Mulliken</th><th class="num">q Löwdin</th>${S?.unrestricted ? '<th class="num">spin</th>' : ''}<th class="num" title="Somma degli ordini di legame di Mayer">valenza</th><th>ibridazione</th><th>VSEPR</th></tr></thead>
          <tbody>${atomRows}</tbody>
        </table></div>
        <p class="desc-muted">Carica formale = e⁻ di valenza − e⁻ non condivisi − ½ e⁻ di legame. Numero di ossidazione: gli elettroni di ogni legame vanno all'atomo più elettronegativo. Le cariche di Mulliken e Löwdin vengono invece dalla funzione d'onda: sono frazionarie e dipendono dalla base.</p>
      </article>
      <article class="wide">
        <h3>Legami: Lewis contro meccanica quantistica</h3>
        <div class="table-wrap"><table class="data-table">
          <thead><tr><th>Legame</th><th>Lewis</th><th class="num">σ / π</th><th class="num" title="Ordine di legame di Mayer, dalla matrice densità">ordine di Mayer</th><th class="num">lunghezza calc. (Å)</th><th class="num">sperim. (Å)</th><th class="num" title="Somma dei raggi covalenti di Pyykkö">raggi covalenti (Å)</th><th class="num">Δχ</th><th class="num" title="Pauling: 1 − exp(−Δχ²/4)">carattere ionico</th><th>tipo</th></tr></thead>
          <tbody>${bondRows}</tbody>
        </table></div>
        <p class="desc-muted">L'ordine di Mayer B<sub>AB</sub> = Σ<sub>μ∈A,ν∈B</sub>(PS)<sub>μν</sub>(PS)<sub>νμ</sub> misura quante coppie di elettroni sono condivise: vale circa 1, 2, 3 per legami singoli, doppi e tripli, ed è minore per i legami polari. Carattere ionico secondo Pauling: 1 − e<sup>−(Δχ)²/4</sup>.</p>
      </article>
      <article>
        <h3>Angoli di legame</h3>
        ${angleRows.length ? `<div class="table-wrap"><table class="data-table"><thead><tr><th>Angolo</th><th class="num">calcolato</th><th class="num">sperim.</th><th>VSEPR</th></tr></thead><tbody>${angleRows.join('')}</tbody></table></div>` : '<p class="desc-muted">Molecola biatomica: nessun angolo.</p>'}
      </article>
      <article>
        <h3>Termodinamica statistica</h3>
        ${thermoHTML}
      </article>
    </div>`;
  const Tin = $('mol-T');
  if (Tin) {
    Tin.addEventListener('input', () => { Tin.nextElementSibling.textContent = `${Tin.value} K`; });
    Tin.addEventListener('change', () => { M.view.T = +Tin.value; renderAnalysis(); });
  }
}

// ---------------------------------------------------------------------------
// Editor: modifica della molecola nel modello 3D
// ---------------------------------------------------------------------------

function highlightSelection() {
  if (!M.mol3d) return;
  M.mol3d.atomMeshes.forEach((m, i) => {
    const on = M.sel.includes(i);
    m.material.emissive = new THREE.Color(on ? cssVar('--accent') : '#000000');
    m.material.emissiveIntensity = on ? 0.55 : 0;
  });
}

const EDIT_ELEMENTS = [1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35];

function renderEditor() {
  const box = $('mol-editor');
  if (!box) return;
  const sel = M.sel;
  const btn = (act, label, extra = '') => `<button type="button" class="mini-btn" data-act="${act}" ${extra}>${label}</button>`;
  let html = '<p class="lbl">Modifica la molecola</p>';
  if (!sel.length) {
    html += '<p class="hint">Clicca un atomo nel modello 3D per selezionarlo. Con due atomi selezionati puoi creare o cambiare il legame tra loro.</p>';
  } else if (sel.length === 1) {
    const a = M.graph.atoms[sel[0]];
    html += `<p class="hint">Selezionato: <b>${sym(a.Z)}${sel[0] + 1}</b>${a.charge ? ` (carica ${sgn(a.charge, 0)})` : ''}</p>
      <div class="edit-row"><span>Sostituisci con</span>${EDIT_ELEMENTS.map(Z => btn(`replace:${Z}`, sym(Z), Z === a.Z ? 'disabled' : '')).join('')}</div>
      <div class="edit-row"><span>Aggiungi legato</span>${[6, 7, 8, 9, 17].map(Z => btn(`add:${Z}`, sym(Z))).join('')}${btn('add:6:2', 'C=')}${btn('add:8:2', '=O')}${btn('add:7:3', '≡N')}</div>
      <div class="edit-row"><span>Carica</span>${btn('charge:-1', '−1')}${btn('charge:1', '+1')}${btn('delete', 'Elimina atomo')}</div>`;
  } else {
    const [i, j] = sel;
    const b = M.graph.bonds.find(x => (x.a === i && x.b === j) || (x.a === j && x.b === i));
    const cur = b ? b.order : 0;
    html += `<p class="hint">Legame tra <b>${sym(M.graph.atoms[i].Z)}${i + 1}</b> e <b>${sym(M.graph.atoms[j].Z)}${j + 1}</b>: ${['nessuno', 'singolo', 'doppio', 'triplo'][cur]}</p>
      <div class="edit-row"><span>Legame</span>${[0, 1, 2, 3].map(o => btn(`bond:${o}`, ['nessuno', 'singolo', 'doppio', 'triplo'][o], o === cur ? 'disabled' : '')).join('')}</div>`;
  }
  box.innerHTML = html;
  box.querySelectorAll('button[data-act]').forEach(b => b.addEventListener('click', () => applyEdit(b.dataset.act)));
}

function applyEdit(act) {
  const atoms = M.graph.atoms.map(a => ({ Z: a.Z, charge: a.charge ?? 0, aromatic: false }));
  let bonds = M.graph.bonds.map(b => ({ a: b.a, b: b.b, order: b.order, aromatic: false }));
  const [cmd, p1, p2] = act.split(':');
  const i = M.sel[0];
  if (cmd === 'replace') {
    atoms[i].Z = +p1;
    atoms[i].charge = 0;
    // un H sostituito diventa un atomo pesante con un solo legame
  } else if (cmd === 'add') {
    atoms.push({ Z: +p1, charge: 0, aromatic: false });
    bonds.push({ a: i, b: atoms.length - 1, order: +(p2 ?? 1), aromatic: false });
    if (atoms[i].Z === 1) atoms[i].Z = 6; // non si lega nulla a un idrogeno: diventa carbonio
  } else if (cmd === 'charge') {
    atoms[i].charge += +p1;
  } else if (cmd === 'delete') {
    const drop = new Set([i]);
    bonds.forEach(b => {
      if (b.a === i && atoms[b.b].Z === 1) drop.add(b.b);
      if (b.b === i && atoms[b.a].Z === 1) drop.add(b.a);
    });
    const map = [];
    const kept = [];
    atoms.forEach((a, k) => { if (!drop.has(k)) { map[k] = kept.length; kept.push(a); } });
    bonds = bonds.filter(b => !drop.has(b.a) && !drop.has(b.b)).map(b => ({ ...b, a: map[b.a], b: map[b.b] }));
    atoms.length = 0;
    atoms.push(...kept);
    if (!atoms.length) return;
  } else if (cmd === 'bond') {
    const [a, b] = M.sel;
    const o = +p1;
    const k = bonds.findIndex(x => (x.a === a && x.b === b) || (x.a === b && x.b === a));
    if (k >= 0) { if (o === 0) bonds.splice(k, 1); else bonds[k].order = o; } else if (o > 0) bonds.push({ a, b, order: o, aromatic: false });
  }
  // gli idrogeni si ricalcolano con le valenze standard
  const filled = fillHydrogens({ atoms, bonds });
  const smiles = writeSmiles(filled);
  if (!smiles) return;
  M.entry = null;
  M.name = 'Molecola personalizzata';
  loadSmiles(smiles);
}
