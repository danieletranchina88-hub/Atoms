// Modalità "Sandbox": dinamica molecolare reattiva in tempo reale. Si aggiungono atomi e molecole in
// una scatola, si controllano temperatura, volume e luce, e si osservano legami che si formano e si
// rompono, specie chimiche, reazioni, pressione e distribuzione di Maxwell–Boltzmann.

import * as THREE from 'three';
import { PRESETS, SANDBOX_MOLECULES, SANDBOX_IONS, SANDBOX_ELEMENTS, SPECIES_NAMES } from '../chem/sandboxData.js';
import { ATOM_PARAMS, KB_EV } from '../chem/reactiveData.js';
import { cpkColor, covalentRadius, PAULING, valenceElectrons } from '../chem/elementData.js';
import { parseSmiles, hillFormula, writeSmiles } from '../chem/smiles.js';
import { ELEMENTS } from '../physics/elements.js';
import { drawXY, wavelengthColor } from './chemCharts.js';
import { MV2 } from '../chem/md.js';

const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const sym = (Z) => ELEMENTS[Z - 1].symbol;
const nf = (v, d = 2) => (v === null || v === undefined || !Number.isFinite(v) ? '—' : Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d, useGrouping: Math.abs(v) >= 1e4 }));
const sgn = (v, d = 2) => (v > 0 ? '+' : v < 0 ? '−' : '') + nf(Math.abs(v), d);
const KJ_PER_EV = 96.48533212;
const NA = 6.02214076e23;

let viewer = null;
let worker = null;
let active = false;
let callback = null;
let openInMolecule = null;

const SB = {
  preset: PRESETS.find(p => p.id === 'empty') ?? PRESETS[0],
  fidelity: 'auto',
  orbital: 'homo',
  cloudNote: '',
  frame: null,
  fresh: false,
  census: null,
  history: [],
  events: [],
  pulses: [],
  phase: null,
  phaseSince: 0,
  bondBorn: new Map(),
  selected: -1,
  tool: 'select',
  add: { kind: 'mol', id: 'H2O', symbol: 'C' },
  count: 5,
  color: 'element',
  style: 'ball',
  chart1: 'energy',
  lambda: 400,
  info: '',
  lastPanel: 0,
  lastCharts: 0,
  userPaused: false,
  photonFlash: null,
  phys: null,
};

// ---------------------------------------------------------------------------
// Worker
// ---------------------------------------------------------------------------

function post(msg) {
  if (!worker) startWorker();
  worker.postMessage(msg);
}

function startWorker() {
  worker = new Worker(new URL('../chem/sandboxWorker.js', import.meta.url), { type: 'module' });
  worker.onmessage = (ev) => {
    const m = ev.data;
    if (m.type === 'frame') {
      const oldMode = SB.frame?.stats.forceField;
      const oldPaused = SB.userPaused;
      const elec = (st) => st ? `${st.charge}|${st.field}|${st.multiplicity}|${st.fidelity}` : '';
      const oldElec = elec(SB.frame?.stats);
      SB.frame = m;
      SB.userPaused = m.stats.paused;
      if (active && (oldMode !== m.stats.forceField || oldPaused !== SB.userPaused || oldElec !== elec(m.stats))) renderControls();
      SB.fresh = true;
      if (SB.tour) advanceTour(m.stats);
      if (m.stats?.phase) notePhase(m.stats.phase);
      if (m.phys) {
        const key = (p) => `${p?.elements?.join(',')}|${p?.rdf?.pair?.join('-')}|${p?.barostat?.on}|${p?.barostat?.P0}`;
        const changed = key(m.phys) !== key(SB.phys);
        SB.phys = m.phys;
        if (changed && active) renderControls();
      }
      if (m.census) {
        SB.census = m.census;
        SB.history = m.census.history;
        renderHud();
        if (m.census.events.length) {
          SB.events.push(...m.census.events);
          if (SB.events.length > 500) SB.events.splice(0, SB.events.length - 500);
          spawnReactions(m.census.events);
        }
      }
    } else if (m.type === 'info') {
      SB.info = m.text;
      renderSide();
    } else if (m.type === 'forcefield') {
      SB.forceField = m.kind;
      SB.fidelity = m.fidelity ?? SB.fidelity;
      if (m.why) SB.info = m.why;
      if (active) { renderControls(); renderSide(); }
    } else if (m.type === 'wave') {
      onWave(m);
    } else if (m.type === 'snapshot') {
      saveSnapshot(m.name, m.data);
    } else if (m.type === 'photon') {
      SB.photonFlash = { i: m.i, j: m.j, t: performance.now(), lambda: m.lambda };
    } else if (m.type === 'error') {
      if (m.paused) { SB.userPaused = true; renderControls(); }
      SB.info = `Errore: ${m.text}`;
      renderSide();
    }
  };
}

// ---------------------------------------------------------------------------
// Attivazione
// ---------------------------------------------------------------------------

export function initSandboxMode(v, { openMolecule }) {
  viewer = v;
  openInMolecule = openMolecule;
  const canvas = viewer.renderer.domElement;
  canvas.addEventListener('pointerdown', onPointerDown);
  canvas.addEventListener('pointermove', onPointerMove);
  canvas.addEventListener('pointerup', onPointerUp);
  $('sb-chart1').addEventListener('change', (e) => { SB.chart1 = e.target.value; drawCharts(); });
}

export function activateSandbox() {
  active = true;
  document.body.dataset.mode = 'sandbox';
  viewer.clear();
  viewer.setAxesVisible(false);
  buildScene();
  if (!worker) { startWorker(); loadPreset(SB.preset); }
  post({ type: 'set', paused: SB.userPaused });
  renderSide();
  renderControls();
  renderAnalysis();
  $('busy').hidden = true;
}

export function deactivateSandbox() {
  active = false;
  clearCloud();
  if (worker) worker.postMessage({ type: 'set', paused: true });
  if (callback) viewer.frameCallbacks.delete(callback);
  callback = null;
  viewer.controls.enabled = true;
}

export function sandboxRedraw() {
  if (!active) return;
  scene.dirtyTheme = true;
  drawCharts();
}


const PHASE_HEX = { gas: '#7ec8ff', liquido: '#f0b429', solido: '#8b93ff', misto: '#e07a9a', vuota: '#9aa3b2' };

function ensureHud() {
  const vp = $('viewport');
  if (!vp || $('sb-hud')) return;
  const hud = document.createElement('div');
  hud.id = 'sb-hud';
  hud.className = 'sb-hud';
  hud.innerHTML = '<div id="sb-phase" class="sb-phase"></div><div id="sb-species" class="sb-species"></div><div id="sb-reactions" class="sb-reactions"></div>';
  vp.appendChild(hud);
}

function notePhase(phase) {
  if (!phase) return;
  if (!SB.phase || SB.phase.label !== phase.label) SB.phaseSince = performance.now();
  SB.phase = phase;
  renderHud();
}

function chargeText(q) {
  if (!q) return '0';
  return `${q > 0 ? '+' : '−'}${Math.abs(q)}`;
}

function formulaText(list) {
  return (list && list.length ? list.join(' + ') : '—').replace(/(\d+)/g, '<sub>$1</sub>');
}

function spawnReactions(events) {
  const now = performance.now();
  for (const e of events) {
    if (e.exchange) continue;
    const born = (e.products?.length ?? 0) > (e.reactants?.length ?? 0);
    SB.pulses.push({
      at: e.at ?? [0, 0, 0],
      t0: now,
      color: born ? '#3dde7a' : '#ff6b57',
      text: `${(e.reactants ?? []).join(' + ') || '·'} → ${(e.products ?? []).join(' + ') || '·'}`,
    });
  }
  if (SB.pulses.length > 24) SB.pulses.splice(0, SB.pulses.length - 24);
  renderHud();
}

function renderHud() {
  ensureHud();
  const phaseEl = $('sb-phase');
  const specEl = $('sb-species');
  const rxEl = $('sb-reactions');
  if (!phaseEl) return;
  const ph = SB.frame?.stats?.phase ?? SB.phase;
  if (ph) {
    phaseEl.hidden = false;
    phaseEl.dataset.phase = ph.label;
    phaseEl.innerHTML = `<b>${ph.title}</b><span>${ph.label === 'pochi' ? '' : `${Math.round((ph.frac ?? 0) * 100)}% coordinati · `}${ph.note ?? ''}</span>`;
  } else phaseEl.hidden = true;
  const species = SB.census?.species ?? [];
  specEl.innerHTML = species.slice(0, 6).map(([f, n]) => `<span>${f} <b>${n}</b></span>`).join('') || '<span>nessuna specie</span>';
  const recent = SB.events.filter(e => !e.exchange).slice(-4).reverse();
  rxEl.innerHTML = recent.map(e => `<div class="sb-rx"><b>${formulaText(e.reactants)} → ${formulaText(e.products)}</b><span>${nf(e.t / 1000, 2)} ps</span></div>`).join('') || '<div class="sb-rx dim">In attesa di una trasformazione…</div>';
}

// Giro delle fasi dell'argon (Lennard–Jones, ε/k = 119,8 K): le soglie sono in tempo simulato, non in secondi
// dell'orologio, così il sistema ha lo stesso tempo per equilibrarsi su qualunque computer.
// 300 K è sopra il punto critico dell'argon (150,7 K): fluido; 90 K: goccia di liquido col suo vapore; 40 K: solido.
const TOUR = [[300, 0], [90, 15], [40, 45]];   // [T in K, da t in ps]

function startPhaseTour() {
  const p = PRESETS.find(x => x.id === 'argon-liquid');
  if (!p) return;
  loadPreset(p);
  SB.userPaused = false;
  SB.tour = { step: 0 };
  post({ type: 'set', paused: false, thermostat: true, T: TOUR[0][0] });
  SB.info = 'Giro delle fasi: argon a 300 K (fluido sopra il punto critico di 151 K), da 15 ps a 90 K (goccia di liquido e vapore), da 45 ps a 40 K (solido). La fase è letta dalla struttura (q₆), le temperature del modello LJ non coincidono esattamente con quelle dell\'argon reale.';
  renderSide();
  renderControls();
}

function advanceTour(stats) {
  const tour = SB.tour;
  if (!tour || SB.preset.id !== 'argon-liquid') { SB.tour = null; return; }
  const next = TOUR[tour.step + 1];
  if (next && stats.t / 1000 >= next[1]) {
    tour.step++;
    post({ type: 'set', T: next[0], thermostat: true });
    if (tour.step === TOUR.length - 1) SB.tour = null;
  }
}

function phaseColor(i, pos, Z, N, out) {
  const cutoff = Z[i] === 1 ? 1.7 : 4.0;
  const c2 = cutoff * cutoff;
  let n = 0;
  const x = pos[3 * i], y = pos[3 * i + 1], z = pos[3 * i + 2];
  const step = N > 180 ? 2 : 1;
  for (let j = 0; j < N; j += step) {
    if (j === i) continue;
    const dx = x - pos[3 * j], dy = y - pos[3 * j + 1], dz = z - pos[3 * j + 2];
    if (dx * dx + dy * dy + dz * dz < c2) n++;
  }
  if (n >= 5) out.set('#7d8cff');
  else if (n >= 2) out.set('#f0b429');
  else out.set('#9fd7ff');
  return out;
}

function pruneBonds(bonds, NB) {
  const live = new Set();
  for (let k = 0; k < NB; k++) {
    const i = bonds[4 * k], j = bonds[4 * k + 1];
    live.add(i < j ? i + '-' + j : j + '-' + i);
  }
  if (SB.bondBorn.size > live.size + 40) {
    for (const key of SB.bondBorn.keys()) if (!live.has(key)) SB.bondBorn.delete(key);
  }
}

function pulseUpdate() {
  if (!scene.pulses) return;
  const now = performance.now();
  SB.pulses = SB.pulses.filter(p => now - p.t0 < 1400);
  scene.pulses.forEach((mesh, i) => {
    const p = SB.pulses[i];
    if (!p) { mesh.visible = false; return; }
    const u = (now - p.t0) / 1400;
    mesh.visible = true;
    mesh.position.set(p.at[0], p.at[1], p.at[2]);
    const s = 0.35 + u * 1.8;
    mesh.scale.setScalar(s);
    mesh.material.color.set(p.color);
    mesh.material.opacity = 0.55 * (1 - u);
  });
}

function loadPreset(p) {
  SB.preset = p;
  SB.tour = null;
  SB.frame = null;
  SB.events = [];
  SB.history = [];
  SB.census = null;
  SB.selected = -1;
  SB.info = '';
  if (p.color) SB.color = p.color;
  SB.style = p.style ?? (SB.style === 'cloud' || SB.style === 'orbital' ? 'ball' : SB.style);
  clearCloud();
  if (p.light) SB.lambda = p.light.lambda;
  post({ type: 'preset', preset: p });
  const quantum = p.forceField === 'hf' || p.forceField === 'mindo3' || p.quantum;
  const empty = !((p.add && p.add.length) || p.atoms);
  post({ type: 'set', paused: quantum, T: p.T });
  SB.userPaused = quantum;
  if (empty) SB.info = 'Scatola libera. Versa molecole o clicca nella scena: il modello si sceglie dalla composizione.';
  SB.pulses = [];
  SB.bondBorn = new Map();
  SB.phase = null;
  scene.framed = false;
  renderHud();
  if (active) { renderSide(); renderControls(); renderAnalysis(); }
}

// ---------------------------------------------------------------------------
// Scena 3D: atomi e legami come "instanced mesh" (migliaia di oggetti in una sola chiamata di disegno)
// ---------------------------------------------------------------------------

const scene = { group: null, atoms: null, bonds: null, box: null, sel: null, grabLine: null, cap: 0, bcap: 0, framed: false, boxSize: 0, dirtyTheme: false };
const tmpM = new THREE.Matrix4(), tmpQ = new THREE.Quaternion(), tmpV = new THREE.Vector3(), tmpS = new THREE.Vector3(), tmpC = new THREE.Color();
const UP = new THREE.Vector3(0, 1, 0);

function buildScene() {
  scene.group = new THREE.Group();
  viewer.add(scene.group);
  scene.cap = 0; scene.bcap = 0; scene.boxSize = 0;
  scene.atoms = null; scene.bonds = null;
  const boxGeo = new THREE.EdgesGeometry(new THREE.BoxGeometry(1, 1, 1));
  scene.box = new THREE.LineSegments(boxGeo, new THREE.LineBasicMaterial({ color: cssVar('--scene-axis'), transparent: true, opacity: 0.9 }));
  scene.group.add(scene.box);
  scene.sel = new THREE.Mesh(new THREE.SphereGeometry(1, 24, 16), new THREE.MeshBasicMaterial({ color: cssVar('--accent'), wireframe: true, transparent: true, opacity: 0.8 }));
  scene.sel.visible = false;
  scene.group.add(scene.sel);
  scene.grabLine = new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3()]), new THREE.LineBasicMaterial({ color: cssVar('--accent') }));
  scene.grabLine.visible = false;
  scene.group.add(scene.grabLine);
  scene.hb = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineDashedMaterial({ color: cssVar('--accent'), dashSize: 0.18, gapSize: 0.12, transparent: true, opacity: 0.85 }));
  scene.hb.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(6 * 2048), 3));
  scene.hb.frustumCulled = false;
  scene.group.add(scene.hb);
  scene.flash = new THREE.Mesh(new THREE.SphereGeometry(1, 16, 12), new THREE.MeshBasicMaterial({ color: '#ffffff', transparent: true, opacity: 0.6, depthWrite: false }));
  scene.flash.visible = false;
  scene.group.add(scene.flash);
  scene.pulses = [];
  for (let i = 0; i < 16; i++) {
    const mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 14, 10), new THREE.MeshBasicMaterial({ color: '#3dde7a', transparent: true, opacity: 0.45, depthWrite: false }));
    mesh.visible = false;
    scene.group.add(mesh);
    scene.pulses.push(mesh);
  }
  ensureHud();
  scene.framed = false;
  callback = () => updateScene();
  viewer.frameCallbacks.add(callback);
}

function ensureCapacity(N, NB) {
  if (!scene.atoms || scene.cap < N) {
    if (scene.atoms) { scene.group.remove(scene.atoms); scene.atoms.geometry.dispose(); scene.atoms.material.dispose(); }
    scene.cap = Math.max(256, 2 * N);
    scene.atoms = new THREE.InstancedMesh(new THREE.SphereGeometry(1, 20, 14), new THREE.MeshStandardMaterial({ roughness: 0.4, metalness: 0.05 }), scene.cap);
    scene.atoms.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(3 * scene.cap), 3);
    scene.atoms.frustumCulled = false;
    scene.group.add(scene.atoms);
  }
  if (!scene.bonds || scene.bcap < NB) {
    if (scene.bonds) { scene.group.remove(scene.bonds); scene.bonds.geometry.dispose(); scene.bonds.material.dispose(); }
    scene.bcap = Math.max(512, 2 * NB);
    const light = cssVar('--scene-mode') === 'light';
    scene.bonds = new THREE.InstancedMesh(new THREE.CylinderGeometry(1, 1, 1, 10, 1, false), new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.04 }), scene.bcap);
    scene.bonds.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(3 * scene.bcap), 3);
    scene.bonds.frustumCulled = false;
    scene.group.add(scene.bonds);
  }
}

function atomRadius(Z) {
  if (SB.style === 'vdw') return 0.5 * ATOM_PARAMS[Z][5] * 0.82;
  // nella nuvola si vedono solo i nuclei (puntiformi alla scala degli elettroni)
  if (SB.style === 'cloud' || SB.style === 'orbital') return Z === 1 ? 0.07 : 0.11;
  const cov = ATOM_PARAMS[Z][0] === 0 ? 0.5 * ATOM_PARAMS[Z][5] * 100 * 0.5 : covalentRadius(Z);
  return 0.12 + 0.0034 * cov;
}

function heatColor(t, out) {
  // blu (freddo) → bianco → rosso (caldo)
  const x = Math.max(0, Math.min(1, t));
  if (x < 0.5) out.setRGB(0.25 + 1.5 * x * 0.75, 0.45 + x * 1.1, 1);
  else out.setRGB(1, 1 - (x - 0.5) * 1.6, 1 - (x - 0.5) * 1.9);
  return out;
}

function chargeRGB(q, out) {
  const t = Math.max(-1, Math.min(1, q / 0.7));
  if (t < 0) out.setRGB(0.88, 0.88 + 0.7 * t, 0.88 + 0.78 * t);
  else out.setRGB(0.88 - 0.72 * t, 0.88 - 0.42 * t, 0.95);
  return out;
}

function updateScene() {
  if (!active) return;
  const f = SB.frame;
  if (scene.dirtyTheme) {
    scene.dirtyTheme = false;
    scene.box.material.color.set(cssVar('--scene-axis'));
    scene.hb.material.color.set(cssVar('--accent'));
    if (scene.bonds) scene.bonds.material.color.set(cssVar('--scene-mode') === 'light' ? '#8a93a3' : '#b7bfcc');
    SB.fresh = true;
  }
  if (!f) return;
  const L = f.stats.box;
  if (scene.boxSize !== L) { scene.box.scale.setScalar(L); scene.boxSize = L; }
  // nuovo esperimento: inquadra la scatola quando arriva il primo fotogramma con le sue dimensioni
  if (!scene.framed) {
    viewer.frame(L * 0.62);
    viewer.camera.position.set(0.55, -0.8, 0.45).normalize().multiplyScalar(L * 0.62 * 3.1);
    viewer.controls.target.set(0, 0, 0);
    scene.framed = true;
    renderControls();
  }
  flashUpdate();
  pulseUpdate();
  if (!SB.fresh) return;
  SB.fresh = false;
  const { N, Z, pos, q, ke, bonds } = f;
  const NB = bonds.length / 4;
  ensureCapacity(N, NB * 3);
  const Tref = Math.max(f.stats.T, 50);
  for (let i = 0; i < N; i++) {
    const r = atomRadius(Z[i]);
    tmpM.makeScale(r, r, r).setPosition(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
    scene.atoms.setMatrixAt(i, tmpM);
    if (SB.color === 'charge') chargeRGB(q[i], tmpC);
    else if (SB.color === 'ke') heatColor(ke[i] / (3 * KB_EV * Tref), tmpC);
    else if (SB.color === 'phase') phaseColor(i, pos, Z, N, tmpC);
    else tmpC.set(cpkColor(Z[i]));
    scene.atoms.setColorAt(i, tmpC);
  }
  scene.atoms.count = N;
  scene.atoms.instanceMatrix.needsUpdate = true;
  scene.atoms.instanceColor.needsUpdate = true;
  scene.atoms.boundingSphere = null;
  // legami: 1, 2 o 3 cilindri secondo l'ordine; gli spostamenti sono perpendicolari alla vista
  let nb = 0;
  const view = new THREE.Vector3().subVectors(viewer.camera.position, viewer.controls.target).normalize();
  const a = new THREE.Vector3(), b = new THREE.Vector3(), d = new THREE.Vector3(), perp = new THREE.Vector3();
  const showBonds = SB.style === 'ball';
  for (let k = 0; showBonds && k < NB; k++) {
    const i = bonds[4 * k], j = bonds[4 * k + 1], n = bonds[4 * k + 2], s = bonds[4 * k + 3];
    a.set(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
    b.set(pos[3 * j], pos[3 * j + 1], pos[3 * j + 2]);
    d.subVectors(b, a);
    const len = d.length();
    if (len < 1e-6) continue;
    d.divideScalar(len);
    perp.crossVectors(d, view);
    if (perp.lengthSq() < 1e-6) perp.set(1, 0, 0).cross(d);
    perp.normalize();
    const strength = Math.min(1, s / Math.max(n, 0.3));
    let offs, rads;
    if (n >= 2.5) { offs = [0, 0.17, -0.17]; rads = [0.055, 0.045, 0.045]; }
    else if (n >= 1.65) { offs = [0.085, -0.085]; rads = [0.06, 0.06]; }
    else if (n >= 1.3) { offs = [0.07, -0.08]; rads = [0.065, 0.03]; }
    else { offs = [0]; rads = [0.085]; }
    tmpQ.setFromUnitVectors(UP, d);
    for (let c = 0; c < offs.length; c++) {
      const rad = rads[c] * (0.35 + 0.65 * strength);
      tmpV.copy(a).addScaledVector(d, len / 2).addScaledVector(perp, offs[c]);
      tmpS.set(rad, len, rad);
      tmpM.compose(tmpV, tmpQ, tmpS);
      scene.bonds.setMatrixAt(nb, tmpM);
      const key = i < j ? i + '-' + j : j + '-' + i;
      if (!SB.bondBorn.has(key)) SB.bondBorn.set(key, performance.now());
      const age = performance.now() - SB.bondBorn.get(key);
      if (scene.bonds.instanceColor) {
        if (age < 700) tmpC.set('#3dde7a');
        else if (strength < 0.55) tmpC.set('#f0b429');
        else tmpC.set(cssVar('--scene-mode') === 'light' ? '#8a93a3' : '#d5dbe6');
        scene.bonds.setColorAt(nb, tmpC);
      }
      nb++;
    }
  }
  pruneBonds(bonds, NB);
  scene.bonds.count = nb;
  scene.bonds.instanceMatrix.needsUpdate = true;
  if (scene.bonds.instanceColor) scene.bonds.instanceColor.needsUpdate = true;
  // legami a idrogeno: linee tratteggiate H···A
  const hb = f.hbonds ?? [];
  const hpos = scene.hb.geometry.attributes.position;
  const nh = Math.min(hb.length / 2, 2048);
  for (let k = 0; k < nh; k++) {
    const h = hb[2 * k], ac = hb[2 * k + 1];
    hpos.setXYZ(2 * k, pos[3 * h], pos[3 * h + 1], pos[3 * h + 2]);
    hpos.setXYZ(2 * k + 1, pos[3 * ac], pos[3 * ac + 1], pos[3 * ac + 2]);
  }
  hpos.needsUpdate = true;
  scene.hb.geometry.setDrawRange(0, 2 * nh);
  scene.hb.computeLineDistances();
  scene.hb.visible = SB.style === 'ball' && nh > 0;
  // selezione
  if (SB.selected >= 0 && SB.selected < N) {
    scene.sel.visible = true;
    scene.sel.position.set(pos[3 * SB.selected], pos[3 * SB.selected + 1], pos[3 * SB.selected + 2]);
    scene.sel.scale.setScalar(atomRadius(Z[SB.selected]) * 1.6);
  } else scene.sel.visible = false;
  if (drag.active && drag.i < N) {
    const g = scene.grabLine.geometry.attributes.position;
    g.setXYZ(0, pos[3 * drag.i], pos[3 * drag.i + 1], pos[3 * drag.i + 2]);
    g.setXYZ(1, drag.target.x, drag.target.y, drag.target.z);
    g.needsUpdate = true;
    scene.grabLine.visible = true;
  } else scene.grabLine.visible = false;
  const now = performance.now();
  pumpCloud(now);
  if (now - SB.lastPanel > 200) { SB.lastPanel = now; renderPanelBody(); }
  if (now - SB.lastCharts > 500) { SB.lastCharts = now; drawCharts(); renderLiveAnalysis(); }
}

// ---------------------------------------------------------------------------
// Nuvola elettronica e orbitali: la funzione d'onda corrente arriva dal worker della dinamica, la griglia
// 3D si calcola in un secondo worker (così la dinamica non si ferma), le isosuperfici qui
// ---------------------------------------------------------------------------

const cloud = { worker: null, busy: false, last: 0, id: 0, group: null, pos: null };
const BOHR_A = 0.52917721090;

function clearCloud() {
  if (cloud.group) {
    cloud.group.traverse(o => { o.geometry?.dispose(); o.material?.dispose(); });
    scene.group?.remove(cloud.group);
    cloud.group = null;
  }
  SB.cloudNote = '';
  cloud.id++;
  cloud.busy = false;
}

function pumpCloud(now) {
  if (SB.style !== 'cloud' && SB.style !== 'orbital') { if (cloud.group) clearCloud(); return; }
  if (cloud.busy || now - cloud.last < 120 || !SB.frame?.N) return;
  cloud.busy = true;
  cloud.last = now;
  post({ type: 'wave', reqId: cloud.id });
}

function onWave(m) {
  if (m.reqId !== cloud.id || !active) return;
  const wantOrb = SB.style === 'orbital';
  const w = m.wave;
  if (wantOrb && !w) {
    cloud.busy = false;
    if (cloud.group) { scene.group.remove(cloud.group); cloud.group = null; }
    SB.cloudNote = '<b>Gli orbitali esistono solo nel calcolo quantistico:</b> scegli "quantistico MINDO/3" o "Hartree–Fock ab initio".';
    return;
  }
  if (!cloud.worker) {
    cloud.worker = new Worker(new URL('../chem/densityWorker.js', import.meta.url), { type: 'module' });
    cloud.worker.onmessage = (ev) => onGrid(ev.data);
  }
  const box = m.box;
  const half = box / 2 + 1.2;
  const res = Math.max(28, Math.min(72, Math.round(2 * half / 0.24) + 1));
  const orbital = wantOrb && w ? (SB.orbital === 'lumo' ? w.lumo : w.homo) : null;
  const req = { id: m.reqId, what: wantOrb ? 'orbital' : 'density', box, res };
  if (!w) Object.assign(req, { mode: 'promolecular', Z: m.Z, pos: m.pos });
  else if (w.kind === 'sto') Object.assign(req, { mode: 'sto', Z: w.Z, pos: w.pos, first: w.first, P: w.P, orb: orbital?.c });
  else Object.assign(req, { mode: 'gauss', Z: w.atoms.map(a => a.Z), pos: Float64Array.from(w.atoms.flatMap(a => a.xyz.map(v => v * BOHR_A))), atoms: w.atoms, basisName: w.basisName, P: w.P, orb: orbital?.c });
  req.surfaces = wantOrb
    ? [{ iso: 0.05, sign: 1 }, { iso: 0.05, sign: -1 }]
    // tre superfici di densità costante: confine di van der Waals (0,002 e/bohr³), regione dei legami, gusci interni
    : [{ iso: 0.002, sign: 1, colorByAtom: true }, { iso: 0.05, sign: 1, colorByAtom: true }, { iso: 0.3, sign: 1, colorByAtom: true }];
  if (!wantOrb) req.colors = Float32Array.from(req.Z.flatMap(z => { const c = new THREE.Color(cpkColor(z)); return [0.35 + 0.65 * c.r, 0.35 + 0.65 * c.g, 0.35 + 0.65 * c.b]; }));
  cloud.meta = { mode: req.mode, orbital, eUnit: w?.eUnit ?? 1, ff: m.forceField };
  cloud.worker.postMessage(req);
}

function onGrid(g) {
  cloud.busy = false;
  if (g.id !== cloud.id || !active || !scene.group) return;
  if (g.error) { SB.cloudNote = `Densità non calcolata: ${escapeHtml(g.error)}`; return; }
  const meta = cloud.meta;
  const grp = new THREE.Group();
  const OPAC = { 0.002: 0.13, 0.05: 0.24, 0.3: 0.5 };
  const pos = cssVar('--phase-pos'), neg = cssVar('--phase-neg');
  for (const sf of g.meshes) {
    if (!sf.positions.length) continue;
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(sf.positions, 3));
    geo.setAttribute('normal', new THREE.BufferAttribute(sf.normals, 3));
    if (sf.colors) geo.setAttribute('color', new THREE.BufferAttribute(sf.colors, 3));
    const orb = g.what !== 'density';
    const mat = new THREE.MeshStandardMaterial({
      color: sf.colors ? '#ffffff' : sf.sign > 0 ? pos : neg, vertexColors: !!sf.colors, roughness: 0.35, metalness: 0,
      transparent: true, opacity: orb ? 0.8 : OPAC[sf.iso] ?? 0.3, side: THREE.DoubleSide, depthWrite: false,
    });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.renderOrder = 2;
    grp.add(mesh);
  }
  if (g.what === 'density') {
    // tre superfici di densità costante: il confine di van der Waals (0,002 e/bohr³), la regione dei legami, i gusci interni
    SB.cloudNote = meta.mode === 'promolecular'
      ? 'Nuvola elettronica <b>promolecolare</b>: somma delle densità degli atomi isolati (calcolate con la DFT del sito). Con il campo classico gli elettroni non si ridistribuiscono nei legami: per la densità vera scegli il motore quantistico.'
      : `Densità elettronica ρ(r) calcolata dalla funzione d'onda ${meta.mode === 'sto' ? 'MINDO/3 (valenza) più il core atomico' : 'Hartree–Fock (tutti gli elettroni)'}: superfici a 0,002 e/bohr³ (confine di van der Waals), 0,05 (legami) e 0,3 (vicino ai nuclei).`;
  } else {
    const o = meta.orbital;
    SB.cloudNote = o ? `${SB.orbital === 'lumo' ? 'LUMO' : 'HOMO'} (spin ${o.spin}), ε = ${nf(o.e * meta.eUnit, 2)} eV: isosuperficie |ψ| = 0,05 bohr<sup>−3/2</sup>, in colore il segno della funzione d'onda.` : 'Orbitale non disponibile.';
  }
  if (cloud.group) { cloud.group.traverse(x => { x.geometry?.dispose(); x.material?.dispose(); }); scene.group.remove(cloud.group); }
  cloud.group = grp;
  if (g.what === 'density') renderLegend();
  scene.group.add(grp);
}

function flashUpdate() {
  const fl = SB.photonFlash;
  const f = SB.frame;
  if (!fl || !f) { scene.flash.visible = false; return; }
  const age = (performance.now() - fl.t) / 600;
  if (age > 1 || (fl.i >= 0 && fl.i >= f.N)) { scene.flash.visible = false; SB.photonFlash = null; return; }
  const p = f.pos;
  scene.flash.visible = true;
  if (fl.at) scene.flash.position.set(...fl.at);
  else scene.flash.position.set((p[3 * fl.i] + p[3 * fl.j]) / 2, (p[3 * fl.i + 1] + p[3 * fl.j + 1]) / 2, (p[3 * fl.i + 2] + p[3 * fl.j + 2]) / 2);
  scene.flash.scale.setScalar((fl.size ?? 1) * (0.4 + 1.6 * age));
  scene.flash.material.color.set(lambdaCss(fl.lambda));
  scene.flash.material.opacity = 0.7 * (1 - age);
}

// ---------------------------------------------------------------------------
// Interazione con il mouse
// ---------------------------------------------------------------------------

const drag = { active: false, i: -1, plane: new THREE.Plane(), target: new THREE.Vector3(), down: null };

function rayFrom(e) {
  const rect = viewer.renderer.domElement.getBoundingClientRect();
  const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
  const ray = new THREE.Raycaster();
  ray.setFromCamera(ndc, viewer.camera);
  return ray;
}

function pickAtom(e) {
  if (!scene.atoms || !SB.frame) return -1;
  const hits = rayFrom(e).intersectObject(scene.atoms, false);
  return hits.length ? hits[0].instanceId : -1;
}

/** Punto sul piano perpendicolare alla vista che passa per il centro della scatola. */
function pointInBox(e) {
  const ray = rayFrom(e).ray;
  const n = new THREE.Vector3().subVectors(viewer.camera.position, viewer.controls.target).normalize();
  const plane = new THREE.Plane().setFromNormalAndCoplanarPoint(n, new THREE.Vector3(0, 0, 0));
  const p = new THREE.Vector3();
  if (!ray.intersectPlane(plane, p)) return null;
  const h = (SB.frame?.stats.box ?? 20) / 2 - 1;
  p.clampScalar(-h, h);
  return p;
}

function onPointerDown(e) {
  if (!active) return;
  drag.down = [e.clientX, e.clientY];
  if (SB.tool !== 'grab') return;
  const i = pickAtom(e);
  if (i < 0) return;
  const f = SB.frame;
  const p = new THREE.Vector3(f.pos[3 * i], f.pos[3 * i + 1], f.pos[3 * i + 2]);
  const n = new THREE.Vector3().subVectors(viewer.camera.position, viewer.controls.target).normalize();
  drag.plane.setFromNormalAndCoplanarPoint(n, p);
  drag.active = true;
  drag.i = i;
  drag.target.copy(p);
  SB.selected = i;
  viewer.controls.enabled = false;
  viewer.renderer.domElement.setPointerCapture?.(e.pointerId);
  post({ type: 'grab', i, target: [p.x, p.y, p.z] });
}

function onPointerMove(e) {
  if (!active || !drag.active) return;
  const p = new THREE.Vector3();
  if (rayFrom(e).ray.intersectPlane(drag.plane, p)) {
    drag.target.copy(p);
    post({ type: 'grab', i: drag.i, target: [p.x, p.y, p.z] });
  }
}

function onPointerUp(e) {
  if (!active) return;
  if (drag.active) {
    drag.active = false;
    post({ type: 'grab', i: null });
    viewer.controls.enabled = true;
    renderLiveAnalysis();
    return;
  }
  if (!drag.down) return;
  const moved = Math.hypot(e.clientX - drag.down[0], e.clientY - drag.down[1]);
  drag.down = null;
  if (moved > 5) return;
  const i = pickAtom(e);
  if (SB.tool === 'select') {
    SB.selected = i;
    SB.fresh = true;
    renderLiveAnalysis();
  } else if (SB.tool === 'spark') {
    const p = i >= 0 ? [SB.frame.pos[3 * i], SB.frame.pos[3 * i + 1], SB.frame.pos[3 * i + 2]] : pointInBox(e)?.toArray();
    if (p) {
      post({ type: 'spark', center: p, radius: 4.5, T: 9000 });
      SB.photonFlash = { at: p, i: -1, j: -1, t: performance.now(), lambda: 600, size: 3 };
    }
  } else if (SB.tool === 'add') {
    const p = pointInBox(e);
    if (p) post({ type: 'add', ...addPayload(1), at: p.toArray() });
  } else if (SB.tool === 'delete' && i >= 0) {
    post({ type: 'remove', indices: [i] });
    SB.selected = -1;
  }
}

function addPayload(count) {
  if (SB.add.kind === 'atom') return { symbol: SB.add.symbol, count };
  if (SB.add.kind === 'smiles') return { smiles: SB.add.smiles, count };
  const m = [...SANDBOX_MOLECULES, ...SANDBOX_IONS].find(x => x.id === SB.add.id);
  return { smiles: m.smiles, count };
}

// ---------------------------------------------------------------------------
// Pannelli
// ---------------------------------------------------------------------------

function renderSide() {
  if (!active) return;
  const p = SB.preset;
  const molBtns = SANDBOX_MOLECULES.map(m => {
    const f = hillFormula(parseSmiles(m.smiles).atoms);
    return `<button type="button" class="mol-chip ${SB.add.kind === 'mol' && SB.add.id === m.id ? 'active' : ''}" data-mol="${m.id}" title="${m.name}">${f}</button>`;
  }).join('');
  const ionBtns = SANDBOX_IONS.map(m => {
    const g = parseSmiles(m.smiles);
    const f = hillFormula(g.atoms, g.atoms.reduce((a, b) => a + b.charge, 0));
    return `<button type="button" class="mol-chip ${SB.add.kind === 'mol' && SB.add.id === m.id ? 'active' : ''}" data-mol="${m.id}" title="${m.name}">${f}</button>`;
  }).join('');
  const elBtns = SANDBOX_ELEMENTS.map(Z => `<button type="button" class="mol-chip ${SB.add.kind === 'atom' && SB.add.symbol === sym(Z) ? 'active' : ''}" data-el="${sym(Z)}" title="${ELEMENTS[Z - 1].name}">${sym(Z)}</button>`).join('');
  $('element-card').innerHTML = `<div class="btn-row"><a class="btn" href="#fasi">Fasi · dinamica MD</a><a class="btn" href="#becher">Soluzioni · cristalli</a></div>
    <h3 class="side-h">Sandbox chimica</h3>
    <label class="lbl" for="sb-preset">Esperimento</label>
    <select id="sb-preset">
      <optgroup label="Quantistica: elettroni calcolati a ogni passo (MINDO/3)">${PRESETS.filter(x => x.quantum).map(x => `<option value="${x.id}" ${x.id === p.id ? 'selected' : ''}>${x.name}</option>`).join('')}</optgroup>
      <optgroup label="Campo di forze reattivo (classico, veloce)">${PRESETS.filter(x => !x.atoms && !x.quantum).map(x => `<option value="${x.id}" ${x.id === p.id ? 'selected' : ''}>${x.name}</option>`).join('')}</optgroup>
      <optgroup label="Ab initio (Hartree–Fock)">${PRESETS.filter(x => x.atoms && !x.quantum).map(x => `<option value="${x.id}" ${x.id === p.id ? 'selected' : ''}>${x.name}</option>`).join('')}</optgroup>
    </select>
    <details><summary>Obiettivo e istruzioni</summary><p class="mol-note">${p.text}</p>
    ${p.tips?.length ? `<ul class="sb-tips">${p.tips.map(t => `<li>${t}</li>`).join('')}</ul>` : ''}</details>
    <div class="btn-row"><button type="button" class="btn" id="sb-restart">Ricomincia</button><button type="button" class="btn" id="sb-clear">Svuota la scatola</button></div>
    <label class="lbl">Versa nella scatola, senza cancellare</label>
    <div class="btn-row">
      <button type="button" class="btn" id="sb-pour-fire">H₂ + O₂</button>
      <button type="button" class="btn" id="sb-pour-water">acqua</button>
      <button type="button" class="btn" id="sb-pour-ar">argon</button>
      <button type="button" class="btn" id="sb-pour-salt">Na + Cl₂</button>
    </div>
    <label class="lbl">Aggiungi molecole</label>
    <div class="mol-list">${molBtns}</div>
    <label class="lbl">Aggiungi ioni (calcolati con MINDO/3 o Hartree–Fock)</label>
    <div class="mol-list">${ionBtns}</div>
    <label class="lbl">Aggiungi atomi (radicali, gas nobili)</label>
    <div class="mol-list">${elBtns}</div>
    <label class="lbl" for="sb-smiles">Oppure una molecola da SMILES</label>
    <div class="smiles-row"><input id="sb-smiles" placeholder="es. CC(=O)O oppure [NH4+]" value="${SB.add.kind === 'smiles' ? escapeHtml(SB.add.smiles) : ''}" spellcheck="false"><button type="button" class="btn" id="sb-smiles-ok">Usa</button></div>
    <div class="ctl" style="margin-top:8px"><label class="lbl" for="sb-count">Quantità: <span id="sb-count-out">${SB.count}</span></label>
      <div class="range-row"><input type="range" id="sb-count" min="1" max="40" step="1" value="${SB.count}"><button type="button" class="btn" id="sb-add">Aggiungi ${escapeHtml(addLabel())}</button></div></div>
    <p class="hint">Con lo strumento <b>Aggiungi</b> (a destra) puoi anche cliccare nella scatola per mettere una molecola dove vuoi.</p>
    ${SB.info ? `<p class="hint warn">${escapeHtml(SB.info)}</p>` : ''}
    <label class="lbl" for="sb-save-name">I tuoi esperimenti (salvati in questo browser)</label>
    <div class="smiles-row"><input id="sb-save-name" placeholder="nome dell'esperimento" spellcheck="false"><button type="button" class="btn" id="sb-save">Salva</button></div>
    ${savedList().length ? `<ul class="sb-saves">${savedList().map((x, k) => `<li><button type="button" class="linkish" data-load="${k}" title="Riprendi da questo stato">${escapeHtml(x.name)}</button> <span class="desc-muted">${x.data.Z.length} atomi · ${new Date(x.date).toLocaleString('it-IT', { dateStyle: 'short', timeStyle: 'short' })}</span> <button type="button" class="linkish" data-del="${k}" title="Elimina">✕</button></li>`).join('')}</ul>` : ''}
    <div class="btn-row"><button type="button" class="btn" id="sb-xyz" title="Coordinate atomiche nel formato XYZ, leggibile da Avogadro, VMD, Jmol…">Esporta .xyz</button></div>`;
  $('sb-save').addEventListener('click', () => post({ type: 'snapshot', name: ($('sb-save-name').value.trim() || `esperimento ${savedList().length + 1}`).slice(0, 60) }));
  $('sb-xyz').addEventListener('click', () => post({ type: 'snapshot', name: '\u0000xyz' }));
  $('element-card').querySelectorAll('[data-load]').forEach(b => b.addEventListener('click', () => {
    const x = savedList()[+b.dataset.load];
    if (!x) return;
    SB.events = []; SB.history = []; SB.census = null; SB.selected = -1; SB.phys = null;
    SB.preset = { ...SB.preset, name: x.name, text: `Esperimento salvato il ${new Date(x.date).toLocaleString('it-IT')}.`, tips: [] };
    post({ type: 'restore', data: x.data });
    post({ type: 'set', paused: false, T: x.data.T });
    SB.userPaused = false; scene.framed = false;
    renderSide(); renderControls(); renderAnalysis();
  }));
  $('element-card').querySelectorAll('[data-del]').forEach(b => b.addEventListener('click', () => {
    const list = savedList();
    list.splice(+b.dataset.del, 1);
    writeSaved(list);
    renderSide();
  }));
  $('sb-preset').addEventListener('change', (e) => loadPreset(PRESETS.find(x => x.id === e.target.value)));
  $('sb-restart').addEventListener('click', () => loadPreset(SB.preset));
  $('sb-clear').addEventListener('click', () => { post({ type: 'clear' }); SB.events = []; SB.history = []; SB.census = null; SB.selected = -1; SB.userPaused = false; post({ type: 'set', paused: false, forceField: 'auto' }); renderAnalysis(); });
  $('sb-pour-fire').addEventListener('click', () => pour([['[H][H]', 6], ['O=O', 3]]));
  $('sb-pour-water').addEventListener('click', () => pour([['O', 8]]));
  $('sb-pour-ar').addEventListener('click', () => pour([['[Ar]', 40]]));
  $('sb-pour-salt').addEventListener('click', () => pour([['[Na]', 8], ['ClCl', 4]]));
  $('element-card').querySelectorAll('[data-mol]').forEach(b => b.addEventListener('click', () => { SB.add = { kind: 'mol', id: b.dataset.mol }; renderSide(); }));
  $('element-card').querySelectorAll('[data-el]').forEach(b => b.addEventListener('click', () => { SB.add = { kind: 'atom', symbol: b.dataset.el }; renderSide(); }));
  $('sb-smiles-ok').addEventListener('click', () => {
    const s = $('sb-smiles').value.trim();
    try {
      const g = parseSmiles(s);
      if (g.atoms.some(a => !ATOM_PARAMS[a.Z])) throw new Error('elemento non disponibile nella sandbox');
      SB.add = { kind: 'smiles', smiles: s };
      SB.info = '';
    } catch (e) { SB.info = `SMILES non valido: ${e.message}`; }
    renderSide();
  });
  $('sb-count').addEventListener('input', (e) => { SB.count = +e.target.value; $('sb-count-out').textContent = SB.count; });
  $('sb-add').addEventListener('click', () => post({ type: 'add', ...addPayload(SB.count) }));
}

const SAVE_KEY = 'atlante-sandbox-esperimenti';
function savedList() {
  try { const v = JSON.parse(localStorage.getItem(SAVE_KEY) ?? '[]'); return Array.isArray(v) ? v : []; } catch { return []; }
}
function writeSaved(list) {
  try { localStorage.setItem(SAVE_KEY, JSON.stringify(list)); return true; } catch { return false; }
}

/** Risposta del worker a "snapshot": salva nel browser oppure esporta le coordinate in formato XYZ. */
function saveSnapshot(name, data) {
  if (name === '\u0000xyz') {
    const lines = [String(data.Z.length), `Atlante Orbitale, sandbox: scatola ${data.box.toFixed(2)} A, T = ${data.T.toFixed(1)} K`];
    data.Z.forEach((z, i) => lines.push(`${sym(z).padEnd(2)} ${data.pos.slice(3 * i, 3 * i + 3).map(v => v.toFixed(5).padStart(11)).join(' ')}`));
    const url = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], { type: 'chemical/x-xyz' }));
    const a = document.createElement('a');
    a.href = url; a.download = 'sandbox.xyz';
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    return;
  }
  const r = (v) => Math.round(v * 1e5) / 1e5;
  const compact = { ...data, pos: data.pos.map(r), vel: data.vel.map(v => Math.round(v * 1e7) / 1e7) };
  const list = savedList();
  list.unshift({ name, date: Date.now(), data: compact });
  if (list.length > 12) list.length = 12;
  SB.info = writeSaved(list) ? `Salvato "${name}".` : 'Impossibile salvare: la memoria del browser è piena o non disponibile.';
  renderSide();
}

function pour(items) {
  SB.userPaused = false;
  post({ type: 'set', paused: false });
  for (const [smiles, count] of items) post({ type: 'add', smiles, count });
}

function addLabel() {
  if (SB.add.kind === 'atom') return SB.add.symbol;
  if (SB.add.kind === 'smiles') return SB.add.smiles;
  const m = [...SANDBOX_MOLECULES, ...SANDBOX_IONS].find(x => x.id === SB.add.id);
  return hillFormula(parseSmiles(m.smiles).atoms);
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

const TOOLS = [
  { id: 'select', label: 'Seleziona', title: 'Clic su un atomo: distanze, ordini di legame, carica parziale. Trascina per ruotare la vista.' },
  { id: 'grab', label: 'Afferra', title: 'Trascina un atomo con una "pinzetta" elastica: puoi rompere legami o spingere molecole a reagire.' },
  { id: 'spark', label: 'Scintilla', title: 'Clic: gli atomi vicini vengono portati a 9000 K per un istante.' },
  { id: 'add', label: 'Aggiungi', title: 'Clic nella scatola: aggiunge la molecola scelta a sinistra.' },
  { id: 'delete', label: 'Elimina', title: 'Clic su un atomo per toglierlo.' },
];

function logSlider(id, label, min, max, value, fmt) {
  return `<div class="ctl"><label class="lbl" for="${id}">${label}</label>
    <div class="range-row"><input type="range" id="${id}" min="${Math.log10(min)}" max="${Math.log10(max)}" step="0.01" value="${Math.log10(value)}"><output id="${id}-out">${fmt(value)}</output></div></div>`;
}

function renderControls() {
  if (!active) return;
  const st = SB.frame?.stats;
  const T = st?.Ttarget ?? SB.preset.T;
  const box = st?.box ?? SB.preset.box;
  const thermo = st?.thermostat ?? true;
  const spf = st?.stepsPerFrame ?? 40;
  const lightOn = st?.light?.on ?? false;
  const photonKJ = 119626.566 / SB.lambda;
  const baro = SB.phys?.barostat ?? { on: false, P0: 1 };
  $('controls').innerHTML = `
    <div class="ctl"><span class="lbl">Simulazione</span>
      <div class="btn-row" style="margin-top:0">
        <button type="button" class="btn" id="sb-play">${SB.userPaused ? '▶ Avvia' : '❚❚ Pausa'}</button>
        <button type="button" class="btn" id="sb-step" ${SB.userPaused ? '' : 'disabled'}>1 passo</button>
        <button type="button" class="btn" id="sb-phases" title="Argon denso: caldo, poi ambiente, poi freddo">Giro delle fasi</button>
      </div></div>
    <div class="ctl"><label class="lbl" for="sb-dt">Passo Δt (fs)</label>
      <input id="sb-dt" type="number" min="0.001" max="2" step="0.001" value="${st?.dt ?? SB.preset.dt ?? 0.2}"></div>
    <div class="btn-row"><button class="btn" id="sb-zero">Azzera misure</button><button class="btn" id="sb-export">Esporta CSV</button></div>
    <div class="ctl"><span class="lbl">Modello delle forze</span></div>
    <div class="seg" id="sb-ff">
      <button type="button" data-v="auto" aria-pressed="${(st?.fidelity ?? SB.fidelity) === 'auto'}" title="Sceglie il modello più fedele che può ancora girare: LJ pubblicato per i nobili, MINDO/3 fino a 36 atomi, Hartree–Fock fino a 6, altrimenti classico">fedeltà automatica</button>
      <button type="button" data-v="reactive" aria-pressed="${(st?.fidelity ?? SB.fidelity) === 'reactive'}" title="Potenziale empirico del progetto: non parametrizzato per prevedere reazioni generali">classico</button>
      <button type="button" data-v="lj" aria-pressed="${(st?.fidelity ?? SB.fidelity) === 'lj'}" title="Lennard–Jones con σ e ε pubblicati. Solo He, Ne, Ar, Kr, Xe">LJ nobili</button>
      <button type="button" data-v="mindo3" aria-pressed="${(st?.fidelity ?? SB.fidelity) === 'mindo3'}" title="MINDO/3 a ogni passo, fino a 90 atomi di H, B, C, N, O, F, P, S, Cl">MINDO/3</button>
      <button type="button" data-v="hf" aria-pressed="${(st?.fidelity ?? SB.fidelity) === 'hf'}" title="Hartree–Fock/STO-3G a ogni passo, fino a 8 atomi">Hartree–Fock</button>
    </div>
    <p class="hint">${escapeHtml(st?.modelWhy ?? SB.info ?? '')}</p>
    <div class="ctl"><span class="lbl">Carica e spin · carica totale ${chargeText(st?.charge ?? 0)}</span>
      <div class="seg" id="sb-mult">${[[0, 'spin libero'], [1, 'singoletto'], [2, 'doppietto'], [3, 'tripletto'], [4, 'quartetto']].map(([v, t]) => `<button type="button" data-v="${v}" aria-pressed="${(st?.multiplicity ?? 0) === v || (v === 0 && !st?.multiplicity)}" title="${v ? `Molteplicità 2S+1 = ${v} fissata` : 'Lo spin esce dal calcolo (livello di Fermi comune per i due spin)'}">${t}</button>`).join('')}</div>
      <p class="hint">Le specie cariche si aggiungono dalla lista degli ioni o da SMILES con la carica, per esempio [NH4+] o C[O-]. La carica totale entra nella funzione d'onda; la carica di ogni molecola si legge dal calcolo.</p></div>
    <div class="ctl"><span class="lbl">Campo elettrico uniforme (V/Å) ${st?.field ? `· acceso, |E| = ${nf(Math.hypot(...st.field), 2)} V/Å` : '· spento'}</span>
      <div class="range-row field-row">${['x', 'y', 'z'].map((c, k) => `<label for="sb-E${c}">E<sub>${c}</sub></label><input id="sb-E${c}" type="text" inputmode="decimal" value="${nf(st?.field?.[k] ?? 0, 2)}" style="width:58px">`).join('')}</div>
      <div class="btn-row" style="margin-top:4px"><button type="button" class="btn" id="sb-field-on">Applica il campo</button><button type="button" class="btn" id="sb-field-off" ${st?.field ? '' : 'disabled'}>Spegni</button></div>
      <p class="hint">Calcolato con MINDO/3 nella funzione d'onda: polarizza le molecole, orienta i dipoli e spinge gli ioni (forza qE). 1 V/Å = 10¹⁰ V/m, come vicino alla punta di un microscopio a effetto di campo. Accendere o cambiare il campo è lavoro sul sistema, contato nel bilancio. Con base minima la ionizzazione per effetto tunnel non è descritta.</p></div>
    ${logSlider('sb-T', `Temperatura del termostato`, 10, 8000, T, v => `${nf(v, 0)} K`)}
    <div class="seg" id="sb-thermo">
      <button type="button" data-v="1" aria-pressed="${thermo}" title="Termostato di Bussi: scambia calore con un bagno a temperatura costante (insieme canonico NVT)">Termostato</button>
      <button type="button" data-v="0" aria-pressed="${!thermo}" title="Nessuno scambio di calore: l'energia totale si conserva (NVE). Le reazioni esotermiche scaldano il sistema.">Isolato</button>
    </div>
    <div class="ctl"><label class="lbl" for="sb-box">Lato della scatola (volume)</label>
      <div class="range-row"><input type="range" id="sb-box" min="10" max="60" step="0.5" value="${box}"><output id="sb-box-out">${nf(box, 1)} Å</output></div></div>
    <div class="ctl"><span class="lbl">Pressione</span>
      <div class="seg" id="sb-baro">
        <button type="button" data-v="0" aria-pressed="${!baro.on}" title="Volume fisso (NVT o NVE)">volume fisso</button>
        <button type="button" data-v="1" aria-pressed="${baro.on}" title="Il volume si adatta per mantenere la pressione scelta (barostato di Berendsen)">pressione costante</button>
      </div>
      <div class="range-row" style="margin-top:4px"><label for="sb-p0">P₀ (bar)</label><input id="sb-p0" type="text" inputmode="decimal" value="${baro.P0}" style="width:80px"></div>
      <p class="hint">Barostato di Berendsen (τ<sub>P</sub> = 10 ps): la scatola si allarga o si stringe finché la pressione sulle pareti vale P₀.</p></div>
    <div class="ctl"><span class="lbl">Analisi</span>
      <div class="range-row"><label for="sb-rdf">g(r) fra</label><select id="sb-rdf">${rdfOptions()}</select><button type="button" class="btn" id="sb-reset-an" title="Azzera g(r), spostamento quadratico medio e fluttuazioni di energia">Azzera</button></div></div>
    <details><summary>Interazioni e aspetto</summary>
    <div class="ctl"><label class="lbl" for="sb-speed">Velocità (passi per fotogramma)</label>
      <div class="range-row"><input type="range" id="sb-speed" min="1" max="200" step="1" value="${spf}"><output id="sb-speed-out">${spf}</output></div></div>
    <div class="ctl"><span class="lbl">Strumento</span></div>
    <div class="seg" id="sb-tool">${TOOLS.map(t => `<button type="button" data-v="${t.id}" aria-pressed="${SB.tool === t.id}" title="${t.title}">${t.label}</button>`).join('')}</div>
    <p class="hint" id="sb-tool-hint">${TOOLS.find(t => t.id === SB.tool).title}</p>
    <div class="ctl"><label class="lbl" for="sb-lambda">Impulso equivalente: λ = <span id="sb-lambda-out">${SB.lambda} nm</span> <i class="swatch" id="sb-lambda-sw" style="background: ${lambdaCss(SB.lambda)}"></i></label>
      <div class="range-row"><input type="range" id="sb-lambda" min="100" max="800" step="5" value="${SB.lambda}"><output id="sb-ephoton">${nf(photonKJ, 0)} kJ/mol</output></div>
      <div class="btn-row" style="margin-top:4px"><button type="button" class="btn" id="sb-light" aria-pressed="${lightOn}">${lightOn ? 'Ferma impulsi' : 'Impulsi ripetuti'}</button><button type="button" class="btn" id="sb-flash">Deposita hc/λ</button></div>
      <p class="hint">Deposito meccanico di hc/λ. Assorbimento, stati eccitati e rese fotochimiche non sono calcolati.</p></div>
    <div class="ctl"><span class="lbl">Colore degli atomi</span></div>
    <div class="seg" id="sb-color">
      <button type="button" data-v="element" aria-pressed="${SB.color === 'element'}">elemento</button>
      <button type="button" data-v="charge" aria-pressed="${SB.color === 'charge'}">carica parziale</button>
      <button type="button" data-v="ke" aria-pressed="${SB.color === 'ke'}">energia cinetica</button>
      <button type="button" data-v="phase" aria-pressed="${SB.color === 'phase'}">fase locale</button>
    </div>
    <div class="seg" id="sb-style" style="margin-top:4px">
      <button type="button" data-v="ball" aria-pressed="${SB.style === 'ball'}">sfere e bastoncini</button>
      <button type="button" data-v="vdw" aria-pressed="${SB.style === 'vdw'}">van der Waals</button>
      <button type="button" data-v="cloud" aria-pressed="${SB.style === 'cloud'}" title="La densità degli elettroni ρ(r): come appaiono davvero atomi e molecole">nuvola elettronica</button>
      <button type="button" data-v="orbital" aria-pressed="${SB.style === 'orbital'}" title="Orbitali di frontiera calcolati: HOMO (l'elettrone più alto) e LUMO (il primo posto libero)">orbitali</button>
    </div>
    ${SB.style === 'orbital' ? `<div class="seg" id="sb-orb" style="margin-top:4px">
      <button type="button" data-v="homo" aria-pressed="${SB.orbital === 'homo'}">HOMO</button>
      <button type="button" data-v="lumo" aria-pressed="${SB.orbital === 'lumo'}">LUMO</button>
    </div>` : ''}
    </details>`;
  const seg = (id, fn) => $(id).querySelectorAll('button').forEach(b => b.addEventListener('click', () => {
    fn(b.dataset.v);
    $(id).querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', String(x === b)));
  }));
  $('sb-dt').addEventListener('change', e => { if (e.target.checkValidity()) post({ type: 'set', dt: +e.target.value }); });
  $('sb-zero').addEventListener('click', () => { SB.history = []; SB.events = []; post({ type: 'reset-measurements' }); });
  $('sb-export').addEventListener('click', () => {
    const rows = ['time_fs,temperature_K,total_eV,potential_eV,conserved_eV,pressure_bar',
      ...SB.history.map(h => [h.t, h.T, h.E, h.Ep, h.Ec, h.P].join(','))];
    const url = URL.createObjectURL(new Blob([rows.join('\n')], { type: 'text/csv;charset=utf-8' }));
    const a = document.createElement('a'); a.href = url; a.download = 'atoms-sandbox.csv'; a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
  $('sb-play').addEventListener('click', () => { SB.userPaused = !SB.userPaused; post({ type: 'set', paused: SB.userPaused }); renderControls(); });
  $('sb-phases').addEventListener('click', startPhaseTour);
  $('sb-step').addEventListener('click', () => post({ type: 'step', n: 1 }));
  $('sb-T').addEventListener('input', (e) => { const v = Math.pow(10, +e.target.value); $('sb-T-out').textContent = `${nf(v, 0)} K`; post({ type: 'set', T: v }); });
  seg('sb-thermo', (v) => post({ type: 'set', thermostat: v === '1' }));
  seg('sb-ff', (v) => post({ type: 'set', forceField: v }));
  seg('sb-mult', (v) => post({ type: 'set', multiplicity: +v || null }));
  $('sb-field-on').addEventListener('click', () => {
    const f = ['x', 'y', 'z'].map(c => parseFloat(String($(`sb-E${c}`).value).replace(',', '.')) || 0);
    post({ type: 'set', field: f });
    setTimeout(renderControls, 150);
  });
  $('sb-field-off').addEventListener('click', () => { post({ type: 'set', field: null }); setTimeout(renderControls, 150); });
  $('sb-box').addEventListener('input', (e) => { const v = +e.target.value; $('sb-box-out').textContent = `${nf(v, 1)} Å`; post({ type: 'set', box: v }); });
  seg('sb-baro', (v) => post({ type: 'set', barostat: { on: v === '1' } }));
  $('sb-p0').addEventListener('change', (e) => { const v = parseFloat(String(e.target.value).replace(',', '.')); if (Number.isFinite(v) && v > 0) post({ type: 'set', barostat: { P0: v } }); });
  $('sb-rdf').addEventListener('change', (e) => post({ type: 'set', rdfPair: e.target.value.split('-').map(Number) }));
  $('sb-reset-an').addEventListener('click', () => post({ type: 'resetAnalysis' }));
  $('sb-speed').addEventListener('input', (e) => { $('sb-speed-out').textContent = e.target.value; post({ type: 'set', stepsPerFrame: +e.target.value }); });
  seg('sb-tool', (v) => { SB.tool = v; $('sb-tool-hint').textContent = TOOLS.find(t => t.id === v).title; });
  $('sb-lambda').addEventListener('input', (e) => {
    SB.lambda = +e.target.value;
    $('sb-lambda-out').textContent = `${SB.lambda} nm`;
    $('sb-ephoton').textContent = `${nf(119626.566 / SB.lambda, 0)} kJ/mol`;
    $('sb-lambda-sw').style.background = lambdaCss(SB.lambda);
    post({ type: 'set', light: { lambda: SB.lambda } });
  });
  $('sb-light').addEventListener('click', () => { post({ type: 'set', light: { on: !lightOn, lambda: SB.lambda } }); setTimeout(renderControls, 80); });
  $('sb-flash').addEventListener('click', () => post({ type: 'photon', lambda: SB.lambda }));
  seg('sb-color', (v) => { SB.color = v; SB.fresh = true; renderLegend(); });
  seg('sb-style', (v) => { SB.style = v; SB.fresh = true; clearCloud(); renderControls(); renderLegend(); });
  if ($('sb-orb')) seg('sb-orb', (v) => { SB.orbital = v; clearCloud(); });
  renderLegend();
}

function rdfOptions() {
  const els = SB.phys?.elements ?? [];
  const cur = SB.phys?.rdf?.pair?.join('-');
  const opts = [];
  for (let a = 0; a < els.length; a++) for (let b = a; b < els.length; b++) {
    const v = `${els[a]}-${els[b]}`;
    opts.push(`<option value="${v}" ${v === cur || `${els[b]}-${els[a]}` === cur ? 'selected' : ''}>${sym(els[a])}–${sym(els[b])}</option>`);
  }
  return opts.join('') || '<option>—</option>';
}

/** Colore di una lunghezza d'onda; ultravioletto e infrarosso con colori convenzionali. */
function lambdaCss(nm) {
  if (nm < 380) return '#a78bfa';
  if (nm > 750) return '#8b1e1e';
  return wavelengthColor(nm);
}

function renderLegend() {
  const leg = $('viewport-legend');
  if (SB.style === 'cloud') {
    const els = [...new Set(Array.from(SB.frame?.Z ?? []))].sort((a, b) => a - b);
    leg.innerHTML = `${els.map(z => `<span><i class="swatch" style="background:${cpkColor(z)}"></i>${sym(z)}</span>`).join(' ')} · punti = nuclei; la nuvola prende il colore dell'atomo più vicino`;
    return;
  }
  if (SB.style === 'orbital') { leg.innerHTML = `<span><i class="swatch" style="background:${cssVar('--phase-pos')}"></i>ψ > 0</span> <span><i class="swatch" style="background:${cssVar('--phase-neg')}"></i>ψ < 0</span> · punti = nuclei`; return; }
  if (SB.color === 'charge') leg.innerHTML = '<span><i class="swatch" style="background:#e0402a"></i>δ− (negativa)</span> <span><i class="swatch" style="background:#dcdce0"></i>neutra</span> <span><i class="swatch" style="background:#3f7fe8"></i>δ+ (positiva)</span>';
  else if (SB.color === 'ke') leg.innerHTML = '<span><i class="swatch" style="background:#4070ff"></i>lento</span> <span><i class="swatch" style="background:#ffffff"></i>≈ 3/2 kT</span> <span><i class="swatch" style="background:#ff3020"></i>veloce</span>';
  else if (SB.color === 'phase') leg.innerHTML = '<span><i class="swatch" style="background:#9fd7ff"></i>isolato · gas</span> <span><i class="swatch" style="background:#f0b429"></i>aggregato</span> <span><i class="swatch" style="background:#7d8cff"></i>coordinato</span>';
  else leg.innerHTML = 'Colori CPK · verde = legame appena nato · ambra = legame che si forma · tratteggio = legame a idrogeno';
}

function renderPanelBody() {
  if (!active || !SB.frame) return;
  const s = SB.frame.stats;
  const N = SB.frame.N;
  const V = s.box ** 3;                                   // Å³
  const nMol = s.nMol || 0;
  const Pid = nMol * KB_EV * s.T / V * 1.602176634e6;     // bar
  const Z = Pid > 0 ? s.P / Pid : NaN;
  const kj = (e) => nf(e * KJ_PER_EV, 1);
  const ph = s.phase;
  $('viewport-title').innerHTML = `${SB.preset.name}<small>t = ${nf(s.t / 1000, 2)} ps · ${N} atomi · ${nMol} molecole${ph ? ` · ${ph.title}` : ''}</small>`;
  $('viewport-note').innerHTML = `${SB.cloudNote ? `${SB.cloudNote} ` : ''}${s.forceField === 'hf' ? 'Dinamica ab initio (Hartree–Fock a ogni passo)' : s.forceField === 'mindo3' ? 'Dinamica quantistica (SCF MINDO/3 a ogni passo)' : s.forceField === 'lj' ? 'Lennard–Jones dei gas nobili' : 'Potenziale classico qualitativo'}: passo Δt = ${nf(s.dt, 2)} fs, ${s.paused ? '<b>in pausa</b>' : `${s.stepsPerFrame} passi per fotogramma`}. Trascina per ruotare, rotellina per ingrandire.`;
  $('live-metrics').innerHTML = `<span>T cinetica <b>${nf(s.T, 0)} K</b></span><span>Δt <b>${nf(s.dt, 3)} fs</b></span><span>Deriva energetica <b>${sgn((s.diagnostics?.drift ?? 0) * KJ_PER_EV, 4)} kJ/mol</b></span><span>Modello <b>${s.forceField === 'hf' ? 'UHF / STO-3G' : s.forceField === 'mindo3' ? 'MINDO/3 quantistico' : s.forceField === 'lj' ? 'Lennard–Jones' : 'classico qualitativo'}</b></span>`;
  const dtIn = $('sb-dt');
  if (dtIn && document.activeElement !== dtIn && Math.abs(+dtIn.value - s.dt) > 1e-9) dtIn.value = +s.dt.toPrecision(4);
  $('panel-body').innerHTML = `
    <div><h3>Stato termodinamico</h3>
    <dl class="info-list">
      <dt>Temperatura istantanea</dt><dd>${nf(s.T, 0)} K</dd>
      <dt>Temperatura del bagno</dt><dd>${s.thermostat ? `${nf(s.Ttarget, 0)} K` : 'isolato (NVE)'}</dd>
      <dt>Pressione sulle pareti</dt><dd>${nf(s.P, s.P < 100 ? 1 : 0)} bar</dd>
      <dt>Gas ideale nkT/V</dt><dd>${nf(Pid, Pid < 100 ? 1 : 0)} bar</dd>
      <dt>Z = PV/nRT</dt><dd>${Number.isFinite(Z) ? nf(Z, 2) : '—'}</dd>
      <dt>Volume</dt><dd>${nf(V / 1000, 2)} nm³</dd>
      <dt>Densità numerica</dt><dd>${nf(nMol / V * 1e27 / NA, 2)} mol/L</dd>
    </dl></div>
    ${s.forceField === 'mindo3' && s.hf?.method ? `<div><h3>Struttura elettronica a ogni passo</h3>
    <dl class="info-list">
      <dt>Metodo</dt><dd>MINDO/3 (UHF, SCF)</dd>
      <dt>Orbitali di valenza</dt><dd>${s.hf.nbf}</dd>
      <dt>Spin totale S<sub>z</sub></dt><dd>${nf(Math.abs(s.hf.Sz), 1)}${Math.abs(s.hf.Sz) > 0.25 ? ' (elettroni spaiati)' : ''}</dd>
      <dt>HOMO / LUMO</dt><dd>${nf(s.hf.homo, 2)} / ${nf(s.hf.lumo, 2)} eV</dd>
      <dt>Gap HOMO–LUMO</dt><dd>${nf(s.hf.gap, 2)} eV</dd>
      <dt>Carica totale</dt><dd>${chargeText(s.charge ?? 0)}</dd>
      ${s.dipole ? `<dt>Momento di dipolo</dt><dd>${nf(Math.hypot(...s.dipole) * 4.80320, 2)} D${s.charge ? ' (dipende dall\'origine: c\'è una carica netta)' : ''}</dd>` : ''}
      ${s.field ? `<dt>Campo elettrico</dt><dd>(${s.field.map(v => nf(v, 2)).join('; ')}) V/Å</dd>` : ''}
      <dt>Calore di formazione ΔfH</dt><dd>${nf(s.hf.Hf * 4.184, 0)} kJ/mol</dd>
      <dt>Iterazioni SCF</dt><dd>${s.hf.iterations}${s.hf.converged ? '' : ' (non convergente)'}</dd>
    </dl>
    <p class="hint">A ogni passo si risolvono le equazioni di Roothaan–Hall per tutti gli elettroni di valenza; le forze sono il gradiente analitico dell'energia. Lo spin non è imposto: radicali e O₂ tripletto escono dal calcolo. Cariche e ordini di legame dalla matrice densità. ΔfH è la somma dei calori di formazione di tutte le molecole della scatola.</p></div>` : ''}
    ${s.forceField === 'hf' && s.hf ? `<div><h3>Calcolo quantistico a ogni passo</h3>
    <dl class="info-list">
      <dt>Metodo</dt><dd>UHF/${s.hf.basis}</dd>
      <dt>Funzioni di base</dt><dd>${s.hf.nbf}</dd>
      <dt>Molteplicità di spin 2S+1</dt><dd>${s.hf.multiplicity} (${['', 'singoletto', 'doppietto', 'tripletto', 'quartetto'][s.hf.multiplicity] ?? ''})</dd>
      <dt>⟨S²⟩ (esatto ${nf((s.hf.multiplicity - 1) / 2 * ((s.hf.multiplicity - 1) / 2 + 1), 2)})</dt><dd>${nf(s.hf.S2, 3)}</dd>
      <dt>Iterazioni SCF</dt><dd>${s.hf.iterations}${s.hf.converged ? '' : ' (non convergente)'}</dd>
    </dl>
    <p class="hint">Cariche di Mulliken e ordini di legame di Mayer dalla funzione d'onda. Energia totale elettronica + nucleare.</p></div>` : ''}
    <div><h3>Energia (kJ/mol di scatola)</h3>
    <dl class="info-list">
      <dt>Cinetica</dt><dd>${kj(s.Ekin)}</dd>
      <dt>Potenziale</dt><dd>${kj(s.Epot + s.Ewall)}</dd>
      ${s.parts && s.parts.hf !== undefined ? `<dt class="sub">· energia Hartree–Fock</dt><dd>${kj(s.parts.hf)}</dd>` : ''}
      ${s.parts && s.parts.bond !== undefined ? `<dt class="sub">· legami covalenti</dt><dd>${kj(s.parts.bond)}</dd>
      <dt class="sub">· angoli</dt><dd>${kj(s.parts.angle)}</dd>
      <dt class="sub">· van der Waals</dt><dd>${kj(s.parts.vdw)}</dd>
      <dt class="sub">· elettrostatica</dt><dd>${kj(s.parts.es)}</dd>
      <dt class="sub">· legami a idrogeno</dt><dd>${kj(s.parts.hbond ?? 0)}</dd>` : ''}
      <dt>Totale</dt><dd>${kj(s.Etot)}</dd>
      <dt>Calore dal termostato</dt><dd>${sgn(s.heatBath * KJ_PER_EV, 1)}</dd>
      <dt>Lavoro esterno (pareti, impulsi, pinzetta)</dt><dd>${sgn(s.work * KJ_PER_EV, 1)}</dd>
      <dt>Scambio per aggiunta/rimozione</dt><dd>${kj(s.matterExchange)}</dd>
      <dt>Energia della pinzetta</dt><dd>${kj(s.Egrab)}</dd>
      <dt>Deriva Δ(U − Q − W − E materia)</dt><dd>${sgn((s.diagnostics?.drift ?? 0) * KJ_PER_EV, 4)} kJ/mol</dd>
      <dt>Somma cariche parziali</dt><dd>${sgn(s.diagnostics?.totalCharge ?? 0, 6)} e</dd>
    </dl>
    <p class="hint">La deriva misura l’errore numerico dall’azzeramento. Confronta Δt e Δt/2 in NVE. 1 eV per scatola = 96,485 kJ per mole di copie della scatola.</p></div>
    ${physHtml(s)}`;
}

/** Grandezze di chimica fisica misurate sulla traiettoria. */
function physHtml(s) {
  const ph = SB.phys;
  if (!ph) return '';
  const R = 8.314462618;
  const cv = ph.cv, nMol = ph.nMol || 1;
  const cvMol = cv ? cv.Cv / nMol * KJ_PER_EV * 1000 : null;       // J/(mol K) per mole di molecole
  const cvErr = cv ? cv.err / nMol * KJ_PER_EV * 1000 : null;
  const monatomic = ph.N === nMol;
  return `<div><h3>Chimica fisica</h3>
    <dl class="info-list">
      <dt>C<sub>V</sub> dalle fluttuazioni</dt><dd>${cv && !cv.stationary ? 'non definita: il sistema non è in equilibrio (reazioni o riscaldamento in corso)' : cv ? `${nf(cvMol, 1)} ± ${nf(cvErr, 1)} J/(mol K)` : (s.thermostat && !ph.barostat?.on ? 'in accumulo…' : 'serve il termostato a volume fisso')}</dd>
      ${cv?.stationary ? `<dt class="sub">in unità di R per molecola</dt><dd>${nf(cvMol / R, 2)} R${monatomic ? ' (gas monoatomico ideale: 1,50 R)' : ''}</dd><dt class="sub">campioni</dt><dd>${cv.samples}</dd>` : ''}
      <dt>Coefficiente di diffusione D</dt><dd>${ph.msd?.D !== null && ph.msd?.D !== undefined ? `${nf(ph.msd.D * 1e5, 3)}·10⁻⁵ cm²/s` : '—'}</dd>
      <dt class="sub">esponente α di MSD ∝ t^α</dt><dd>${ph.msd?.alpha !== null && ph.msd?.alpha !== undefined ? nf(ph.msd.alpha, 2) : '—'}</dd>
      ${ph.rdf ? `<dt>Primo picco di g(r) ${sym(ph.rdf.pair[0])}–${sym(ph.rdf.pair[1])}</dt><dd>${nf(Math.max(...ph.rdf.g), 2)} a ${nf(ph.rdf.r[ph.rdf.g.indexOf(Math.max(...ph.rdf.g))], 2)} Å</dd>` : ''}
    </dl>
    <p class="hint">C<sub>V</sub> = (⟨E²⟩ − ⟨E⟩²)/(k<sub>B</sub>T²) nell'insieme canonico (errore dalla media a blocchi). È classica: ogni vibrazione conta k<sub>B</sub> anche quando nella realtà è "congelata" dalla quantizzazione. D dalla relazione di Einstein, MSD = 6Dt.</p></div>`;
}

// ---------------------------------------------------------------------------
// Grafici
// ---------------------------------------------------------------------------

const SPECIES_COLORS = ['#f0b429', '#3fa7ff', '#ff6b57', '#62d48f', '#c77dff', '#4dd6c9', '#e879a6', '#a3a3a3'];

function drawCharts() {
  if (!active) return;
  const H = SB.history;
  const muted = cssVar('--muted');
  // 1. energia / temperatura / composizione nel tempo
  $('radial-title').textContent = { rdf: 'Struttura: funzione di distribuzione radiale', msd: 'Diffusione: spostamento quadratico medio' }[SB.chart1] ?? 'Andamento nel tempo';
  if (SB.chart1 === 'rdf') drawRDF();
  else if (SB.chart1 === 'msd') drawMSD();
  else if (H.length > 1) {
    const ts = H.map(h => h.t / 1000);
    const xmin = ts[0], xmax = Math.max(ts[ts.length - 1], xmin + 0.1);
    if (SB.chart1 === 'energy') {
      const E0 = H[0].E, C0 = H[0].Ec ?? H[0].E;
      const ek = H.map(h => (h.E - h.Ep) * KJ_PER_EV), ep = H.map(h => (h.Ep - H[0].Ep) * KJ_PER_EV), et = H.map(h => (h.E - E0) * KJ_PER_EV);
      const ec = H.map(h => ((h.Ec ?? h.E) - C0) * KJ_PER_EV);
      const all = [...ek, ...ep, ...et, ...ec];
      const lo = Math.min(...all), hi = Math.max(...all);
      drawXY($('chart-radial'), {
        series: [
          { xs: ts, ys: ek, color: cssVar('--phase-pos'), label: 'cinetica' },
          { xs: ts, ys: ep, color: cssVar('--phase-neg'), label: 'potenziale (Δ)' },
          { xs: ts, ys: et, color: cssVar('--accent'), label: 'totale U (Δ)', width: 2 },
          { xs: ts, ys: ec, color: cssVar('--text'), label: 'U − Q − W − E materia', dash: [5, 4], width: 1.6 },
        ],
        xmin, xmax, ymin: lo - 0.05 * (hi - lo + 1), ymax: hi + 0.05 * (hi - lo + 1), xlabel: 't (ps)', ylabel: 'energia (kJ/mol)',
      });
      $('radial-note').innerHTML = 'Quando si formano legami l\'energia potenziale scende e quella cinetica (la temperatura) sale. La linea tratteggiata sottrae calore, lavoro e scambi di materia: una deriva segnala errore numerico.';
    } else {
      const Ts = H.map(h => h.T);
      const hi = Math.max(...Ts, SB.frame?.stats.Ttarget ?? 0) * 1.1;
      drawXY($('chart-radial'), {
        series: [{ xs: ts, ys: Ts, color: cssVar('--phase-pos'), label: 'T istantanea', width: 2 }],
        hlines: SB.frame?.stats.thermostat ? [{ y: SB.frame.stats.Ttarget, label: 'bagno termico', color: muted }] : [],
        xmin, xmax, ymin: 0, ymax: hi || 1, xlabel: 't (ps)', ylabel: 'T (K)',
      });
      $('radial-note').innerHTML = 'T = 2E<sub>cin</sub>/(3Nk<sub>B</sub>). Le fluttuazioni si riducono come 1/√N: sono il segno che la temperatura è una grandezza statistica.';
    }
  } else clearCanvas($('chart-radial'));
  // 2. Maxwell–Boltzmann
  const mb = SB.frame?.mb;
  $('levels-title').textContent = mb ? `Velocità delle molecole di ${mb.species}` : 'Distribuzione delle velocità';
  if (mb && mb.samples > 0) {
    const nb = mb.bins.length, dv = mb.vmax / nb;
    const toMS = 1e5; // Å/fs → m/s
    const xs = [], ys = [];
    mb.bins.forEach((c, k) => {
      const y = c / (mb.samples * dv * toMS);
      xs.push(k * dv * toMS, (k + 1) * dv * toMS); ys.push(y, y);
    });
    const theo = { xs: [], ys: [] };
    const a = mb.M * MV2 / (2 * KB_EV * mb.T);
    for (let k = 0; k <= 120; k++) {
      const v = k / 120 * mb.vmax;
      const f = 4 * Math.PI * Math.pow(a / Math.PI, 1.5) * v * v * Math.exp(-a * v * v);
      theo.xs.push(v * toMS); theo.ys.push(f / toMS);
    }
    const ymax = Math.max(...ys, ...theo.ys) * 1.15;
    drawXY($('chart-levels'), {
      series: [
        { xs, ys, color: cssVar('--phase-neg'), fill: `color-mix(in srgb, ${cssVar('--phase-neg')} 25%, transparent)`, label: 'simulazione' },
        { xs: theo.xs, ys: theo.ys, color: cssVar('--accent'), width: 2.2, label: `Maxwell–Boltzmann ${nf(mb.T, 0)} K` },
      ],
      xmin: 0, xmax: mb.vmax * toMS, ymin: 0, ymax, xlabel: 'v (m/s)', ylabel: 'f(v) (s/m)',
    });
    const vmean = Math.sqrt(8 * KB_EV * mb.T / (Math.PI * mb.M * MV2)) * toMS;
    $('levels-note').innerHTML = `f(v) = 4π (m/2πk<sub>B</sub>T)<sup>3/2</sup> v² e<sup>−mv²/2k<sub>B</sub>T</sup>, con m = ${nf(mb.M, 2)} u; T del bagno in NVT, cinetica istantanea in NVE. Il confronto presuppone equilibrio termico. Velocità media teorica √(8k<sub>B</sub>T/πm) = ${nf(vmean, 0)} m/s. ${mb.count} molecole, media mobile su molti fotogrammi.`;
  } else { clearCanvas($('chart-levels')); $('levels-note').textContent = ''; }
  // 3. specie nel tempo
  $('slice-title').textContent = 'Composizione nel tempo';
  if (H.length > 1) {
    const last = H[H.length - 1].counts;
    const keys = new Set();
    H.forEach(h => Object.keys(h.counts).forEach(k => keys.add(k)));
    const ranked = [...keys].sort((a, b) => Math.max(...H.map(h => h.counts[b] ?? 0)) - Math.max(...H.map(h => h.counts[a] ?? 0))).slice(0, 7);
    const ts = H.map(h => h.t / 1000);
    const series = ranked.map((k, c) => ({ xs: ts, ys: H.map(h => h.counts[k] ?? 0), color: SPECIES_COLORS[c], label: `${k} (${last[k] ?? 0})`, width: 1.8 }));
    const ymax = Math.max(1, ...series.flatMap(s => s.ys)) * 1.1;
    drawXY($('chart-slice'), { series, xmin: ts[0], xmax: Math.max(ts[ts.length - 1], ts[0] + 0.1), ymin: 0, ymax, xlabel: 't (ps)', ylabel: 'numero di molecole' });
    $('slice-note').innerHTML = 'Frammenti per connettività: soglie 0,55 / 0,35. Formule uguali possono indicare isomeri diversi; le curve non sono costanti cinetiche sperimentali.';
  } else { clearCanvas($('chart-slice')); $('slice-note').textContent = ''; }
}

function drawRDF() {
  const r = SB.phys?.rdf;
  if (!r) { clearCanvas($('chart-radial')); $('radial-note').textContent = 'g(r) si accumula durante la simulazione.'; return; }
  const [za, zb] = r.pair;
  const gmax = Math.max(1.5, ...r.g) * 1.1;
  drawXY($('chart-radial'), {
    series: [{ xs: r.r, ys: r.g, color: cssVar('--accent'), width: 2.2, label: `g(r) ${sym(za)}–${sym(zb)}`, fill: `color-mix(in srgb, ${cssVar('--accent')} 15%, transparent)` }],
    hlines: [{ y: 1, label: 'gas ideale', color: cssVar('--muted') }],
    xmin: 0, xmax: r.r[r.r.length - 1] + 0.05, ymin: 0, ymax: gmax, xlabel: 'r (Å)', ylabel: 'g(r)',
  });
  const k = r.g.indexOf(Math.max(...r.g));
  $('radial-note').innerHTML = `Funzione di distribuzione radiale fra ${sym(za)} e ${sym(zb)} di molecole diverse, media su ${r.samples} configurazioni. È normalizzata sulla distribuzione delle distanze di punti uniformi nel cubo di lato ${nf(r.L, 1)} Å, quindi g = 1 significa "nessuna struttura". Primo picco: g = ${nf(r.g[k], 2)} a ${nf(r.r[k], 2)} Å (primo guscio di vicini). Nei liquidi compaiono altri gusci smorzati; nel gas g ≈ 1 oltre il diametro atomico.`;
}

function drawMSD() {
  const m = SB.phys?.msd;
  if (!m || m.points.length < 3) { clearCanvas($('chart-radial')); $('radial-note').textContent = 'Lo spostamento quadratico medio si accumula durante la simulazione.'; return; }
  const xs = m.points.map(p => p[0] / 1000), ys = m.points.map(p => p[1]);
  const series = [{ xs, ys, color: cssVar('--accent'), width: 2.2, label: 'MSD simulato' }];
  if (m.fit) series.push({ xs: [m.fit.t0 / 1000, m.fit.t1 / 1000], ys: [m.fit.a + m.fit.b * m.fit.t0, m.fit.a + m.fit.b * m.fit.t1], color: cssVar('--phase-neg'), width: 2, dash: [5, 4], label: 'retta di Einstein' });
  const ymax = Math.max(...ys, 1) * 1.1;
  drawXY($('chart-radial'), {
    series, hlines: m.saturation < ymax * 1.5 ? [{ y: m.saturation, label: 'pareti: L²/2', color: cssVar('--muted') }] : [],
    xmin: 0, xmax: Math.max(xs[xs.length - 1], 0.01), ymin: 0, ymax: Math.max(ymax, m.saturation < ymax * 1.5 ? m.saturation * 1.08 : 0), xlabel: 't (ps)', ylabel: 'MSD (Å²)',
  });
  const regime = m.alpha === null ? '' : m.alpha > 1.6 ? 'moto balistico (α ≈ 2: gli atomi volano liberi fra un urto e l\'altro, come in un gas rarefatto)' : m.alpha > 1.25 ? 'transizione fra moto balistico e diffusivo' : m.alpha > 0.75 ? 'regime diffusivo (α ≈ 1)' : 'moto confinato (α < 1: gabbia dei vicini, solido o pareti)';
  $('radial-note').innerHTML = `⟨|r(t) − r(0)|²⟩ degli atomi. Einstein: MSD = 6 D t nel regime diffusivo. Esponente α = d ln MSD/d ln t = ${m.alpha === null ? '—' : nf(m.alpha, 2)}: ${regime}. ${m.D !== null ? `<b>D = ${nf(m.D * 1e5, 3)}·10⁻⁵ cm²/s</b>.` : 'D non si riporta finché il moto non è diffusivo.'} Le pareti limitano MSD a L²/2.`;
}

function clearCanvas(c) {
  const g = c.getContext('2d');
  g.clearRect(0, 0, c.width, c.height);
}

// ---------------------------------------------------------------------------
// Analisi: composizione, reazioni, atomo selezionato
// ---------------------------------------------------------------------------

function renderAnalysis() {
  if (!active) return;
  $('analysis').innerHTML = `
    <h2>Che cosa succede nella scatola</h2>
    <div class="analysis-grid">
      <div id="sb-selected"></div>
      <div id="sb-species"></div>
      <div id="sb-events" class="wide"></div>
      <details class="wide"><summary>Modello, equazioni e limiti</summary>
        <p>Velocity Verlet: F = −∇U, x(t+Δt) = x + vΔt + ½aΔt². NVE: energia approssimativamente conservata con errore dipendente dal passo. NVT: termostato CSVR di Bussi.</p>
        <p>MINDO/3: nuclei classici, elettroni di valenza con SCF UHF semiempirica (Bingham, Dewar, Lo 1975) a ogni passo; fino a 90 atomi di H, B, C, N, O, F, P, S, Cl. Spin libero (livello di Fermi comune), parametri verificati contro PySCF. È un metodo semiempirico: errori tipici sui calori di formazione di circa 11 kcal/mol, legami a idrogeno sottostimati.</p>
        <p>Hartree–Fock: nuclei classici, elettroni UHF/STO-3G; fino a 8 atomi. SCF non convergente: arresto. Lo spin iniziale è scelto fra le due molteplicità più basse, non fra tutti gli stati possibili.</p>
        <p>Campo classico: potenziale empirico specifico del progetto, ispirato a forme pubblicate. Barriere e reazioni non validate in generale; non è un’implementazione parametrizzata di ReaxFF o REBO. Assenti solvente, fotofisica e cinetica elettronica.</p>
        <p>Le cariche parziali non sono numeri di ossidazione. Gli ordini frazionari dinamici non determinano univocamente Lewis, ibridazione o specie chimica.</p>
      </details>
    </div>`;
  renderLiveAnalysis();
}

function renderLiveAnalysis() {
  if (!active || !$('sb-species')) return;
  const c = SB.census;
  // composizione
  if (c) {
    const total = c.species.reduce((s, [, n]) => s + n, 0);
    $('sb-species').innerHTML = `<h3>Composizione (${total} molecole)</h3>
      <div class="table-scroll short"><table class="data-table"><thead><tr><th>Specie</th><th>Nome</th><th class="num">n</th><th class="num">frazione molare</th><th></th></tr></thead><tbody>
      ${c.species.map(([f, n]) => `<tr><td>${f}</td><td>${SPECIES_NAMES[f] ? `compatibile con ${SPECIES_NAMES[f]}` : 'formula, identità non assegnata'}</td><td class="num">${n}</td><td class="num">${nf(n / total, 3)}</td>
        <td><button type="button" class="linkish" data-mb="${f}" title="Mostra la distribuzione delle velocità di questa specie">v</button></td></tr>`).join('')}
      </tbody></table></div>`;
    $('sb-species').querySelectorAll('[data-mb]').forEach(b => b.addEventListener('click', () => post({ type: 'set', mbSpecies: b.dataset.mb })));
  }
  // reazioni
  // con pochi atomi (ab initio) anche lo scambio di atomi uguali è la reazione da osservare
  const showExchange = (SB.frame?.stats.forceField === 'hf') || (SB.frame?.N ?? 99) <= 8;
  const real = SB.events.filter(e => showExchange || !e.exchange);
  const agg = new Map();
  for (const e of real) {
    const k = `${e.reactants.join(' + ')} → ${e.products.join(' + ')}`;
    agg.set(k, (agg.get(k) ?? 0) + 1);
  }
  const top = [...agg].sort((a, b) => b[1] - a[1]).slice(0, 12);
  const recent = real.slice(-14).reverse();
  $('sb-events').innerHTML = `<h3>Cambi di connettività osservati</h3>
    ${real.length ? `<div class="analysis-grid">
      <div><p class="lbl">Ultimi eventi</p><ul class="sb-log">${recent.map(e => `<li><span class="t">${nf(e.t / 1000, 3)} ps</span> ${e.reactants.join(' + ')} → ${e.products.join(' + ')}${e.exchange ? ' <span class="desc-muted">(scambio di atomi)</span>' : ''}</li>`).join('')}</ul></div>
      <div><p class="lbl">Reazioni elementari più frequenti</p><ul class="sb-log">${top.map(([k, n]) => `<li><span class="t">${n}×</span> ${k}</li>`).join('')}</ul></div>
    </div>` : '<p class="desc-muted">Nessun cambio di connettività registrato. Avvia o avanza di un passo per osservare la traiettoria.</p>'}
    <p class="desc-muted">Una reazione è registrata quando cambia la connettività: un legame si forma (ordine > 0,55) o si rompe (< 0,35). Le reazioni che scambiano solo atomi uguali (H + H₂ → H₂ + H) non sono elencate.</p>`;
  renderSelected();
}

function renderSelected() {
  const box = $('sb-selected');
  if (!box) return;
  const f = SB.frame;
  const i = SB.selected;
  if (!f || i < 0 || i >= f.N) {
    box.innerHTML = '<h3>Atomo selezionato</h3><p class="desc-muted">Clicca un atomo nella scatola (strumento "Seleziona") per vederne distanze, ordini di legame e carica parziale.</p>';
    return;
  }
  const Z = f.Z[i];
  const partners = [];
  for (let k = 0; k < f.bonds.length / 4; k++) {
    const a = f.bonds[4 * k], b = f.bonds[4 * k + 1], n = f.bonds[4 * k + 2], s = f.bonds[4 * k + 3];
    if (s < 0.3) continue;
    if (a === i) partners.push({ j: b, n, s }); else if (b === i) partners.push({ j: a, n, s });
  }
  const sumN = partners.reduce((sum, p) => sum + p.n, 0);
  const frag = SB.census?.frags?.find(fr => fr.includes(i));
  const fragFormula = frag ? hillFormula(frag.map(k => ({ Z: f.Z[k] }))) : sym(Z);
  const orderName = (n) => (n > 2.5 ? 'triplo' : n > 1.75 ? 'doppio' : n > 1.3 ? 'frazionario' : n > 0.75 ? 'singolo' : 'parziale');
  box.innerHTML = `<h3>Atomo selezionato: ${sym(Z)}${i + 1} (${ELEMENTS[Z - 1].name})</h3>
    <dl class="info-list">
      <dt>Carica parziale</dt><dd>${sgn(f.q[i], 2)} e</dd>
      <dt>Elettronegatività di Pauling</dt><dd>${nf(PAULING[Z], 2)}</dd>
      <dt>Numero di ossidazione</dt><dd>Richiede struttura di Lewis assegnata</dd>
      <dt>Somma degli ordini di legame</dt><dd>${nf(sumN, 2)} (valenza ${ATOM_PARAMS[Z][0]})</dd>

      <dt>Energia cinetica</dt><dd>${nf(f.ke[i] * KJ_PER_EV, 1)} kJ/mol</dd>
      <dt>Molecola</dt><dd>${fragFormula}${SPECIES_NAMES[fragFormula] ? ` (${SPECIES_NAMES[fragFormula]})` : ''}</dd>
    </dl>
    ${partners.length ? `<table class="data-table" style="margin-top:8px"><thead><tr><th>Legato a</th><th class="num">ordine</th><th>tipo</th><th class="num">distanza</th><th class="num">q vicino</th></tr></thead><tbody>
      ${partners.map(p => {
        const d = Math.hypot(f.pos[3 * i] - f.pos[3 * p.j], f.pos[3 * i + 1] - f.pos[3 * p.j + 1], f.pos[3 * i + 2] - f.pos[3 * p.j + 2]);
        return `<tr><td>${sym(f.Z[p.j])}${p.j + 1}</td><td class="num">${nf(p.n, 2)}</td><td>${orderName(p.n)}</td><td class="num">${nf(d * 100, 0)} pm</td><td class="num">${sgn(f.q[p.j], 2)}</td></tr>`;
      }).join('')}</tbody></table>` : '<p class="desc-muted">Atomo libero (nessun legame).</p>'}
    ${frag && frag.length > 1 && frag.length <= 30 ? '<div class="btn-row"><button type="button" class="btn" id="sb-to-mol">Studia questa molecola con Hartree–Fock</button></div>' : ''}`;
  $('sb-to-mol')?.addEventListener('click', () => {
    const idx = new Map(frag.map((k, n) => [k, n]));
    const atoms = frag.map(k => ({ Z: f.Z[k], charge: 0, aromatic: false }));
    const bonds = [];
    for (let k = 0; k < f.bonds.length / 4; k++) {
      const a = f.bonds[4 * k], b = f.bonds[4 * k + 1], n = f.bonds[4 * k + 2], s = f.bonds[4 * k + 3];
      if (s < 0.4 || !idx.has(a) || !idx.has(b)) continue;
      bonds.push({ a: idx.get(a), b: idx.get(b), order: Math.max(1, Math.min(3, Math.round(n))), aromatic: false });
    }
    const smiles = writeSmiles({ atoms, bonds });
    if (smiles && openInMolecule) openInMolecule(smiles, `${fragFormula} dalla sandbox`);
  });
}
