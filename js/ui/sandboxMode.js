// Modalità "Sandbox": dinamica molecolare reattiva in tempo reale. Si aggiungono atomi e molecole in
// una scatola, si controllano temperatura, volume e luce, e si osservano legami che si formano e si
// rompono, specie chimiche, reazioni, pressione e distribuzione di Maxwell–Boltzmann.

import * as THREE from 'three';
import { PRESETS, SANDBOX_MOLECULES, SANDBOX_ELEMENTS, SPECIES_NAMES } from '../chem/sandboxData.js';
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
  preset: PRESETS[3],
  frame: null,
  fresh: false,
  census: null,
  history: [],
  events: [],
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
      SB.frame = m;
      SB.fresh = true;
      if (m.census) {
        SB.census = m.census;
        SB.history = m.census.history;
        if (m.census.events.length) {
          SB.events.push(...m.census.events);
          if (SB.events.length > 500) SB.events.splice(0, SB.events.length - 500);
        }
      }
    } else if (m.type === 'info') {
      SB.info = m.text;
      renderSide();
    } else if (m.type === 'forcefield') {
      SB.forceField = m.kind;
      if (active) renderControls();
    } else if (m.type === 'photon') {
      SB.photonFlash = { i: m.i, j: m.j, t: performance.now(), lambda: m.lambda };
    } else if (m.type === 'error') {
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

function loadPreset(p) {
  SB.preset = p;
  SB.events = [];
  SB.history = [];
  SB.census = null;
  SB.selected = -1;
  SB.info = '';
  if (p.color) SB.color = p.color;
  if (p.light) SB.lambda = p.light.lambda;
  post({ type: 'preset', preset: p });
  post({ type: 'set', paused: false, T: p.T });
  SB.userPaused = false;
  scene.framed = false;
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
  scene.cap = 0; scene.bcap = 0;
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
    scene.bonds = new THREE.InstancedMesh(new THREE.CylinderGeometry(1, 1, 1, 10, 1, false), new THREE.MeshStandardMaterial({ color: light ? '#8a93a3' : '#b7bfcc', roughness: 0.5 }), scene.bcap);
    scene.bonds.frustumCulled = false;
    scene.group.add(scene.bonds);
  }
}

function atomRadius(Z) {
  if (SB.style === 'vdw') return 0.5 * ATOM_PARAMS[Z][5] * 0.82;
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
  if (!scene.framed && Math.abs(L - SB.preset.box) < 1e-6) {
    viewer.frame(L * 0.62);
    viewer.camera.position.set(0.55, -0.8, 0.45).normalize().multiplyScalar(L * 0.62 * 3.1);
    viewer.controls.target.set(0, 0, 0);
    scene.framed = true;
    renderControls();
  }
  flashUpdate();
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
  const showBonds = SB.style !== 'vdw';
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
      scene.bonds.setMatrixAt(nb++, tmpM);
    }
  }
  scene.bonds.count = nb;
  scene.bonds.instanceMatrix.needsUpdate = true;
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
  scene.hb.visible = SB.style !== 'vdw' && nh > 0;
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
  if (now - SB.lastPanel > 200) { SB.lastPanel = now; renderPanelBody(); }
  if (now - SB.lastCharts > 500) { SB.lastCharts = now; drawCharts(); renderLiveAnalysis(); }
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
  const m = SANDBOX_MOLECULES.find(x => x.id === SB.add.id);
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
  const elBtns = SANDBOX_ELEMENTS.map(Z => `<button type="button" class="mol-chip ${SB.add.kind === 'atom' && SB.add.symbol === sym(Z) ? 'active' : ''}" data-el="${sym(Z)}" title="${ELEMENTS[Z - 1].name}">${sym(Z)}</button>`).join('');
  $('element-card').innerHTML = `
    <h3 class="side-h">Sandbox chimica</h3>
    <label class="lbl" for="sb-preset">Esperimento</label>
    <select id="sb-preset">
      <optgroup label="Campo di forze reattivo">${PRESETS.filter(x => !x.atoms).map(x => `<option value="${x.id}" ${x.id === p.id ? 'selected' : ''}>${x.name}</option>`).join('')}</optgroup>
      <optgroup label="Ab initio (Hartree–Fock)">${PRESETS.filter(x => x.atoms).map(x => `<option value="${x.id}" ${x.id === p.id ? 'selected' : ''}>${x.name}</option>`).join('')}</optgroup>
    </select>
    <p class="mol-note">${p.text}</p>
    ${p.tips?.length ? `<ul class="sb-tips">${p.tips.map(t => `<li>${t}</li>`).join('')}</ul>` : ''}
    <div class="btn-row"><button type="button" class="btn" id="sb-restart">Ricomincia</button><button type="button" class="btn" id="sb-clear">Svuota la scatola</button></div>
    <label class="lbl">Aggiungi molecole</label>
    <div class="mol-list">${molBtns}</div>
    <label class="lbl">Aggiungi atomi (radicali, gas nobili)</label>
    <div class="mol-list">${elBtns}</div>
    <label class="lbl" for="sb-smiles">Oppure una molecola da SMILES</label>
    <div class="smiles-row"><input id="sb-smiles" placeholder="es. CC(=O)O" value="${SB.add.kind === 'smiles' ? SB.add.smiles : ''}" spellcheck="false"><button type="button" class="btn" id="sb-smiles-ok">Usa</button></div>
    <div class="ctl" style="margin-top:8px"><label class="lbl" for="sb-count">Quantità: <span id="sb-count-out">${SB.count}</span></label>
      <div class="range-row"><input type="range" id="sb-count" min="1" max="40" step="1" value="${SB.count}"><button type="button" class="btn" id="sb-add">Aggiungi ${escapeHtml(addLabel())}</button></div></div>
    <p class="hint">Con lo strumento <b>Aggiungi</b> (a destra) puoi anche cliccare nella scatola per mettere una molecola dove vuoi.</p>
    ${SB.info ? `<p class="hint warn">${SB.info}</p>` : ''}`;
  $('sb-preset').addEventListener('change', (e) => loadPreset(PRESETS.find(x => x.id === e.target.value)));
  $('sb-restart').addEventListener('click', () => loadPreset(SB.preset));
  $('sb-clear').addEventListener('click', () => { post({ type: 'clear' }); SB.events = []; SB.history = []; SB.selected = -1; renderAnalysis(); });
  $('element-card').querySelectorAll('[data-mol]').forEach(b => b.addEventListener('click', () => { SB.add = { kind: 'mol', id: b.dataset.mol }; renderSide(); }));
  $('element-card').querySelectorAll('[data-el]').forEach(b => b.addEventListener('click', () => { SB.add = { kind: 'atom', symbol: b.dataset.el }; renderSide(); }));
  $('sb-smiles-ok').addEventListener('click', () => {
    const s = $('sb-smiles').value.trim();
    try {
      const g = parseSmiles(s);
      if (g.atoms.some(a => !ATOM_PARAMS[a.Z])) throw new Error('elemento non disponibile nella sandbox');
      if (g.atoms.some(a => a.charge)) throw new Error('la sandbox contiene solo specie neutre');
      SB.add = { kind: 'smiles', smiles: s };
      SB.info = '';
    } catch (e) { SB.info = `SMILES non valido: ${e.message}`; }
    renderSide();
  });
  $('sb-count').addEventListener('input', (e) => { SB.count = +e.target.value; $('sb-count-out').textContent = SB.count; });
  $('sb-add').addEventListener('click', () => post({ type: 'add', ...addPayload(SB.count) }));
}

function addLabel() {
  if (SB.add.kind === 'atom') return SB.add.symbol;
  if (SB.add.kind === 'smiles') return SB.add.smiles;
  const m = SANDBOX_MOLECULES.find(x => x.id === SB.add.id);
  return hillFormula(parseSmiles(m.smiles).atoms);
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

const TOOLS = [
  { id: 'select', label: 'Seleziona', title: 'Clic su un atomo: legami, carica, ossidazione. Trascina per ruotare la vista.' },
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
  $('controls').innerHTML = `
    <div class="ctl"><span class="lbl">Simulazione</span>
      <div class="btn-row" style="margin-top:0">
        <button type="button" class="btn" id="sb-play">${SB.userPaused ? '▶ Avvia' : '❚❚ Pausa'}</button>
        <button type="button" class="btn" id="sb-step" ${SB.userPaused ? '' : 'disabled'}>+10 fs</button>
      </div></div>
    <div class="ctl"><span class="lbl">Forze sugli atomi</span></div>
    <div class="seg" id="sb-ff">
      <button type="button" data-v="reactive" aria-pressed="${(st?.forceField ?? SB.forceField ?? 'reactive') === 'reactive'}" title="Campo di forze reattivo: centinaia di atomi in tempo reale">campo reattivo</button>
      <button type="button" data-v="hf" aria-pressed="${(st?.forceField ?? SB.forceField) === 'hf'}" title="Hartree–Fock/STO-3G a ogni passo: forze dalla meccanica quantistica, fino a 8 atomi">Hartree–Fock ab initio</button>
    </div>
    ${logSlider('sb-T', `Temperatura del termostato`, 10, 8000, T, v => `${nf(v, 0)} K`)}
    <div class="seg" id="sb-thermo">
      <button type="button" data-v="1" aria-pressed="${thermo}" title="Termostato di Bussi: scambia calore con un bagno a temperatura costante (insieme canonico NVT)">Termostato</button>
      <button type="button" data-v="0" aria-pressed="${!thermo}" title="Nessuno scambio di calore: l'energia totale si conserva (NVE). Le reazioni esotermiche scaldano il sistema.">Isolato</button>
    </div>
    <div class="ctl"><label class="lbl" for="sb-box">Lato della scatola (volume)</label>
      <div class="range-row"><input type="range" id="sb-box" min="10" max="60" step="0.5" value="${box}"><output id="sb-box-out">${nf(box, 1)} Å</output></div></div>
    <div class="ctl"><label class="lbl" for="sb-speed">Velocità (passi per fotogramma)</label>
      <div class="range-row"><input type="range" id="sb-speed" min="1" max="200" step="1" value="${spf}"><output id="sb-speed-out">${spf}</output></div></div>
    <div class="ctl"><span class="lbl">Strumento</span></div>
    <div class="seg" id="sb-tool">${TOOLS.map(t => `<button type="button" data-v="${t.id}" aria-pressed="${SB.tool === t.id}" title="${t.title}">${t.label}</button>`).join('')}</div>
    <p class="hint" id="sb-tool-hint">${TOOLS.find(t => t.id === SB.tool).title}</p>
    <div class="ctl"><label class="lbl" for="sb-lambda">Luce: λ = <span id="sb-lambda-out">${SB.lambda} nm</span> <i class="swatch" id="sb-lambda-sw" style="background: ${lambdaCss(SB.lambda)}"></i></label>
      <div class="range-row"><input type="range" id="sb-lambda" min="100" max="800" step="5" value="${SB.lambda}"><output id="sb-ephoton">${nf(photonKJ, 0)} kJ/mol</output></div>
      <div class="btn-row" style="margin-top:4px"><button type="button" class="btn" id="sb-light" aria-pressed="${lightOn}">${lightOn ? 'Spegni la luce' : 'Accendi la luce'}</button><button type="button" class="btn" id="sb-flash">Un fotone</button></div>
      <p class="hint">E = hc/λ: ogni fotone assorbito deposita la sua energia in un legame.</p></div>
    <div class="ctl"><span class="lbl">Colore degli atomi</span></div>
    <div class="seg" id="sb-color">
      <button type="button" data-v="element" aria-pressed="${SB.color === 'element'}">elemento</button>
      <button type="button" data-v="charge" aria-pressed="${SB.color === 'charge'}">carica parziale</button>
      <button type="button" data-v="ke" aria-pressed="${SB.color === 'ke'}">energia cinetica</button>
    </div>
    <div class="seg" id="sb-style" style="margin-top:4px">
      <button type="button" data-v="ball" aria-pressed="${SB.style === 'ball'}">sfere e bastoncini</button>
      <button type="button" data-v="vdw" aria-pressed="${SB.style === 'vdw'}">van der Waals</button>
    </div>`;
  const seg = (id, fn) => $(id).querySelectorAll('button').forEach(b => b.addEventListener('click', () => {
    fn(b.dataset.v);
    $(id).querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', String(x === b)));
  }));
  $('sb-play').addEventListener('click', () => { SB.userPaused = !SB.userPaused; post({ type: 'set', paused: SB.userPaused }); renderControls(); });
  $('sb-step').addEventListener('click', () => post({ type: 'step', n: Math.round(10 / (st?.dt ?? 0.4)) }));
  $('sb-T').addEventListener('input', (e) => { const v = Math.pow(10, +e.target.value); $('sb-T-out').textContent = `${nf(v, 0)} K`; post({ type: 'set', T: v }); });
  seg('sb-thermo', (v) => post({ type: 'set', thermostat: v === '1' }));
  seg('sb-ff', (v) => post({ type: 'set', forceField: v }));
  $('sb-box').addEventListener('input', (e) => { const v = +e.target.value; $('sb-box-out').textContent = `${nf(v, 1)} Å`; post({ type: 'set', box: v }); });
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
  seg('sb-style', (v) => { SB.style = v; SB.fresh = true; });
  renderLegend();
}

/** Colore di una lunghezza d'onda; ultravioletto e infrarosso con colori convenzionali. */
function lambdaCss(nm) {
  if (nm < 380) return '#a78bfa';
  if (nm > 750) return '#8b1e1e';
  return wavelengthColor(nm);
}

function renderLegend() {
  const leg = $('viewport-legend');
  if (SB.color === 'charge') leg.innerHTML = '<span><i class="swatch" style="background:#e0402a"></i>δ− (negativa)</span> <span><i class="swatch" style="background:#dcdce0"></i>neutra</span> <span><i class="swatch" style="background:#3f7fe8"></i>δ+ (positiva)</span>';
  else if (SB.color === 'ke') leg.innerHTML = '<span><i class="swatch" style="background:#4070ff"></i>lento</span> <span><i class="swatch" style="background:#ffffff"></i>≈ 3/2 kT</span> <span><i class="swatch" style="background:#ff3020"></i>veloce</span>';
  else leg.innerHTML = 'Colori CPK · legami: 1, 2 o 3 cilindri secondo l\'ordine di legame; sottili = legami che si formano o si rompono; tratteggio = legame a idrogeno';
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
  $('viewport-title').innerHTML = `${SB.preset.name}<small>t = ${nf(s.t / 1000, 2)} ps · ${N} atomi · ${nMol} molecole</small>`;
  $('viewport-note').innerHTML = `${s.forceField === 'hf' ? 'Dinamica ab initio (Hartree–Fock a ogni passo)' : 'Dinamica molecolare reattiva'}: passo Δt = ${nf(s.dt, 2)} fs, ${s.paused ? '<b>in pausa</b>' : `${s.stepsPerFrame} passi per fotogramma`}. Trascina per ruotare, rotellina per ingrandire.`;
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
      <dt>Lavoro esterno (scintille, luce)</dt><dd>${sgn(s.work * KJ_PER_EV, 1)}</dd>
      <dt>Totale − scambi</dt><dd>${kj(s.Etot - s.heatBath - s.work)}</dd>
    </dl>
    <p class="hint">"Totale − scambi" è costante: è il primo principio della termodinamica (ΔU = Q + W) verificato passo per passo.</p></div>`;
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
  $('radial-title').textContent = 'Andamento nel tempo';
  if (H.length > 1) {
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
          { xs: ts, ys: ec, color: cssVar('--text'), label: 'U − Q − W', dash: [5, 4], width: 1.6 },
        ],
        xmin, xmax, ymin: lo - 0.05 * (hi - lo + 1), ymax: hi + 0.05 * (hi - lo + 1), xlabel: 't (ps)', ylabel: 'energia (kJ/mol)',
      });
      $('radial-note').innerHTML = 'Quando si formano legami l\'energia potenziale scende e quella cinetica (la temperatura) sale. La linea tratteggiata U − Q − W resta piatta: è il primo principio, ΔU = Q + W.';
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
    $('levels-note').innerHTML = `f(v) = 4π (m/2πk<sub>B</sub>T)<sup>3/2</sup> v² e<sup>−mv²/2k<sub>B</sub>T</sup>, con m = ${nf(mb.M, 2)} u e T del termostato. Velocità media teorica √(8k<sub>B</sub>T/πm) = ${nf(vmean, 0)} m/s. ${mb.count} molecole, media mobile su molti fotogrammi.`;
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
    $('slice-note').innerHTML = 'Le specie sono i gruppi di atomi uniti da legami (ordine di legame > 0,5). Le curve sono la cinetica chimica che emerge dagli urti.';
  } else { clearCanvas($('chart-slice')); $('slice-note').textContent = ''; }
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
      <article class="wide"><h3>Il modello fisico</h3>
        <div class="analysis-grid">
          <div>
            <p class="desc"><b>Equazioni del moto.</b> Ogni atomo segue la legge di Newton F = ma. Le equazioni sono integrate con l'algoritmo velocity Verlet (passo 0,4 fs), che conserva l'energia: in modalità "isolato" l'energia totale resta costante.</p>
            <p class="eq">x(t+Δt) = x + vΔt + ½aΔt²,&nbsp; v(t+Δt) = v + ½[a(t) + a(t+Δt)]Δt</p>
            <p class="desc"><b>Legami che si formano e si rompono.</b> L'energia di legame segue il potenziale di Abell–Tersoff–Brenner. L'ordine di legame n nasce dalla valenza libera di ciascun atomo (ottetto): C ne ha 4, N 3, O 2, H e gli alogeni 1. Lunghezza ed energia seguono la relazione di Pauling:</p>
            <p class="eq">r(n) = r₁ − c ln n,&nbsp; D(n) = D₁ nᵖ</p>
            <p class="desc">D₁, p e c vengono dalle energie medie di legame e dalle lunghezze tabulate (C–C, C=C, C≡C…). Così 2 H₂ + O₂ → 2 H₂O libera circa 500 kJ/mol, come la legge di Hess applicata alle energie di legame.</p>
          </div>
          <div>
            <p class="desc"><b>Elettronegatività e cariche.</b> Le cariche parziali si ricalcolano a ogni passo con l'equalizzazione dell'elettronegatività (principio di Sanderson): gli elettroni fluiscono lungo i legami verso l'atomo più elettronegativo, finché il potenziale chimico è uguale ovunque.</p>
            <p class="eq">E<sub>Q</sub> = Σ (χ<sub>i</sub>q<sub>i</sub> + ½J<sub>i</sub>q<sub>i</sub>²) + Σ q<sub>i</sub>q<sub>j</sub>γ<sub>ij</sub>(r)</p>
            <p class="desc"><b>Angoli.</b> Dal teorema di Coulson, cos θ = −1/λ per gli ibridi spᵏ (180°, 120°, 109,5°). Ogni coppia solitaria chiude l'angolo di circa 2,5° (VSEPR: NH₃ 107°, H₂O 104,5°). <b>Forze di van der Waals</b> dal campo di forze UFF. <b>Legami a idrogeno</b> D–H···A (D, A = N, O, F) con il termine del campo DREIDING, E = D<sub>hb</sub>[5(R₀/R)¹² − 6(R₀/R)¹⁰] cos⁴θ: il dimero d'acqua risulta legato di circa 15 kJ/mol a O···O = 2,9 Å.</p>
            <p class="desc"><b>Temperatura e pressione.</b> Il termostato di Bussi–Donadio–Parrinello riscala le velocità in modo stocastico e riproduce l'insieme canonico. La pressione è la forza media degli urti sulle pareti divisa per l'area: la teoria cinetica dei gas, misurata.</p>
            <p class="desc-muted">Limiti: è un modello classico. Gli elettroni non sono trattati esplicitamente: non ci sono stati di spin, né ioni in soluzione, né ipervalenza (SF₆, SO₃). Le barriere di reazione sono qualitative (H + H₂: circa 25 kJ/mol contro 40 misurati). Per la chimica quantistica esatta di una singola molecola usa la modalità Molecole.</p>
          </div>
        </div>
      </article>
    </div>`;
  renderLiveAnalysis();
}

function renderLiveAnalysis() {
  if (!active || !$('sb-species')) return;
  const c = SB.census;
  // composizione
  if (c) {
    const total = c.species.reduce((s, [, n]) => s + n, 0) || 1;
    $('sb-species').innerHTML = `<h3>Composizione (${total} molecole)</h3>
      <div class="table-scroll short"><table class="data-table"><thead><tr><th>Specie</th><th>Nome</th><th class="num">n</th><th class="num">frazione molare</th><th></th></tr></thead><tbody>
      ${c.species.map(([f, n]) => `<tr><td>${f}</td><td>${SPECIES_NAMES[f] ?? ''}</td><td class="num">${n}</td><td class="num">${nf(n / total, 3)}</td>
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
  $('sb-events').innerHTML = `<h3>Registro delle reazioni</h3>
    ${real.length ? `<div class="analysis-grid">
      <div><p class="lbl">Ultimi eventi</p><ul class="sb-log">${recent.map(e => `<li><span class="t">${nf(e.t / 1000, 3)} ps</span> ${e.reactants.join(' + ')} → ${e.products.join(' + ')}${e.exchange ? ' <span class="desc-muted">(scambio di atomi)</span>' : ''}</li>`).join('')}</ul></div>
      <div><p class="lbl">Reazioni elementari più frequenti</p><ul class="sb-log">${top.map(([k, n]) => `<li><span class="t">${n}×</span> ${k}</li>`).join('')}</ul></div>
    </div>` : '<p class="desc-muted">Ancora nessuna reazione. Alza la temperatura, usa una scintilla o accendi la luce.</p>'}
    <p class="desc-muted">Una reazione è registrata quando cambia la connettività: un legame si forma (ordine > 0,55) o si rompe (< 0,35). Le reazioni che scambiano solo atomi uguali (H + H₂ → H₂ + H) non sono elencate.</p>`;
  renderSelected();
}

function renderSelected() {
  const box = $('sb-selected');
  if (!box) return;
  const f = SB.frame;
  const i = SB.selected;
  if (!f || i < 0 || i >= f.N) {
    box.innerHTML = '<h3>Atomo selezionato</h3><p class="desc-muted">Clicca un atomo nella scatola (strumento "Seleziona") per vederne legami, carica parziale, ibridazione e numero di ossidazione.</p>';
    return;
  }
  const Z = f.Z[i];
  const partners = [];
  for (let k = 0; k < f.bonds.length / 4; k++) {
    const a = f.bonds[4 * k], b = f.bonds[4 * k + 1], n = f.bonds[4 * k + 2], s = f.bonds[4 * k + 3];
    if (s < 0.3) continue;
    if (a === i) partners.push({ j: b, n, s }); else if (b === i) partners.push({ j: a, n, s });
  }
  const chi = (z) => PAULING[z] ?? 2;
  let ox = 0;
  for (const p of partners) {
    const o = Math.round(p.n);
    if (f.Z[p.j] === Z || Math.abs(chi(f.Z[p.j]) - chi(Z)) < 1e-9) continue;
    ox += chi(f.Z[p.j]) > chi(Z) ? o : -o;
  }
  const sumN = partners.reduce((s, p) => s + p.n, 0);
  const lp = Math.max(0, (valenceElectrons(Z) - sumN) / 2);
  const steric = partners.length + Math.round(lp);
  const hyb = { 2: 'sp', 3: 'sp²', 4: 'sp³' }[steric] ?? '—';
  const frag = SB.census?.frags?.find(fr => fr.includes(i));
  const fragFormula = frag ? hillFormula(frag.map(k => ({ Z: f.Z[k] }))) : sym(Z);
  const orderName = (n) => (n > 2.5 ? 'triplo' : n > 1.75 ? 'doppio' : n > 1.3 ? 'aromatico / delocalizzato' : n > 0.75 ? 'singolo' : 'parziale');
  box.innerHTML = `<h3>Atomo selezionato: ${sym(Z)}${i + 1} (${ELEMENTS[Z - 1].name})</h3>
    <dl class="info-list">
      <dt>Carica parziale</dt><dd>${sgn(f.q[i], 2)} e</dd>
      <dt>Elettronegatività di Pauling</dt><dd>${nf(PAULING[Z], 2)}</dd>
      <dt>Numero di ossidazione</dt><dd>${ox > 0 ? '+' : ox < 0 ? '−' : ''}${Math.abs(ox)}</dd>
      <dt>Somma degli ordini di legame</dt><dd>${nf(sumN, 2)} (valenza ${ATOM_PARAMS[Z][0]})</dd>
      <dt>Coppie solitarie (stima)</dt><dd>${nf(lp, 1)}</dd>
      <dt>Ibridazione (VSEPR)</dt><dd>${partners.length >= 2 ? hyb : '—'}</dd>
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
