// Worker della sandbox: fa girare la dinamica molecolare e invia al disegno 3D un fotogramma
// alla volta (posizioni, legami, cariche, grandezze termodinamiche, specie e reazioni).

import { Simulation, MV2 } from './md.js';
import { ReactiveFF } from './reactive.js';
import { parseSmiles } from './smiles.js';
import { embedMolecule } from './embed.js';
import { analyzeStructure } from './structure.js';
import { KB_EV, ATOM_PARAMS } from './reactiveData.js';
import { makeHFProvider, aimdFeasible } from './aimd.js';

const BOHR_ANG = 0.52917721090;
let sim = new Simulation({ box: 20, T: 300 });
let paused = true;
let stepsPerFrame = 40;
let light = { on: false, lambda: 400, rate: 2 };   // fotoni per picosecondo
let lastLightTime = 0;
let msTarget = 14;
let mbSpecies = null;
let mbHist = null;
let lastEventSent = 0;
let censusSent = -1;
const templates = new Map();

/** Geometria di una molecola dal suo SMILES: modello VSEPR, poi minimizzazione con il campo reattivo. */
function template(smiles) {
  if (templates.has(smiles)) return templates.get(smiles);
  const g = parseSmiles(smiles);
  if (!g.atoms.length || g.atoms.some(a => !ATOM_PARAMS[a.Z] || a.charge))
    throw new Error('Sandbox: solo specie neutre con elementi supportati.');
  const an = analyzeStructure({ atoms: g.atoms, bonds: g.bonds, charge: 0 });
  const at = embedMolecule(g, an.atoms.map(a => a.lonePairs + a.radical), { seeds: 2, iterations: 1500 });
  const Z = at.map(a => a.Z);
  const pos = Float64Array.from(at.flatMap(a => a.xyz.map(v => v * BOHR_ANG)));
  if (Z.length > 1) {
    const ff = new ReactiveFF();
    const F = new Float64Array(pos.length), v = new Float64Array(pos.length);
    let dt = 0.05, alpha = 0.1, npos = 0;
    for (let s = 0; s < 2500; s++) {
      ff.compute(Z, pos, F);
      let P = 0, vn = 0, fn = 0;
      for (let k = 0; k < pos.length; k++) { P += F[k] * v[k]; vn += v[k] * v[k]; fn += F[k] * F[k]; }
      if (Math.sqrt(fn) < 1e-4) break;
      vn = Math.sqrt(vn); fn = Math.sqrt(fn);
      for (let k = 0; k < pos.length; k++) v[k] = (1 - alpha) * v[k] + alpha * F[k] / fn * vn;
      if (P > 0) { if (++npos > 5) { dt = Math.min(dt * 1.1, 0.3); alpha *= 0.99; } } else { npos = 0; dt *= 0.5; v.fill(0); alpha = 0.1; }
      for (let k = 0; k < pos.length; k++) { v[k] += dt * F[k]; pos[k] += dt * v[k]; }
    }
  }
  const t = { Z, pos };
  templates.set(smiles, t);
  return t;
}

function atomSmiles(sym) { return `[${sym}]`; }

function add({ smiles, symbol, count = 1, at = null, T }) {
  const s = smiles ?? atomSmiles(symbol);
  const t = template(s);
  // con atomi che formano legami serve un passo corto (vibrazioni di 10 fs); i gas nobili tollerano 2 fs
  if (t.Z.some(z => ATOM_PARAMS[z][0] > 0)) sim.dt = Math.min(sim.dt, 0.4);
  if (forceField === 'hf') {
    const f = aimdFeasible([...sim.Z, ...Array.from({ length: count }, () => t.Z).flat()], hfBasis);
    if (!f.ok) throw new Error(f.reason);
  }
  mbHist = null;
  const placed = sim.editInventory(() => sim.addMolecule(t, count, T ?? sim.T, at));
  return placed;
}

let forceField = 'reactive';
let hfBasis = 'STO-3G';

/** Sceglie il modello delle forze; con Hartree–Fock verifica che il sistema sia abbastanza piccolo. */
function setForceField(kind, basis = hfBasis) {
  hfBasis = basis;
  if (kind === 'hf') {
    const f = aimdFeasible(sim.Z, basis);
    if (!f.ok) throw new Error(f.reason);
  }
  if (!['hf', 'reactive'].includes(kind)) throw new Error('Modello non valido.');
  forceField = kind;
  sim.provider = kind === 'hf' ? makeHFProvider({ basis }) : null;
  sim.res = null;
  if (kind === 'hf') sim.dt = Math.min(sim.dt, 0.25);
  sim.resetMeasurements(); censusSent = -1; lastEventSent = 0; mbHist = null;
  postMessage({ type: 'forcefield', kind });
}

function loadPreset(p) {
  sim = new Simulation({ box: p.box, T: p.T, dt: p.dt ?? 0.2, seed: p.seed ?? 2024 });
  censusSent = -1; lastEventSent = 0; mbHist = null;
  stepsPerFrame = p.speed ?? 40;
  if (p.atoms) {
    sim.addAtoms(p.atoms);
    if (p.thermalize) for (let i = 0; i < sim.N; i++) sim.thermalize(i, p.thermalize);
    sim.thermostat = p.thermostat ?? false;
    setForceField(p.forceField ?? 'reactive');
    sim.forces();
    sim.resetMeasurements();
    light = { ...light, on: false };
    lastEventSent = 0;
    mbHist = null;
    return { placed: p.atoms.length, wanted: p.atoms.length };
  }
  setForceField('reactive');
  let placed = 0, wanted = 0;
  for (const [s, n] of p.add) {
    wanted += n;
    placed += add({ smiles: s, count: n });
  }
  if (p.photons) {
    sim.forces();
    const [za, zb] = p.photons.pair ?? [];
    for (let k = 0; k < p.photons.count; k++) sim.photon(p.photons.lambda, za ? (a, b) => (a === za && b === zb) || (a === zb && b === za) : null);
  }
  if (p.spark) sim.spark([0, 0, 0], p.spark.radius, p.spark.T);
  light = p.light ? { ...light, ...p.light } : { ...light, on: false };
  lastLightTime = 0;
  lastEventSent = 0;
  mbHist = null;
  sim.resetMeasurements();
  return { placed, wanted };
}

/** Istogramma delle velocità del centro di massa delle molecole di una specie. */
function speedHistogram() {
  if (!sim.frags?.length) return null;
  const counts = new Map();
  for (const f of sim.frags) counts.set(f.formula, (counts.get(f.formula) ?? 0) + 1);
  let species = mbSpecies && counts.has(mbSpecies) ? mbSpecies : null;
  if (!species) species = [...counts.entries()].sort((a, b) => b[1] - a[1])[0][0];
  const frags = sim.frags.filter(f => f.formula === species);
  const M = frags[0].atoms.reduce((s, i) => s + sim.mass[i], 0);
  const referenceT = Math.max(1e-6, sim.thermostat ? sim.T : sim.temperature());
  const vmp = Math.sqrt(2 * KB_EV * referenceT / (M * MV2)); // velocità più probabile (Å/fs)
  const vmax = 3.2 * vmp;
  const nb = 30;
  if (!mbHist || mbHist.species !== species || Math.abs(mbHist.vmax - vmax) / vmax > 0.02) {
    mbHist = { species, vmax, M, bins: new Float64Array(nb), samples: 0 };
  }
  if (mbHist.step === sim.stepCount) return { species, M, vmax, bins: Array.from(mbHist.bins), samples: mbHist.samples, count: frags.length, T: referenceT };
  mbHist.step = sim.stepCount;
  // media mobile esponenziale: il grafico segue i cambi di temperatura
  const decay = 0.97;
  for (let k = 0; k < nb; k++) mbHist.bins[k] *= decay;
  mbHist.samples = mbHist.samples * decay;
  for (const f of frags) {
    let vx = 0, vy = 0, vz = 0;
    for (const i of f.atoms) { vx += sim.mass[i] * sim.vel[3 * i]; vy += sim.mass[i] * sim.vel[3 * i + 1]; vz += sim.mass[i] * sim.vel[3 * i + 2]; }
    const v = Math.hypot(vx, vy, vz) / M;
    const k = Math.floor(v / vmax * nb);
    if (k >= 0 && k < nb) mbHist.bins[k] += 1;
    mbHist.samples += 1;
  }
  return { species, M, vmax, bins: Array.from(mbHist.bins), samples: mbHist.samples, count: frags.length, T: referenceT };
}

function frame() {
  if (!sim.res) sim.forces();
  const N = sim.N;
  const res = sim.res;
  const pos = Float32Array.from(sim.pos);
  const q = res ? Float32Array.from(res.q) : new Float32Array(N);
  const bonds = [];
  for (const b of res?.bonds ?? []) { const s = b.n * b.w; if (s > 0.12) bonds.push(b.i, b.j, b.n, s); }
  const ke = new Float32Array(N);
  for (let i = 0; i < N; i++) {
    const vx = sim.vel[3 * i], vy = sim.vel[3 * i + 1], vz = sim.vel[3 * i + 2];
    ke[i] = 0.5 * sim.mass[i] * (vx * vx + vy * vy + vz * vz) * MV2;
  }
  const Ekin = sim.kinetic();
  const msg = {
    type: 'frame',
    N, Z: Int8Array.from(sim.Z), pos, q, ke, bonds: Float32Array.from(bonds), hbonds: Int32Array.from(res?.hbonds ?? []),
    stats: {
      t: sim.time, T: sim.temperature(), Ttarget: sim.T, thermostat: sim.thermostat, box: sim.box,
      Ekin, Epot: res?.E ?? 0, parts: res?.parts ?? null, Ewall: sim.Ewall, Egrab: sim.Egrab, Etot: sim.totalEnergy(),
      heatBath: sim.heatBath, work: sim.work, matterExchange: sim.matterExchange, diagnostics: sim.diagnostics(), P: sim.measurePressure(), dt: sim.dt,
      nMol: sim.frags?.length ?? 0, paused, stepsPerFrame, light, clamped: sim.clamped,
      forceField, hf: sim.provider?.info ?? null,
    },
    mb: speedHistogram(),
  };
  if (censusSent !== sim.censusSerial) {
    censusSent = sim.censusSerial;
    msg.census = {
      species: [...sim.species.entries()].sort((a, b) => b[1] - a[1]),
      history: sim.history.slice(-600),
      events: sim.events.filter(e => e.serial > lastEventSent),
      frags: sim.frags?.map(f => f.atoms) ?? [],
    };
    lastEventSent = sim.eventSerial;
  }
  postMessage(msg, [pos.buffer, q.buffer, ke.buffer]);
}

function loop() {
  const t0 = performance.now();
  try {
    if (!paused && sim.N) {
    let n = 0;
    while (n < stepsPerFrame && performance.now() - t0 < msTarget) {
      sim.step();
      n++;
      if (light.on) {
        // Impulsi Poisson, frequenza imposta: non una sezione di assorbimento.
        const p = light.rate * sim.dt / 1000;
        if (sim.rng.uni() < p) {
          const hit = sim.photon(light.lambda);
          if (hit) postMessage({ type: 'photon', ...hit, lambda: light.lambda });
        }
      }
    }
  }
  frame();
  } catch (e) {
    paused = true;
    postMessage({ type: 'error', text: e.message, paused: true });
  }
  setTimeout(loop, paused ? 60 : 0);
}

onmessage = (ev) => {
  const m = ev.data;
  try {
    switch (m.type) {
      case 'preset': {
        const r = loadPreset(m.preset);
        postMessage({ type: 'info', text: r.placed < r.wanted ? `Inserite ${r.placed} molecole su ${r.wanted}: la scatola è piena.` : '' });
        break;
      }
      case 'add': {
        const placed = add(m);
        postMessage({ type: 'info', text: placed < (m.count ?? 1) ? `Inserite ${placed} su ${m.count}: non c'è spazio libero.` : '' });
        break;
      }
      case 'clear': sim.clear(); forceField = 'reactive'; light.on = false; paused = true; censusSent = -1; lastEventSent = 0; mbHist = null; postMessage({ type: 'forcefield', kind: forceField }); break;
      case 'remove': sim.editInventory(() => sim.remove(m.indices)); lastEventSent = sim.eventSerial; break;
      case 'set':
        if (m.T !== undefined) { if (!Number.isFinite(m.T) || m.T < 0) throw new Error('Temperatura non valida.'); sim.T = m.T; }
        if (m.thermostat !== undefined) sim.thermostat = m.thermostat;
        if (m.tau !== undefined) { if (!Number.isFinite(m.tau) || m.tau <= 0) throw new Error('Tempo del termostato non valido.'); sim.tau = m.tau; }
        if (m.box !== undefined) sim.changeBox(m.box);
        if (m.dt !== undefined) { if (!Number.isFinite(m.dt) || m.dt < 0.001 || m.dt > 2) throw new Error('Δt ammesso: 0,001–2 fs.'); sim.dt = m.dt; }
        if (m.paused !== undefined) paused = m.paused;
        if (m.stepsPerFrame !== undefined) stepsPerFrame = m.stepsPerFrame;
        if (m.light) light = { ...light, ...m.light };
        if (m.mbSpecies !== undefined) { mbSpecies = m.mbSpecies; mbHist = null; }
        if (m.forceField !== undefined) setForceField(m.forceField, m.basis ?? hfBasis);
        break;
      case 'step': paused = true; for (let k = 0; k < Math.min(1000, m.n ?? 1); k++) sim.step(); sim.census(); break;
      case 'grab': sim.setGrab(m.i === null ? null : { i: m.i, target: m.target }); break;
      case 'spark': sim.spark(m.center, m.radius ?? 4, m.T ?? 8000); break;
      case 'photon': {
        const hit = sim.photon(m.lambda);
        if (hit) postMessage({ type: 'photon', ...hit, lambda: m.lambda });
        break;
      }
      case 'thermalize': sim.rethermalize(); break;
      case 'reset-measurements': sim.resetMeasurements(); censusSent = -1; lastEventSent = 0; break;
      default: break;
    }
  } catch (e) {
    paused = true;
    postMessage({ type: 'error', text: e.message, paused: true });
  }
};

loop();
