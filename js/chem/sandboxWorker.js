// Worker della sandbox: fa girare la dinamica molecolare e invia al disegno 3D un fotogramma
// alla volta (posizioni, legami, cariche, grandezze termodinamiche, specie e reazioni).

import { Simulation, MV2 } from './md.js';
import { ReactiveFF } from './reactive.js';
import { parseSmiles } from './smiles.js';
import { embedMolecule } from './embed.js';
import { analyzeStructure } from './structure.js';
import { KB_EV, ATOM_PARAMS } from './reactiveData.js';
import { makeHFProvider, aimdFeasible } from './aimd.js';
import { RDF, MSD, HeatCapacity } from './mdAnalysis.js';

const BOHR_ANG = 0.52917721090;
let sim = new Simulation({ box: 20, T: 300 });
let paused = false;
let stepsPerFrame = 40;
let light = { on: false, lambda: 400, rate: 2 };   // fotoni per picosecondo
let lastLightTime = 0;
let msTarget = 14;
let mbSpecies = null;
let mbHist = null;
let lastEventSent = 0;
let censusSent = -1;
const templates = new Map();
// chimica fisica: g(r), spostamento quadratico medio, capacità termica
const rdf = new RDF(10, 100);
const msd = new MSD();
const heatCap = new HeatCapacity();
let rdfPairUser = null;
let physSent = -1;

function resetAnalysis() {
  rdf.reset(rdfPairUser);
  msd.reset(sim);
  heatCap.reset();
}

/** Coppia predefinita per g(r): l'elemento più abbondante diverso dall'idrogeno, con sé stesso. */
function defaultPair() {
  const counts = new Map();
  for (const z of sim.Z) counts.set(z, (counts.get(z) ?? 0) + 1);
  const ranked = [...counts].sort((a, b) => (a[0] === 1) - (b[0] === 1) || b[1] - a[1]);
  return ranked.length ? [ranked[0][0], ranked[0][0]] : null;
}

function analysisSample() {
  const present = new Set(sim.Z);
  if (!rdf.pair || !present.has(rdf.pair[0]) || !present.has(rdf.pair[1])) rdf.reset(rdfPairUser && present.has(rdfPairUser[0]) && present.has(rdfPairUser[1]) ? rdfPairUser : defaultPair());
  rdf.accumulate(sim);
  msd.sample(sim);
  // C_V dalle fluttuazioni vale nell'insieme canonico: termostato acceso, volume e composizione fissi
  if (sim.thermostat && !sim.barostat.on) heatCap.sample(sim, `${sim.N}|${sim.T}|${sim.box}|${forceField}`);
  else heatCap.reset();
}

function physics() {
  const counts = new Map();
  for (const z of sim.Z) counts.set(z, (counts.get(z) ?? 0) + 1);
  return {
    rdf: rdf.result(),
    msd: msd.result(sim.box),
    cv: sim.thermostat && !sim.barostat.on ? heatCap.result(sim.T) : null,
    elements: [...counts].sort((a, b) => b[1] - a[1]).map(([z]) => z),
    nMol: sim.frags?.length ?? 0, N: sim.N, T: sim.T, box: sim.box,
    barostat: sim.barostat,
  };
}

function snapshot() {
  return { Z: sim.Z.slice(), pos: Array.from(sim.pos), vel: Array.from(sim.vel), box: sim.box, T: sim.T, thermostat: sim.thermostat, dt: sim.dt, forceField, barostat: { ...sim.barostat } };
}

function restore(d) {
  sim = new Simulation({ box: d.box, T: d.T, dt: d.dt ?? 0.4 });
  sim.addAtoms(d.Z.map((z, i) => ({ Z: z, pos: d.pos.slice(3 * i, 3 * i + 3), vel: d.vel.slice(3 * i, 3 * i + 3) })));
  sim.thermostat = d.thermostat;
  if (d.barostat) sim.barostat = { ...d.barostat };
  setForceField(d.forceField ?? 'reactive');
  sim.forces();
  sim.census();
  lastEventSent = 0; mbHist = null;
  resetAnalysis();
}

/** Geometria di una molecola dal suo SMILES: modello VSEPR, poi minimizzazione con il campo reattivo. */
function template(smiles) {
  if (templates.has(smiles)) return templates.get(smiles);
  const g = parseSmiles(smiles);
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
  const placed = sim.addMolecule(t, count, T ?? sim.T, at);
  sim.res = null;
  if (forceField === 'hf') setForceField('hf');
  sim.census();
  return placed;
}

let forceField = 'reactive';
let hfBasis = 'STO-3G';

/** Sceglie il modello delle forze; con Hartree–Fock verifica che il sistema sia abbastanza piccolo. */
function setForceField(kind, basis = hfBasis) {
  hfBasis = basis;
  if (kind === 'hf') {
    const f = aimdFeasible(sim.Z, basis);
    if (!f.ok) { postMessage({ type: 'info', text: f.reason }); kind = 'reactive'; }
  }
  forceField = kind;
  sim.provider = kind === 'hf' ? makeHFProvider({ basis }) : null;
  sim.res = null;
  if (kind === 'hf') sim.dt = Math.min(sim.dt, 0.25);
  postMessage({ type: 'forcefield', kind });
}

function loadPreset(p) {
  sim = new Simulation({ box: p.box, T: p.T, dt: p.dt ?? 0.4 });
  stepsPerFrame = p.speed ?? 40;
  if (p.atoms) {
    sim.addAtoms(p.atoms);
    if (p.thermalize) for (let i = 0; i < sim.N; i++) sim.thermalize(i, p.thermalize);
    sim.thermostat = p.thermostat ?? false;
    setForceField(p.forceField ?? 'reactive');
    sim.forces();
    sim.census();
    light = { ...light, on: false };
    lastEventSent = 0;
    mbHist = null;
    resetAnalysis();
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
  sim.census();
  resetAnalysis();
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
  const vmp = Math.sqrt(2 * KB_EV * sim.T / (M * MV2)); // velocità più probabile (Å/fs)
  const vmax = 3.2 * vmp;
  const nb = 30;
  if (!mbHist || mbHist.species !== species || Math.abs(mbHist.vmax - vmax) / vmax > 0.02) {
    mbHist = { species, vmax, M, bins: new Float64Array(nb), samples: 0 };
  }
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
  return { species, M, vmax, bins: Array.from(mbHist.bins), samples: mbHist.samples, count: frags.length, T: sim.T };
}

function frame() {
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
      Ekin, Epot: res?.E ?? 0, parts: res?.parts ?? null, Ewall: sim.Ewall, Etot: sim.totalEnergy(),
      heatBath: sim.heatBath, work: sim.work, P: sim.measurePressure(), dt: sim.dt,
      nMol: sim.frags?.length ?? 0, paused, stepsPerFrame, light, clamped: sim.clamped,
      forceField, hf: sim.provider?.info ?? null,
    },
    mb: speedHistogram(),
  };
  if (censusSent !== sim.history.length) {
    censusSent = sim.history.length;
    msg.census = {
      species: [...sim.species.entries()].sort((a, b) => b[1] - a[1]),
      history: sim.history.slice(-600),
      events: sim.events.slice(lastEventSent),
      frags: sim.frags?.map(f => f.atoms) ?? [],
    };
    lastEventSent = sim.events.length;
    // le analisi si ricalcolano ogni 5 censimenti (la normalizzazione di g(r) costa qualche millisecondo)
    if (physSent < 0 || sim.history.length - physSent >= 5 || sim.history.length < physSent) { msg.phys = physics(); physSent = sim.history.length; }
  }
  postMessage(msg, [pos.buffer, q.buffer, ke.buffer]);
}

function loop() {
  const t0 = performance.now();
  if (!paused && sim.N) {
    let n = 0;
    while (n < stepsPerFrame && performance.now() - t0 < msTarget) {
      sim.step();
      n++;
      if (sim.stepCount % sim.censusEvery === 0) analysisSample();
      if (light.on) {
        // fotoni: processo di Poisson con `rate` assorbimenti per picosecondo
        const p = light.rate * sim.dt / 1000;
        if (sim.rng.uni() < p) {
          const hit = sim.photon(light.lambda);
          if (hit) postMessage({ type: 'photon', ...hit, lambda: light.lambda });
        }
      }
    }
  }
  frame();
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
      case 'clear': sim.clear(); lastEventSent = 0; mbHist = null; resetAnalysis(); break;
      case 'remove': sim.remove(m.indices); sim.provider?.reset(); sim.census(); lastEventSent = sim.events.length; break;
      case 'set':
        if (m.T !== undefined) sim.T = m.T;
        if (m.thermostat !== undefined) sim.thermostat = m.thermostat;
        if (m.tau !== undefined) sim.tau = m.tau;
        if (m.box !== undefined) sim.box = m.box;
        if (m.dt !== undefined) sim.dt = m.dt;
        if (m.paused !== undefined) paused = m.paused;
        if (m.stepsPerFrame !== undefined) stepsPerFrame = m.stepsPerFrame;
        if (m.light) light = { ...light, ...m.light };
        if (m.mbSpecies !== undefined) { mbSpecies = m.mbSpecies; mbHist = null; }
        if (m.forceField !== undefined) setForceField(m.forceField, m.basis ?? hfBasis);
        if (m.barostat) sim.barostat = { ...sim.barostat, ...m.barostat };
        if (m.rdfPair !== undefined) { rdfPairUser = m.rdfPair; rdf.reset(m.rdfPair); physSent = -1; }
        break;
      case 'resetAnalysis': resetAnalysis(); physSent = -1; break;
      case 'snapshot': postMessage({ type: 'snapshot', name: m.name, data: snapshot() }); break;
      case 'restore': restore(m.data); physSent = -1; break;
      case 'step': for (let k = 0; k < (m.n ?? 1); k++) sim.step(); break;
      case 'grab': sim.grab = m.i === null ? null : { i: m.i, target: m.target }; break;
      case 'spark': sim.spark(m.center, m.radius ?? 4, m.T ?? 8000); break;
      case 'photon': {
        const hit = sim.photon(m.lambda);
        if (hit) postMessage({ type: 'photon', ...hit, lambda: m.lambda });
        break;
      }
      case 'thermalize': for (let i = 0; i < sim.N; i++) sim.thermalize(i, sim.T); break;
      default: break;
    }
  } catch (e) {
    postMessage({ type: 'error', text: e.message });
  }
};

loop();
