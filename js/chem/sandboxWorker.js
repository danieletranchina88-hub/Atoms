// Worker della sandbox: fa girare la dinamica molecolare e invia al disegno 3D un fotogramma
// alla volta (posizioni, legami, cariche, grandezze termodinamiche, specie e reazioni).

import { prepareCollision } from './collision.js';
import { BondMonitor, TrajectoryRecorder } from './reactionTrace.js';
import { hillFormula } from './smiles.js';
import { Simulation, MV2, ACC } from './md.js';
import { ReactiveFF } from './reactive.js';
import { parseSmiles } from './smiles.js';
import { embedMolecule } from './embed.js';
import { analyzeStructure } from './structure.js';
import { KB_EV, ATOM_PARAMS } from './reactiveData.js';
import { makeHFProvider, aimdFeasible } from './aimd.js';
import { makeMindo3Provider, mindo3Supports } from './mindo3.js';
import { makeGFN2Provider, gfn2Feasible, GFN2_MAX_ATOMS, gfn2Electrons as gfn2El } from './xtb/provider.js';
import { GFN2xTB } from './xtb/gfn2.js';
import { RDF, MSD, HeatCapacity } from './mdAnalysis.js';
import { makeNobleLJProvider, nobleOnly, NOBLE_LJ } from './ljNoble.js';
import { localOrder, phaseVerdict, Mobility } from './structureOrder.js';

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
  return { Z: sim.Z.slice(), formal: sim.formal.slice(), pos: Array.from(sim.pos), vel: Array.from(sim.vel), box: sim.box, T: sim.T, thermostat: sim.thermostat, dt: sim.dt, forceField, barostat: { ...sim.barostat }, field: fieldVec, multiplicity: spinMult, solvent: solventName };
}

function restore(d) {
  sim = new Simulation({ box: d.box, T: d.T, dt: d.dt ?? 0.4 });
  sim.addAtoms(d.Z.map((z, i) => ({ Z: z, formal: d.formal?.[i] ?? 0, pos: d.pos.slice(3 * i, 3 * i + 3), vel: d.vel.slice(3 * i, 3 * i + 3) })));
  fieldVec = d.field ?? null; spinMult = d.multiplicity ?? null; solventName = d.solvent ?? null;
  sim.thermostat = d.thermostat;
  if (d.barostat) sim.barostat = { ...d.barostat };
  setForceField(d.forceField ?? 'reactive');
  sim.forces();
  sim.census();
  censusSent = -1; lastEventSent = 0; mbHist = null;
  resetAnalysis();
}

/**
 * Geometria di una molecola dal suo SMILES: modello VSEPR, poi minimizzazione. Le specie cariche (ioni, zwitterioni)
 * sono ammesse: la carica formale resta sugli atomi che la portano e la somma è la carica della scatola.
 * Gli ioni e le molecole con elementi fuori dal campo classico si rilassano con GFN2-xTB alla loro carica.
 */
function template(smiles, quantumGeometry = false) {
  const cacheKey = `${smiles}|${quantumGeometry}`;
  if (templates.has(cacheKey)) return templates.get(cacheKey);
  const g = parseSmiles(smiles);
  if (!g.atoms.length || g.atoms.some(a => !(a.Z >= 1 && a.Z <= 86)))
    throw new Error('Sandbox: elemento non supportato nella stringa SMILES (ammessi H–Rn).');
  const charge = g.atoms.reduce((s, a) => s + a.charge, 0);
  const an = analyzeStructure({ atoms: g.atoms, bonds: g.bonds, charge });
  const at = embedMolecule(g, an.atoms.map(a => a.lonePairs + a.radical), { seeds: 2, iterations: 1500 });
  const Z = at.map(a => a.Z);
  const formal = at.map((a, i) => (i < g.atoms.length && g.atoms[i].Z === a.Z ? g.atoms[i].charge : 0));
  if (formal.reduce((s, q) => s + q, 0) !== charge) { formal.fill(0); formal[0] = charge; }
  const pos = Float64Array.from(at.flatMap(a => a.xyz.map(v => v * BOHR_ANG)));
  if (Z.length > 1) {
    const quantum = quantumGeometry || charge !== 0 || !classicalOK(Z);
    const ff = quantum ? null : new ReactiveFF();
    const xtb = quantum ? new GFN2xTB({ etol: 1e-8, ptol: 1e-6 }) : null;
    const nel = quantum ? gfn2El(Z) - charge : 0;
    const F = new Float64Array(pos.length), v = new Float64Array(pos.length);
    let dt = 0.05, alpha = 0.1, npos = 0;
    for (let s = 0; s < (quantum ? 400 : 2500); s++) {
      if (quantum) {
        const r = xtb.compute(Z, Float64Array.from(pos, x => x / BOHR_ANG), { charge, uhf: ((nel % 2) + 2) % 2 });
        if (!r.converged || !Number.isFinite(r.energy) || !r.gradient.every(Number.isFinite)) throw new Error('GFN2-xTB: geometria iniziale non rilassabile con SCF convergente.');
        for (let k = 0; k < pos.length; k++) F[k] = -r.gradient[k] * 27.211386245988 / BOHR_ANG;
      } else ff.compute(Z, pos, F);
      let P = 0, vn = 0, fn = 0;
      for (let k = 0; k < pos.length; k++) { P += F[k] * v[k]; vn += v[k] * v[k]; fn += F[k] * F[k]; }
      if (Math.sqrt(fn) < (quantum ? 1e-3 : 1e-4)) break;
      vn = Math.sqrt(vn); fn = Math.sqrt(fn);
      for (let k = 0; k < pos.length; k++) v[k] = (1 - alpha) * v[k] + alpha * F[k] / fn * vn;
      if (P > 0) { if (++npos > 5) { dt = Math.min(dt * 1.1, quantum ? 0.1 : 0.3); alpha *= 0.99; } } else { npos = 0; dt *= 0.5; v.fill(0); alpha = 0.1; }
      for (let k = 0; k < pos.length; k++) { v[k] += dt * F[k]; pos[k] += dt * v[k]; }
    }
  }
  const t = { Z, pos, formal, charge };
  templates.set(cacheKey, t);
  return t;
}

function atomSmiles(sym) { return `[${sym}]`; }

/** Il campo classico reattivo ha parametri solo per alcuni elementi (H–Ca, Br, Kr, I, Xe). */
const classicalOK = (Z) => Z.every(z => !!ATOM_PARAMS[z]);

function add({ smiles, symbol, count = 1, at = null, T }) {
  const s = smiles ?? atomSmiles(symbol);
  const t = template(s, fidelity === 'auto' || forceField === 'gfn2');
  // con atomi che formano legami serve un passo corto (vibrazioni di 10 fs); i gas nobili tollerano 2 fs
  if (t.Z.some(z => (ATOM_PARAMS[z]?.[0] ?? 1) > 0)) sim.dt = Math.min(sim.dt, 0.4);
  if (forceField === 'hf') {
    const f = aimdFeasible([...sim.Z, ...Array.from({ length: count }, () => t.Z).flat()], hfBasis);
    if (!f.ok) throw new Error(f.reason);
  }
  const nextZ = [...sim.Z, ...Array.from({ length: count }, () => t.Z).flat()];
  if (forceField === 'mindo3') mindo3Check(nextZ);
  if (forceField === 'gfn2') gfn2Check(nextZ);
  // una carica netta, un campo elettrico o un elemento senza parametri classici richiedono un motore quantistico
  const nextCharge = sim.netCharge + count * t.charge;
  if (nextCharge !== 0 || t.charge !== 0 || fieldVec || !classicalOK(t.Z)) {
    const why = quantumFor(nextZ, nextCharge);
    if (!why.ok) throw new Error(why.reason);
    if (fidelity !== 'auto' && (forceField === 'reactive' || forceField === 'lj'))
      throw new Error(`${classicalOK(t.Z) ? 'Specie cariche e campo elettrico richiedono' : 'Gli elementi fuori dal campo classico richiedono'} un motore quantistico: scegli "fedeltà automatica" o GFN2-xTB.`);
    if (fidelity !== 'auto' && !classicalOK(nextZ) && forceField !== 'gfn2')
      throw new Error('Con elementi fuori da H–Ca, Br, Kr, I, Xe serve GFN2-xTB (copre tutti gli elementi fino al radon).');
  }
  mbHist = null;
  // il modello corrente deve poter calcolare le forze sulla scatola con i nuovi atomi
  if (fidelity === 'auto') {
    const pick = bestModel(nextZ, nextCharge);
    modelWhy = pick.why;
    if (pick.kind !== forceField) applyModel(pick.kind);
  }
  const placed = sim.editInventory(() => sim.addMolecule(t, count, T ?? sim.T, at));
  return placed;
}

let forceField = 'reactive';
let fidelity = 'auto';
let fieldVec = null;   // campo elettrico uniforme [Ex, Ey, Ez] in V/Å (GFN2-xTB o MINDO/3)
let spinMult = null;   // molteplicità di spin fissata (2S+1) oppure null = libera
let solventName = null; // solvente implicito ALPB attorno alle molecole ('water') oppure null (vuoto)
let hfBasis = 'STO-3G';
let modelWhy = 'Scatola vuota: aggiungi atomi o molecole.';

/** Sceglie il modello delle forze; con Hartree–Fock verifica che il sistema sia abbastanza piccolo. */
const MINDO3_MAX_ATOMS = 90;
const SYM = { 1: 'H', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl' };

/** MINDO/3 si usa solo se tutti gli elementi e le coppie hanno parametri pubblicati: altrimenti errore esplicito. */
function mindo3Check(Z) {
  if (!Z.length) return;
  const s = mindo3Supports(Z);
  if (!s.ok) {
    const names = s.missing.map(k => k.split('-').map(z => SYM[z] ?? `Z=${z}`).join('–')).join(', ');
    throw new Error(`MINDO/3 non ha parametri per ${names}: usa il campo classico. Il metodo copre H, B, C, N, O, F, P, S, Cl (non tutte le coppie).`);
  }
  if (Z.length > MINDO3_MAX_ATOMS) throw new Error(`Il calcolo quantistico MINDO/3 è limitato a ${MINDO3_MAX_ATOMS} atomi (ce ne sarebbero ${Z.length}): ogni passo richiede la diagonalizzazione dell'hamiltoniana.`);
}

/** GFN2-xTB copre H–Rn; il limite di atomi tiene la dinamica interattiva. */
function gfn2Check(Z) {
  if (!Z.length) return;
  const f = gfn2Feasible(Z);
  if (!f.ok) throw new Error(f.reason);
}

/** Esiste un motore quantistico per questa composizione e carica? (necessario per ioni, campo elettrico, metalli…) */
function quantumFor(Z, charge) {
  const g = gfn2Feasible(Z);
  if (g.ok) return { ok: true, kind: 'gfn2' };
  const what = charge ? `Carica totale ${charge > 0 ? '+' : ''}${charge}` : fieldVec ? 'Campo elettrico' : 'Elementi senza parametri classici';
  return { ok: false, reason: `${what}: serve un calcolo quantistico, ma ${g.reason}` };
}

const GFN2_AUTO_MAX = GFN2_MAX_ATOMS;
const GFN2_WHY = 'GFN2-xTB a ogni passo: tight binding quantistico con elettrostatica fino ai quadrupoli atomici e dispersione D4, validato su tblite.';

function bestModel(Z, charge = sim.netCharge) {
  if (Z.length && (charge !== 0 || fieldVec || solventName || !classicalOK(Z))) {
    const q = quantumFor(Z, charge);
    const what = charge ? `Carica totale ${charge > 0 ? '+' : ''}${charge}` : fieldVec ? 'Campo elettrico' : solventName ? 'Solvente implicito' : 'Elementi oltre il campo classico';
    if (q.ok) return { kind: q.kind, why: `${what}: ${GFN2_WHY}` };
    throw new Error(q.reason);
  }
  if (!Z.length) return { kind: 'reactive', why: 'Scatola vuota. Versa atomi o molecole: il modello si sceglie da solo.' };
  if (nobleOnly(Z)) {
    const names = [...new Set(Z)].map(z => NOBLE_LJ[z].symbol).join(', ');
    return { kind: 'lj', why: `${names}: Lennard–Jones con σ e ε pubblicati, taglio 2,5σ. Nessun legame chimico.` };
  }
  if (Z.length <= GFN2_AUTO_MAX) return { kind: 'gfn2', why: GFN2_WHY };
  throw new Error(`Fedeltà automatica: GFN2-xTB arriva a ${GFN2_MAX_ATOMS} atomi. Riduci il numero di molecole, oppure scegli esplicitamente il campo classico qualitativo.`);
}

function applyModel(kind, basis = hfBasis) {
  hfBasis = basis;
  if (kind === 'hf') {
    const f = aimdFeasible(sim.Z, basis);
    if (!f.ok) throw new Error(f.reason);
  }
  if (kind === 'mindo3') mindo3Check(sim.Z);
  if (kind === 'gfn2') gfn2Check(sim.Z);
  if (kind === 'lj' && sim.Z.length && !nobleOnly(sim.Z)) throw new Error('Il Lennard–Jones pubblicato vale solo per He, Ne, Ar, Kr e Xe.');
  if (!['hf', 'reactive', 'mindo3', 'lj', 'gfn2'].includes(kind)) throw new Error('Modello non valido.');
  if (kind !== 'gfn2' && !classicalOK(sim.Z)) throw new Error('Nella scatola ci sono elementi coperti solo da GFN2-xTB (oltre H–Ca, Br, Kr, I, Xe).');
  if ((kind === 'reactive' || kind === 'lj') && sim.netCharge !== 0) throw new Error(`La scatola ha carica ${sim.netCharge > 0 ? '+' : ''}${sim.netCharge}: il campo classico descrive solo frammenti neutri. Usa GFN2-xTB, MINDO/3 o Hartree–Fock.`);
  if (kind !== 'gfn2' && solventName) throw new Error('Il solvente implicito ALPB è disponibile con GFN2-xTB: toglilo o scegli GFN2-xTB.');
  if (kind !== 'mindo3' && kind !== 'gfn2' && fieldVec) throw new Error('Il campo elettrico esterno è calcolato con GFN2-xTB o MINDO/3: spegnilo o scegli uno dei due.');
  const same = kind === forceField && (kind === 'reactive' || sim.provider);
  forceField = kind;
  if (!same) {
    sim.provider = kind === 'hf' ? makeHFProvider({ basis, fixedMultiplicity: spinMult })
      : kind === 'mindo3' ? makeMindo3Provider({ multiplicity: spinMult, field: fieldVec })
      : kind === 'gfn2' ? makeGFN2Provider({ multiplicity: spinMult, field: fieldVec, solvent: solventName })
      : kind === 'lj' ? makeNobleLJProvider()
      : null;
  }
  sim.res = null;
  if (kind === 'hf') sim.dt = Math.min(sim.dt, 0.25);
  if (kind === 'mindo3' || kind === 'gfn2') sim.dt = Math.min(sim.dt, 0.4);
  if (kind === 'lj') sim.dt = Math.min(Math.max(sim.dt, 0.5), 2);
  sim.resetMeasurements(); censusSent = -1; lastEventSent = 0; mbHist = null;
  postMessage({ type: 'forcefield', kind, fidelity, why: modelWhy });
}

function setForceField(kind, basis = hfBasis) {
  if (kind === 'auto') {
    const pick = bestModel(sim.Z);
    fidelity = 'auto';
    modelWhy = pick.why;
    applyModel(pick.kind, basis);
    return;
  }
  fidelity = kind;
  modelWhy = kind === 'lj' ? 'Lennard–Jones pubblicato, scelto a mano.'
    : kind === 'mindo3' ? 'MINDO/3 scelto a mano.'
    : kind === 'gfn2' ? 'GFN2-xTB scelto a mano.'
    : kind === 'hf' ? 'Hartree–Fock scelto a mano.'
    : 'Campo classico scelto a mano: utile per scatole grandi, non quantitativo sulle barriere.';
  applyModel(kind, basis);
}

/**
 * Accende, cambia o spegne il campo elettrico uniforme (V/Å). Cambiare il campo cambia l'energia potenziale:
 * la differenza è lavoro esterno sul sistema, contabilizzato come tale (non deriva numerica).
 */
function setField(f) {
  const v = f && f.some(x => x) ? f.map(Number) : null;
  if (v && (v.length !== 3 || !v.every(Number.isFinite) || Math.hypot(...v) > 10)) throw new Error('Campo elettrico: tre componenti finite, modulo al massimo 10 V/Å.');
  const old = fieldVec;
  if (!sim.res && sim.N) sim.forces();
  const before = sim.res?.E ?? 0;
  fieldVec = v;
  try {
    if (fidelity === 'auto') refreshFidelity();
    else if (v && forceField !== 'mindo3' && forceField !== 'gfn2') throw new Error('Il campo elettrico esterno è calcolato con GFN2-xTB o MINDO/3: sceglili o usa la fedeltà automatica.');
    sim.provider?.setField?.(v);
    sim.res = null;
    if (sim.N) { sim.forces(); sim.work += sim.res.E - before; }
  } catch (e) {
    fieldVec = old; sim.provider?.setField?.(old); sim.res = null;
    throw e;
  }
}

/**
 * Mette le molecole in un solvente implicito (ALPB, acqua) o le riporta nel vuoto. Il cambio di energia è lavoro esterno,
 * come per il campo elettrico. Richiede GFN2-xTB.
 */
function setSolvent(name) {
  const v = name ? String(name) : null;
  if (v && v !== 'water') throw new Error('Solvente implicito disponibile: acqua (ALPB).');
  const old = solventName;
  if (!sim.res && sim.N) sim.forces();
  const before = sim.res?.E ?? 0;
  solventName = v;
  try {
    if (fidelity === 'auto') refreshFidelity();
    else if (v && forceField !== 'gfn2') throw new Error('Il solvente implicito ALPB è calcolato con GFN2-xTB: sceglilo o usa la fedeltà automatica.');
    sim.provider?.setSolvent?.(v);
    sim.res = null;
    if (sim.N) { sim.forces(); sim.work += sim.res.E - before; }
  } catch (e) {
    solventName = old; sim.provider?.setSolvent?.(old); sim.res = null;
    throw e;
  }
}

/** Fissa la molteplicità di spin 2S+1 (null = spin libero, livello di Fermi comune). Il cambio è lavoro esterno. */
function setMultiplicity(mult) {
  const v = mult === null || mult === 0 ? null : Number(mult);
  if (v !== null && (!Number.isInteger(v) || v < 1 || v > 7)) throw new Error('Molteplicità ammessa: 1–7.');
  if (!sim.res && sim.N) sim.forces();
  const before = sim.res?.E ?? 0, old = spinMult;
  spinMult = v;
  try {
    if (sim.provider && 'multiplicity' in sim.provider) { sim.provider.multiplicity = v; sim.provider.reset(); }
    if (sim.provider && 'fixedMultiplicity' in sim.provider) { sim.provider.fixedMultiplicity = v; sim.provider.reset(); }
    sim.res = null;
    if (sim.N) { sim.forces(); sim.work += sim.res.E - before; }
  } catch (e) {
    spinMult = old;
    if (sim.provider && 'multiplicity' in sim.provider) sim.provider.multiplicity = old;
    if (sim.provider && 'fixedMultiplicity' in sim.provider) sim.provider.fixedMultiplicity = old;
    sim.provider?.reset(); sim.res = null;
    throw e;
  }
}

function refreshFidelity() {
  if (fidelity !== 'auto') return;
  const pick = bestModel(sim.Z);
  modelWhy = pick.why;
  if (pick.kind !== forceField) applyModel(pick.kind);
  else postMessage({ type: 'forcefield', kind: forceField, fidelity, why: modelWhy });
}

function loadPreset(p) {
  sim = new Simulation({ box: p.box, T: p.T, dt: p.dt ?? 0.2, seed: p.seed ?? 2024 });
  fieldVec = p.field ?? null; spinMult = p.multiplicity ?? null; solventName = p.solvent ?? null;
  censusSent = -1; lastEventSent = 0; mbHist = null;
  stepsPerFrame = p.speed ?? 40;
  if (p.atoms) {
    sim.addAtoms(p.atoms);
    if (p.thermalize) for (let i = 0; i < sim.N; i++) sim.thermalize(i, p.thermalize);
    sim.thermostat = p.thermostat ?? false;
    if (p.forceField) setForceField(p.forceField);
    else refreshFidelity();
    sim.forces();
    sim.resetMeasurements();
    light = { ...light, on: false };
    lastEventSent = 0;
    mbHist = null;
    resetAnalysis();
    return { placed: p.atoms.length, wanted: p.atoms.length };
  }
  if (p.forceField) setForceField(p.forceField);
  else refreshFidelity();
  if (p.thermostat !== undefined) sim.thermostat = p.thermostat;
  let placed = 0, wanted = 0;
  for (const [s, n] of p.add) {
    wanted += n;
    placed += add({ smiles: s, count: n });
  }
  if (p.forceField) setForceField(p.forceField);
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
  resetAnalysis();
  refreshFidelity();
  return { placed, wanted };
}

/** Controlled binary encounter: only initial geometry and COM velocities are prescribed. */
function loadCollision(m) {
  const kind = m.model ?? 'gfn2';
  if (!['gfn2', 'hf'].includes(kind)) throw new Error('Urti controllati: scegli GFN2-xTB o Hartree–Fock.');
  const left = template(m.left, true), right = template(m.right, true);
  const Z = [...left.Z, ...right.Z];
  if (kind === 'hf') { const test = aimdFeasible(Z); if (!test.ok) throw new Error(test.reason); }
  else gfn2Check(Z);
  const collision = prepareCollision(left, right, m);
  const saved = { sim, forceField, fidelity, fieldVec, spinMult, solventName, modelWhy };
  try { loadPreset({ atoms: collision.atoms, box: collision.box, T: 300, dt: 0.1, speed: 2, forceField: kind, thermostat: false, multiplicity: m.multiplicity ?? null, solvent: m.solvent ?? null }); } catch (e) {
    ({ sim, forceField, fidelity, fieldVec, spinMult, solventName, modelWhy } = saved);
    censusSent = -1; phaseCache = null;
    postMessage({ type: 'forcefield', kind: forceField, fidelity, why: modelWhy });
    throw e;
  }
  paused = true;
  postMessage({ type: 'collision-ready', collision: { energy: collision.energy, distance: collision.distance, impact: collision.impact, speed: collision.speed, reducedMass: collision.reducedMass, leftCount: left.Z.length, rightCount: right.Z.length } });
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


/** Indicazione strutturale, non un diagramma di fase sperimentale: vicini entro 4,2 Å. */
function classifyPhase() {
  const N = sim.N;
  if (!N) return { label: 'vuota', title: 'Scatola vuota', frac: 0, clusters: 0, T: 0, note: 'Aggiungi atomi o molecole.' };
  if (forceField === 'lj' && nobleOnly(sim.Z)) {
    // struttura locale (q₆ di Steinhardt) entro 1,5σ, il primo minimo di g(r): non soglie su T* e ρ*,
    // perché con pareti e coesistenza la densità media della scatola non è quella delle fasi
    const sig = sim.Z.reduce((s, z) => s + NOBLE_LJ[z].sigma, 0) / N;
    const eps = sim.Z.reduce((s, z) => s + NOBLE_LJ[z].epsilonK, 0) / N;
    const o = localOrder(sim.pos, { cutoff: 1.5 * sig });
    // tempo caratteristico τ = σ √(m/ε) in fs (m in u, ε in eV): per l'argon 2,16 ps
    const m = sim.mass.reduce((a, b) => a + b, 0) / N;
    const tau = sig * Math.sqrt(m / (eps * KB_EV) / ACC);
    const v = phaseVerdict(o.fraction, mobility.sample(sim.pos, sim.time, tau, sig));
    const label = { solido: 'solido', liquido: 'liquido', gas: 'gas', misto: 'misto' }[v.id];
    return { label, title: v.name, frac: 1 - o.fraction.vapor, clusters: 0, T: sim.temperature(), note: `${v.why}. T* = kT/ε = ${(sim.temperature() / eps).toFixed(2)}${mobility.lostFraction != null ? ` · primi vicini persi in 3τ: ${Math.round(100 * mobility.lostFraction)}%` : ''} · struttura locale q₆ (ten Wolde–Frenkel).` };
  }
  // una fase (gas, liquido, solido) è una proprietà collettiva: con poche molecole non ha senso
  const nFrag = sim.frags?.length ?? 0;
  if (nFrag && nFrag < 8) return { label: 'pochi', title: `${nFrag} ${nFrag === 1 ? 'molecola' : 'molecole'}: nessuna fase`, frac: 0, clusters: nFrag, T: sim.temperature(), note: 'Gas, liquido e solido descrivono molte particelle: qui si seguono le singole molecole.' };
  const pos = sim.pos;
  const heavy = [];
  for (let i = 0; i < N; i++) if (sim.Z[i] !== 1) heavy.push(i);
  const use = heavy.length ? heavy : Array.from({ length: N }, (_, i) => i);
  const c2 = 4.2 * 4.2;
  const neigh = new Int16Array(N);
  for (let a = 0; a < use.length; a++) {
    const i = use[a];
    for (let b = a + 1; b < use.length; b++) {
      const j = use[b];
      const dx = pos[3 * i] - pos[3 * j], dy = pos[3 * i + 1] - pos[3 * j + 1], dz = pos[3 * i + 2] - pos[3 * j + 2];
      if (dx * dx + dy * dy + dz * dz < c2) { neigh[i]++; neigh[j]++; }
    }
  }
  let coordinated = 0;
  for (const i of use) if (neigh[i] >= 4) coordinated++;
  const frac = coordinated / use.length;
  const big = (sim.frags ?? []).filter(f => f.atoms.some(i => sim.Z[i] !== 1) && f.atoms.length >= 4).length;
  const T = sim.temperature();
  let label = 'gas', title = 'Gas';
  if (frac >= 0.55 && T < 180) { label = 'solido'; title = 'Condensato freddo'; }
  else if (frac >= 0.38) { label = 'liquido'; title = 'Liquido o aggregato'; }
  else if (frac >= 0.12 || big >= 2) { label = 'misto'; title = 'Fase mista'; }
  return { label, title, frac, clusters: big, T, note: 'Vicini entro 4,2 Å nel modello. Non è la fase sperimentale.' };
}

// la lettura della fase costa O(N²): si aggiorna ogni 10 fotogrammi o quando cambia il numero di atomi
let phaseCache = null, phaseFrame = 0;
const mobility = new Mobility();
function cachedPhase() {
  if (!phaseCache || phaseCache.N !== sim.N || ++phaseFrame % 10 === 0) phaseCache = { N: sim.N, value: classifyPhase() };
  return phaseCache.value;
}

function stateFrame() {
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
    N, Z: Int16Array.from(sim.Z), pos, q, ke, vel: Float32Array.from(sim.vel), forces: Float32Array.from(sim.F), spin: res?.spin ? Float32Array.from(res.spin) : null, bonds: Float32Array.from(bonds), hbonds: Int32Array.from(res?.hbonds ?? []),
    stats: {
      t: sim.time, T: sim.temperature(), Ttarget: sim.T, thermostat: sim.thermostat, box: sim.box,
      Ekin, Epot: res?.E ?? 0, parts: res?.parts ?? null, Ewall: sim.Ewall, Egrab: sim.Egrab, Etot: sim.totalEnergy(),
      heatBath: sim.heatBath, work: sim.work, matterExchange: sim.matterExchange, diagnostics: sim.diagnostics(), P: sim.measurePressure(), dt: sim.dt,
      nMol: sim.frags?.length ?? 0, paused, stepsPerFrame, light, clamped: sim.clamped,
      forceField, fidelity, modelWhy, hf: sim.provider?.info ?? null,
      charge: sim.netCharge, field: fieldVec, multiplicity: spinMult, solvent: solventName, dipole: sim.provider?.info?.dipole ?? null,
    },
  };
  return msg;
}

function frame() {
  recordTrace();
  const msg = stateFrame();
  msg.stats.phase = cachedPhase();
  msg.mb = speedHistogram();
  msg.trace = { ...recorder.meta(), config: { ...traceConfig }, Z: sim.Z.slice(), samples: pendingSamples, events: traceEvents.slice(-100) };
  pendingSamples = [];
  if (censusSent !== sim.censusSerial) {
    censusSent = sim.censusSerial;
    msg.census = {
      species: [...sim.species.entries()].sort((a, b) => b[1] - a[1]),
      history: sim.history.slice(-600),
      events: sim.events.filter(e => e.serial > lastEventSent),
      frags: sim.frags?.map(f => f.atoms) ?? [],
    };
    lastEventSent = sim.eventSerial;
    // le analisi si ricalcolano ogni 5 censimenti (la normalizzazione di g(r) costa qualche millisecondo)
    if (physSent < 0 || sim.history.length - physSent >= 5 || sim.history.length < physSent) { msg.phys = physics(); physSent = sim.history.length; }
  }
  postMessage(msg, [msg.pos.buffer, msg.q.buffer, msg.ke.buffer, msg.vel.buffer, msg.forces.buffer]);
}

// All bond changes are measured on accepted steps, including changes within a connected fragment.
const recorder = new TrajectoryRecorder();
const bondMonitor = new BondMonitor();
const traceConfig = { enabled: true, stride: null, waves: true, pauseOnBond: false };
let inventoryVersion = 0, traceKey = '', pendingSamples = [], traceEvents = [], lastTraceStep = -1, lastStepChanged = false;

function traceCensus(f) {
  const parent = Array.from({ length: f.N }, (_, i) => i);
  const find = i => { while (parent[i] !== i) i = parent[i]; return i; };
  for (const key of bondMonitor.connected) {
    const [i, j] = key.split(':').map(Number); parent[find(i)] = find(j);
  }
  const groups = new Map();
  for (let i = 0; i < f.N; i++) { const root = find(i); if (!groups.has(root)) groups.set(root, []); groups.get(root).push(i); }
  const counts = new Map();
  for (const ids of groups.values()) {
    const q = ids.reduce((v, i) => v + f.q[i], 0), rounded = Math.round(q);
    const formula = hillFormula(ids.map(i => ({ Z: f.Z[i] })), Math.abs(q - rounded) < 0.35 ? rounded : 0);
    counts.set(formula, (counts.get(formula) ?? 0) + 1);
  }
  return { frags: [...groups.values()], species: [...counts], events: [], history: [] };
}

function recordTrace() {
  const key = `${inventoryVersion}|${forceField}|${fieldVec}|${solventName}|${spinMult}|${sim.netCharge}|${sim.Z.join(',')}`;
  if (key !== traceKey) {
    traceKey = key; recorder.reset('Nuovo inventario o modello elettronico'); bondMonitor.reset();
    pendingSamples = []; traceEvents = []; lastTraceStep = -1; phaseCache = null; physSent = -1;
  }
  if (lastTraceStep === sim.stepCount) return;
  if (!sim.res) sim.forces();
  const f = stateFrame();
  const changes = bondMonitor.sample({ t: sim.time, pos: f.pos, q: f.q, bonds: f.bonds });
  lastStepChanged = changes.length > 0;
  traceEvents.push(...changes);
  if (traceEvents.length > 1000) traceEvents.splice(0, traceEvents.length - 1000);
  if (changes.length && traceConfig.pauseOnBond) paused = true;
  const stride = traceConfig.stride ?? (sim.N <= 8 ? 1 : sim.N <= 40 ? 5 : 20);
  if (traceConfig.enabled && (lastTraceStep < 0 || sim.stepCount % stride === 0 || changes.length)) {
    f.stats.phase = cachedPhase();
    f.census = traceCensus(f); f.stats.nMol = f.census.frags.length;
    const wave = traceConfig.waves && sim.N <= 40 ? sim.provider?.wavefunction?.() ?? null : null;
    const sample = recorder.add(f, wave);
    const { wave: unused, bytes, ...small } = sample;
    pendingSamples.push(small);
    if (pendingSamples.length > 600) pendingSamples.shift();
  }
  lastTraceStep = sim.stepCount;
}

/**
 * Un passo di velocity Verlet con rigetto: se lo spostamento previsto è instabile il passo non viene eseguito,
 * Δt si dimezza e si riprova. Δt non viene mai aumentato da solo: crescere oltre il valore scelto (per esempio
 * fino a 1–2 fs con atomi di idrogeno) rompe la conservazione dell'energia.
 */
const DT_MIN = 0.01;
function safeStep() {
  recordTrace();
  for (;;) {
    try { sim.step(); recordTrace(); return; } catch (e) {
      if (!/instabile/.test(e.message) || sim.dt / 2 < DT_MIN) throw e;
      sim.dt /= 2;
      postMessage({ type: 'info', text: `Passo instabile: Δt dimezzato a ${sim.dt.toLocaleString('it-IT', { maximumSignificantDigits: 3 })} fs e passo ripetuto. Nessuna velocità è stata tagliata.` });
    }
  }
}

function loop() {
  const t0 = performance.now();
  try {
    if (!paused && sim.N) {
    let n = 0;
    while (n < stepsPerFrame && performance.now() - t0 < msTarget) {
      safeStep();
      n++;
      if (sim.stepCount % sim.censusEvery === 0) analysisSample();
      if (paused) break;
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
      case 'collision': loadCollision(m); inventoryVersion++; break;
      case 'preset': {
        inventoryVersion++;
        const r = loadPreset(m.preset);
        postMessage({ type: 'info', text: r.placed < r.wanted ? `Inserite ${r.placed} molecole su ${r.wanted}: la scatola è piena.` : '' });
        break;
      }
      case 'add': {
        const placed = add(m);
        if (placed) inventoryVersion++;
        refreshFidelity();
        postMessage({ type: 'info', text: placed < (m.count ?? 1) ? `Inserite ${placed} su ${m.count}: non c'è spazio libero. ${modelWhy}` : modelWhy });
        break;
      }
      case 'clear': inventoryVersion++; sim.clear(); fieldVec = null; spinMult = null; solventName = null; light.on = false; paused = false; censusSent = -1; lastEventSent = 0; mbHist = null; resetAnalysis(); fidelity = 'auto'; refreshFidelity(); break;
      case 'remove': sim.editInventory(() => sim.remove(m.indices)); inventoryVersion++; lastEventSent = sim.eventSerial; refreshFidelity(); break;
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
        if (m.field !== undefined) setField(m.field);
        if (m.solvent !== undefined) setSolvent(m.solvent);
        if (m.multiplicity !== undefined) setMultiplicity(m.multiplicity);
        if (m.barostat) sim.barostat = { ...sim.barostat, ...m.barostat };
        if (m.rdfPair !== undefined) { rdfPairUser = m.rdfPair; rdf.reset(m.rdfPair); physSent = -1; }
        break;
      case 'resetAnalysis': resetAnalysis(); physSent = -1; break;
      case 'wave': {
        // funzione d'onda corrente per il disegno della densità e degli orbitali
        const sample = m.traceId !== undefined ? recorder.get(m.traceId) : null;
        if (m.traceId !== undefined && !sample) throw new Error('Fotogramma non più nel buffer.');
        const w = sample ? sample.wave : sim.provider?.wavefunction?.() ?? null;
        postMessage({ type: 'wave', reqId: m.reqId, purpose: m.purpose, epoch: recorder.epoch, t: sample?.stats.t ?? sim.time, box: sample?.stats.box ?? sim.box, Z: sample?.Z ?? sim.Z.slice(), pos: sample?.pos ?? Float64Array.from(sim.pos), wave: w, forceField: sample?.stats.forceField ?? forceField });
        break;
      }
      case 'trace-config': {
        if (m.stride !== undefined && m.stride !== null && (!Number.isInteger(m.stride) || m.stride < 1 || m.stride > 100)) throw new Error('Campionamento: 1–100 passi o automatico.');
        for (const key of ['enabled', 'waves', 'pauseOnBond', 'stride']) if (key in m) traceConfig[key] = m[key];
        break;
      }
      case 'trace-reset': traceKey = ''; break;
      case 'trace-inspect': {
        if (m.epoch !== recorder.epoch) break;
        const sample = recorder.get(m.id);
        if (!sample) throw new Error('Fotogramma non più nel buffer.');
        paused = true;
        postMessage({ type: 'trace-frame', reqId: m.reqId, frame: { ...sample, stats: { ...sample.stats, paused: true }, replay: true } });
        break;
      }
      case 'trace-export': postMessage({ type: 'trace-export', data: { version: 1, kind: 'accepted-md-trajectory', units: { pos: 'angstrom', vel: 'angstrom/fs', forces: 'eV/angstrom', t: 'fs', energy: 'eV', charge: 'e', spin: 'n_alpha-n_beta' }, meta: recorder.meta(), config: traceConfig, events: traceEvents, frames: recorder.frames.map(({ wave, bytes, ...f }) => f) } }); break;
      case 'snapshot': postMessage({ type: 'snapshot', name: m.name, data: snapshot() }); break;
      case 'restore': inventoryVersion++; restore(m.data); physSent = -1; break;
      case 'step': paused = true; for (let k = 0; k < Math.min(1000, m.n ?? 1); k++) { safeStep(); if (traceConfig.pauseOnBond && lastStepChanged) break; } sim.census(); break;
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
