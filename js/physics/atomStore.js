// Gestione dei calcoli: esecuzione nel Web Worker (con ripiego sul thread principale) e cache per Z.

import { makeGrid, solveRadial, radialStats, enclosingRadius } from './scf.js';
import { RadialFunction } from './wavefunction.js';

const cache = new Map();
const pending = new Map();
let worker = null;
let nextId = 1;

function getWorker() {
  if (worker !== null) return worker;
  try {
    worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = (ev) => {
      const { id, atom, error } = ev.data;
      const p = pending.get(id);
      if (!p) return;
      pending.delete(id);
      if (error) p.reject(new Error(error)); else p.resolve(atom);
    };
    worker.onerror = () => {
      // Worker non disponibile: si ripiega sul calcolo nel thread principale.
      for (const p of pending.values()) p.fallback();
      pending.clear();
      worker = false;
    };
  } catch {
    worker = false;
  }
  return worker;
}

async function computeInMainThread(Z) {
  const { computeAtom } = await import('./atom.js');
  // lascia il tempo al browser di mostrare l'indicatore di calcolo
  await new Promise(r => setTimeout(r, 30));
  return computeAtom(Z);
}

function decorate(raw) {
  const grid = makeGrid(raw.Z);
  const atom = { ...raw, fullGrid: grid, virtual: new Map() };
  for (const o of atom.orbitals) {
    o.radial = new RadialFunction(grid, o.u, o.l);
    o.occupied = true;
  }
  // Potenziale di Latter per gli orbitali non occupati: coda −(Z−N+1)/r
  const Nel = raw.config.reduce((s, c) => s + c.occ, 0);
  const zion = raw.Z - Nel + 1;
  atom.latter = new Float64Array(grid.N);
  for (let i = 0; i < grid.N; i++) atom.latter[i] = Math.min(raw.V[i], -zion / grid.r[i]);
  return atom;
}

/** Restituisce (con cache) l'atomo calcolato per Z. */
export function loadAtom(Z) {
  if (cache.has(Z)) return cache.get(Z);
  const promise = new Promise((resolve, reject) => {
    const w = getWorker();
    const fallback = () => computeInMainThread(Z).then(resolve, reject);
    if (!w) { fallback(); return; }
    const id = nextId++;
    pending.set(id, { resolve, reject, fallback });
    w.postMessage({ id, Z });
  }).then(decorate);
  cache.set(Z, promise);
  promise.catch(() => cache.delete(Z));
  return promise;
}

/**
 * Orbitale (n, l) dell'atomo: quello autoconsistente se occupato, altrimenti uno stato "virtuale"
 * calcolato nel potenziale dell'atomo con la corretta coda coulombiana.
 */
export function getOrbital(atom, n, l) {
  const occ = atom.orbitals.find(o => o.n === n && o.l === l);
  if (occ) return occ;
  const key = `${n},${l}`;
  if (atom.virtual.has(key)) return atom.virtual.get(key);
  const grid = atom.fullGrid;
  const res = solveRadial(grid, atom.latter, n, l, undefined, atom.Z);
  const stats = radialStats(grid, res.u);
  const orb = {
    n, l, occ: 0, e: res.e, u: res.u,
    label: `${n}${'spdfghi'[l]}`,
    rAvg: stats.rAvg, rMaxProb: stats.rMaxProb, nodes: stats.nodes,
    r90: enclosingRadius(grid, res.u, 0.9),
    radial: new RadialFunction(grid, res.u, l),
    occupied: false,
    converged: res.converged,
  };
  atom.virtual.set(key, orb);
  return orb;
}
