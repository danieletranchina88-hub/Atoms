// Dinamica molecolare ab initio (Born–Oppenheimer): a ogni passo si risolvono le equazioni di
// Hartree–Fock per gli elettroni e le forze sui nuclei vengono dal gradiente analitico dell'energia,
// F = −∂E/∂R. Nessun parametro empirico: solo l'equazione di Schrödinger (nell'approssimazione HF).

import { runHF, HARTREE_EV } from './hf.js';
import { hfGradient } from './gradient.js';
import { populationAnalysis } from './properties.js';
import { buildBasis } from './integrals.js';

const BOHR_ANG = 0.52917721090;
export const AIMD_MAX_ATOMS = 8;

/**
 * Fornitore di forze Hartree–Fock (UHF) per la classe Simulation.
 * Lo stato di spin è scelto una volta, all'inizio: il più stabile tra le due molteplicità più basse
 * (per O₂ il tripletto, come prevede la regola di Hund).
 */
export function makeHFProvider({ basis = 'STO-3G', fixedMultiplicity = null } = {}) {
  let guess = null;
  let multiplicity = null;
  let charge = 0;
  let key = '';
  let step = 0;
  const provider = {
    kind: 'hf',
    basis,
    info: {},
    fixedMultiplicity,
    reset() { guess = null; multiplicity = null; key = ''; step = 0; },
    setField(f) { if (f && f.some(v => v)) throw new Error('Il campo elettrico esterno è disponibile con MINDO/3, non con Hartree–Fock.'); },
    compute(Z, pos, F, ctx = {}) {
      if ((ctx.charge ?? 0) !== charge) { charge = ctx.charge ?? 0; guess = null; multiplicity = null; }
      const atoms = Z.map((z, i) => ({ Z: z, xyz: [pos[3 * i] / BOHR_ANG, pos[3 * i + 1] / BOHR_ANG, pos[3 * i + 2] / BOHR_ANG] }));
      const k = Array.from(Z).join(',');
      if (k !== key) { guess = null; multiplicity = null; key = k; }
      const nel = Z.reduce((s, z) => s + z, 0) - charge;
      const run = (mult, g, bs) => runHF(atoms, { basis, charge, unrestricted: true, multiplicity: mult, guess: g, breakSymmetry: bs, conv: 1e-9, maxIter: 150 });
      let res;
      if (provider.fixedMultiplicity && multiplicity === null) {
        const m = provider.fixedMultiplicity;
        if ((nel - (m - 1)) % 2 !== 0 || m - 1 > nel) throw new Error(`Molteplicità ${m} impossibile con ${nel} elettroni.`);
        res = run(m, null, m === 1);
        multiplicity = m;
      } else if (multiplicity === null) {
        // stato fondamentale di spin: si confrontano le due molteplicità più basse
        const low = nel % 2 === 0 ? 1 : 2;
        const cands = [];
        for (const m of [low, low + 2]) {
          if (m - 1 > nel) continue;
          try { const r = run(m, null, m === 1); if (r.converged && Number.isFinite(r.energy)) cands.push(r); } catch { /* molteplicità impossibile */ }
        }
        if (!cands.length) throw new Error('SCF non convergente: nessuna forza HF utilizzabile.');
        res = cands.reduce((a, b) => (b.energy < a.energy ? b : a));
        multiplicity = res.multiplicity;
      } else {
        res = run(multiplicity, guess, false);
        // ogni pochi passi (e se la SCF non converge) si prova una soluzione a simmetria di spin rotta:
        // serve per descrivere correttamente la rottura di un legame (per esempio H₂ → 2 H)
        if (!res.converged || (multiplicity === 1 && step % (Z.length <= 3 ? 1 : 10) === 0)) {
          try {
            const bs = run(multiplicity, null, true);
            if (!res.converged || bs.energy < res.energy - 1e-7) res = bs;
          } catch { /* resta la soluzione precedente */ }
        }
      }
      if (!res.converged || !Number.isFinite(res.energy)) throw new Error('SCF non convergente: dinamica sospesa.');
      step++;
      guess = { Pa: res.Pa, Pb: res.Pb };
      const grad = hfGradient(res);
      const f = HARTREE_EV / BOHR_ANG;
      for (let i = 0; i < Z.length; i++) for (let c = 0; c < 3; c++) F[3 * i + c] = -grad[i][c] * f;
      const pop = populationAnalysis(res);
      const bonds = [];
      for (let i = 0; i < Z.length; i++) {
        for (let j = i + 1; j < Z.length; j++) {
          const n = pop.bondOrder[i][j];
          if (n > 0.05) bonds.push({ i, j, n, w: 1 });
        }
      }
      provider.info = { multiplicity, S2: res.S2, converged: res.converged, iterations: res.iterations, nbf: res.n, basis, charge };
      provider.last = res;
      return {
        E: res.energy * HARTREE_EV,
        parts: { hf: res.energy * HARTREE_EV },
        q: Float64Array.from(pop.mulliken),
        spin: Float64Array.from(pop.spin),
        bonds,
        hbonds: [],
      };
    },
  };
  /** Funzione d'onda corrente: densità totale e orbitali di frontiera nella base gaussiana (per il disegno). */
  provider.wavefunction = () => {
    const r = provider.last;
    if (!r) return null;
    const n = r.n;
    const P = new Float64Array(n * n), Ps = new Float64Array(n * n);
    for (let k = 0; k < n * n; k++) { P[k] = r.Pa[k] + r.Pb[k]; Ps[k] = r.Pa[k] - r.Pb[k]; }
    const pick = (C, k, e, spin) => (k >= 0 && k < n ? { e, spin, c: Float64Array.from({ length: n }, (_, i) => C[i * n + k]) } : null);
    const hA = r.nalpha - 1, hB = r.nbeta - 1;
    const homo = hB >= 0 && r.epsB[hB] > r.epsA[hA] ? pick(r.Cb, hB, r.epsB[hB], 'β') : pick(r.Ca, hA, r.epsA[hA], 'α');
    const lumo = r.nbeta < n && r.epsB[r.nbeta] < (r.epsA[r.nalpha] ?? Infinity) ? pick(r.Cb, r.nbeta, r.epsB[r.nbeta], 'β') : pick(r.Ca, r.nalpha, r.epsA[r.nalpha], 'α');
    const orbitals = [];
    for (const [C, eps, occ, spin] of [[r.Ca, r.epsA, r.nalpha, 'α'], [r.Cb, r.epsB, r.nbeta, 'β']])
      for (let k = Math.max(0, occ - 5); k < Math.min(n, occ + 5); k++) orbitals.push({ ...pick(C, k, eps[k], spin), id: `${spin}:${k}`, index: k, occ: k < occ ? 1 : 0 });
    return { kind: 'gauss', Ps, orbitals, basisName: basis, atoms: r.atoms.map(a => ({ Z: a.Z, xyz: a.xyz.slice() })), n, P, homo, lumo, eUnit: HARTREE_EV };
  };
  return provider;
}

/** Controlla che il sistema sia trattabile con la dinamica ab initio. */
export function aimdFeasible(Z, basis = 'STO-3G') {
  if (!Z.length) return { ok: false, reason: 'La scatola è vuota.' };
  if (Z.length > AIMD_MAX_ATOMS) return { ok: false, reason: `La dinamica ab initio è possibile fino a ${AIMD_MAX_ATOMS} atomi (ce ne sono ${Z.length}): ogni passo richiede un calcolo Hartree–Fock completo.` };
  try {
    const b = buildBasis(Z.map(z => ({ Z: z, xyz: [0, 0, 0] })), basis);
    if (b.nbf > 40) return { ok: false, reason: `Troppe funzioni di base (${b.nbf}).` };
  } catch {
    return { ok: false, reason: 'Un elemento non è disponibile nella base scelta.' };
  }
  return { ok: true };
}
