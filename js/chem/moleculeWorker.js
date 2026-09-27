// Web Worker per la chimica quantistica molecolare: calcoli pesanti fuori dal thread dell'interfaccia.
// Messaggi: { id, type, ...payload } → risposte { id, result } | { id, error } | { id, progress }.

import { runHF, HARTREE_EV } from './hf.js';
import { hfGradient } from './gradient.js';
import { optimizeGeometry } from './optimize.js';
import { harmonicFrequencies, thermochemistry, ISOTOPE_MASS } from './vibrations.js';
import { populationAnalysis, dipoleMoment, mp2Energy, orbitalSummary, orbitalComposition, electrostaticPotential, moCharacter } from './properties.js';
import { connectivityFromGeometry } from './structure.js';
import { inertia } from './symmetry.js';
import { analyzeSymmetry } from './symmetry.js';
import { basisValues } from './integrals.js';

let current = null; // ultimo risultato SCF (con base e integrali), per i calcoli su griglia
let currentId = null;
const sendProgress = (p) => self.postMessage({ id: currentId, progress: p });

function summarize(res, withMP2) {
  const pop = populationAnalysis(res);
  const dip = dipoleMoment(res);
  const orb = orbitalSummary(res);
  let mp2 = null;
  if (withMP2) {
    try { mp2 = mp2Energy(res).energy; } catch { mp2 = null; }
  }
  const sym = res.atoms.length > 1 ? analyzeSymmetry(res.atoms, res.atoms.map(a => ISOTOPE_MASS[a.Z] ?? 2 * a.Z)) : { pointGroup: 'K_h', sigma: 1 };
  // composizione di ciascun orbitale molecolare (per la tabella degli OM)
  const compositions = [];
  const nShow = Math.min(res.n, Math.max(res.nalpha, res.nbeta) + 8);
  // geometria: lineare (asse = momento d'inerzia minimo) o planare (normale = momento massimo)
  const bonds = connectivityFromGeometry(res.atoms);
  let geom = null;
  if (res.atoms.length > 1) {
    const { moments, axes } = inertia(res.atoms, res.atoms.map(a => ISOTOPE_MASS[a.Z] ?? 2 * a.Z));
    const planar = Math.abs(moments[2] - moments[0] - moments[1]) < 1e-3 * Math.max(moments[2], 1);
    if (moments[0] < 1e-3 * moments[2]) geom = { kind: 'linear', axis: axes[0] };
    else if (planar) geom = { kind: 'planar', axis: axes[2] };
  }
  for (let k = 0; k < nShow; k++) {
    compositions.push({
      alpha: { ...orbitalComposition(res, k, 'α'), ...moCharacter(res, k, 'α', bonds, geom) },
      beta: res.unrestricted ? { ...orbitalComposition(res, k, 'β'), ...moCharacter(res, k, 'β', bonds, geom) } : null,
    });
  }
  return {
    energy: res.energy, Enuc: res.Enuc, converged: res.converged, iterations: res.iterations,
    nbf: res.n, basis: res.basis.name, nel: res.nel, nalpha: res.nalpha, nbeta: res.nbeta,
    multiplicity: res.multiplicity, charge: res.charge, unrestricted: res.unrestricted, S2: res.S2,
    epsA: Array.from(res.epsA), epsB: Array.from(res.epsB),
    mulliken: pop.mulliken, lowdin: pop.lowdin, spin: pop.spin,
    bondOrder: pop.bondOrder.map(r => Array.from(r)), valence: pop.valence,
    dipole: { vector: dip.vector, debye: dip.debye },
    homo: orb.homo, lumo: orb.lumo, gapEV: orb.gapEV, koopmansIE: orb.koopmansIE,
    mp2,
    pointGroup: sym.pointGroup, sigma: sym.sigma,
    compositions,
    time: res.time,
    history: res.history,
  };
}

/** Griglia cubica: valori di un OM o della densità elettronica. */
function gridField({ kind, index, spin, res: N, half, center }) {
  const r = current;
  const n = r.n;
  const vals = new Float32Array(N * N * N);
  const phi = new Float64Array(n);
  const step = 2 * half / (N - 1);
  const C = spin === 'β' ? r.Cb : r.Ca;
  // orbitali occupati per la densità
  const occ = [];
  if (kind === 'density') {
    for (let k = 0; k < r.nalpha; k++) occ.push([r.Ca, k]);
    for (let k = 0; k < r.nbeta; k++) occ.push([r.Cb, k]);
  }
  let idx = 0;
  for (let kz = 0; kz < N; kz++) {
    const z = center[2] - half + kz * step;
    for (let ky = 0; ky < N; ky++) {
      const y = center[1] - half + ky * step;
      for (let kx = 0; kx < N; kx++) {
        const x = center[0] - half + kx * step;
        basisValues(r.basis, x, y, z, phi);
        if (kind === 'mo') {
          let s = 0;
          for (let m = 0; m < n; m++) s += C[m * n + index] * phi[m];
          vals[idx++] = s;
        } else {
          let rho = 0;
          for (const [Cm, k] of occ) {
            let s = 0;
            for (let m = 0; m < n; m++) s += Cm[m * n + k] * phi[m];
            rho += s * s;
          }
          vals[idx++] = rho;
        }
      }
    }
    if ((kz & 7) === 0) sendProgress(kz / N);
  }
  return { values: vals, res: N, half, step, center };
}

self.onmessage = (ev) => {
  const { id, type } = ev.data;
  currentId = id;
  const progress = (p, extra) => self.postMessage({ id, progress: p, ...extra });
  try {
    let result;
    if (type === 'scf') {
      const { atoms, opts } = ev.data;
      const res = runHF(atoms, opts);
      current = res;
      result = summarize(res, opts.mp2 !== false && res.n <= 110);
    } else if (type === 'optimize') {
      const { atoms, opts } = ev.data;
      const r = optimizeGeometry(atoms, {
        ...opts,
        onStep: (s) => progress(null, { step: s.step, energy: s.energy, gmax: s.gmax, atoms: s.atoms }),
      });
      current = r.result;
      result = { atoms: r.atoms, converged: r.converged, history: r.history, summary: summarize(r.result, opts.mp2 !== false && r.result.n <= 110) };
    } else if (type === 'frequencies') {
      const { atoms, opts, T } = ev.data;
      const f = harmonicFrequencies(atoms, { ...opts, onProgress: (p) => progress(p) });
      const th = thermochemistry(atoms, f.modes, f.reference.energy, { T: T ?? 298.15, multiplicity: f.reference.multiplicity });
      result = {
        modes: f.modes.map(m => ({ freq: m.freq, ir: m.ir, mu: m.reducedMass, d: Array.from(m.displacement) })),
        thermo: th,
      };
    } else if (type === 'grid') {
      if (!current) throw new Error('Nessun calcolo disponibile');
      const g = gridField(ev.data);
      result = g;
    } else if (type === 'esp') {
      if (!current) throw new Error('Nessun calcolo disponibile');
      result = { values: electrostaticPotential(current, ev.data.points, sendProgress) };
    } else if (type === 'scan') {
      // curva di energia potenziale di una biatomica: RHF e UHF a confronto
      const { Z1, Z2, distances, opts } = ev.data;
      const out = [];
      let guessR = null, guessU = null;
      distances.forEach((d, k) => {
        const atoms = [{ Z: Z1, xyz: [0, 0, 0] }, { Z: Z2, xyz: [0, 0, d] }];
        const rh = runHF(atoms, { ...opts, guess: guessR ?? undefined });
        guessR = { Pa: rh.Pa, Pb: rh.Pb };
        let uh = null;
        try {
          // a ogni distanza si parte da una soluzione a simmetria rotta (miscela HOMO–LUMO):
          // oltre il punto di Coulson–Fischer UHF scende sotto RHF e dissocia correttamente
          uh = runHF(atoms, { ...opts, unrestricted: true, breakSymmetry: opts.multiplicity === 1 || !opts.multiplicity });
          if (guessU) {
            const alt = runHF(atoms, { ...opts, unrestricted: true, guess: guessU });
            if (alt.energy < uh.energy) uh = alt;
          }
          guessU = { Pa: uh.Pa, Pb: uh.Pb };
        } catch { uh = null; }
        out.push({ d, rhf: rh.energy, uhf: uh ? uh.energy : null, S2: uh ? uh.S2 : null });
        progress((k + 1) / distances.length);
      });
      result = { points: out };
    } else if (type === 'gradient') {
      if (!current) throw new Error('Nessun calcolo disponibile');
      result = { gradient: hfGradient(current) };
    }
    self.postMessage({ id, result });
  } catch (err) {
    self.postMessage({ id, error: String(err?.message ?? err) });
  }
};

export { HARTREE_EV };
