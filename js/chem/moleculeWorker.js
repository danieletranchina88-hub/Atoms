// Web Worker per la chimica quantistica molecolare: calcoli pesanti fuori dal thread dell'interfaccia.
// Messaggi: { id, type, ...payload } → risposte { id, result } | { id, error } | { id, progress }.

import { runHF, HARTREE_EV } from './hf.js';
import { hfGradient } from './gradient.js';
import { optimizeGeometry } from './optimize.js';
import { harmonicFrequencies, thermochemistry, ISOTOPE_MASS } from './vibrations.js';
import { populationAnalysis, dipoleMoment, mp2Energy, orbitalSummary, orbitalComposition } from './properties.js';
import { analyzeSymmetry } from './symmetry.js';
import { basisValues, shellPair, hermiteR } from './integrals.js';

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
  for (let k = 0; k < nShow; k++) {
    compositions.push({ alpha: orbitalComposition(res, k, 'α'), beta: res.unrestricted ? orbitalComposition(res, k, 'β') : null });
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

/**
 * Potenziale elettrostatico molecolare esatto nei punti dati:
 *   V(r) = Σ_A Z_A/|r − R_A| − Σ_μν P_μν ∫ φ_μ φ_ν / |r − r'| dr'
 */
function electrostaticPotential(points) {
  const r = current;
  const { basis, atoms, n } = r;
  const P = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) P[i] = r.Pa[i] + r.Pb[i];
  const pairs = [];
  const { shells } = basis;
  for (let A = 0; A < shells.length; A++) {
    for (let B = 0; B <= A; B++) {
      const pr = shellPair(shells[A], shells[B]);
      const sa = shells[A], sb = shells[B];
      // pesi di densità per le componenti della coppia (fattore 2 fuori diagonale)
      const w = [];
      let wmax = 0;
      sa.comps.forEach((ca, ia) => sb.comps.forEach((cb, ib) => {
        const mu = sa.offset + ia, nu = sb.offset + ib;
        if (A === B && nu > mu) return;
        const f = (A === B && mu === nu) ? 1 : 2;
        const v = f * P[mu * n + nu] * sa.compNorm[ia] * sb.compNorm[ib];
        w.push([ca, cb, v]);
        wmax = Math.max(wmax, Math.abs(v));
      }));
      if (wmax < 1e-10) continue;
      pairs.push({ pr, w, L: sa.l + sb.l });
    }
  }
  const out = new Float32Array(points.length / 3);
  for (let q = 0; q < out.length; q++) {
    const cx = points[3 * q], cy = points[3 * q + 1], cz = points[3 * q + 2];
    let v = 0;
    for (const a of atoms) {
      const d = Math.hypot(cx - a.xyz[0], cy - a.xyz[1], cz - a.xyz[2]);
      v += a.Z / Math.max(d, 1e-6);
    }
    let el = 0;
    for (const { pr, w, L } of pairs) {
      const D = L + 1;
      const lb = pr.lb, T = pr.T;
      for (const p of pr.prims) {
        if (p.K < 1e-14) continue;
        const R = hermiteR(L, p.p, p.P[0] - cx, p.P[1] - cy, p.P[2] - cz);
        let s = 0;
        for (const [ca, cb, wv] of w) {
          let t0 = 0;
          for (let t = 0; t <= ca[0] + cb[0]; t++) {
            const ex = p.Ex[(ca[0] * (lb + 1) + cb[0]) * T + t];
            if (ex === 0) continue;
            for (let u = 0; u <= ca[1] + cb[1]; u++) {
              const ey = p.Ey[(ca[1] * (lb + 1) + cb[1]) * T + u];
              if (ey === 0) continue;
              for (let vv = 0; vv <= ca[2] + cb[2]; vv++) {
                const ez = p.Ez[(ca[2] * (lb + 1) + cb[2]) * T + vv];
                if (ez === 0) continue;
                t0 += ex * ey * ez * R[(t * D + u) * D + vv];
              }
            }
          }
          s += wv * t0;
        }
        el += p.c * 2 * Math.PI / p.p * s;
      }
    }
    out[q] = v - el;
    if ((q & 255) === 0) sendProgress(q / out.length);
  }
  return out;
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
      result = { values: electrostaticPotential(ev.data.points) };
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
          uh = runHF(atoms, { ...opts, unrestricted: true, breakSymmetry: !guessU, guess: guessU ?? undefined });
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
