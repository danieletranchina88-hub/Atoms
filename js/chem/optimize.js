// Ottimizzazione della geometria: minimo dell'energia sulla superficie di Born–Oppenheimer.
// Quasi-Newton in coordinate cartesiane con:
//  • Hessiana iniziale modello di Lindh (Chem. Phys. Lett. 241, 423, 1995), costruita da tutte le
//    distanze, gli angoli e i diedri pesati con ρ_ij = exp[α_ij (r_ref² − r_ij²)];
//  • aggiornamento BFGS dell'Hessiana;
//  • passo RFO (rational function optimization) limitato da un raggio di fiducia.

import { runHF } from './hf.js';
import { hfGradient } from './gradient.js';
import { eigh } from './linalg.js';

const rowOf = (Z) => (Z <= 2 ? 0 : Z <= 10 ? 1 : 2);
const ALPHA = [[1.0, 0.3949, 0.3949], [0.3949, 0.28, 0.28], [0.3949, 0.28, 0.28]];
const RREF = [[1.35, 2.10, 2.53], [2.10, 2.87, 3.40], [2.53, 3.40, 3.40]];

function dist(x, i, j) {
  return Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]);
}

export function angle(x, i, j, k) {
  const u = [x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]];
  const v = [x[3 * k] - x[3 * j], x[3 * k + 1] - x[3 * j + 1], x[3 * k + 2] - x[3 * j + 2]];
  const c = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (Math.hypot(...u) * Math.hypot(...v));
  return Math.acos(Math.max(-1, Math.min(1, c)));
}

export function dihedral(x, i, j, k, l) {
  const p = (a) => [x[3 * a], x[3 * a + 1], x[3 * a + 2]];
  const [a, b, c, d] = [p(i), p(j), p(k), p(l)];
  const sub = (u, v) => [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
  const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const b1 = sub(b, a), b2 = sub(c, b), b3 = sub(d, c);
  const n1 = cross(b1, b2), n2 = cross(b2, b3);
  const m1 = cross(n1, b2);
  const b2n = Math.hypot(...b2);
  return Math.atan2(dot(m1, n2) / b2n, dot(n1, n2));
}

/** Hessiana modello di Lindh in coordinate cartesiane (3N × 3N). */
export function lindhHessian(atoms, x) {
  const N = atoms.length;
  const n3 = 3 * N;
  const H = new Float64Array(n3 * n3);
  const rho = Array.from({ length: N }, () => new Float64Array(N));
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      if (i === j) continue;
      const a = rowOf(atoms[i].Z), b = rowOf(atoms[j].Z);
      const r = dist(x, i, j);
      rho[i][j] = Math.exp(ALPHA[a][b] * (RREF[a][b] ** 2 - r * r));
    }
  }
  const addTerm = (k, atomsIdx, fn) => {
    // riga di Wilson B per differenze finite della coordinata interna
    const h = 1e-5;
    const b = new Float64Array(n3);
    for (const a of atomsIdx) {
      for (let c = 0; c < 3; c++) {
        const xp = Float64Array.from(x); xp[3 * a + c] += h;
        const xm = Float64Array.from(x); xm[3 * a + c] -= h;
        let d = fn(xp) - fn(xm);
        if (d > Math.PI) d -= 2 * Math.PI;
        if (d < -Math.PI) d += 2 * Math.PI;
        b[3 * a + c] = d / (2 * h);
      }
    }
    for (const a of atomsIdx) for (let c = 0; c < 3; c++) {
      const bi = b[3 * a + c];
      if (bi === 0) continue;
      for (const a2 of atomsIdx) for (let c2 = 0; c2 < 3; c2++) H[(3 * a + c) * n3 + 3 * a2 + c2] += k * bi * b[3 * a2 + c2];
    }
  };
  for (let i = 0; i < N; i++) {
    for (let j = i + 1; j < N; j++) {
      const k = 0.45 * rho[i][j];
      if (k < 1e-4) continue;
      addTerm(k, [i, j], (xx) => dist(xx, i, j));
    }
  }
  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      for (let k = i + 1; k < N; k++) {
        if (i === j || k === j) continue;
        const kf = 0.15 * rho[i][j] * rho[j][k];
        if (kf < 1e-4) continue;
        if (Math.sin(angle(x, i, j, k)) < 0.02) continue; // angoli lineari: coordinate mal definite
        addTerm(kf, [i, j, k], (xx) => angle(xx, i, j, k));
      }
    }
  }
  for (let j = 0; j < N; j++) {
    for (let k = 0; k < N; k++) {
      if (k === j) continue;
      if (rho[j][k] < 1e-2) continue;
      for (let i = 0; i < N; i++) {
        if (i === j || i === k) continue;
        for (let l = 0; l < N; l++) {
          if (l === i || l === j || l === k || l < i) continue;
          const kt = 0.005 * rho[i][j] * rho[j][k] * rho[k][l];
          if (kt < 1e-4) continue;
          if (Math.sin(angle(x, i, j, k)) < 0.05 || Math.sin(angle(x, j, k, l)) < 0.05) continue;
          addTerm(kt, [i, j, k, l], (xx) => dihedral(xx, i, j, k, l));
        }
      }
    }
  }
  return H;
}

/** Passo RFO: minimizza il modello quadratico nello spazio degli autovettori con autovalore non nullo. */
function rfoStep(H, g, n, trust) {
  const { values, vectors } = eigh(H, n);
  // componenti del gradiente lungo gli autovettori (si scartano traslazioni e rotazioni, λ ≈ 0)
  const modes = [];
  for (let k = 0; k < n; k++) {
    if (Math.abs(values[k]) < 1e-7) continue;
    let gk = 0;
    for (let i = 0; i < n; i++) gk += vectors[i * n + k] * g[i];
    modes.push({ k, lam: values[k], g: gk });
  }
  // autovalore RFO: λ tale che λ = Σ g_k² / (λ − h_k), cercato sotto il minimo h_k
  const hmin = Math.min(...modes.map(m => m.lam), 0);
  let lo = hmin - 1e4, hi = hmin - 1e-12;
  const f = (lam) => lam - modes.reduce((s, m) => s + m.g * m.g / (lam - m.lam), 0);
  for (let it = 0; it < 200; it++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) > 0) hi = mid; else lo = mid;
  }
  const shift = 0.5 * (lo + hi);
  const step = new Float64Array(n);
  for (const m of modes) {
    const c = -m.g / (m.lam - shift);
    for (let i = 0; i < n; i++) step[i] += c * vectors[i * n + m.k];
  }
  // raggio di fiducia sullo spostamento massimo di un atomo
  let maxDisp = 0;
  for (let a = 0; a < n / 3; a++) maxDisp = Math.max(maxDisp, Math.hypot(step[3 * a], step[3 * a + 1], step[3 * a + 2]));
  if (maxDisp > trust) for (let i = 0; i < n; i++) step[i] *= trust / maxDisp;
  return step;
}

/**
 * Ottimizza la geometria.
 * @param atoms [{Z, xyz}] (bohr)
 * @param opts opzioni HF + { maxSteps, onStep(info) }
 */
export function optimizeGeometry(atoms, opts = {}) {
  const N = atoms.length;
  const n3 = 3 * N;
  let x = new Float64Array(n3);
  atoms.forEach((a, i) => { x[3 * i] = a.xyz[0]; x[3 * i + 1] = a.xyz[1]; x[3 * i + 2] = a.xyz[2]; });
  const make = (xx) => atoms.map((a, i) => ({ ...a, xyz: [xx[3 * i], xx[3 * i + 1], xx[3 * i + 2]] }));
  let H = lindhHessian(atoms, x);
  let trust = 0.3;
  let guess = null;
  let res = runHF(make(x), { ...opts, guess });
  let E = res.energy;
  let g = Float64Array.from(hfGradient(res).flat());
  const history = [];
  const maxSteps = opts.maxSteps ?? 60;
  let converged = false;
  for (let step = 0; step <= maxSteps; step++) {
    const gmax = Math.max(...g.map(Math.abs));
    const grms = Math.sqrt(g.reduce((s, v) => s + v * v, 0) / n3);
    history.push({ step, energy: E, gmax, grms });
    opts.onStep?.({ step, energy: E, gmax, grms, atoms: make(x) });
    if (gmax < (opts.gmax ?? 4.5e-4) && grms < (opts.grms ?? 3e-4)) { converged = true; break; }
    if (step === maxSteps) break;
    const dx = rfoStep(H, g, n3, trust);
    // energia prevista dal modello quadratico
    let pred = 0;
    for (let i = 0; i < n3; i++) {
      pred += g[i] * dx[i];
      let hs = 0;
      for (let j = 0; j < n3; j++) hs += H[i * n3 + j] * dx[j];
      pred += 0.5 * dx[i] * hs;
    }
    const xn = new Float64Array(n3);
    for (let i = 0; i < n3; i++) xn[i] = x[i] + dx[i];
    guess = { Pa: res.Pa, Pb: res.Pb };
    const resN = runHF(make(xn), { ...opts, guess });
    const En = resN.energy;
    const gn = Float64Array.from(hfGradient(resN).flat());
    // aggiornamento del raggio di fiducia
    const ratio = pred !== 0 ? (En - E) / pred : 1;
    if (ratio < 0.25) trust = Math.max(0.03, trust * 0.5);
    else if (ratio > 0.75) trust = Math.min(0.5, trust * 1.5);
    if (En > E + 1e-6 && trust > 0.03) {
      // passo rifiutato: si resta nel punto precedente con raggio ridotto
      trust *= 0.5;
      continue;
    }
    // BFGS: H ← H + y yᵀ/(yᵀs) − (Hs)(Hs)ᵀ/(sᵀHs)
    const y = new Float64Array(n3);
    for (let i = 0; i < n3; i++) y[i] = gn[i] - g[i];
    const Hs = new Float64Array(n3);
    for (let i = 0; i < n3; i++) { let s = 0; for (let j = 0; j < n3; j++) s += H[i * n3 + j] * dx[j]; Hs[i] = s; }
    const ys = y.reduce((s, v, i) => s + v * dx[i], 0);
    const sHs = Hs.reduce((s, v, i) => s + v * dx[i], 0);
    if (ys > 1e-10 && sHs > 1e-10) {
      for (let i = 0; i < n3; i++) for (let j = 0; j < n3; j++) H[i * n3 + j] += y[i] * y[j] / ys - Hs[i] * Hs[j] / sHs;
    }
    x = xn; E = En; g = gn; res = resN;
  }
  return { atoms: make(x), energy: E, result: res, converged, history };
}
