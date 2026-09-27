// Vibrazioni molecolari e termodinamica statistica (gas ideale, rotore rigido, oscillatore armonico).
//
// Hessiana: H_ij = ∂²E/∂x_i∂x_j per differenze centrali dei gradienti analitici.
// Modi normali: autovettori della Hessiana pesata con le masse, F_ij = H_ij / √(m_i m_j),
// dopo aver proiettato via traslazioni e rotazioni (condizioni di Eckart).
//   ν̃_k = √λ_k / (2πc)
// Intensità IR (approssimazione del doppio armonico):  I_k ∝ |∂μ/∂Q_k|²

import { runHF } from './hf.js';
import { hfGradient } from './gradient.js';
import { dipoleMoment } from './properties.js';
import { eigh } from './linalg.js';
import { analyzeSymmetry, inertia } from './symmetry.js';

// Masse dell'isotopo più abbondante (u), Z = 1…36 (AME 2016)
export const ISOTOPE_MASS = [0,
  1.00782503, 4.00260325, 7.01600344, 9.0121831, 11.0093054, 12.0, 14.0030740, 15.9949146, 18.9984032, 19.9924402,
  22.9897693, 23.9850417, 26.9815385, 27.9769265, 30.9737620, 31.9720712, 34.9688527, 39.9623831, 38.9637065, 39.9625909,
  44.9559083, 47.9479409, 50.9439570, 51.9405063, 54.9380439, 55.9349363, 58.9331943, 57.9353424, 62.9295977, 63.9291420,
  68.9255735, 73.9211778, 74.9215946, 79.9165218, 78.9183376, 83.9114977,
];

// Costanti fisiche CODATA 2018
const H_PLANCK = 6.62607015e-34;
const K_B = 1.380649e-23;
const C_LIGHT = 2.99792458e10;     // cm/s
const N_A = 6.02214076e23;
const R_GAS = 8.314462618;         // J/(mol K)
const AMU = 1.66053906660e-27;
const BOHR_M = 0.52917721090e-10;
const HARTREE_J = 4.3597447222071e-18;
export const AU_TO_CM1 = Math.sqrt(HARTREE_J / (BOHR_M * BOHR_M * AMU)) / (2 * Math.PI * C_LIGHT); // ≈ 5140,49
const IR_KM_MOL = 974.8801;        // (e²/amu) → km/mol

/**
 * Calcolo delle frequenze armoniche.
 * @param atoms geometria (bohr), preferibilmente ottimizzata
 * @param opts opzioni HF + { step (bohr), onProgress }
 */
export function harmonicFrequencies(atoms, opts = {}) {
  const N = atoms.length;
  const n3 = 3 * N;
  const h = opts.step ?? 0.005;
  const masses = atoms.map(a => ISOTOPE_MASS[a.Z] ?? 2 * a.Z);
  const ref = runHF(atoms, opts);
  const guess = { Pa: ref.Pa, Pb: ref.Pb };
  const Hs = new Float64Array(n3 * n3);
  const dMu = Array.from({ length: n3 }, () => [0, 0, 0]); // ∂μ/∂x_i
  let done = 0;
  for (let i = 0; i < n3; i++) {
    const g = [];
    const mu = [];
    for (const s of [1, -1]) {
      const disp = atoms.map((a, k) => ({ ...a, xyz: a.xyz.map((v, c) => (3 * k + c === i ? v + s * h : v)) }));
      const r = runHF(disp, { ...opts, guess });
      g.push(hfGradient(r).flat());
      mu.push(dipoleMoment(r).vector);
      opts.onProgress?.(++done / (2 * n3));
    }
    for (let j = 0; j < n3; j++) Hs[i * n3 + j] = (g[0][j] - g[1][j]) / (2 * h);
    for (let c = 0; c < 3; c++) dMu[i][c] = (mu[0][c] - mu[1][c]) / (2 * h);
  }
  // simmetrizza
  for (let i = 0; i < n3; i++) for (let j = 0; j < i; j++) {
    const v = 0.5 * (Hs[i * n3 + j] + Hs[j * n3 + i]);
    Hs[i * n3 + j] = Hs[j * n3 + i] = v;
  }
  return { ...normalModes(atoms, masses, Hs, dMu), reference: ref, hessian: Hs, masses };
}

/** Modi normali da Hessiana cartesiana, con proiezione di traslazioni e rotazioni. */
export function normalModes(atoms, masses, Hs, dMu) {
  const N = atoms.length;
  const n3 = 3 * N;
  const F = new Float64Array(n3 * n3);
  const sm = [];
  for (let a = 0; a < N; a++) for (let c = 0; c < 3; c++) sm.push(Math.sqrt(masses[a]));
  for (let i = 0; i < n3; i++) for (let j = 0; j < n3; j++) F[i * n3 + j] = Hs[i * n3 + j] / (sm[i] * sm[j]);

  // vettori di traslazione e rotazione (pesati con le masse), ortonormalizzati
  const { center } = inertia(atoms, masses);
  const trv = [];
  for (let c = 0; c < 3; c++) {
    const v = new Float64Array(n3);
    for (let a = 0; a < N; a++) v[3 * a + c] = sm[3 * a];
    trv.push(v);
  }
  for (let c = 0; c < 3; c++) {
    const v = new Float64Array(n3);
    for (let a = 0; a < N; a++) {
      const r = [0, 1, 2].map(k => atoms[a].xyz[k] - center[k]);
      const e = [0, 0, 0]; e[c] = 1;
      const cr = [e[1] * r[2] - e[2] * r[1], e[2] * r[0] - e[0] * r[2], e[0] * r[1] - e[1] * r[0]];
      for (let k = 0; k < 3; k++) v[3 * a + k] = cr[k] * sm[3 * a];
    }
    trv.push(v);
  }
  const basis = [];
  for (const v of trv) {
    const w = Float64Array.from(v);
    for (const b of basis) { const d = w.reduce((s, x, i) => s + x * b[i], 0); for (let i = 0; i < n3; i++) w[i] -= d * b[i]; }
    const l = Math.hypot(...w);
    if (l > 1e-6) basis.push(w.map(x => x / l));
  }
  const nExt = basis.length; // 5 per molecole lineari, 6 altrimenti
  // proiettore P = 1 − Σ b bᵀ,  F' = P F P
  const P = new Float64Array(n3 * n3);
  for (let i = 0; i < n3; i++) P[i * n3 + i] = 1;
  for (const b of basis) for (let i = 0; i < n3; i++) for (let j = 0; j < n3; j++) P[i * n3 + j] -= b[i] * b[j];
  const PF = new Float64Array(n3 * n3);
  for (let i = 0; i < n3; i++) for (let k = 0; k < n3; k++) {
    const p = P[i * n3 + k];
    if (p === 0) continue;
    for (let j = 0; j < n3; j++) PF[i * n3 + j] += p * F[k * n3 + j];
  }
  const Fp = new Float64Array(n3 * n3);
  for (let i = 0; i < n3; i++) for (let k = 0; k < n3; k++) {
    const v = PF[i * n3 + k];
    if (v === 0) continue;
    for (let j = 0; j < n3; j++) Fp[i * n3 + j] += v * P[k * n3 + j];
  }
  const { values, vectors } = eigh(Fp, n3);
  // scarta i 5/6 autovalori più vicini a zero (moti esterni)
  const order = [...values.keys()].sort((a, b) => Math.abs(values[a]) - Math.abs(values[b]));
  const external = new Set(order.slice(0, nExt));
  const modes = [];
  for (let k = 0; k < n3; k++) {
    if (external.has(k)) continue;
    const lam = values[k];
    const freq = Math.sign(lam) * Math.sqrt(Math.abs(lam)) * AU_TO_CM1;
    // spostamenti cartesiani (non pesati) e intensità IR
    const disp = new Float64Array(n3);
    let redMass = 0;
    for (let i = 0; i < n3; i++) { disp[i] = vectors[i * n3 + k] / sm[i]; redMass += disp[i] * disp[i]; }
    redMass = 1 / redMass;
    const dq = [0, 0, 0];
    if (dMu) for (let i = 0; i < n3; i++) for (let c = 0; c < 3; c++) dq[c] += dMu[i][c] * vectors[i * n3 + k] / sm[i];
    const ir = dMu ? IR_KM_MOL * (dq[0] ** 2 + dq[1] ** 2 + dq[2] ** 2) : null;
    // normalizza gli spostamenti per l'animazione (massimo spostamento = 1)
    const maxd = Math.max(...[...Array(n3 / 3).keys()].map(a => Math.hypot(disp[3 * a], disp[3 * a + 1], disp[3 * a + 2])));
    modes.push({ freq, reducedMass: redMass, ir, displacement: disp.map(v => v / maxd) });
  }
  modes.sort((a, b) => a.freq - b.freq);
  return { modes, linear: nExt === 5 };
}

/**
 * Funzioni termodinamiche a temperatura T e pressione p (gas ideale).
 * energyHartree: energia elettronica; restituisce contributi in kJ/mol e J/(mol K).
 */
export function thermochemistry(atoms, modes, energyHartree, { T = 298.15, p = 101325, multiplicity = 1, symmetry } = {}) {
  const masses = atoms.map(a => ISOTOPE_MASS[a.Z] ?? 2 * a.Z);
  const M = masses.reduce((s, m) => s + m, 0) * AMU;
  const kT = K_B * T;
  const sym = symmetry ?? analyzeSymmetry(atoms, masses);
  // traslazione (Sackur–Tetrode)
  const qt = Math.pow(2 * Math.PI * M * kT / (H_PLANCK * H_PLANCK), 1.5) * kT / p;
  const St = R_GAS * (Math.log(qt) + 2.5);
  const Et = 1.5 * R_GAS * T;
  // rotazione (rotore rigido classico)
  let Sr = 0, Er = 0, rotConst = [];
  if (atoms.length > 1) {
    const { moments } = inertia(atoms, masses);
    const I = moments.map(m => m * AMU * BOHR_M * BOHR_M); // kg m²
    const theta = I.map(i => (i > 1e-50 ? H_PLANCK * H_PLANCK / (8 * Math.PI * Math.PI * i * K_B) : Infinity));
    rotConst = I.map(i => (i > 1e-50 ? H_PLANCK / (8 * Math.PI * Math.PI * i * C_LIGHT) : null)); // cm⁻¹
    if (sym.linear) {
      const qr = T / (sym.sigma * theta[2]);
      Sr = R_GAS * (Math.log(qr) + 1);
      Er = R_GAS * T;
    } else {
      const qr = Math.sqrt(Math.PI) / sym.sigma * Math.sqrt(T ** 3 / (theta[0] * theta[1] * theta[2]));
      Sr = R_GAS * (Math.log(qr) + 1.5);
      Er = 1.5 * R_GAS * T;
    }
  }
  // vibrazione (oscillatori armonici quantistici)
  let zpe = 0, Ev = 0, Sv = 0, Cv = 0;
  for (const m of modes) {
    if (m.freq <= 0) continue;
    const th = H_PLANCK * C_LIGHT * m.freq / K_B;
    const x = th / T;
    zpe += 0.5 * R_GAS * th;
    Ev += R_GAS * th * (0.5 + 1 / Math.expm1(x));
    Sv += R_GAS * (x / Math.expm1(x) - Math.log(1 - Math.exp(-x)));
    Cv += R_GAS * x * x * Math.exp(x) / Math.expm1(x) ** 2;
  }
  const Se = R_GAS * Math.log(multiplicity);
  const Eel = energyHartree * HARTREE_J * N_A; // J/mol
  const U = Eel + Et + Er + Ev;
  const H = U + R_GAS * T;
  const S = St + Sr + Sv + Se;
  const G = H - T * S;
  const kJ = (v) => v / 1000;
  return {
    T, p, sigma: sym.sigma, pointGroup: sym.pointGroup, rotationalConstants: rotConst,
    zpe: kJ(zpe), Etherm: kJ(Et + Er + Ev), H: kJ(H), G: kJ(G), S, Cv: Cv + 1.5 * R_GAS + (sym.linear ? R_GAS : 1.5 * R_GAS),
    parts: { St, Sr, Sv, Se, Et: kJ(Et), Er: kJ(Er), Ev: kJ(Ev) },
    // correzioni rispetto all'energia elettronica (hartree), come nei programmi di chimica quantistica
    Hcorr: (H - Eel) / (HARTREE_J * N_A),
    Gcorr: (G - Eel) / (HARTREE_J * N_A),
  };
}
