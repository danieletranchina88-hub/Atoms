// MINDO/3: metodo quantistico semiempirico SCF (Bingham, Dewar, Lo, J. Am. Chem. Soc. 97, 1285, 1975).
//
// Gli elettroni di valenza sono descritti da orbitali di Slater (s, p) espansi in 6 gaussiane (Stewart,
// J. Chem. Phys. 52, 431, 1970); l'hamiltoniana è quella INDO con gli integrali a uno e due centri e i parametri
// di legame β_AB e α_AB pubblicati. Le equazioni di Roothaan–Hall non ristrette (UHF) si risolvono a ogni passo
// della dinamica: le forze sui nuclei sono il gradiente analitico dell'energia.
//
// Parametri e forma delle equazioni: implementazione di riferimento di PySCF (pyscf-semiempirical 0.1.1,
// derivata da PyQuante 1.6), con cui il codice è verificato (tests/mindo3.mjs).
//
// Unità interne: Å ed eV (E2 = e²/(4πε₀) = 14,399 eV·Å, come in MOPAC).

import { eigh, symFunction } from './linalg.js';

const E2 = 14.399;
const EV2KCAL = 23.061;
const KB_EV = 8.617333262e-5;
const BOHR = 0.52917721090;

// --- parametri per numero atomico (eV) ----------------------------------------------------------
const P = {
  //      Uss      Upp      Gss     Gpp     Gsp     Gp2     Hsp    Hp2    F03     Vs       Vp       zs        zp        Eheat   Eisol
  1: [-12.505, 0, 12.848, 0, 0, 0, 0, 0, 12.848, -13.605, 0, 1.30, 0, 52.102, -12.505],
  5: [-33.61, -25.11, 10.59, 8.86, 9.56, 7.86, 1.81, 0.50, 8.958, -15.160, -8.520, 1.211156, 0.972826, 135.7, -61.70],
  6: [-51.79, -39.18, 12.23, 11.08, 11.47, 9.84, 2.43, 0.62, 10.833, -21.340, -11.540, 1.739391, 1.709645, 170.89, -119.47],
  7: [-66.06, -56.40, 13.59, 12.98, 12.66, 11.59, 3.14, 0.70, 12.377, -27.510, -14.340, 2.704546, 1.870839, 113.0, -187.51],
  8: [-91.73, -78.80, 15.42, 14.52, 14.48, 12.98, 3.94, 0.77, 13.985, -35.300, -17.910, 3.640575, 2.168448, 59.559, -307.07],
  9: [-129.86, -105.93, 16.92, 16.71, 17.25, 14.91, 4.83, 0.90, 16.250, -43.700, -20.890, 3.111270, 1.419860, 18.86, -475.00],
  15: [-56.23, -42.31, 11.56, 8.64, 10.08, 7.68, 1.92, 0.48, 9.00, -21.100, -10.290, 1.926108, 1.590665, 79.8, -150.81],
  16: [-73.39, -57.25, 12.88, 9.90, 11.26, 8.83, 2.26, 0.54, 10.20, -23.840, -12.410, 1.719480, 1.403205, 65.65, -229.15],
  17: [-98.99, -76.43, 15.03, 11.30, 13.16, 9.97, 2.42, 0.67, 11.73, -25.260, -15.090, 3.430887, 1.627017, 28.95, -345.93],
};
// carica del "core" (elettroni di valenza)
const CORE = { 1: 1, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 15: 5, 16: 6, 17: 7 };

// parametri di coppia β_AB (risonanza) e α_AB (repulsione dei core), Bingham–Dewar–Lo
const PAIR = {};
const pair = (a, b, beta, alpha) => { PAIR[a < b ? `${a}-${b}` : `${b}-${a}`] = [beta, alpha]; };
pair(1, 1, 0.244770, 1.489450);
pair(1, 5, 0.185347, 2.090352); pair(5, 5, 0.151324, 2.280544);
pair(1, 6, 0.315011, 1.475836); pair(5, 6, 0.250031, 2.138291); pair(6, 6, 0.419907, 1.371208);
pair(1, 7, 0.360776, 0.589380); pair(5, 7, 0.310959, 1.909763); pair(6, 7, 0.410886, 1.635259); pair(7, 7, 0.377342, 2.029618);
pair(1, 8, 0.417759, 0.478901); pair(5, 8, 0.349745, 2.484827); pair(6, 8, 0.464514, 1.820975); pair(7, 8, 0.458110, 1.873859); pair(8, 8, 0.659407, 1.537190);
pair(1, 9, 0.195242, 3.771362); pair(5, 9, 0.219591, 2.862183); pair(6, 9, 0.247494, 2.725913); pair(7, 9, 0.205347, 2.861667); pair(8, 9, 0.334044, 2.266949); pair(9, 9, 0.197464, 3.864997);
pair(1, 15, 0.320118, 0.923170); pair(6, 15, 0.457816, 1.029693); pair(8, 15, 0.470000, 1.662500); pair(9, 15, 0.300000, 1.750000); pair(15, 15, 0.311790, 1.186652);
pair(1, 16, 0.220654, 1.700698); pair(6, 16, 0.284620, 1.761370); pair(7, 16, 0.313170, 1.878176); pair(8, 16, 0.422890, 2.077240); pair(16, 16, 0.202489, 1.751617);
pair(1, 17, 0.231653, 2.089404); pair(6, 17, 0.315480, 1.676222); pair(7, 17, 0.302298, 1.817064); pair(15, 17, 0.277322, 1.543720); pair(16, 17, 0.221764, 1.950318); pair(17, 17, 0.258969, 1.792125);

export const MINDO3_ELEMENTS = Object.keys(P).map(Number);

/** Il sistema è trattabile con MINDO/3 se ogni elemento e ogni coppia di elementi presenti ha i parametri. */
export function mindo3Supports(Z) {
  const els = [...new Set(Z)];
  const missingEl = els.filter(z => !P[z]);
  if (missingEl.length) return { ok: false, missing: missingEl.map(String), reason: 'elementi senza parametri MINDO/3' };
  const missing = [];
  for (let a = 0; a < els.length; a++) for (let b = a; b < els.length; b++) {
    const k = els[a] < els[b] ? `${els[a]}-${els[b]}` : `${els[b]}-${els[a]}`;
    if (!PAIR[k]) missing.push(k);
  }
  return missing.length ? { ok: false, missing, reason: 'coppie di elementi senza parametri MINDO/3' } : { ok: true };
}

// --- base STO-6G (Stewart 1970): esponenti e coefficienti per ζ = 1 -------------------------------
const STO6G = {
  '1,0': [[2.310303149e01, 4.235915534e00, 1.185056519e00, 4.070988982e-01, 1.580884151e-01, 6.510953954e-02],
    [9.163596280e-03, 4.936149294e-02, 1.685383049e-01, 3.705627997e-01, 4.164915298e-01, 1.303340841e-01]],
  '2,0': [[2.768496241e01, 5.077140627e00, 1.426786050e00, 2.040335729e-01, 9.260298399e-02, 4.416183978e-02],
    [-4.151277819e-03, -2.067024148e-02, -5.150303337e-02, 3.346271174e-01, 5.621061301e-01, 1.712994697e-01]],
  '2,1': [[5.868285913e00, 1.530329631e00, 5.475665231e-01, 2.288932733e-01, 1.046655969e-01, 4.948220127e-02],
    [7.924233646e-03, 5.144104825e-02, 1.898400060e-01, 4.049863191e-01, 4.012362861e-01, 1.051855189e-01]],
  '3,0': [[3.273031938e00, 9.200611311e-01, 3.593349765e-01, 8.636686991e-02, 4.797373812e-02, 2.724741144e-02],
    [-6.775596947e-03, -5.639325779e-02, -1.587856086e-01, 5.534527651e-01, 5.015351020e-01, 7.223633674e-02]],
  '3,1': [[5.077973607e00, 1.340786940e00, 2.248434849e-01, 1.131741848e-01, 6.076408893e-02, 3.315424265e-02],
    [-3.329929840e-03, -1.419488340e-02, 1.639395770e-01, 4.485358256e-01, 3.908813050e-01, 7.411456232e-02]],
};
const principalN = (z) => (z < 3 ? 1 : z < 10 ? 2 : 3);

// sovrapposizione 1D fra primitive (Obara–Saika), fino a l = 2 su A e l = 1 su B
function os1d(XPA, XPB, p, out) {
  // out[i*2 + j], i = 0..2, j = 0..1
  const s00 = Math.sqrt(Math.PI / p), o = 1 / (2 * p);
  const s10 = XPA * s00, s01 = XPB * s00;
  const s11 = XPB * s10 + o * s00;
  const s20 = XPA * s10 + o * s00;
  const s21 = XPB * s20 + o * 2 * s10;
  out[0] = s00; out[1] = s01; out[2] = s10; out[3] = s11; out[4] = s20; out[5] = s21;
}

/** Funzioni di base di un elemento: shell s e (se Z > 2) p, esponenti in bohr⁻², coefficienti normalizzati. */
const shellCache = new Map();
function shells(z) {
  if (shellCache.has(z)) return shellCache.get(z);
  const n = principalN(z);
  const make = (l, zeta) => {
    const [ex, co] = STO6G[`${n},${l}`];
    const e = ex.map(a => a * zeta * zeta);
    // normalizzazione delle primitive e poi della contrazione (come PySCF)
    let c = co.map((cc, k) => cc * Math.pow(2 * e[k] / Math.PI, 0.75) * (l ? 2 * Math.sqrt(e[k]) : 1));
    let s = 0;
    for (let a = 0; a < 6; a++) for (let b = 0; b < 6; b++) {
      const p = e[a] + e[b];
      s += c[a] * c[b] * Math.pow(Math.PI / p, 1.5) * (l ? 1 / (2 * p) : 1);
    }
    c = c.map(x => x / Math.sqrt(s));
    return { l, e, c };
  };
  const out = [make(0, P[z][11])];
  if (z > 2) out.push(make(1, P[z][12]));
  shellCache.set(z, out);
  return out;
}

/**
 * Sovrapposizione fra le funzioni di valenza di due atomi (A in a, B in b, in Å) e sua derivata rispetto
 * alla posizione di A (in Å⁻¹). Ordine delle funzioni: s, px, py, pz.
 */
function pairOverlap(za, zb, a, b) {
  const sa = shells(za), sb = shells(zb);
  const na = za > 2 ? 4 : 1, nb = zb > 2 ? 4 : 1;
  const S = new Float64Array(na * nb), dS = new Float64Array(3 * na * nb);
  const A = [a[0] / BOHR, a[1] / BOHR, a[2] / BOHR], B = [b[0] / BOHR, b[1] / BOHR, b[2] / BOHR];
  const R2 = (A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2;
  const ox = new Float64Array(6), oy = new Float64Array(6), oz = new Float64Array(6);
  // componenti cartesiane: indice della funzione → [lx, ly, lz]
  const comps = (l, off) => (l === 0 ? [[off, 0, 0, 0]] : [[off, 1, 0, 0], [off + 1, 0, 1, 0], [off + 2, 0, 0, 1]]);
  for (const shA of sa) for (const shB of sb) {
    const ca = comps(shA.l, shA.l ? 1 : 0), cb = comps(shB.l, shB.l ? 1 : 0);
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
      const ea = shA.e[i], eb = shB.e[j], p = ea + eb;
      const mu = ea * eb / p;
      const K = Math.exp(-mu * R2);
      if (K < 1e-18) continue;
      const Px = (ea * A[0] + eb * B[0]) / p, Py = (ea * A[1] + eb * B[1]) / p, Pz = (ea * A[2] + eb * B[2]) / p;
      os1d(Px - A[0], Px - B[0], p, ox);
      os1d(Py - A[1], Py - B[1], p, oy);
      os1d(Pz - A[2], Pz - B[2], p, oz);
      const cc = shA.c[i] * shB.c[j] * K;
      for (const [ia, lx, ly, lz] of ca) {
        for (const [ib, mx, my, mz] of cb) {
          const sx = ox[lx * 2 + mx], sy = oy[ly * 2 + my], sz = oz[lz * 2 + mz];
          S[ia * nb + ib] += cc * sx * sy * sz;
          // ∂/∂A_x di una gaussiana cartesiana: 2a·g(l+1) − l·g(l−1)
          const dx = 2 * ea * ox[(lx + 1) * 2 + mx] - (lx ? lx * ox[(lx - 1) * 2 + mx] : 0);
          const dy = 2 * ea * oy[(ly + 1) * 2 + my] - (ly ? ly * oy[(ly - 1) * 2 + my] : 0);
          const dz = 2 * ea * oz[(lz + 1) * 2 + mz] - (lz ? lz * oz[(lz - 1) * 2 + mz] : 0);
          const k = ia * nb + ib;
          dS[k] += cc * dx * sy * sz / BOHR;
          dS[na * nb + k] += cc * sx * dy * sz / BOHR;
          dS[2 * na * nb + k] += cc * sx * sy * dz / BOHR;
        }
      }
    }
  }
  return { S, dS, na, nb };
}

/** Valori delle funzioni di base di un atomo in un punto (per disegnare densità e orbitali). r in Å. */
export function basisOnAtom(z, dx, dy, dz, out, off) {
  const X = dx / BOHR, Y = dy / BOHR, Zc = dz / BOHR;
  const r2 = X * X + Y * Y + Zc * Zc;
  const sh = shells(z);
  let s = 0;
  const e0 = sh[0].e, c0 = sh[0].c;
  for (let k = 0; k < 6; k++) s += c0[k] * Math.exp(-e0[k] * r2);
  out[off] = s;
  if (z > 2) {
    const e1 = sh[1].e, c1 = sh[1].c;
    let p = 0;
    for (let k = 0; k < 6; k++) p += c1[k] * Math.exp(-e1[k] * r2);
    out[off + 1] = p * X; out[off + 2] = p * Y; out[off + 3] = p * Zc;
  }
}
export const nBasis = (z) => (z > 2 ? 4 : 1);

// --- calcolo SCF ------------------------------------------------------------------------------------

/**
 * Stato persistente fra un passo e l'altro della dinamica (densità di partenza, DIIS).
 * opts.Tel: temperatura elettronica (K) per l'occupazione di Fermi; 0 = aufbau.
 * opts.charge: carica totale; opts.spinPolarize: rompe la simmetria di spin nel primo tentativo.
 */
export class Mindo3 {
  constructor(opts = {}) {
    this.Tel = opts.Tel ?? 300;
    this.charge = opts.charge ?? 0;
    this.maxIter = opts.maxIter ?? 200;
    this.conv = opts.conv ?? 1e-7;
    this.nalpha = opts.nalpha ?? null; // se dati, occupazioni fisse per spin (spin vincolato)
    this.nbeta = opts.nbeta ?? null;
    this.spinGuess = opts.spinGuess ?? null; // polarizzazione di spin iniziale per atomo (−0,5 … 0,5)
    this.field = opts.field ?? null;          // campo elettrico esterno uniforme [Ex, Ey, Ez] in V/Å
    this.key = '';
    this.Pa = null; this.Pb = null;
    this.seed = 12345;
    this.last = null;
  }

  reset() { this.key = ''; this.Pa = null; this.Pb = null; this.Pa1 = null; this.Pb1 = null; this.good = null; }

  rand() { this.seed = (this.seed * 1103515245 + 12345) % 2147483648; return this.seed / 2147483648; }

  /**
   * Gruppi di atomi connessi da sovrapposizioni (entro RCUT) e loro numero di elettroni, oppure null se c'è un solo
   * gruppo o la molteplicità è fissata globalmente. Elettroni per gruppo: popolazioni del passo precedente
   * arrotondate; all'inizio Σ(Z_core − carica formale) se le cariche formali sono note (this.formal).
   */
  coupledBlocks(Z, N, first, zc, nel, overlaps, fresh) {
    if (this.nalpha !== null && this.nbeta !== null) return null;
    const parent = Int32Array.from({ length: N }, (_, i) => i);
    const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    for (const o of overlaps) parent[find(o.A)] = find(o.B);
    const groups = new Map();
    for (let A = 0; A < N; A++) { const r = find(A); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(A); }
    if (groups.size < 2) return null;
    const n = first[N];
    let pop = null;
    if (!fresh && this.Pa && this.Pa.length === n * n) {
      pop = new Float64Array(N);
      for (let A = 0; A < N; A++) for (let m = first[A]; m < first[A + 1]; m++) pop[A] += this.Pa[m * n + m] + this.Pb[m * n + m];
    } else if (this.formal && this.formal.length === N) {
      pop = Float64Array.from(zc, (z, A) => z - this.formal[A]);
    } else return null;   // primo calcolo senza cariche formali: livello di Fermi comune, poi gruppi dal passo dopo
    const list = [...groups.values()].map(atoms => {
      const exact = atoms.reduce((s, A) => s + pop[A], 0);
      const idx = [];
      for (const A of atoms) for (let m = first[A]; m < first[A + 1]; m++) idx.push(m);
      return { atoms, idx, exact, nel: Math.round(exact) };
    });
    // l'arrotondamento deve conservare il numero totale di elettroni
    let diff = nel - list.reduce((s, b) => s + b.nel, 0);
    const byResidual = [...list].sort((a, b) => (b.exact - b.nel) - (a.exact - a.nel));
    for (let k = 0; diff !== 0 && k < 4 * list.length; k++) {
      const b = diff > 0 ? byResidual[k % list.length] : byResidual[list.length - 1 - (k % list.length)];
      if (diff > 0 && b.nel < 2 * b.idx.length) { b.nel++; diff--; } else if (diff < 0 && b.nel > 0) { b.nel--; diff++; }
    }
    if (diff !== 0 || list.some(b => b.nel < 0 || b.nel > 2 * b.idx.length)) return null;
    return list;
  }

  /**
   * Energia (eV), forze (eV/Å, in F) e analisi della funzione d'onda.
   * @param Z numeri atomici; pos Float64Array 3N in Å; F Float64Array 3N (uscita, può essere null)
   */
  compute(Z, pos, F = null) {
    const N = Z.length;
    const sup = mindo3Supports(Z);
    if (!sup.ok) throw new Error(`MINDO/3: ${sup.reason} (${sup.missing.join(', ')})`);
    // indici delle funzioni di base
    const first = new Int32Array(N + 1);
    for (let A = 0; A < N; A++) first[A + 1] = first[A] + nBasis(Z[A]);
    const n = first[N];
    const zc = Z.map(z => CORE[z]);
    const nel = zc.reduce((a, b) => a + b, 0) - this.charge;
    const key = `${Array.from(Z).join(',')}|${this.charge}`;
    const fresh = key !== this.key || !this.Pa;
    this.key = key;

    // --- integrali a due centri: γ_AB, sovrapposizioni, hamiltoniana di core -------------------
    const rho = Z.map(z => E2 / P[z][8]);
    const gam = new Float64Array(N * N), dgam = new Float64Array(N * N), Rm = new Float64Array(N * N);
    for (let A = 0; A < N; A++) for (let B = A + 1; B < N; B++) {
      const dx = pos[3 * A] - pos[3 * B], dy = pos[3 * A + 1] - pos[3 * B + 1], dz = pos[3 * A + 2] - pos[3 * B + 2];
      const R2 = dx * dx + dy * dy + dz * dz, R = Math.sqrt(R2);
      const c2 = 0.25 * (rho[A] + rho[B]) ** 2;
      const g = E2 / Math.sqrt(R2 + c2);
      gam[A * N + B] = gam[B * N + A] = g;
      dgam[A * N + B] = dgam[B * N + A] = -E2 * R / Math.pow(R2 + c2, 1.5);
      Rm[A * N + B] = Rm[B * N + A] = R;
    }
    const H = new Float64Array(n * n);
    const Ival = new Float64Array(n);
    for (let A = 0; A < N; A++) {
      const p = P[Z[A]];
      Ival[first[A]] = p[9];
      for (let k = 1; k < nBasis(Z[A]); k++) Ival[first[A] + k] = p[10];
      let vnuc = 0;
      for (let B = 0; B < N; B++) if (B !== A) vnuc += zc[B] * gam[A * N + B];
      H[first[A] * n + first[A]] = p[0] - vnuc;
      for (let k = 1; k < nBasis(Z[A]); k++) { const m = first[A] + k; H[m * n + m] = p[1] - vnuc; }
    }
    const overlaps = []; // coppie vicine: { A, B, S, dS, beta }
    const RCUT = 7; // Å: oltre, β·S è sotto 10⁻³ eV anche per H–H
    for (let A = 0; A < N; A++) for (let B = A + 1; B < N; B++) {
      if (Rm[A * N + B] > RCUT) continue;
      const [beta] = PAIR[Z[A] < Z[B] ? `${Z[A]}-${Z[B]}` : `${Z[B]}-${Z[A]}`];
      const o = pairOverlap(Z[A], Z[B], [pos[3 * A], pos[3 * A + 1], pos[3 * A + 2]], [pos[3 * B], pos[3 * B + 1], pos[3 * B + 2]]);
      overlaps.push({ A, B, ...o, beta });
      for (let i = 0; i < o.na; i++) for (let j = 0; j < o.nb; j++) {
        const mu = first[A] + i, nu = first[B] + j;
        const h = o.S[i * o.nb + j] * (Ival[mu] + Ival[nu]) * beta;
        H[mu * n + nu] = H[nu * n + mu] = h;
      }
    }

    // --- blocchi di atomi accoppiati -------------------------------------------------------------
    // La Fock fra atomi senza sovrapposizione (oltre RCUT) è nulla: un elettrone non può passare da un gruppo
    // di molecole a un altro lontano. Ogni gruppo connesso conserva il suo numero di elettroni (preso dal passo
    // precedente o, all'inizio, dalle cariche formali delle specie inserite); il livello di Fermi è per gruppo.
    const blocks = this.coupledBlocks(Z, N, first, zc, nel, overlaps, fresh);

    // --- campo elettrico esterno uniforme -------------------------------------------------------
    // energia di un elettrone (carica −e) nel potenziale φ = −E·r: +E·r. In V/Å e Å il risultato è in eV.
    // Elementi di matrice nella base di valenza: ⟨μ|r|μ⟩ = R_A e il dipolo atomico ⟨s|r_c|p_c⟩ = D_A.
    const field = this.field && this.field.some(v => v) ? this.field : null;
    if (field) {
      for (let A = 0; A < N; A++) {
        const dot = field[0] * pos[3 * A] + field[1] * pos[3 * A + 1] + field[2] * pos[3 * A + 2];
        for (let m = first[A]; m < first[A + 1]; m++) H[m * n + m] += dot;
        if (nBasis(Z[A]) === 4) {
          const D = spDipole(Z[A]), sIdx = first[A];
          for (let c = 0; c < 3; c++) { const k = sIdx * n + sIdx + 1 + c; H[k] += field[c] * D; H[(sIdx + 1 + c) * n + sIdx] += field[c] * D; }
        }
      }
    }

    // --- densità di partenza ------------------------------------------------------------------
    let Pa, Pb;
    if (fresh) {
      Pa = new Float64Array(n * n); Pb = new Float64Array(n * n);
      // occupazione media degli orbitali di ogni atomo, con una piccola polarizzazione di spin casuale:
      // permette alla SCF di trovare stati aperti (O₂ tripletto, radicali, legami rotti in modo omolitico)
      for (let A = 0; A < N; A++) {
        const nb = nBasis(Z[A]);
        const pol = this.spinGuess ? (this.spinGuess[A] ?? 0) : (this.rand() - 0.5) * 0.4;
        for (let k = 0; k < nb; k++) {
          const m = first[A] + k;
          Pa[m * n + m] = zc[A] / nb * (0.5 + pol);
          Pb[m * n + m] = zc[A] / nb * (0.5 - pol);
        }
      }
    } else if (this.Pa1 && this.Pa1.length === this.Pa.length) {
      // estrapolazione lineare dai due passi precedenti: la SCF parte già vicina alla soluzione
      Pa = new Float64Array(n * n); Pb = new Float64Array(n * n);
      for (let k = 0; k < n * n; k++) { Pa[k] = 2 * this.Pa[k] - this.Pa1[k]; Pb[k] = 2 * this.Pb[k] - this.Pb1[k]; }
    } else { Pa = Float64Array.from(this.Pa); Pb = Float64Array.from(this.Pb); }

    // --- ciclo SCF con DIIS ---------------------------------------------------------------------
    // Se non converge si riprova dall'ultima densità convergente (senza estrapolazione) con uno smorzamento
    // forte: forze da una SCF non convergente non sono il gradiente di nessuna energia.
    const oneC = Z.map(z => oneCenter(z));
    const Fa = new Float64Array(n * n), Fb = new Float64Array(n * n);
    let E = 0, converged = false, it = 0, res = null;
    const kT = KB_EV * this.Tel;
    for (let attempt = 0; attempt < 2 && !converged; attempt++) {
      const damp = attempt === 0 ? (fresh ? 4 : 0) : 40;
      if (attempt === 1) {
        if (this.good && this.good.key === key && this.good.Pa.length === n * n) { Pa = Float64Array.from(this.good.Pa); Pb = Float64Array.from(this.good.Pb); }
        res = null;
      }
      const diis = [];
      let Eold = Infinity;
      const maxIt = attempt === 0 ? this.maxIter : 2 * this.maxIter;
      for (it = 0; it < maxIt; it++) {
        buildFock(Z, first, N, n, H, gam, oneC, Pa, Pb, Fa, Fb);
        // energia elettronica: ½ Σ_σ Pσ (H + Fσ)
        let Eel = 0;
        for (let k = 0; k < n * n; k++) Eel += 0.5 * (Pa[k] * (H[k] + Fa[k]) + Pb[k] * (H[k] + Fb[k]));
        // errore DIIS: FP − PF per ogni spin (a guscio chiuso i due spin coincidono)
        let closed = true;
        for (let k = 0; k < n * n; k++) if (Math.abs(Pa[k] - Pb[k]) > 1e-9) { closed = false; break; }
        const err = new Float64Array(2 * n * n);
        let emax = 0;
        commutator(Fa, Pa, n, err, 0);
        if (closed) err.copyWithin(n * n, 0, n * n); else commutator(Fb, Pb, n, err, n * n);
        for (let k = 0; k < err.length; k++) emax = Math.max(emax, Math.abs(err[k]));
        E = Eel - (res ? kT * res.entropy : 0);
        if (it > 0 && Math.abs(E - Eold) < this.conv && emax < Math.sqrt(this.conv) * 10) { converged = true; break; }
        Eold = E;
        const Fcat = new Float64Array(2 * n * n);
        Fcat.set(Fa); Fcat.set(Fb, n * n);
        diis.push({ F: Fcat, e: err });
        if (diis.length > 8) diis.shift();
        let Fx = null;
        if (diis.length >= 2 && it >= damp / 2) Fx = diisExtrapolate(diis);
        const Fua = Fx ? Fx.subarray(0, n * n) : Fa, Fub = Fx ? Fx.subarray(n * n) : Fb;
        res = blocks ? occupyBlocks(Fua, Fub, n, blocks, kT, closed) : occupy(Fua, Fub, n, nel, kT, this.nalpha, this.nbeta, closed);
        // miscelazione con la densità precedente nelle prime iterazioni (evita oscillazioni)
        const mix = it < damp ? 0.5 : 0;
        for (let k = 0; k < n * n; k++) {
          Pa[k] = mix * Pa[k] + (1 - mix) * res.Pa[k];
          Pb[k] = mix * Pb[k] + (1 - mix) * res.Pb[k];
        }
      }
    }
    if (converged) this.good = { key, Pa: Float64Array.from(Pa), Pb: Float64Array.from(Pb) };
    if (!res) res = blocks ? occupyBlocks(Fa, Fb, n, blocks, kT) : occupy(Fa, Fb, n, nel, kT, this.nalpha, this.nbeta);
    this.Pa1 = fresh ? null : this.Pa; this.Pb1 = fresh ? null : this.Pb;
    this.Pa = Float64Array.from(Pa); this.Pb = Float64Array.from(Pb);

    // --- energia totale ---------------------------------------------------------------------------
    let Enuc = 0;
    for (let A = 0; A < N; A++) for (let B = A + 1; B < N; B++) {
      const R = Rm[A * N + B], g = gam[A * N + B];
      Enuc += zc[A] * zc[B] * (g + (E2 / R - g) * coreScale(Z[A], Z[B], R).f);
    }
    // energia dei core nel campo: carica +Z_core nel potenziale φ = −E·R
    let Efield = 0;
    if (field) for (let A = 0; A < N; A++) Efield -= zc[A] * (field[0] * pos[3 * A] + field[1] * pos[3 * A + 1] + field[2] * pos[3 * A + 2]);
    const Etot = E + Enuc + Efield;
    let Hf = Etot * EV2KCAL;
    for (const z of Z) Hf += P[z][13] - P[z][14] * EV2KCAL;

    // --- gradiente analitico a densità fissa (la SCF è variazionale e la base è ortonormale) --------
    if (F) {
      F.fill(0);
      const Pt = new Float64Array(n * n);
      for (let k = 0; k < n * n; k++) Pt[k] = Pa[k] + Pb[k];
      const popA = new Float64Array(N);
      for (let A = 0; A < N; A++) for (let m = first[A]; m < first[A + 1]; m++) popA[A] += Pt[m * n + m];
      for (let A = 0; A < N; A++) for (let B = A + 1; B < N; B++) {
        const R = Rm[A * N + B], g = gam[A * N + B], dg = dgam[A * N + B];
        let x2 = 0;
        for (let m = first[A]; m < first[A + 1]; m++) for (let q = first[B]; q < first[B + 1]; q++) x2 += Pa[m * n + q] ** 2 + Pb[m * n + q] ** 2;
        const { f, df } = coreScale(Z[A], Z[B], R);
        const zz = zc[A] * zc[B];
        let dEdR = dg * (popA[A] * popA[B] - x2 - zc[B] * popA[A] - zc[A] * popA[B]);
        dEdR += zz * (dg * (1 - f) - g * df - E2 / (R * R) * f + E2 / R * df);
        for (let c = 0; c < 3; c++) {
          const u = (pos[3 * A + c] - pos[3 * B + c]) / R;
          F[3 * A + c] -= dEdR * u;
          F[3 * B + c] += dEdR * u;
        }
      }
      // campo uniforme: forza q_A·E su ogni atomo (Z_core meno la popolazione elettronica)
      if (field) for (let A = 0; A < N; A++) for (let c = 0; c < 3; c++) F[3 * A + c] += (zc[A] - popA[A]) * field[c];
      for (const o of overlaps) {
        const { A, B, na, nb, dS, beta } = o;
        const gx = [0, 0, 0];
        for (let i = 0; i < na; i++) for (let j = 0; j < nb; j++) {
          const mu = first[A] + i, nu = first[B] + j;
          const w = 2 * Pt[mu * n + nu] * (Ival[mu] + Ival[nu]) * beta;
          for (let c = 0; c < 3; c++) gx[c] += w * dS[c * na * nb + i * nb + j];
        }
        for (let c = 0; c < 3; c++) { F[3 * A + c] -= gx[c]; F[3 * B + c] += gx[c]; }
      }
    }

    // --- analisi: cariche, ordini di legame, spin ------------------------------------------------
    const q = new Float64Array(N), spin = new Float64Array(N);
    for (let A = 0; A < N; A++) {
      let pa = 0, pb = 0;
      for (let m = first[A]; m < first[A + 1]; m++) { pa += Pa[m * n + m]; pb += Pb[m * n + m]; }
      q[A] = zc[A] - pa - pb;
      spin[A] = pa - pb;
    }
    const bonds = [];
    for (let A = 0; A < N; A++) for (let B = A + 1; B < N; B++) {
      if (Rm[A * N + B] > 4) continue;
      let s = 0;
      for (let m = first[A]; m < first[A + 1]; m++) for (let q2 = first[B]; q2 < first[B + 1]; q2++) s += 2 * (Pa[m * n + q2] ** 2 + Pb[m * n + q2] ** 2);
      if (s > 0.05) bonds.push({ i: A, j: B, n: s, w: 1 });
    }
    // momento di dipolo (e·Å): Σ q_A R_A − 2 Σ D_A P(s, p_c). Per uno ione dipende dall'origine.
    const dipole = [0, 0, 0];
    for (let A = 0; A < N; A++) {
      for (let c = 0; c < 3; c++) dipole[c] += q[A] * pos[3 * A + c];
      if (nBasis(Z[A]) === 4) {
        const D = spDipole(Z[A]), sIdx = first[A];
        for (let c = 0; c < 3; c++) dipole[c] -= 2 * D * (Pa[sIdx * n + sIdx + 1 + c] + Pb[sIdx * n + sIdx + 1 + c]);
      }
    }
    const Sz = 0.5 * (res.na - res.nb);
    this.last = {
      Z: Array.from(Z), first: Array.from(first), n, pos: Float64Array.from(pos), Pa: Float64Array.from(Pa), Pb: Float64Array.from(Pb),
      Ca: res.Ca, Cb: res.Cb, ea: res.ea, eb: res.eb, fa: res.fa, fb: res.fb,
    };
    return {
      E: Etot, Eel: E, Enuc, Efield, Hf, converged, iterations: it + 1, q, spin, bonds, dipole,
      nalpha: res.na, nbeta: res.nb, Sz, gap: res.gap, homo: res.homo, lumo: res.lumo, nbf: n,
    };
  }
}

/**
 * Integrale di dipolo atomico ⟨ns|x|np_x⟩ fra orbitali di Slater con lo stesso n (in Å):
 * D = (2n+1)/√3 · (4ζ_s ζ_p)^(n+½) / (ζ_s+ζ_p)^(2n+2) a₀.
 */
function spDipole(z) {
  const n = z <= 10 ? 2 : 3, zs = P[z][11], zp = P[z][12];
  return (2 * n + 1) / Math.sqrt(3) * Math.pow(4 * zs * zp, n + 0.5) / Math.pow(zs + zp, 2 * n + 2) * 0.52917721090;
}

/** Fattore della repulsione dei core: e^(−αR), oppure α e^(−R) per le coppie N–H e O–H. Con la derivata in R. */
function coreScale(za, zb, R) {
  const [, alpha] = PAIR[za < zb ? `${za}-${zb}` : `${zb}-${za}`];
  const nh = (za === 1 && (zb === 7 || zb === 8)) || (zb === 1 && (za === 7 || za === 8));
  if (nh) { const f = alpha * Math.exp(-R); return { f, df: -f }; }
  const f = Math.exp(-alpha * R);
  return { f, df: -alpha * f };
}

/** Integrali a un centro: matrici J e K (4×4 o 1×1) come in MINDO/3. */
function oneCenter(z) {
  const p = P[z];
  if (z < 3) return { n: 1, J: [p[2]], K: [0] };
  const [, , Gss, Gpp, Gsp, Gp2, Hsp, Hp2] = p;
  const J = new Float64Array(16), K = new Float64Array(16);
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
    if (i === 0 && j === 0) J[0] = Gss;
    else if (i === 0 || j === 0) { J[i * 4 + j] = Gsp; K[i * 4 + j] = Hsp; }
    else if (i === j) J[i * 4 + j] = Gpp;
    else { J[i * 4 + j] = Gp2; K[i * 4 + j] = Hp2; }
  }
  return { n: 4, J, K };
}

/** Matrici di Fock UHF: Fσ = H + J(P) − K(Pσ). */
function buildFock(Z, first, N, n, H, gam, oneC, Pa, Pb, Fa, Fb) {
  Fa.set(H); Fb.set(H);
  const popA = new Float64Array(N);
  for (let A = 0; A < N; A++) for (let m = first[A]; m < first[A + 1]; m++) popA[A] += Pa[m * n + m] + Pb[m * n + m];
  for (let A = 0; A < N; A++) {
    const o = oneC[A], p0 = first[A], nb = o.n;
    for (let i = 0; i < nb; i++) {
      const mi = p0 + i;
      // parte coulombiana e di scambio a un centro
      let jd = 0, kda = 0, kdb = 0;
      for (let j = 0; j < nb; j++) {
        const mj = p0 + j;
        const pt = Pa[mj * n + mj] + Pb[mj * n + mj];
        jd += o.J[i * nb + j] * pt;
        kda += o.K[i * nb + j] * Pa[mj * n + mj];
        kdb += o.K[i * nb + j] * Pb[mj * n + mj];
      }
      Fa[mi * n + mi] += jd - kda; Fb[mi * n + mi] += jd - kdb;
      for (let j = 0; j < nb; j++) {
        const mj = p0 + j;
        const pa = Pa[mi * n + mj], pb = Pb[mi * n + mj];
        const kk = o.K[i * nb + j], jj = o.J[i * nb + j];
        // J[i,j] += 2 (ij|ij) P_ij (i ≠ j); K[i,j] += ((ii|jj) + (ij|ij)) Pσ_ij
        if (i !== j) { Fa[mi * n + mj] += 2 * kk * (pa + pb); Fb[mi * n + mj] += 2 * kk * (pa + pb); }
        Fa[mi * n + mj] -= (jj + kk) * pa; Fb[mi * n + mj] -= (jj + kk) * pb;
      }
    }
    // due centri
    for (let B = 0; B < N; B++) {
      if (B === A) continue;
      const g = gam[A * N + B];
      for (let i = p0; i < first[A + 1]; i++) {
        Fa[i * n + i] += g * popA[B]; Fb[i * n + i] += g * popA[B];
        for (let j = first[B]; j < first[B + 1]; j++) { Fa[i * n + j] -= g * Pa[i * n + j]; Fb[i * n + j] -= g * Pb[i * n + j]; }
      }
    }
  }
}

/** FP − PF = M − Mᵀ con M = FP (F e P simmetriche). */
function commutator(F, D, n, out, off) {
  const M = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    const row = i * n;
    for (let k = 0; k < n; k++) {
      const f = F[row + k];
      if (f === 0) continue;
      const kr = k * n;
      for (let j = 0; j < n; j++) M[row + j] += f * D[kr + j];
    }
  }
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) out[off + i * n + j] = M[i * n + j] - M[j * n + i];
}

function diisExtrapolate(list) {
  const m = list.length;
  const B = Array.from({ length: m + 1 }, () => new Float64Array(m + 1));
  for (let i = 0; i < m; i++) for (let j = 0; j <= i; j++) {
    let s = 0;
    const a = list[i].e, b = list[j].e;
    for (let k = 0; k < a.length; k++) s += a[k] * b[k];
    B[i][j] = B[j][i] = s;
  }
  for (let i = 0; i < m; i++) { B[i][m] = B[m][i] = -1; }
  B[m][m] = 0;
  const rhs = new Float64Array(m + 1); rhs[m] = -1;
  // eliminazione di Gauss con pivot
  const M = B.map((r, i) => [...r, rhs[i]]);
  for (let k = 0; k <= m; k++) {
    let piv = k;
    for (let i = k + 1; i <= m; i++) if (Math.abs(M[i][k]) > Math.abs(M[piv][k])) piv = i;
    if (Math.abs(M[piv][k]) < 1e-14) return null;
    [M[k], M[piv]] = [M[piv], M[k]];
    for (let i = k + 1; i <= m; i++) { const f = M[i][k] / M[k][k]; for (let j = k; j <= m + 1; j++) M[i][j] -= f * M[k][j]; }
  }
  const c = new Float64Array(m + 1);
  for (let i = m; i >= 0; i--) { let s = M[i][m + 1]; for (let j = i + 1; j <= m; j++) s -= M[i][j] * c[j]; c[i] = s / M[i][i]; }
  const F = new Float64Array(list[0].F.length);
  for (let i = 0; i < m; i++) { const Fi = list[i].F; for (let k = 0; k < F.length; k++) F[k] += c[i] * Fi[k]; }
  return F;
}

/**
 * Occupazioni: livello di Fermi comune ai due spin (la magnetizzazione è libera, come nella teoria dello
 * spin non ristretta), distribuzione di Fermi–Dirac a temperatura elettronica kT (0 = aufbau).
 */
function occupy(Fa, Fb, n, nel, kT, fixA, fixB, closed = false) {
  const a = eigh(Fa, n);
  let b = a;
  if (!closed) { for (let k = 0; k < n * n; k++) if (Math.abs(Fa[k] - Fb[k]) > 1e-10) { b = eigh(Fb, n); break; } }
  const levels = [];
  for (let k = 0; k < n; k++) { levels.push([a.values[k], 0, k]); levels.push([b.values[k], 1, k]); }
  const fa = new Float64Array(n), fb = new Float64Array(n);
  let entropy = 0;
  if (fixA !== null && fixB !== null) {
    for (let k = 0; k < Math.min(fixA, n); k++) fa[k] = 1;
    for (let k = 0; k < Math.min(fixB, n); k++) fb[k] = 1;
  } else if (kT <= 0) {
    levels.sort((x, y) => x[0] - y[0]);
    for (let k = 0; k < Math.min(nel, levels.length); k++) (levels[k][1] ? fb : fa)[levels[k][2]] = 1;
  } else {
    const count = (mu) => levels.reduce((s, [e]) => s + fermi(e, mu, kT), 0);
    let lo = Math.min(...levels.map(l => l[0])) - 20, hi = Math.max(...levels.map(l => l[0])) + 20;
    for (let it = 0; it < 64 && hi - lo > 1e-10; it++) { const mid = 0.5 * (lo + hi); if (count(mid) < nel) lo = mid; else hi = mid; }
    const mu = 0.5 * (lo + hi);
    for (const [e, s, k] of levels) {
      const f = fermi(e, mu, kT);
      (s ? fb : fa)[k] = f;
      if (f > 1e-12 && f < 1 - 1e-12) entropy -= f * Math.log(f) + (1 - f) * Math.log(1 - f);
    }
  }
  const dens = (V, f) => {
    // D = C f Cᵀ sulle sole colonne occupate (copiate in righe contigue per l'accesso in memoria)
    const occ = [];
    for (let k = 0; k < n; k++) if (f[k] > 1e-10) occ.push(k);
    const m = occ.length;
    const Ct = new Float64Array(m * n);
    occ.forEach((k, q) => { for (let i = 0; i < n; i++) Ct[q * n + i] = V[i * n + k]; });
    const D = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let q = 0; q < m; q++) {
        const ci = Ct[q * n + i] * f[occ[q]];
        if (ci === 0) continue;
        const cq = q * n, row = i * n;
        for (let j = i; j < n; j++) D[row + j] += ci * Ct[cq + j];
      }
    }
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) D[i * n + j] = D[j * n + i];
    return D;
  };
  const na = fa.reduce((s, x) => s + x, 0), nb = fb.reduce((s, x) => s + x, 0);
  // HOMO e LUMO sull'insieme dei due spin
  const occ = levels.map(([e, s, k]) => [e, (s ? fb : fa)[k]]);
  const homo = Math.max(...occ.filter(x => x[1] > 0.5).map(x => x[0]), -Infinity);
  const lumo = Math.min(...occ.filter(x => x[1] < 0.5).map(x => x[0]), Infinity);
  const Pa = dens(a.vectors, fa);
  let same = b === a;
  if (same) for (let k = 0; k < n; k++) if (fa[k] !== fb[k]) { same = false; break; }
  return {
    Pa, Pb: same ? Pa.slice() : dens(b.vectors, fb), Ca: a.vectors, Cb: b.vectors, ea: a.values, eb: b.values,
    fa, fb, na, nb, entropy, homo, lumo, gap: lumo - homo,
  };
}
/** Occupazione separata per gruppi di orbitali senza accoppiamento: ogni gruppo ha il suo numero di elettroni. */
function occupyBlocks(Fa, Fb, n, blocks, kT, closed = false) {
  const Pa = new Float64Array(n * n), Pb = new Float64Array(n * n);
  const Ca = new Float64Array(n * n), Cb = new Float64Array(n * n);
  const ea = new Float64Array(n), eb = new Float64Array(n), fa = new Float64Array(n), fb = new Float64Array(n);
  let col = 0, entropy = 0, na = 0, nb = 0, homo = -Infinity, lumo = Infinity;
  for (const b of blocks) {
    const m = b.idx.length;
    const sa = new Float64Array(m * m), sb = new Float64Array(m * m);
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) { sa[i * m + j] = Fa[b.idx[i] * n + b.idx[j]]; sb[i * m + j] = Fb[b.idx[i] * n + b.idx[j]]; }
    const r = occupy(sa, sb, m, b.nel, kT, null, null, closed);
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) {
      Pa[b.idx[i] * n + b.idx[j]] = r.Pa[i * m + j];
      Pb[b.idx[i] * n + b.idx[j]] = r.Pb[i * m + j];
    }
    for (let k = 0; k < m; k++) {
      for (let i = 0; i < m; i++) { Ca[b.idx[i] * n + col + k] = r.Ca[i * m + k]; Cb[b.idx[i] * n + col + k] = r.Cb[i * m + k]; }
      ea[col + k] = r.ea[k]; eb[col + k] = r.eb[k]; fa[col + k] = r.fa[k]; fb[col + k] = r.fb[k];
    }
    col += m; entropy += r.entropy; na += r.na; nb += r.nb;
    homo = Math.max(homo, r.homo); lumo = Math.min(lumo, r.lumo);
  }
  return { Pa, Pb, Ca, Cb, ea, eb, fa, fb, na, nb, entropy, homo, lumo, gap: lumo - homo };
}
const fermi = (e, mu, kT) => { const x = (e - mu) / kT; return x > 40 ? 0 : x < -40 ? 1 : 1 / (1 + Math.exp(x)); };

/**
 * Fornitore di forze per la dinamica molecolare della sandbox (stessa interfaccia del campo reattivo):
 * a ogni passo una SCF MINDO/3 completa. Espone l'ultima funzione d'onda per disegnare densità e orbitali.
 */
/** Elettroni di valenza (carica dei core) di MINDO/3 per una lista di numeri atomici. */
export const valenceElectrons = (Z) => Z.reduce((s, z) => s + (CORE[z] ?? 0), 0);

/**
 * Forze MINDO/3 per la dinamica. La carica totale arriva dalla simulazione (ctx.charge, ioni inseriti);
 * multiplicity: null = spin libero (livello di Fermi comune), altrimenti 2S+1 fissato;
 * field: campo elettrico uniforme [Ex, Ey, Ez] in V/Å.
 */
export function makeMindo3Provider({ Tel = 300, conv = 1e-5, spinGuess = null, multiplicity = null, field = null } = {}) {
  // tolleranza SCF di 10⁻⁵ eV: in NVE l'energia totale deriva di circa 0,5 meV per atomo e per ps
  const m = new Mindo3({ Tel, conv, spinGuess, field });
  const provider = {
    kind: 'mindo3',
    info: {},
    multiplicity,
    reset() { m.reset(); },
    setField(f) { m.field = f && f.some(v => v) ? f.slice() : null; },
    get field() { return m.field; },
    compute(Z, pos, F, ctx = {}) {
      const charge = ctx.charge ?? 0;
      if (charge !== m.charge) { m.charge = charge; m.reset(); }
      m.formal = ctx.formal ?? null;
      if (provider.multiplicity) {
        const nel = valenceElectrons(Z) - charge, unp = provider.multiplicity - 1;
        if (nel < 0 || (nel - unp) % 2 !== 0 || unp > nel) throw new Error(`Molteplicità ${provider.multiplicity} impossibile con ${nel} elettroni di valenza: scegli ${nel % 2 ? 'doppietto, quartetto…' : 'singoletto, tripletto…'}`);
        m.nalpha = (nel + unp) / 2; m.nbeta = (nel - unp) / 2;
      } else { m.nalpha = null; m.nbeta = null; }
      const r = m.compute(Z, pos, F);
      provider.info = { method: 'MINDO/3', nbf: r.nbf, Sz: r.Sz, gap: r.gap, homo: r.homo, lumo: r.lumo, iterations: r.iterations, converged: r.converged, Hf: r.Hf, Tel, charge, multiplicity: provider.multiplicity, dipole: r.dipole, field: m.field };
      return { E: r.E, parts: { mindo3: r.E, Hf: r.Hf, field: r.Efield }, q: r.q, bonds: r.bonds, hbonds: [], spin: r.spin };
    },
    /** Funzione d'onda corrente: densità totale e orbitali di frontiera (per il disegno). */
    wavefunction() {
      const w = m.last;
      if (!w) return null;
      const n = w.n, N = w.Z.length;
      // la base di MINDO/3 è intesa ortogonalizzata alla Löwdin (φ' = φ S^(−1/2)): per disegnare densità e
      // orbitali con gli orbitali di Slater veri si trasformano P' = S^(−1/2) P S^(−1/2) e c' = S^(−1/2) c
      const S = new Float64Array(n * n);
      for (let k = 0; k < n; k++) S[k * n + k] = 1;
      for (let A = 0; A < N; A++) for (let B = A + 1; B < N; B++) {
        const a = [w.pos[3 * A], w.pos[3 * A + 1], w.pos[3 * A + 2]], b = [w.pos[3 * B], w.pos[3 * B + 1], w.pos[3 * B + 2]];
        if (Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]) > 9) continue;
        const o = pairOverlap(w.Z[A], w.Z[B], a, b);
        for (let i = 0; i < o.na; i++) for (let j = 0; j < o.nb; j++) {
          const mu = w.first[A] + i, nu = w.first[B] + j;
          S[mu * n + nu] = S[nu * n + mu] = o.S[i * o.nb + j];
        }
      }
      const X = symFunction(S, n, (x) => 1 / Math.sqrt(Math.max(x, 1e-8)));
      const Pt = new Float64Array(n * n);
      for (let k = 0; k < n * n; k++) Pt[k] = w.Pa[k] + w.Pb[k];
      const XP = new Float64Array(n * n), P = new Float64Array(n * n);
      for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const x = X[i * n + k]; if (x) for (let j = 0; j < n; j++) XP[i * n + j] += x * Pt[k * n + j]; }
      for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const x = XP[i * n + k]; if (x) for (let j = 0; j < n; j++) P[i * n + j] += x * X[k * n + j]; }
      // HOMO e LUMO sull'insieme dei due spin (livello di Fermi comune)
      const levels = [];
      for (let k = 0; k < n; k++) { levels.push({ e: w.ea[k], f: w.fa[k], C: w.Ca, k, spin: 'α' }); levels.push({ e: w.eb[k], f: w.fb[k], C: w.Cb, k, spin: 'β' }); }
      const occ = levels.filter(l => l.f > 0.5).sort((a, b) => b.e - a.e);
      const vir = levels.filter(l => l.f <= 0.5).sort((a, b) => a.e - b.e);
      const vec = (l) => {
        if (!l) return null;
        const c = new Float64Array(n);
        for (let i = 0; i < n; i++) { let s2 = 0; for (let k = 0; k < n; k++) s2 += X[i * n + k] * l.C[k * n + l.k]; c[i] = s2; }
        return { e: l.e, spin: l.spin, c };
      };
      const St = new Float64Array(n * n), XS = new Float64Array(n * n), Ps = new Float64Array(n * n);
      for (let k = 0; k < n * n; k++) St[k] = w.Pa[k] - w.Pb[k];
      for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const x = X[i * n + k]; if (x) for (let j = 0; j < n; j++) XS[i * n + j] += x * St[k * n + j]; }
      for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const x = XS[i * n + k]; if (x) for (let j = 0; j < n; j++) Ps[i * n + j] += x * X[k * n + j]; }
      const orbitals = [...occ.slice(0, 6), ...vir.slice(0, 6)].map(l => ({ ...vec(l), id: `${l.spin}:${l.k}`, index: l.k, occ: l.f }));
      return { kind: 'sto', Z: w.Z, pos: w.pos, first: w.first, n, P, Ps, orbitals, homo: vec(occ[0]), lumo: vec(vir[0]), eUnit: 1 };
    },
  };
  return provider;
}
