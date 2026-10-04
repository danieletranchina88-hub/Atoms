// GFN2-xTB: tight binding autoconsistente con elettrostatica anisotropa fino ai quadrupoli atomici,
// terzo ordine diagonale per shell, dispersione D4 autoconsistente e termine a tre corpi ATM.
// Riferimento: C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 15, 1652 (2019).
// L'implementazione segue tblite (github.com/tblite/tblite) ed è validata su di esso in tests/gfn2.mjs.
// Unità interne: bohr e hartree.

import { GFN2 as PAR } from './gfn2Data.js';
import { buildBasis, moleculeIntegrals, nSph, QIDX, shellPair, pairBuffers } from './gto.js';
import { eigh } from '../linalg.js';
import { ALPBWater } from './alpb.js';

const KB = 3.166808578545117e-6; // costante di Boltzmann in hartree/K
const EL = PAR.elements;
const SHELL3 = [PAR.thirdorder.s, PAR.thirdorder.p, PAR.thirdorder.d];
const LCHAR = 'spdf';
const MP = PAR.multipole, DSP = PAR.dispersion, D4 = PAR.d4;
const QSCALE = [1, 2, 1, 2, 2, 1];
const CN_CUTOFF = 25, REP_CUTOFF = 25, D4CN_CUTOFF = 30, DISP2_CUTOFF = 50, DISP3_CUTOFF = 25;

export const GFN2_MAX_Z = 86;
export const gfn2Supports = (Z) => Z.every(z => z >= 1 && z <= GFN2_MAX_Z);

function kshell(l1, l2) {
  const k = PAR.hamiltonian.kshell;
  const a = LCHAR[Math.min(l1, l2)] + LCHAR[Math.max(l1, l2)];
  if (k[a] !== undefined) return k[a];
  return 0.5 * (k[LCHAR[l1] + LCHAR[l1]] + k[LCHAR[l2] + LCHAR[l2]]);
}

/** erf(x) con la serie a termini positivi erf(x) = 2x/√π e^{−x²} Σ (2x²)ⁿ/(2n+1)!!, esatta a ~1e-16. */
export function erf(x) {
  const ax = Math.abs(x);
  if (ax > 6) return Math.sign(x);
  const x2 = ax * ax;
  let term = 1, sum = 1;
  for (let n = 1; n < 200; n++) { term *= 2 * x2 / (2 * n + 1); sum += term; if (term < 1e-17 * sum) break; }
  return Math.sign(x) * 2 / Math.sqrt(Math.PI) * ax * Math.exp(-x2) * sum;
}

/* ------------------------------------------------------------------ */
/* numeri di coordinazione                                             */
/* ------------------------------------------------------------------ */

/** CN di GFN2: doppia funzione esponenziale (k = 10 e 20, raggio spostato di 2 bohr), raggi covalenti D3. */
function cnGFN2(Z, pos) {
  const N = Z.length, cn = new Float64Array(N), pairs = [];
  const rc = PAR.covD3;
  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) {
    const dx = pos[3 * i] - pos[3 * j], dy = pos[3 * i + 1] - pos[3 * j + 1], dz = pos[3 * i + 2] - pos[3 * j + 2];
    const r = Math.hypot(dx, dy, dz);
    if (r > CN_CUTOFF) continue;
    const r0 = rc[Z[i] - 1] + rc[Z[j] - 1];
    const e1 = Math.exp(-10 * (r0 / r - 1)), f1 = 1 / (1 + e1);
    const e2 = Math.exp(-20 * ((r0 + 2) / r - 1)), f2 = 1 / (1 + e2);
    const c = f1 * f2;
    cn[i] += c; cn[j] += c;
    const df1 = -f1 * f1 * e1 * 10 * r0 / (r * r), df2 = -f2 * f2 * e2 * 20 * (r0 + 2) / (r * r);
    pairs.push({ i, j, d: (df1 * f2 + f1 * df2) / r, dx, dy, dz }); // d: (dc/dr)/r
  }
  return { cn, pairs };
}

/** CN di D4: funzione errore (k = 7,5) pesata con la differenza di elettronegatività di Pauling. */
function cnD4(Z, pos) {
  const N = Z.length, cn = new Float64Array(N), pairs = [];
  const rc = PAR.covD3, en = PAR.pauling;
  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) {
    const dx = pos[3 * i] - pos[3 * j], dy = pos[3 * i + 1] - pos[3 * j + 1], dz = pos[3 * i + 2] - pos[3 * j + 2];
    const r = Math.hypot(dx, dy, dz);
    if (r > D4CN_CUTOFF) continue;
    const r0 = rc[Z[i] - 1] + rc[Z[j] - 1];
    const w = D4.k4 * Math.exp(-((Math.abs(en[Z[i] - 1] - en[Z[j] - 1]) + D4.k5) ** 2) / D4.k6);
    const t = -D4.kcn * (r / r0 - 1);
    const c = w * 0.5 * (1 + erf(t));
    cn[i] += c; cn[j] += c;
    const dc = w * Math.exp(-t * t) / Math.sqrt(Math.PI) * (-D4.kcn / r0);
    pairs.push({ i, j, d: dc / r, dx, dy, dz });
  }
  return { cn, pairs };
}

/** Accumula Σ_i v_i ∂CN_i/∂R nel gradiente g. */
function cnGradient(pairs, v, g) {
  for (const p of pairs) {
    const s = (v[p.i] + v[p.j]) * p.d;
    g[3 * p.i] += s * p.dx; g[3 * p.i + 1] += s * p.dy; g[3 * p.i + 2] += s * p.dz;
    g[3 * p.j] -= s * p.dx; g[3 * p.j + 1] -= s * p.dy; g[3 * p.j + 2] -= s * p.dz;
  }
}

/* ------------------------------------------------------------------ */
/* DFT-D4                                                             */
/* ------------------------------------------------------------------ */

const TRAPZD = [2.49995e-2, 4.99995e-2, 7.5e-2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.2, 0.2, 0.2, 0.2, 0.35, 0.5, 0.75, 1, 1.75, 2.5, 1.25];
const rc6Cache = new Map();
/** C6 di riferimento fra tutte le coppie di riferimenti degli elementi za e zb (integrazione di Casimir–Polder). */
function refC6(za, zb) {
  const key = za * 100 + zb;
  let m = rc6Cache.get(key);
  if (m) return m;
  const ra = PAR.d4ref[za - 1], rb = PAR.d4ref[zb - 1];
  m = new Float64Array(ra.length * rb.length);
  for (let a = 0; a < ra.length; a++) for (let b = 0; b < rb.length; b++) {
    let s = 0;
    for (let w = 0; w < 23; w++) s += TRAPZD[w] * ra[a].alpha[w] * rb[b].alpha[w];
    m[a * rb.length + b] = 3 / Math.PI * s;
  }
  rc6Cache.set(key, m);
  return m;
}

/** Pesi dei riferimenti D4: gaussiane nel CN (normalizzate) per la funzione di carica ζ(q). Con le derivate. */
function d4Weights(z, cn, q) {
  const refs = PAR.d4ref[z - 1], nr = refs.length;
  const wf = D4.wf, zeff = D4.zeff[z - 1], gam = D4.gam[z - 1] * D4.gc, ga = D4.ga;
  const expw = new Float64Array(nr), dexpw = new Float64Array(nr);
  let norm = 0, dnorm = 0, maxcn = -Infinity;
  for (let a = 0; a < nr; a++) {
    const dcn = cn - refs[a].cn, tmp = Math.exp(-dcn * dcn);
    maxcn = Math.max(maxcn, refs[a].cn);
    let e = 0, de = 0;
    for (let k = 1; k <= refs[a].count; k++) { const t = tmp ** (k * wf); e += t; de += -2 * k * wf * dcn * t; }
    expw[a] = e; dexpw[a] = de; norm += e; dnorm += de;
  }
  const gw = new Float64Array(nr), dgw = new Float64Array(nr);
  const exceptional = !(norm > 1e-300) || !Number.isFinite(1 / norm);
  for (let a = 0; a < nr; a++) {
    if (exceptional) { gw[a] = refs[a].cn === maxcn ? 1 : 0; dgw[a] = 0; continue; }
    gw[a] = expw[a] / norm;
    dgw[a] = (dexpw[a] - expw[a] * dnorm / norm) / norm;
  }
  const w = new Float64Array(nr), dwdcn = new Float64Array(nr), dwdq = new Float64Array(nr);
  const qm = q + zeff;
  for (let a = 0; a < nr; a++) {
    const qr = refs[a].q + zeff;
    let zeta, dzeta;
    if (qm > 0) {
      const scale = Math.exp(gam * (1 - qr / qm));
      zeta = Math.exp(ga * (1 - scale));
      dzeta = -ga * gam * scale * zeta * qr / (qm * qm);
    } else { zeta = Math.exp(ga); dzeta = 0; }
    w[a] = zeta * gw[a]; dwdcn[a] = zeta * dgw[a]; dwdq[a] = dzeta * gw[a];
  }
  return { w, dwdcn, dwdq };
}

/** Parametri di smorzamento razionale (Becke–Johnson) e matrice di dispersione a due corpi senza pesi. */
function d4Setup(Z, pos) {
  const N = Z.length;
  const cnd = cnD4(Z, pos);
  const r4r2 = D4.r4r2;
  const pairs = []; // {i, j, rc6, f (energia per C6 unitario), df/dr / r, dx…}
  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) {
    const dx = pos[3 * i] - pos[3 * j], dy = pos[3 * i + 1] - pos[3 * j + 1], dz = pos[3 * i + 2] - pos[3 * j + 2];
    const r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > DISP2_CUTOFF * DISP2_CUTOFF) continue;
    const qq = 3 * r4r2[Z[i] - 1] * r4r2[Z[j] - 1];
    const r0 = DSP.a1 * Math.sqrt(qq) + DSP.a2;
    const r6 = r2 ** 3, r8 = r6 * r2, t6 = 1 / (r6 + r0 ** 6), t8 = 1 / (r8 + r0 ** 8);
    const f = -(DSP.s6 * t6 + DSP.s8 * qq * t8);
    // d/dr divisa per r: d(1/(r⁶+c))/dr / r = −6 r⁴ t6²
    const df = DSP.s6 * 6 * r2 * r2 * t6 * t6 + DSP.s8 * qq * 8 * r6 * t8 * t8;
    pairs.push({ i, j, rc6: refC6(Z[i], Z[j]), f, df, dx, dy, dz });
  }
  return { cnd, pairs };
}

/** Energia D4 a due corpi con le cariche q: E = ½ Σ w_i·C6ref·w_j·f(r). Restituisce energia e dE/dq. */
function d4Energy(Z, setup, q, wantGrad = false) {
  const N = Z.length;
  const W = [];
  for (let i = 0; i < N; i++) W.push(d4Weights(Z[i], setup.cnd.cn[i], q[i]));
  let E = 0;
  const dEdq = new Float64Array(N), dEdcn = new Float64Array(N);
  const grad = wantGrad ? new Float64Array(3 * N) : null;
  for (const p of setup.pairs) {
    const wi = W[p.i], wj = W[p.j], nb = wj.w.length;
    let c6 = 0, dqi = 0, dqj = 0, dci = 0, dcj = 0;
    for (let a = 0; a < wi.w.length; a++) for (let b = 0; b < nb; b++) {
      const r = p.rc6[a * nb + b];
      c6 += wi.w[a] * r * wj.w[b];
      dqi += wi.dwdq[a] * r * wj.w[b]; dqj += wi.w[a] * r * wj.dwdq[b];
      dci += wi.dwdcn[a] * r * wj.w[b]; dcj += wi.w[a] * r * wj.dwdcn[b];
    }
    E += c6 * p.f;
    dEdq[p.i] += dqi * p.f; dEdq[p.j] += dqj * p.f;
    dEdcn[p.i] += dci * p.f; dEdcn[p.j] += dcj * p.f;
    if (grad) {
      const s = c6 * p.df;
      grad[3 * p.i] += s * p.dx; grad[3 * p.i + 1] += s * p.dy; grad[3 * p.i + 2] += s * p.dz;
      grad[3 * p.j] -= s * p.dx; grad[3 * p.j + 1] -= s * p.dy; grad[3 * p.j + 2] -= s * p.dz;
    }
  }
  if (grad) cnGradient(setup.cnd.pairs, dEdcn, grad);
  return { E, dEdq, grad };
}

/** Termine a tre corpi di Axilrod–Teller–Muto con C9 ≈ √|C6ij C6ik C6jk| (pesi a carica nulla) e smorzamento zero. */
function atmEnergy(Z, pos, cnd, wantGrad) {
  const N = Z.length;
  const W = [];
  for (let i = 0; i < N; i++) W.push(d4Weights(Z[i], cnd.cn[i], 0));
  const C6 = new Float64Array(N * N), dC6 = new Float64Array(N * N); // dC6[i*N+j] = ∂C6_ij/∂CN_i
  for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
    if (i === j) continue;
    const rc6 = refC6(Z[i], Z[j]), wi = W[i], wj = W[j], nb = wj.w.length;
    let c = 0, d = 0;
    for (let a = 0; a < wi.w.length; a++) for (let b = 0; b < nb; b++) { c += wi.w[a] * rc6[a * nb + b] * wj.w[b]; d += wi.dwdcn[a] * rc6[a * nb + b] * wj.w[b]; }
    C6[i * N + j] = c; dC6[i * N + j] = d;
  }
  const r4r2 = D4.r4r2, s9 = DSP.s9, alp = 16 / 3;
  const r0 = (i, j) => DSP.a1 * Math.sqrt(3 * r4r2[Z[i] - 1] * r4r2[Z[j] - 1]) + DSP.a2;
  const vec = (i, j) => [pos[3 * j] - pos[3 * i], pos[3 * j + 1] - pos[3 * i + 1], pos[3 * j + 2] - pos[3 * i + 2]];
  const cut2 = DISP3_CUTOFF * DISP3_CUTOFF;
  let E = 0;
  const grad = wantGrad ? new Float64Array(3 * N) : null;
  const dEdcn = new Float64Array(N);
  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) {
    const vij = vec(i, j), r2ij = vij[0] ** 2 + vij[1] ** 2 + vij[2] ** 2;
    if (r2ij > cut2) continue;
    for (let k = 0; k < j; k++) {
      const vik = vec(i, k), vjk = vec(j, k);
      const r2ik = vik[0] ** 2 + vik[1] ** 2 + vik[2] ** 2, r2jk = vjk[0] ** 2 + vjk[1] ** 2 + vjk[2] ** 2;
      if (r2ik > cut2 || r2jk > cut2) continue;
      const cij = C6[i * N + j], cik = C6[i * N + k], cjk = C6[j * N + k];
      const c9 = Math.sqrt(Math.abs(cij * cik * cjk));
      const R0 = r0(i, j) * r0(i, k) * r0(j, k);
      const r2 = r2ij * r2ik * r2jk, r1 = Math.sqrt(r2), r3 = r1 * r2, r5 = r3 * r2;
      const tt = (R0 / r1) ** alp, fdmp = 1 / (1 + 6 * tt);
      const s = (r2ij + r2jk - r2ik) * (r2ij - r2jk + r2ik) * (-r2ij + r2jk + r2ik);
      const ang = 0.375 * s / r5 + 1 / r3;
      const e = s9 * c9 * ang * fdmp;
      E += e;
      if (!grad) continue;
      // derivate rispetto ai quadrati delle distanze a = r²ij, r²ik, r²jk
      const de = [0, 1, 2].map(w => {
        const [a, b, c] = w === 0 ? [r2ij, r2jk, r2ik] : w === 1 ? [r2ik, r2jk, r2ij] : [r2jk, r2ij, r2ik];
        const u = a + b - c, v = a - b + c, x = -a + b + c;
        const ds = v * x + u * x - u * v;
        const dangda = 0.375 * (ds / r5 - 2.5 * s / (r5 * a)) - 1.5 / (r3 * a);
        // fdmp = 1/(1+6 (R0/r1)^alp), r1 = √(a b c): ∂fdmp/∂a = fdmp² · 6 tt · alp /(2a)
        const dfda = fdmp * fdmp * 6 * tt * alp / (2 * a);
        return s9 * c9 * (dangda * fdmp + ang * dfda);
      });
      const add = (p, q, v, f) => { for (let d = 0; d < 3; d++) { grad[3 * p + d] -= 2 * f * v[d]; grad[3 * q + d] += 2 * f * v[d]; } };
      add(i, j, vij, de[0]); add(i, k, vik, de[1]); add(j, k, vjk, de[2]);
      // dipendenza di C9 dai CN: ∂c9/∂C6ij = c9/(2 C6ij)
      if (c9 > 0) {
        const f = e / 2;
        dEdcn[i] += f * (dC6[i * N + j] / cij + dC6[i * N + k] / cik);
        dEdcn[j] += f * (dC6[j * N + i] / cij + dC6[j * N + k] / cjk);
        dEdcn[k] += f * (dC6[k * N + i] / cik + dC6[k * N + j] / cjk);
      }
    }
  }
  if (grad) cnGradient(cnd.pairs, dEdcn, grad);
  return { E, grad };
}

/* ------------------------------------------------------------------ */
/* repulsione                                                          */
/* ------------------------------------------------------------------ */

function repulsion(Z, pos, wantGrad) {
  const N = Z.length;
  let E = 0;
  const grad = wantGrad ? new Float64Array(3 * N) : null;
  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) {
    const dx = pos[3 * i] - pos[3 * j], dy = pos[3 * i + 1] - pos[3 * j + 1], dz = pos[3 * i + 2] - pos[3 * j + 2];
    const r = Math.hypot(dx, dy, dz);
    if (r > REP_CUTOFF) continue;
    const ei = EL[Z[i] - 1], ej = EL[Z[j] - 1];
    const k = Z[i] <= 2 && Z[j] <= 2 ? PAR.repulsion.klight : PAR.repulsion.kexp;
    const a = Math.sqrt(ei.arep * ej.arep), rk = r ** k;
    const ex = Math.exp(-a * rk), zz = ei.zeff * ej.zeff;
    const e = zz * ex / r;
    E += e;
    if (grad) {
      const de = -e * (a * k * rk + 1) / (r * r); // (dE/dr)/r
      grad[3 * i] += de * dx; grad[3 * i + 1] += de * dy; grad[3 * i + 2] += de * dz;
      grad[3 * j] -= de * dx; grad[3 * j + 1] -= de * dy; grad[3 * j + 2] -= de * dz;
    }
  }
  return { E, grad };
}

/* ------------------------------------------------------------------ */
/* calcolatore                                                         */
/* ------------------------------------------------------------------ */

export class GFN2xTB {
  /** field: campo elettrico uniforme in unità atomiche (hartree/(e·bohr)), oppure null. */
  /** solvent: 'water' per il solvente implicito ALPB, oppure null (fase gassosa). */
  constructor({ Tel = 300, maxIter = 300, etol = 1e-10, ptol = 1e-8, damp = 0.4, field = null, solvent = null } = {}) {
    Object.assign(this, { Tel, maxIter, etol, ptol, damp, field, solvent });
    this.guess = null; // vettore (cariche di shell, dipoli, quadrupoli) dell'ultimo calcolo, per ripartire
  }
  reset() { this.guess = null; }

  /**
   * Z: numeri atomici; pos: coordinate in bohr (3N); charge: carica totale; uhf: elettroni spaiati (Nα − Nβ).
   * Restituisce energia (hartree), gradiente (hartree/bohr) se richiesto, cariche, dipolo e dati della funzione d'onda.
   */
  compute(Z, pos, { charge = 0, uhf = 0, gradient = true } = {}) {
    const N = Z.length;
    if (!gfn2Supports(Z)) throw new Error('GFN2-xTB è parametrizzato solo da H a Rn');
    const key = Z.join(',');
    if (this.basisKey !== key) { this.basis = buildBasis(Z, EL); this.basisKey = key; this.guess = null; }
    const bas = this.basis, n = bas.nao, sh = bas.shells, nsh = sh.length;
    const { S, D, Q, derivs } = moleculeIntegrals(bas, pos, { withGrad: gradient });
    // quadrupoli a traccia nulla: 1,5 r_a r_b − ½ δ_ab r²
    const Qt = Q.map(() => new Float64Array(n * n));
    for (let k = 0; k < n * n; k++) {
      const tr = 0.5 * (Q[0][k] + Q[2][k] + Q[5][k]);
      Qt[0][k] = 1.5 * Q[0][k] - tr; Qt[1][k] = 1.5 * Q[1][k]; Qt[2][k] = 1.5 * Q[2][k] - tr;
      Qt[3][k] = 1.5 * Q[3][k]; Qt[4][k] = 1.5 * Q[4][k]; Qt[5][k] = 1.5 * Q[5][k] - tr;
    }

    // ---- Hamiltoniano di core ----
    const cng = cnGFN2(Z, pos);
    const cn = cng.cn;
    const se = new Float64Array(nsh);
    for (let I = 0; I < nsh; I++) se[I] = sh[I].ref.level - sh[I].ref.kcn * cn[sh[I].atom];
    // H0_μν = S_μν · ½(se_I + se_J) · pk_IJ; pk = π(R)·K fuori dal sito, 1 sullo stesso atomo; dpk = (dpk/dR)/R
    const hfac = new Float64Array(nsh * nsh), pk = new Float64Array(nsh * nsh), dpk = new Float64Array(nsh * nsh);
    const wexp = PAR.hamiltonian.wexp, enscale = PAR.hamiltonian.enscale, rad = PAR.atomicRadii;
    for (let I = 0; I < nsh; I++) for (let J = 0; J <= I; J++) {
      const A = sh[I].atom, B = sh[J].atom;
      let f = 1, df = 0;
      if (A !== B) {
        const r = Math.hypot(pos[3 * A] - pos[3 * B], pos[3 * A + 1] - pos[3 * B + 1], pos[3 * A + 2] - pos[3 * B + 2]);
        const rr = Math.sqrt(r / (rad[Z[A] - 1] + rad[Z[B] - 1]));
        const zi = sh[I].ref.zeta, zj = sh[J].ref.zeta, ki = sh[I].ref.shpoly, kj = sh[J].ref.shpoly;
        const den = EL[Z[A] - 1].en - EL[Z[B] - 1].en;
        const K = kshell(sh[I].l, sh[J].l) * (2 * Math.sqrt(zi * zj) / (zi + zj)) ** wexp * (1 + enscale * den * den);
        f = (1 + ki * rr) * (1 + kj * rr) * K;
        df = (ki * (1 + kj * rr) + kj * (1 + ki * rr)) * 0.5 * rr / r * K / r;
      }
      pk[I * nsh + J] = pk[J * nsh + I] = f;
      dpk[I * nsh + J] = dpk[J * nsh + I] = df;
      hfac[I * nsh + J] = hfac[J * nsh + I] = 0.5 * (se[I] + se[J]) * f;
    }
    const H0 = new Float64Array(n * n);
    for (let mu = 0; mu < n; mu++) for (let nu = 0; nu < n; nu++) H0[mu * n + nu] = S[mu * n + nu] * hfac[bas.aoShell[mu] * nsh + bas.aoShell[nu]];

    // ---- elettrostatica isotropa (secondo e terzo ordine, per shell) ----
    const eta = new Float64Array(nsh), hd = new Float64Array(nsh), n0 = new Float64Array(nsh);
    for (let I = 0; I < nsh; I++) {
      const el = EL[Z[sh[I].atom] - 1];
      eta[I] = el.gam * sh[I].ref.lgam; hd[I] = el.gam3 * SHELL3[sh[I].l]; n0[I] = sh[I].ref.refocc;
    }
    const gamma = new Float64Array(nsh * nsh);
    for (let I = 0; I < nsh; I++) for (let J = 0; J <= I; J++) {
      const A = sh[I].atom, B = sh[J].atom;
      const r2 = (pos[3 * A] - pos[3 * B]) ** 2 + (pos[3 * A + 1] - pos[3 * B + 1]) ** 2 + (pos[3 * A + 2] - pos[3 * B + 2]) ** 2;
      const e = 0.5 * (eta[I] + eta[J]);
      gamma[I * nsh + J] = gamma[J * nsh + I] = 1 / Math.sqrt(r2 + 1 / (e * e));
    }

    // ---- elettrostatica anisotropa: matrici atomiche ----
    const mrad = new Float64Array(N);
    for (let A = 0; A < N; A++) {
      const el = EL[Z[A] - 1];
      mrad[A] = el.mprad + (MP.rmax - el.mprad) / (1 + Math.exp(-MP.kexp * (cn[A] - el.mpvcn - MP.shift)));
    }
    const amSD = new Float64Array(N * N * 3), amDD = new Float64Array(N * N * 9), amSQ = new Float64Array(N * N * 6);
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      if (i === j) continue;
      const rij = [pos[3 * i] - pos[3 * j], pos[3 * i + 1] - pos[3 * j + 1], pos[3 * i + 2] - pos[3 * j + 2]];
      const r = Math.hypot(...rij), g1 = 1 / r, g3 = g1 * g1 * g1, g5 = g3 * g1 * g1;
      const rr = 0.5 * (mrad[i] + mrad[j]) * g1;
      const f3 = 1 / (1 + 6 * rr ** MP.dmp3), f5 = 1 / (1 + 6 * rr ** MP.dmp5);
      const o = i * N + j;
      for (let c = 0; c < 3; c++) amSD[3 * o + c] = rij[c] * g3 * f3;
      for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) amDD[9 * o + 3 * a + b] = (a === b ? g3 * f5 : 0) - 3 * rij[a] * rij[b] * g5 * f5;
      for (let c = 0; c < 6; c++) { const [u, w] = QIDX[c]; amSQ[6 * o + c] = (u === w ? 1 : 2) * rij[u] * rij[w] * g5 * f5; }
    }
    const dk = Z.map(z => EL[z - 1].dkernel), qk = Z.map(z => EL[z - 1].qkernel);

    // ---- dispersione ----
    const d4 = d4Setup(Z, pos);
    const atm = atmEnergy(Z, pos, d4.cnd, gradient);
    const rep = repulsion(Z, pos, gradient);
    let solv = null;
    if (this.solvent) {
      if (this.solvent !== 'water') throw new Error(`Solvente implicito non disponibile: ${this.solvent}`);
      solv = this._alpb ??= new ALPBWater();
      solv.setup(Z, pos);
    }

    // ---- elettroni e occupazioni ----
    let nel = -charge;
    for (let I = 0; I < nsh; I++) nel += n0[I];
    const na = (nel + uhf) / 2, nb = (nel - uhf) / 2;
    if (na < 0 || nb < 0) throw new Error('numero di elettroni e spin incompatibili');
    const kT = Math.max(KB * this.Tel, 1e-12);

    // ortogonalizzazione di Löwdin X = S^{−1/2}
    const es = eigh(S, n);
    const X = new Float64Array(n * n);
    for (let k = 0; k < n; k++) {
      const f = 1 / Math.sqrt(Math.max(es.values[k], 1e-10));
      for (let i = 0; i < n; i++) { const v = es.vectors[i * n + k] * f; if (v) for (let j = 0; j < n; j++) X[i * n + j] += v * es.vectors[j * n + k]; }
    }

    // ---- ciclo autoconsistente ----
    const nvar = nsh + 3 * N + 6 * N;
    let x = this.guess && this.guess.length === nvar ? Float64Array.from(this.guess) : new Float64Array(nvar);
    const mixer = new Broyden(nvar, this.damp);
    const F = new Float64Array(n * n), P = new Float64Array(n * n);
    let Eel = 0, Eold = 0, converged = false, iter = 0, last = null;
    for (iter = 1; iter <= this.maxIter; iter++) {
      const pot = this._potentials(x, { N, nsh, sh, gamma, hd, amSD, amDD, amSQ, dk, qk, Z, d4, pos, solv });
      this._fock(F, H0, S, D, Qt, pot, bas, nsh);
      const orb = solveFock(F, X, n);
      const occ = fermiOccupations(orb.e, na, nb, kT);
      density(P, orb, occ.f);
      const xo = this._moments(P, S, D, Qt, bas, n0, N, nsh);
      Eel = this._electronicEnergy(P, H0, xo, { N, nsh, sh, gamma, hd, amSD, amDD, amSQ, dk, qk, Z, d4, pos, solv }) + occ.ts;
      let diff = 0;
      for (let k = 0; k < nvar; k++) diff = Math.max(diff, Math.abs(xo[k] - x[k]));
      last = { orb, occ, xo };
      if (iter > 1 && diff < this.ptol && Math.abs(Eel - Eold) < this.etol) { converged = true; x = xo; break; }
      Eold = Eel;
      x = mixer.next(x, xo);
    }
    this.guess = Float64Array.from(last.xo);
    const { orb, occ, xo } = last;
    const qsh = xo.subarray(0, nsh);
    const qat = new Float64Array(N);
    for (let I = 0; I < nsh; I++) qat[sh[I].atom] += qsh[I];
    const dpat = xo.subarray(nsh, nsh + 3 * N), qpat = xo.subarray(nsh + 3 * N);
    const dipole = [0, 0, 0];
    for (let A = 0; A < N; A++) for (let c = 0; c < 3; c++) dipole[c] += pos[3 * A + c] * qat[A] + dpat[3 * A + c];
    const energy = Eel + rep.E + atm.E + (solv ? solv.eClassical : 0);
    const res = {
      energy, parts: { electronic: Eel, repulsion: rep.E, atm: atm.E, ts: occ.ts, solvation: solv ? solv.eClassical + solv.energy(qat) : 0 },
      charges: qat, shellCharges: Float64Array.from(qsh), atomicDipoles: Float64Array.from(dpat), atomicQuadrupoles: Float64Array.from(qpat),
      dipole, converged, iterations: iter, orbitalEnergies: orb.e, occupations: occ.f, nao: n,
      homoLumo: frontier(orb.e, occ.f), cn,
    };
    const PS = new Float64Array(n * n);
    for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const a = P[i * n + k]; if (a) for (let j = 0; j < n; j++) PS[i * n + j] += a * S[k * n + j]; }
    this.last = { Z: Z.slice(), pos: Float64Array.from(pos), P: Float64Array.from(P), PS, Ct: orb.Ct, nC: orb.nC, e: orb.e, f: occ.f, basis: bas, n };
    if (gradient) {
      res.gradient = this._gradient({ Z, pos, N, n, nsh, sh, bas, P, orb, occ, xo, cng, se, hfac, pk, dpk, S, D, Qt, derivs, x, eta, gamma, hd, mrad, dk, qk, d4, atm, rep, kT, amSD, amDD, amSQ, solv });
    }
    return res;
  }

  /** Potenziali dalla densità corrente x = (q_shell, dipoli, quadrupoli). */
  _potentials(x, c) {
    const { N, nsh, sh, gamma, hd, amSD, amDD, amSQ, dk, qk, Z, d4 } = c;
    const qsh = x.subarray(0, nsh), dp = x.subarray(nsh, nsh + 3 * N), qp = x.subarray(nsh + 3 * N);
    const vsh = new Float64Array(nsh);
    for (let I = 0; I < nsh; I++) {
      let s = 0;
      for (let J = 0; J < nsh; J++) s += gamma[I * nsh + J] * qsh[J];
      vsh[I] = s + hd[I] * qsh[I] * qsh[I];
    }
    const qat = new Float64Array(N);
    for (let I = 0; I < nsh; I++) qat[sh[I].atom] += qsh[I];
    const vat = new Float64Array(N), vd = new Float64Array(3 * N), vq = new Float64Array(6 * N);
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      if (i === j) continue;
      const o = i * N + j;
      let s = 0;
      for (let c2 = 0; c2 < 3; c2++) s += amSD[3 * o + c2] * dp[3 * j + c2];
      for (let c2 = 0; c2 < 6; c2++) s += amSQ[6 * o + c2] * qp[6 * j + c2];
      vat[i] += s;
      // potenziale sul dipolo e sul quadrupolo di j dovuto a i
      for (let c2 = 0; c2 < 3; c2++) {
        let t = amSD[3 * o + c2] * qat[i];
        for (let a = 0; a < 3; a++) t += amDD[9 * o + 3 * a + c2] * dp[3 * i + a];
        vd[3 * j + c2] += t;
      }
      for (let c2 = 0; c2 < 6; c2++) vq[6 * j + c2] += amSQ[6 * o + c2] * qat[i];
    }
    for (let j = 0; j < N; j++) {
      for (let c2 = 0; c2 < 3; c2++) vd[3 * j + c2] += 2 * dk[j] * dp[3 * j + c2];
      for (let c2 = 0; c2 < 6; c2++) vq[6 * j + c2] += 2 * qk[j] * QSCALE[c2] * qp[6 * j + c2];
    }
    const disp = d4Energy(Z, d4, qat);
    for (let i = 0; i < N; i++) vat[i] += disp.dEdq[i];
    // campo elettrico uniforme: φ(r) = −F·r agisce sulle cariche e sui dipoli atomici
    const fld = this.field;
    if (fld) for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) { vat[i] -= fld[d] * c.pos[3 * i + d]; vd[3 * i + d] -= fld[d]; }
    if (c.solv) c.solv.potential(qat, vat);
    return { vsh, vat, vd, vq };
  }

  _fock(F, H0, S, D, Qt, pot, bas, nsh) {
    const n = bas.nao, aoS = bas.aoShell, aoA = bas.aoAtom;
    const v = new Float64Array(n);
    for (let mu = 0; mu < n; mu++) v[mu] = pot.vsh[aoS[mu]] + pot.vat[aoA[mu]];
    for (let mu = 0; mu < n; mu++) for (let nu = 0; nu < n; nu++) {
      const k = mu * n + nu;
      F[k] = H0[k] - 0.5 * S[k] * (v[mu] + v[nu]);
    }
    // −(T + Tᵀ), T_μν = ½ Σ_c M_c[μ,ν] v_c(atomo di ν)
    for (let mu = 0; mu < n; mu++) for (let nu = 0; nu < n; nu++) {
      const k = mu * n + nu, B = aoA[nu];
      let t = 0;
      for (let c = 0; c < 3; c++) t += D[c][k] * pot.vd[3 * B + c];
      for (let c = 0; c < 6; c++) t += Qt[c][k] * pot.vq[6 * B + c];
      t *= 0.5;
      F[k] -= t; F[nu * n + mu] -= t;
    }
  }

  /** Cariche di shell (Mulliken), dipoli e quadrupoli atomici dalla matrice densità. */
  _moments(P, S, D, Qt, bas, n0, N, nsh) {
    const n = bas.nao, aoS = bas.aoShell, aoA = bas.aoAtom;
    const x = new Float64Array(nsh + 9 * N);
    for (let I = 0; I < nsh; I++) x[I] = n0[I];
    for (let mu = 0; mu < n; mu++) {
      let pop = 0;
      for (let nu = 0; nu < n; nu++) pop += P[mu * n + nu] * S[mu * n + nu];
      x[aoS[mu]] -= pop;
    }
    const off = nsh, offq = nsh + 3 * N;
    for (let mu = 0; mu < n; mu++) for (let nu = 0; nu < n; nu++) {
      const p = P[mu * n + nu], B = aoA[nu], k = mu * n + nu;
      if (!p) continue;
      for (let c = 0; c < 3; c++) x[off + 3 * B + c] -= p * D[c][k];
      for (let c = 0; c < 6; c++) x[offq + 6 * B + c] -= p * Qt[c][k];
    }
    return x;
  }

  _electronicEnergy(P, H0, x, c) {
    const { N, nsh, sh, gamma, hd, amSD, amDD, amSQ, dk, qk, Z, d4 } = c;
    let E = 0;
    for (let k = 0; k < P.length; k++) E += P[k] * H0[k];
    const qsh = x.subarray(0, nsh), dp = x.subarray(nsh, nsh + 3 * N), qp = x.subarray(nsh + 3 * N);
    for (let I = 0; I < nsh; I++) {
      let s = 0;
      for (let J = 0; J < nsh; J++) s += gamma[I * nsh + J] * qsh[J];
      E += 0.5 * qsh[I] * s + hd[I] * qsh[I] ** 3 / 3;
    }
    const qat = new Float64Array(N);
    for (let I = 0; I < nsh; I++) qat[sh[I].atom] += qsh[I];
    for (let j = 0; j < N; j++) {
      for (let c2 = 0; c2 < 3; c2++) E += dk[j] * dp[3 * j + c2] ** 2;
      for (let c2 = 0; c2 < 6; c2++) E += qk[j] * QSCALE[c2] * qp[6 * j + c2] ** 2;
    }
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      if (i === j) continue;
      const o = i * N + j;
      for (let c2 = 0; c2 < 3; c2++) {
        let t = amSD[3 * o + c2] * qat[i];
        for (let a = 0; a < 3; a++) t += 0.5 * amDD[9 * o + 3 * a + c2] * dp[3 * i + a];
        E += t * dp[3 * j + c2];
      }
      for (let c2 = 0; c2 < 6; c2++) E += amSQ[6 * o + c2] * qat[i] * qp[6 * j + c2];
    }
    E += d4Energy(Z, d4, qat).E;
    const fld = this.field;
    if (fld) for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) E -= fld[d] * (qat[i] * c.pos[3 * i + d] + dp[3 * i + d]);
    if (c.solv) E += c.solv.energy(qat);
    return E;
  }

  /** Gradiente analitico dell'energia libera elettronica (formula di Hellmann–Feynman con il termine di Pulay −W·dS). */
  _gradient(c) {
    const { Z, pos, N, n, nsh, sh, bas, P, orb, occ, xo, cng, se, hfac, pk, dpk, S, D, derivs, gamma, mrad, dk, qk, d4, atm, rep } = c;
    const g = new Float64Array(3 * N);
    for (let k = 0; k < 3 * N; k++) g[k] = rep.grad[k] + atm.grad[k];
    const pot = this._potentials(xo, c);
    const qsh = xo.subarray(0, nsh), dp = xo.subarray(nsh, nsh + 3 * N), qp = xo.subarray(nsh + 3 * N);
    const qat = new Float64Array(N);
    for (let I = 0; I < nsh; I++) qat[sh[I].atom] += qsh[I];
    if (c.solv) c.solv.gradient(qat, g);
    // forza del campo esterno sulle cariche atomiche (a cariche fissate)
    if (this.field) for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) g[3 * i + d] -= qat[i] * this.field[d];
    // matrice densità pesata con le energie orbitali
    const W = new Float64Array(n * n);
    const we = new Float64Array(orb.nC);
    for (let k = 0; k < orb.nC; k++) we[k] = occ.f[k] > 1e-14 ? occ.f[k] * orb.e[k] : 0;
    weightedSum(W, orb, we);
    const aoS = bas.aoShell, aoA = bas.aoAtom;
    const v = new Float64Array(n);
    for (let mu = 0; mu < n; mu++) v[mu] = pot.vsh[aoS[mu]] + pot.vat[aoA[mu]];
    // dipendenza delle energie di sito dal CN
    const dEdcn = new Float64Array(N), dEdse = new Float64Array(nsh);
    for (let mu = 0; mu < n; mu++) for (let nu = 0; nu < n; nu++) dEdse[aoS[mu]] += P[mu * n + nu] * S[mu * n + nu] * pk[aoS[mu] * nsh + aoS[nu]];
    for (let I = 0; I < nsh; I++) dEdcn[sh[I].atom] -= sh[I].ref.kcn * dEdse[I];
    // secondo ordine: ∂γ/∂R = −γ³ R
    for (let I = 0; I < nsh; I++) for (let J = 0; J < I; J++) {
      const A = sh[I].atom, B = sh[J].atom;
      if (A === B) continue;
      const gm = gamma[I * nsh + J], f = -qsh[I] * qsh[J] * gm * gm * gm;
      for (let d = 0; d < 3; d++) { const x = f * (pos[3 * A + d] - pos[3 * B + d]); g[3 * A + d] += x; g[3 * B + d] -= x; }
    }
    // elettrostatica anisotropa a momenti fissati, compresa la dipendenza dei raggi di smorzamento dal CN
    const dEdrad = new Float64Array(N);
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      if (i === j) continue;
      const rij = [pos[3 * i] - pos[3 * j], pos[3 * i + 1] - pos[3 * j + 1], pos[3 * i + 2] - pos[3 * j + 2]];
      const r = Math.hypot(...rij), rho = 0.5 * (mrad[i] + mrad[j]);
      const r2 = r * r, r4 = r2 * r2, rho3 = rho * rho * rho, rho4 = rho3 * rho;
      const h3 = 1 / (r2 * r + 6 * rho3), h3r = -3 * r2 * h3 * h3, h3p = -18 * rho * rho * h3 * h3;
      const h5 = 1 / (r4 * r + 6 * rho4 * r), h5r = -(5 * r4 + 6 * rho4) * h5 * h5, h5p = -24 * rho3 * r * h5 * h5;
      const den = r4 + 6 * rho4, k5 = r / den, k5r = (6 * rho4 - 3 * r4) / (den * den), k5p = -24 * rho3 * r / (den * den);
      const di = [dp[3 * i], dp[3 * i + 1], dp[3 * i + 2]], dj = [dp[3 * j], dp[3 * j + 1], dp[3 * j + 2]];
      const q = qp.subarray(6 * j, 6 * j + 6);
      const Qr = [q[0] * rij[0] + q[1] * rij[1] + q[3] * rij[2], q[1] * rij[0] + q[2] * rij[1] + q[4] * rij[2], q[3] * rij[0] + q[4] * rij[1] + q[5] * rij[2]];
      const a = rij[0] * dj[0] + rij[1] * dj[1] + rij[2] * dj[2], b = rij[0] * di[0] + rij[1] * di[1] + rij[2] * di[2];
      const ddij = di[0] * dj[0] + di[1] * dj[1] + di[2] * dj[2];
      const t = rij[0] * Qr[0] + rij[1] * Qr[1] + rij[2] * Qr[2];
      const qi = qat[i];
      const radial = (qi * a * h3r + 0.5 * (ddij * k5r - 3 * b * a * h5r) + qi * t * h5r) / r;
      for (let d = 0; d < 3; d++) {
        const x = qi * dj[d] * h3 - 1.5 * (di[d] * a + dj[d] * b) * h5 + 2 * qi * Qr[d] * h5 + radial * rij[d];
        g[3 * i + d] += x; g[3 * j + d] -= x;
      }
      const drho = qi * a * h3p + 0.5 * (ddij * k5p - 3 * b * a * h5p) + qi * t * h5p;
      dEdrad[i] += 0.5 * drho; dEdrad[j] += 0.5 * drho;
    }
    for (let A = 0; A < N; A++) {
      const el = EL[Z[A] - 1], e = Math.exp(-MP.kexp * (cng.cn[A] - el.mpvcn - MP.shift));
      dEdcn[A] += dEdrad[A] * (MP.rmax - el.mprad) * MP.kexp * e / ((1 + e) ** 2);
    }
    // dispersione D4 a cariche fissate
    const disp = d4Energy(Z, d4, qat, true);
    for (let k = 0; k < 3 * N; k++) g[k] += disp.grad[k];
    // derivate degli integrali: sovrapposizione (H0, cariche di Mulliken, Pulay), dipolo e quadrupolo
    const dQA = new Float64Array(6), tq = new Float64Array(6), tqa = new Float64Array(6);
    const traceless = (src, dst) => {
      const tr = 0.5 * (src[0] + src[2] + src[5]);
      dst[0] = 1.5 * src[0] - tr; dst[1] = 1.5 * src[1]; dst[2] = 1.5 * src[2] - tr;
      dst[3] = 1.5 * src[3]; dst[4] = 1.5 * src[4]; dst[5] = 1.5 * src[5] - tr;
    };
    const rawQ = new Float64Array(6);
    for (const dv of derivs) {
      const { I, J, R } = dv, si = sh[I], sj = sh[J], A = si.atom, B = sj.atom;
      {
        const nb = nSph(sj.l);
        const hIJ = hfac[I * nsh + J], sIJ = 0.5 * (se[I] + se[J]) * dpk[I * nsh + J];
        const vdA = pot.vd.subarray(3 * A, 3 * A + 3), vdB = pot.vd.subarray(3 * B, 3 * B + 3);
        const vqA = pot.vq.subarray(6 * A, 6 * A + 6), vqB = pot.vq.subarray(6 * B, 6 * B + 6);
        for (let i = 0; i < nSph(si.l); i++) for (let j = 0; j < nb; j++) {
          const mu = si.ao + i, nu = sj.ao + j, k = i * nb + j, o = mu * n + nu;
          const p = P[o], s = S[o];
          const ws = 2 * (p * (hIJ - 0.5 * (v[mu] + v[nu])) - W[o]);
          const DB = [D[0][o], D[1][o], D[2][o]];
          for (let gd = 0; gd < 3; gd++) {
            const dS = dv.dS[gd][k];
            let x = ws * dS - 2 * p * s * sIJ * R[gd];
            // dipoli: ⟨μ|(r−B)|ν⟩ pesato con il potenziale di B e ⟨ν|(r−A)|μ⟩ = ⟨μ|(r−B)|ν⟩ + R S con quello di A
            for (let cc = 0; cc < 3; cc++) {
              const dD = dv.dD[gd * 3 + cc][k];
              x -= p * vdB[cc] * dD;
              x -= p * vdA[cc] * (dD - (cc === gd ? s : 0) + R[cc] * dS);
            }
            // quadrupoli (derivate della forma a traccia nulla)
            for (let cc = 0; cc < 6; cc++) {
              const [u, w] = QIDX[cc];
              const dQ = dv.dQ[gd * 6 + cc][k];
              const dDu = dv.dD[gd * 3 + u][k], dDw = dv.dD[gd * 3 + w][k];
              rawQ[cc] = dQ;
              dQA[cc] = dQ - (u === gd ? DB[w] : 0) + R[u] * dDw + dDu * R[w] - (w === gd ? DB[u] : 0)
                - ((u === gd ? R[w] : 0) + (w === gd ? R[u] : 0)) * s + R[u] * R[w] * dS;
            }
            traceless(rawQ, tq); traceless(dQA, tqa);
            for (let cc = 0; cc < 6; cc++) x -= p * (vqB[cc] * tq[cc] + vqA[cc] * tqa[cc]);
            g[3 * A + gd] += x; g[3 * B + gd] -= x;
          }
        }
      }
    }
    cnGradient(cng.pairs, dEdcn, g);
    return g;
  }
}

/* ------------------------------------------------------------------ */
/* algebra della funzione d'onda                                       */
/* ------------------------------------------------------------------ */

function solveFock(F, X, n) {
  // F' = X F X (simmetrica: si calcola solo il triangolo superiore), poi autovalori e autovettori di F'
  const T = new Float64Array(n * n), Fp = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const a = X[i * n + k]; if (a) for (let j = 0; j < n; j++) T[i * n + j] += a * F[k * n + j]; }
  for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const a = T[i * n + k]; if (a) for (let j = i; j < n; j++) Fp[i * n + j] += a * X[k * n + j]; }
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) Fp[i * n + j] = Fp[j * n + i];
  const { values, vectors } = eigh(Fp, n);
  return { e: values, V: vectors, X, n, Ct: null, nC: 0 };
}

/** Coefficienti degli orbitali 0 … m−1 nella base atomica, per righe: Ct[k·n + μ] = C_μk = Σ_m X_μm V_mk. */
function backTransform(orb, m) {
  const { X, V, n } = orb;
  const Ct = new Float64Array(m * n);
  for (let k = 0; k < m; k++) for (let q = 0; q < n; q++) {
    const v = V[q * n + k];
    if (!v) continue;
    const row = k * n, xr = q * n;
    for (let i = 0; i < n; i++) Ct[row + i] += X[xr + i] * v;
  }
  orb.Ct = Ct; orb.nC = m;
}

/** Σ_k w_k c_k c_kᵀ sugli orbitali calcolati (matrice simmetrica). */
function weightedSum(out, orb, w) {
  const n = orb.n, Ct = orb.Ct;
  out.fill(0);
  for (let k = 0; k < orb.nC; k++) {
    const o = w[k];
    if (!o) continue;
    const row = k * n;
    for (let i = 0; i < n; i++) {
      const ci = Ct[row + i] * o;
      if (!ci) continue;
      const oi = i * n;
      for (let j = i; j < n; j++) out[oi + j] += ci * Ct[row + j];
    }
  }
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) out[i * n + j] = out[j * n + i];
}

/** Occupazioni di Fermi a temperatura elettronica kT con livelli di Fermi separati per α e β. */
function fermiOccupations(e, na, nb, kT) {
  const n = e.length, f = new Float64Array(n);
  let ts = 0;
  for (const nel of [na, nb]) {
    const { occ, ts: t } = fermiChannel(e, nel, kT);
    for (let k = 0; k < n; k++) f[k] += occ[k];
    ts += t;
  }
  return { f, ts };
}
function fermiChannel(e, nel, kT) {
  const n = e.length, occ = new Float64Array(n);
  if (nel <= 0) return { occ, ts: 0 };
  // Newton sul potenziale chimico, protetto da un intervallo di bisezione; partenza a metà fra HOMO e LUMO
  const k0 = Math.min(Math.max(Math.ceil(nel) - 1, 0), n - 1);
  let mu = k0 + 1 < n ? 0.5 * (e[k0] + e[k0 + 1]) : e[k0] + 1;
  let lo = e[0] - 1, hi = e[n - 1] + 1;
  for (let it = 0; it < 100; it++) {
    let cnt = 0, der = 0;
    for (let k = 0; k < n; k++) { const o = fermi((e[k] - mu) / kT); cnt += o; der += o * (1 - o); }
    const r = cnt - nel;
    if (Math.abs(r) < 1e-13) break;
    if (r < 0) lo = mu; else hi = mu;
    der /= kT;
    let next = der > 1e-30 ? mu - r / der : 0.5 * (lo + hi);
    if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
    if (Math.abs(next - mu) < 1e-15) break;
    mu = next;
  }
  let ts = 0;
  for (let k = 0; k < n; k++) {
    const o = fermi((e[k] - mu) / kT);
    occ[k] = o;
    if (o > 1e-300 && o < 1 - 1e-16) ts += kT * (o * Math.log(o) + (1 - o) * Math.log(1 - o));
  }
  return { occ, ts, mu };
}
const fermi = (x) => x > 0 ? Math.exp(-x) / (1 + Math.exp(-x)) : 1 / (1 + Math.exp(x));

/** Matrice densità dagli orbitali occupati (più il primo vuoto, utile per il disegno del LUMO). */
function density(P, orb, f) {
  let m = 0;
  for (let k = 0; k < f.length; k++) if (f[k] > 1e-14) m = k + 1;
  backTransform(orb, Math.min(f.length, m + 1));
  const w = new Float64Array(orb.nC);
  for (let k = 0; k < orb.nC; k++) w[k] = f[k] > 1e-14 ? f[k] : 0;
  weightedSum(P, orb, w);
}

function frontier(e, f) {
  let homo = -Infinity, lumo = Infinity;
  for (let k = 0; k < e.length; k++) { if (f[k] > 0.5) homo = Math.max(homo, e[k]); else lumo = Math.min(lumo, e[k]); }
  return { homo, lumo, gap: lumo - homo };
}

/** Miscelamento di Broyden modificato (D. D. Johnson, Phys. Rev. B 38, 12807, 1988). */
class Broyden {
  constructor(n, alpha, memory = 40) { this.n = n; this.alpha = alpha; this.mem = memory; this.dF = []; this.dX = []; this.prev = null; }
  next(xin, xout) {
    const n = this.n, a = this.alpha;
    const Fm = new Float64Array(n);
    for (let k = 0; k < n; k++) Fm[k] = xout[k] - xin[k];
    if (this.prev) {
      const dF = new Float64Array(n), dX = new Float64Array(n);
      let nrm = 0;
      for (let k = 0; k < n; k++) { dF[k] = Fm[k] - this.prev.F[k]; nrm += dF[k] * dF[k]; }
      nrm = Math.sqrt(nrm) || 1;
      for (let k = 0; k < n; k++) { dF[k] /= nrm; dX[k] = (xin[k] - this.prev.x[k]) / nrm; }
      this.dF.push(dF); this.dX.push(dX);
      if (this.dF.length > this.mem) { this.dF.shift(); this.dX.shift(); }
    }
    this.prev = { x: Float64Array.from(xin), F: Fm };
    const m = this.dF.length, out = new Float64Array(n);
    for (let k = 0; k < n; k++) out[k] = xin[k] + a * Fm[k];
    if (!m) return out;
    const w0 = 0.01;
    const A = new Float64Array(m * m), c = new Float64Array(m);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j <= i; j++) {
        let s = 0; for (let k = 0; k < n; k++) s += this.dF[i][k] * this.dF[j][k];
        A[i * m + j] = A[j * m + i] = s + (i === j ? w0 * w0 : 0);
      }
      let s = 0; for (let k = 0; k < n; k++) s += this.dF[i][k] * Fm[k];
      c[i] = s;
    }
    const g = solveSmall(A, c, m);
    for (let i = 0; i < m; i++) for (let k = 0; k < n; k++) out[k] -= g[i] * (a * this.dF[i][k] + this.dX[i][k]);
    return out;
  }
}
function solveSmall(A, b, n) {
  const M = Float64Array.from(A), x = Float64Array.from(b);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r * n + c]) > Math.abs(M[p * n + c])) p = r;
    if (p !== c) { for (let k = 0; k < n; k++) [M[c * n + k], M[p * n + k]] = [M[p * n + k], M[c * n + k]]; [x[c], x[p]] = [x[p], x[c]]; }
    const d = M[c * n + c] || 1e-300;
    for (let r = c + 1; r < n; r++) { const f = M[r * n + c] / d; if (!f) continue; for (let k = c; k < n; k++) M[r * n + k] -= f * M[c * n + k]; x[r] -= f * x[c]; }
  }
  for (let r = n - 1; r >= 0; r--) { let s = x[r]; for (let k = r + 1; k < n; k++) s -= M[r * n + k] * x[k]; x[r] = s / (M[r * n + r] || 1e-300); }
  return x;
}

export { shellPair, pairBuffers, cnGFN2, cnD4, d4Weights, refC6 };
