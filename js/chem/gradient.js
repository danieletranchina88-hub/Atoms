// Gradiente analitico dell'energia Hartree–Fock rispetto alle coordinate nucleari (forze sugli atomi).
//
//   dE/dX = Σ P_μν ∂H_μν/∂X + ½ Σ Γ_μνλσ ∂(μν|λσ)/∂X − Σ W_μν ∂S_μν/∂X + ∂V_NN/∂X
//   Γ_μνλσ = P_μν P_λσ − ½ (P^α_μλ P^α_νσ + P^α_μσ P^α_νλ + P^β_μλ P^β_νσ + P^β_μσ P^β_νλ)
//   W_μν = Σ_σ Σ_i^occ ε_i C_μi C_νi        (densità pesata con l'energia)
//
// La derivata di una gaussiana rispetto al suo centro è una combinazione di gaussiane con
// momento angolare ±1:  ∂/∂A_x [x_A^i e^{−a x_A²}] = 2a x_A^{i+1} e^{−a x_A²} − i x_A^{i−1} e^{−a x_A²}.
// Per invarianza traslazionale la derivata rispetto al quarto centro è −(A + B + C).

import { shellPair, hermiteR } from './integrals.js';

const eIdx = (lb, T, i, j, t) => (i * (lb + 1) + j) * T + t;

function energyWeighted(C, eps, nocc, n) {
  const W = new Float64Array(n * n);
  for (let k = 0; k < nocc; k++) {
    for (let i = 0; i < n; i++) {
      const ci = C[i * n + k] * eps[k];
      for (let j = 0; j < n; j++) W[i * n + j] += ci * C[j * n + k];
    }
  }
  return W;
}

/** Gradiente (array di [gx, gy, gz] per atomo, in hartree/bohr). */
export function hfGradient(res) {
  const { atoms, basis, Pa, Pb, Ca, Cb, epsA, epsB, nalpha, nbeta } = res;
  const n = basis.nbf;
  const nat = atoms.length;
  const grad = Array.from({ length: nat }, () => [0, 0, 0]);
  const Pt = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) Pt[i] = Pa[i] + Pb[i];
  const Wa = energyWeighted(Ca, epsA, nalpha, n);
  const Wb = energyWeighted(Cb, epsB, nbeta, n);
  const W = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) W[i] = Wa[i] + Wb[i];

  // --- repulsione nucleare ---
  for (let i = 0; i < nat; i++) {
    for (let j = 0; j < nat; j++) {
      if (i === j) continue;
      const a = atoms[i].xyz, b = atoms[j].xyz;
      const d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
      const r = Math.hypot(d[0], d[1], d[2]);
      const f = -atoms[i].Z * atoms[j].Z / (r * r * r);
      for (let k = 0; k < 3; k++) grad[i][k] += f * d[k];
    }
  }

  // --- termini monoelettronici ---
  const { shells } = basis;
  for (let A = 0; A < shells.length; A++) {
    for (let B = 0; B <= A; B++) {
      const sa = shells[A];
      const sb = shells[B];
      const pair = shellPair(sa, sb, 1, 3);
      const { lb, T } = pair;
      const L = sa.l + sb.l + 1;
      const D = L + 1;
      const fac = A === B ? 1 : 2;
      const atA = sa.atom, atB = sb.atom;
      for (const pr of pair.prims) {
        const pref = pr.c * Math.pow(Math.PI / pr.p, 1.5);
        const { Ex, Ey, Ez, a, b } = pr;
        const E3 = [Ex, Ey, Ez];
        const s1 = (E, i, j) => (i < 0 || j < 0 ? 0 : E[eIdx(lb, T, i, j, 0)]);
        const k1 = (E, i, j) => (i < 0 ? 0 : -2 * b * b * s1(E, i, j + 2) + b * (2 * j + 1) * s1(E, i, j) - 0.5 * j * (j - 1) * s1(E, i, j - 2));
        const dsA = (E, i, j) => 2 * a * s1(E, i + 1, j) - i * s1(E, i - 1, j);
        const dkA = (E, i, j) => 2 * a * k1(E, i + 1, j) - i * k1(E, i - 1, j);
        const Rs = atoms.map(C => hermiteR(L, pr.p, pr.P[0] - C.xyz[0], pr.P[1] - C.xyz[1], pr.P[2] - C.xyz[2]).slice(0, D * D * D));
        sa.comps.forEach((ca, ia) => {
          sb.comps.forEach((cb, ib) => {
            const mu = sa.offset + ia;
            const nu = sb.offset + ib;
            const nrm = sa.compNorm[ia] * sb.compNorm[ib];
            const p = Pt[mu * n + nu] * fac * nrm;
            const w = W[mu * n + nu] * fac * nrm;
            if (p === 0 && w === 0) return;
            const s = [s1(Ex, ca[0], cb[0]), s1(Ey, ca[1], cb[1]), s1(Ez, ca[2], cb[2])];
            const kk = [k1(Ex, ca[0], cb[0]), k1(Ey, ca[1], cb[1]), k1(Ez, ca[2], cb[2])];
            for (let d = 0; d < 3; d++) {
              const o1 = (d + 1) % 3, o2 = (d + 2) % 3;
              const E = E3[d];
              const dS = dsA(E, ca[d], cb[d]) * s[o1] * s[o2];
              const dT = dkA(E, ca[d], cb[d]) * s[o1] * s[o2] + dsA(E, ca[d], cb[d]) * (kk[o1] * s[o2] + s[o1] * kk[o2]);
              const g = pref * (p * dT - w * dS);
              grad[atA][d] += g;
              grad[atB][d] -= g;
            }
            // attrazione nucleare: derivate rispetto ad A, B e (per invarianza) al nucleo C
            for (let c = 0; c < atoms.length; c++) {
              const R = Rs[c];
              const Zc = atoms[c].Z;
              for (let d = 0; d < 3; d++) {
                let vA = 0, vB = 0;
                // coefficienti di Hermite derivati lungo d
                const ci = [ca[0], ca[1], ca[2]];
                const cj = [cb[0], cb[1], cb[2]];
                for (const [which, shiftSign] of [['A', 1], ['A', -1], ['B', 1], ['B', -1]]) {
                  const ii = [...ci];
                  const jj = [...cj];
                  let coef;
                  if (which === 'A') {
                    ii[d] += shiftSign;
                    coef = shiftSign > 0 ? 2 * a : -ci[d];
                  } else {
                    jj[d] += shiftSign;
                    coef = shiftSign > 0 ? 2 * b : -cj[d];
                  }
                  if (coef === 0 || ii[d] < 0 || jj[d] < 0) continue;
                  let v = 0;
                  for (let t = 0; t <= ii[0] + jj[0]; t++) {
                    const ex = Ex[eIdx(lb, T, ii[0], jj[0], t)];
                    if (ex === 0) continue;
                    for (let u = 0; u <= ii[1] + jj[1]; u++) {
                      const ey = Ey[eIdx(lb, T, ii[1], jj[1], u)];
                      if (ey === 0) continue;
                      for (let w2 = 0; w2 <= ii[2] + jj[2]; w2++) {
                        const ez = Ez[eIdx(lb, T, ii[2], jj[2], w2)];
                        if (ez === 0) continue;
                        v += ex * ey * ez * R[(t * D + u) * D + w2];
                      }
                    }
                  }
                  if (which === 'A') vA += coef * v; else vB += coef * v;
                }
                const f = -Zc * pr.c * 2 * Math.PI / pr.p * p;
                grad[atA][d] += f * vA;
                grad[atB][d] += f * vB;
                grad[c][d] -= f * (vA + vB);
              }
            }
          });
        });
      }
    }
  }

  // --- termini bielettronici ---
  twoElectronGradient(res, Pt, grad);
  return grad;
}

function twoElectronGradient(res, Pt, grad) {
  const { basis, Pa, Pb } = res;
  const n = basis.nbf;
  const { shells } = basis;
  const ns = shells.length;
  const pairs = [];
  for (let a = 0; a < ns; a++) for (let b = 0; b <= a; b++) pairs.push(shellPair(shells[a], shells[b], 1, 1));
  const Q = res.eri.schwarz; // √max|(ab|ab)| per coppia di gusci, stesso ordine
  // massimo di |P| per coppia di gusci, per lo screening con la densità
  const pmax = pairs.map(pr => {
    let m = 0;
    for (let i = 0; i < pr.sa.comps.length; i++) for (let j = 0; j < pr.sb.comps.length; j++) {
      const mu = pr.sa.offset + i, nu = pr.sb.offset + j;
      m = Math.max(m, Math.abs(Pt[mu * n + nu]), Math.abs(Pa[mu * n + nu]), Math.abs(Pb[mu * n + nu]));
    }
    return m;
  });
  const PI52 = 2 * Math.pow(Math.PI, 2.5);
  for (let P = 0; P < pairs.length; P++) {
    const bra = pairs[P];
    for (let K = 0; K <= P; K++) {
      if (Q && Q[P] * Q[K] * Math.max(pmax[P] * pmax[K], 1e-3 * Math.max(pmax[P], pmax[K])) < 1e-11) continue;
      quartetGradient(bra, pairs[K], P === K, n, Pt, Pa, Pb, grad, PI52);
    }
  }
}

// buffer riutilizzati tra i quartetti
let BUF_G0 = new Float64Array(1024);
let BUF_GC = new Float64Array(1024);
let BUF_GAM = new Float64Array(1024);
const BUF_H0 = new Float64Array(512);
const BUF_HC = new Float64Array(1536);
function ensure(buf, size) {
  return buf.length >= size ? buf : new Float64Array(Math.max(size, 2 * buf.length));
}

function quartetGradient(bra, ket, samePair, n, Pt, Pa, Pb, grad, PI52) {
  const { sa, sb } = bra;
  const { sa: sc, sb: sd } = ket;
  const na = sa.comps.length, nb = sb.comps.length, nc = sc.comps.length, nd = sd.comps.length;
  const deg = (sa === sb ? 1 : 2) * (sc === sd ? 1 : 2) * (samePair ? 1 : 2);
  const nk = nc * nd;
  BUF_GAM = ensure(BUF_GAM, na * nb * nk);
  const Gam = BUF_GAM;
  let gmax = 0;
  for (let ia = 0; ia < na; ia++) {
    const i = sa.offset + ia;
    for (let ib = 0; ib < nb; ib++) {
      const j = sb.offset + ib;
      for (let ic = 0; ic < nc; ic++) {
        const k = sc.offset + ic;
        for (let id = 0; id < nd; id++) {
          const l = sd.offset + id;
          const g = Pt[i * n + j] * Pt[k * n + l]
            - 0.5 * (Pa[i * n + k] * Pa[j * n + l] + Pa[i * n + l] * Pa[j * n + k]
              + Pb[i * n + k] * Pb[j * n + l] + Pb[i * n + l] * Pb[j * n + k]);
          const v = g * 0.5 * deg * sa.compNorm[ia] * sb.compNorm[ib] * sc.compNorm[ic] * sd.compNorm[id];
          Gam[((ia * nb + ib) * nc + ic) * nd + id] = v;
          if (Math.abs(v) > gmax) gmax = Math.abs(v);
        }
      }
    }
  }
  if (gmax < 1e-14) return;
  const Lab = sa.l + sb.l;
  const Lcd = sc.l + sd.l;
  const L = Lab + Lcd + 1;
  const D = L + 1;
  const HB = Lab + 1;
  const S1 = HB + 1;
  const nH = S1 * S1 * S1;
  BUF_G0 = ensure(BUF_G0, nk * nH);
  BUF_GC = ensure(BUF_GC, 3 * nk * nH);
  const G0 = BUF_G0, GC = BUF_GC, H0 = BUF_H0, HC = BUF_HC;
  const lbB = bra.lb, TB = bra.T, lbK = ket.lb, TK = ket.T;
  const atA = sa.atom, atB = sb.atom, atC = sc.atom, atD = sd.atom;
  const acc = new Float64Array(9); // dA, dB, dC per x, y, z
  // indici (t,u,v) con t+u+v ≤ m nel cubo di lato S1
  const hLists = [];
  for (let m = 0; m <= HB; m++) {
    const list = [];
    for (let t = 0; t <= m; t++) for (let u = 0; u <= m - t; u++) for (let v = 0; v <= m - t - u; v++) list.push((t * S1 + u) * S1 + v);
    hLists.push(Int32Array.from(list));
  }

  for (const pb of bra.prims) {
    if (pb.K < 1e-18) continue;
    for (const pk of ket.prims) {
      if (pb.K * pk.K * gmax < 1e-16) continue;
      const p = pb.p;
      const q = pk.p;
      const alpha = p * q / (p + q);
      const R = hermiteR(L, alpha, pb.P[0] - pk.P[0], pb.P[1] - pk.P[1], pb.P[2] - pk.P[2]);
      const pref = PI52 / (p * q * Math.sqrt(p + q)) * pb.c * pk.c;
      const c2 = 2 * pk.a;
      G0.fill(0, 0, nk * nH);
      GC.fill(0, 0, 3 * nk * nH);
      let kc = 0;
      for (let ic = 0; ic < nc; ic++) {
        const cc = sc.comps[ic];
        for (let id = 0; id < nd; id++, kc++) {
          const cd = sd.comps[id];
          addKet(G0, kc * nH, pk, lbK, TK, cc[0], cc[1], cc[2], cd, R, D, S1, Lab + 1, 1);
          for (let dir = 0; dir < 3; dir++) {
            const base = (dir * nk + kc) * nH;
            const u0 = cc[0] + (dir === 0 ? 1 : 0), u1 = cc[1] + (dir === 1 ? 1 : 0), u2 = cc[2] + (dir === 2 ? 1 : 0);
            addKet(GC, base, pk, lbK, TK, u0, u1, u2, cd, R, D, S1, Lab, c2);
            if (cc[dir] > 0) {
              const d0 = cc[0] - (dir === 0 ? 1 : 0), d1 = cc[1] - (dir === 1 ? 1 : 0), d2 = cc[2] - (dir === 2 ? 1 : 0);
              addKet(GC, base, pk, lbK, TK, d0, d1, d2, cd, R, D, S1, Lab, -cc[dir]);
            }
          }
        }
      }
      acc.fill(0);
      for (let ia = 0; ia < na; ia++) {
        const ca = sa.comps[ia];
        for (let ib = 0; ib < nb; ib++) {
          const cb = sb.comps[ib];
          const gbase = (ia * nb + ib) * nk;
          H0.fill(0, 0, nH);
          HC.fill(0, 0, 3 * nH);
          let any = false;
          const hl0 = hLists[HB], hl1 = hLists[Lab];
          for (let k = 0; k < nk; k++) {
            const g = Gam[gbase + k];
            if (g === 0) continue;
            any = true;
            const o0 = k * nH;
            for (let q = 0; q < hl0.length; q++) { const h = hl0[q]; H0[h] += g * G0[o0 + h]; }
            for (let dir = 0; dir < 3; dir++) {
              const oc = (dir * nk + k) * nH;
              const oh = dir * nH;
              for (let q = 0; q < hl1.length; q++) { const h = hl1[q]; HC[oh + h] += g * GC[oc + h]; }
            }
          }
          if (!any) continue;
          for (let dir = 0; dir < 3; dir++) {
            const ex = dir === 0 ? 1 : 0, ey = dir === 1 ? 1 : 0, ez = dir === 2 ? 1 : 0;
            let dA = 2 * pb.a * braSum(pb, lbB, TB, ca[0] + ex, ca[1] + ey, ca[2] + ez, cb[0], cb[1], cb[2], H0, 0, S1);
            if (ca[dir] > 0) dA -= ca[dir] * braSum(pb, lbB, TB, ca[0] - ex, ca[1] - ey, ca[2] - ez, cb[0], cb[1], cb[2], H0, 0, S1);
            let dB = 2 * pb.b * braSum(pb, lbB, TB, ca[0], ca[1], ca[2], cb[0] + ex, cb[1] + ey, cb[2] + ez, H0, 0, S1);
            if (cb[dir] > 0) dB -= cb[dir] * braSum(pb, lbB, TB, ca[0], ca[1], ca[2], cb[0] - ex, cb[1] - ey, cb[2] - ez, H0, 0, S1);
            const dC = braSum(pb, lbB, TB, ca[0], ca[1], ca[2], cb[0], cb[1], cb[2], HC, dir * nH, S1);
            acc[dir] += dA;
            acc[3 + dir] += dB;
            acc[6 + dir] += dC;
          }
        }
      }
      for (let dir = 0; dir < 3; dir++) {
        const dA = pref * acc[dir], dB = pref * acc[3 + dir], dC = pref * acc[6 + dir];
        grad[atA][dir] += dA;
        grad[atB][dir] += dB;
        grad[atC][dir] += dC;
        grad[atD][dir] -= dA + dB + dC;
      }
    }
  }
}

/** Σ_tuv E^{ab}_t E_u E_v H[off + h(t,u,v)] */
function braSum(pb, lb, T, a0, a1, a2, b0, b1, b2, H, off, S1) {
  let s = 0;
  const Ex = pb.Ex, Ey = pb.Ey, Ez = pb.Ez;
  const bx = (a0 * (lb + 1) + b0) * T, by = (a1 * (lb + 1) + b1) * T, bz = (a2 * (lb + 1) + b2) * T;
  for (let t = 0; t <= a0 + b0; t++) {
    const ex = Ex[bx + t];
    if (ex === 0) continue;
    for (let u = 0; u <= a1 + b1; u++) {
      const ey = Ey[by + u];
      if (ey === 0) continue;
      const exy = ex * ey;
      const row = off + (t * S1 + u) * S1;
      for (let v = 0; v <= a2 + b2; v++) s += exy * Ez[bz + v] * H[row + v];
    }
  }
  return s;
}

/** G[base + h(t,u,v)] += w Σ_τνφ (−1)^{τ+ν+φ} E^{cd} R_{t+τ,u+ν,v+φ}  per t+u+v ≤ tmax */
function addKet(G, base, pk, lb, T, c0, c1, c2, cd, R, D, S1, tmax, w) {
  const Ex = pk.Ex, Ey = pk.Ey, Ez = pk.Ez;
  const bx = (c0 * (lb + 1) + cd[0]) * T, by = (c1 * (lb + 1) + cd[1]) * T, bz = (c2 * (lb + 1) + cd[2]) * T;
  for (let tau = 0; tau <= c0 + cd[0]; tau++) {
    const ex = Ex[bx + tau];
    if (ex === 0) continue;
    for (let nu = 0; nu <= c1 + cd[1]; nu++) {
      const ey = Ey[by + nu];
      if (ey === 0) continue;
      for (let phi = 0; phi <= c2 + cd[2]; phi++) {
        const ez = Ez[bz + phi];
        if (ez === 0) continue;
        const e = w * ((tau + nu + phi) & 1 ? -1 : 1) * ex * ey * ez;
        for (let t = 0; t <= tmax; t++) {
          for (let u = 0; u <= tmax - t; u++) {
            const rowR = ((t + tau) * D + (u + nu)) * D + phi;
            const rowG = base + (t * S1 + u) * S1;
            for (let v = 0; v <= tmax - t - u; v++) G[rowG + v] += e * R[rowR + v];
          }
        }
      }
    }
  }
}
