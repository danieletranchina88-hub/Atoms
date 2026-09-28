// Integrali molecolari su gaussiane cartesiane contratte (schema di McMurchie–Davidson, 1978).
//
// Funzione primitiva: G = (x−A_x)^i (y−A_y)^j (z−A_z)^k · exp(−α|r−A|²)
// Il prodotto di due gaussiane è una gaussiana centrata in P = (αA + βB)/(α+β), sviluppata
// in gaussiane di Hermite Λ_tuv con coefficienti E_t^{ij} (ricorsione di McMurchie–Davidson):
//   E_0^{00} = exp(−μ X_AB²),  μ = αβ/(α+β)
//   E_t^{i+1,j} = E_{t−1}^{ij}/(2p) + X_PA E_t^{ij} + (t+1) E_{t+1}^{ij}
// Gli integrali coulombiani si riducono agli integrali di Hermite R_tuv, funzioni della
// funzione di Boys F_n(T).

import { BASIS_SETS } from './basisData.js';
import { boys } from './boys.js';

// Componenti cartesiane nell'ordine standard (xx, xy, xz, yy, yz, zz, …)
export function cartesianComponents(l) {
  const out = [];
  for (let x = l; x >= 0; x--) {
    for (let y = l - x; y >= 0; y--) out.push([x, y, l - x - y]);
  }
  return out;
}

function doubleFactorial(n) {
  let r = 1;
  for (let k = n; k > 1; k -= 2) r *= k;
  return r;
}

const L_NAMES = ['s', 'p', 'd', 'f'];
const AXES = ['x', 'y', 'z'];

function componentLabel(c) {
  const l = c[0] + c[1] + c[2];
  if (l === 0) return 's';
  return L_NAMES[l] + AXES[0].repeat(c[0]) + AXES[1].repeat(c[1]) + AXES[2].repeat(c[2]);
}

/**
 * Costruisce la base per una molecola.
 * atoms: [{ Z, xyz: [x,y,z] in bohr }]
 */
export function buildBasis(atoms, basisName) {
  const data = BASIS_SETS[basisName];
  if (!data) throw new Error(`Base ${basisName} non disponibile`);
  const shells = [];
  const functions = [];
  atoms.forEach((atom, ai) => {
    const list = data[atom.Z];
    if (!list) throw new Error(`La base ${basisName} non contiene l'elemento Z=${atom.Z}`);
    for (const [l, exps, coefs] of list) {
      const comps = cartesianComponents(l);
      // coefficienti normalizzati per la componente (l,0,0)
      const norm = exps.map(a => Math.pow(2 * a / Math.PI, 0.75) * Math.pow(4 * a, l / 2) / Math.sqrt(doubleFactorial(2 * l - 1)));
      let s = 0;
      for (let i = 0; i < exps.length; i++) {
        for (let j = 0; j < exps.length; j++) {
          const p = exps[i] + exps[j];
          s += coefs[i] * coefs[j] * norm[i] * norm[j] * Math.pow(Math.PI / p, 1.5) * doubleFactorial(2 * l - 1) / Math.pow(2 * p, l);
        }
      }
      const cnorm = 1 / Math.sqrt(s);
      const shell = {
        l,
        atom: ai,
        center: atom.xyz,
        exps: Float64Array.from(exps),
        coefs: Float64Array.from(coefs.map((c, i) => c * norm[i] * cnorm)),
        comps,
        // fattore per normalizzare ciascuna componente cartesiana
        compNorm: Float64Array.from(comps.map(c => Math.sqrt(doubleFactorial(2 * l - 1) / (doubleFactorial(2 * c[0] - 1) * doubleFactorial(2 * c[1] - 1) * doubleFactorial(2 * c[2] - 1))))),
        offset: functions.length,
      };
      shells.push(shell);
      comps.forEach((c) => functions.push({ shell: shells.length - 1, atom: ai, l, comp: c, label: componentLabel(c) }));
    }
  });
  return { name: basisName, shells, functions, nbf: functions.length };
}

// ---------------------------------------------------------------------------
// Coefficienti di Hermite E_t^{ij} in una dimensione
// ---------------------------------------------------------------------------

/** Tabella E[i][j][t] appiattita: indice (i·(lb+1) + j)·(la+lb+1) + t */
function hermiteE(la, lb, a, b, Ax, Bx) {
  const p = a + b;
  const mu = a * b / p;
  const Px = (a * Ax + b * Bx) / p;
  const XPA = Px - Ax;
  const XPB = Px - Bx;
  const XAB = Ax - Bx;
  const T = la + lb + 1;
  const E = new Float64Array((la + 1) * (lb + 1) * T);
  const idx = (i, j, t) => (i * (lb + 1) + j) * T + t;
  const oo2p = 0.5 / p;
  E[idx(0, 0, 0)] = Math.exp(-mu * XAB * XAB);
  for (let i = 0; i <= la; i++) {
    for (let j = 0; j <= lb; j++) {
      if (i === 0 && j === 0) continue;
      for (let t = 0; t <= i + j; t++) {
        let v = 0;
        if (i > 0) {
          const im = i - 1;
          if (t > 0) v += oo2p * E[idx(im, j, t - 1)];
          if (t <= im + j) v += XPA * E[idx(im, j, t)];
          if (t + 1 <= im + j) v += (t + 1) * E[idx(im, j, t + 1)];
        } else {
          const jm = j - 1;
          if (t > 0) v += oo2p * E[idx(i, jm, t - 1)];
          if (t <= i + jm) v += XPB * E[idx(i, jm, t)];
          if (t + 1 <= i + jm) v += (t + 1) * E[idx(i, jm, t + 1)];
        }
        E[idx(i, j, t)] = v;
      }
    }
  }
  return E;
}

// ---------------------------------------------------------------------------
// Integrali di Hermite R_tuv(α, P−C)
// ---------------------------------------------------------------------------

const Fbuf = new Float64Array(64);
const RMAX_D = 13;
const RbufA = new Float64Array(RMAX_D ** 3);
const RbufB = new Float64Array(RMAX_D ** 3);

/**
 * R^0_tuv per t+u+v ≤ L, indicizzato (t·(L+1) + u)·(L+1) + v.
 * Si calcola R^n per n = L … 0 alternando due buffer:
 *   R^n_{t+1,u,v} = t R^{n+1}_{t−1,u,v} + X_PC R^{n+1}_{t,u,v}   (analoghe per u, v)
 * Il risultato è in un buffer condiviso: va usato prima della chiamata successiva.
 */
export function hermiteR(L, alpha, X, Y, Z) {
  const D = L + 1;
  const T = alpha * (X * X + Y * Y + Z * Z);
  boys(T, L, Fbuf);
  let prev = RbufA;
  let cur = RbufB;
  const m2a = -2 * alpha;
  let f = Math.pow(m2a, L);
  for (let n = L; n >= 0; n--) {
    const top = L - n;
    cur[0] = f * Fbuf[n];
    f /= m2a;
    for (let tot = 1; tot <= top; tot++) {
      for (let t = 0; t <= tot; t++) {
        for (let u = 0; u <= tot - t; u++) {
          const v = tot - t - u;
          let val;
          if (t > 0) {
            val = X * prev[((t - 1) * D + u) * D + v];
            if (t > 1) val += (t - 1) * prev[((t - 2) * D + u) * D + v];
          } else if (u > 0) {
            val = Y * prev[(t * D + u - 1) * D + v];
            if (u > 1) val += (u - 1) * prev[(t * D + u - 2) * D + v];
          } else {
            val = Z * prev[(t * D + u) * D + v - 1];
            if (v > 1) val += (v - 1) * prev[(t * D + u) * D + v - 2];
          }
          cur[(t * D + u) * D + v] = val;
        }
      }
    }
    const tmp = prev; prev = cur; cur = tmp;
  }
  return prev;
}

// ---------------------------------------------------------------------------
// Dati per coppie di gusci
// ---------------------------------------------------------------------------

/** Prepara le primitive di una coppia di gusci (A, B). */
export function shellPair(sa, sb, extraA = 0, extraB = 0) {
  const la = sa.l + extraA;
  const lb = sb.l + extraB;
  const A = sa.center;
  const B = sb.center;
  const prims = [];
  for (let i = 0; i < sa.exps.length; i++) {
    for (let j = 0; j < sb.exps.length; j++) {
      const a = sa.exps[i];
      const b = sb.exps[j];
      const p = a + b;
      prims.push({
        a, b, p,
        P: [(a * A[0] + b * B[0]) / p, (a * A[1] + b * B[1]) / p, (a * A[2] + b * B[2]) / p],
        c: sa.coefs[i] * sb.coefs[j],
        Ex: hermiteE(la, lb, a, b, A[0], B[0]),
        Ey: hermiteE(la, lb, a, b, A[1], B[1]),
        Ez: hermiteE(la, lb, a, b, A[2], B[2]),
      });
      const pr = prims[prims.length - 1];
      // stima della grandezza della coppia (per lo screening delle primitive)
      pr.K = Math.abs(pr.c) * pr.Ex[0] * pr.Ey[0] * pr.Ez[0] * Math.pow(Math.PI / p, 1.25);
    }
  }
  return { sa, sb, la, lb, prims, T: la + lb + 1 };
}

const eIdx = (lb, T, i, j, t) => (i * (lb + 1) + j) * T + t;

// ---------------------------------------------------------------------------
// Integrali monoelettronici: sovrapposizione S, energia cinetica T, attrazione nucleare V, dipolo
// ---------------------------------------------------------------------------

export function oneElectron(basis, atoms) {
  const n = basis.nbf;
  const S = new Float64Array(n * n);
  const Tm = new Float64Array(n * n);
  const V = new Float64Array(n * n);
  const Dx = new Float64Array(n * n);
  const Dy = new Float64Array(n * n);
  const Dz = new Float64Array(n * n);
  const { shells } = basis;
  for (let A = 0; A < shells.length; A++) {
    for (let B = 0; B <= A; B++) {
      const sa = shells[A];
      const sb = shells[B];
      // per l'energia cinetica serve l_B + 2
      const pair = shellPair(sa, sb, 0, 2);
      const { lb, T } = pair;
      const la0 = sa.l;
      const lb0 = sb.l;
      const L = la0 + lb0;
      for (const pr of pair.prims) {
        const pref = pr.c * Math.pow(Math.PI / pr.p, 1.5);
        const { Ex, Ey, Ez, b } = pr;
        const s1 = (E, i, j) => (j < 0 ? 0 : E[eIdx(lb, T, i, j, 0)]);
        const kin1 = (E, i, j) => -2 * b * b * s1(E, i, j + 2) + b * (2 * j + 1) * s1(E, i, j) - 0.5 * j * (j - 1) * s1(E, i, j - 2);
        const mom1 = (E, i, j, XPC) => E[eIdx(lb, T, i, j, 1)] + XPC * E[eIdx(lb, T, i, j, 0)];
        // attrazione nucleare
        const Rs = atoms.map(C => hermiteR(L, pr.p, pr.P[0] - C.xyz[0], pr.P[1] - C.xyz[1], pr.P[2] - C.xyz[2]).slice(0, (L + 1) ** 3));
        const D = L + 1;
        sa.comps.forEach((ca, ia) => {
          sb.comps.forEach((cb, ib) => {
            const mu = sa.offset + ia;
            const nu = sb.offset + ib;
            if (nu > mu) return;
            const nrm = sa.compNorm[ia] * sb.compNorm[ib];
            const sx = s1(Ex, ca[0], cb[0]);
            const sy = s1(Ey, ca[1], cb[1]);
            const sz = s1(Ez, ca[2], cb[2]);
            S[mu * n + nu] += nrm * pref * sx * sy * sz;
            Tm[mu * n + nu] += nrm * pref * (kin1(Ex, ca[0], cb[0]) * sy * sz + sx * kin1(Ey, ca[1], cb[1]) * sz + sx * sy * kin1(Ez, ca[2], cb[2]));
            // dipolo rispetto all'origine
            Dx[mu * n + nu] += nrm * pref * mom1(Ex, ca[0], cb[0], pr.P[0]) * sy * sz;
            Dy[mu * n + nu] += nrm * pref * sx * mom1(Ey, ca[1], cb[1], pr.P[1]) * sz;
            Dz[mu * n + nu] += nrm * pref * sx * sy * mom1(Ez, ca[2], cb[2], pr.P[2]);
            let v = 0;
            for (let t = 0; t <= ca[0] + cb[0]; t++) {
              const ex = Ex[eIdx(lb, T, ca[0], cb[0], t)];
              if (ex === 0) continue;
              for (let u = 0; u <= ca[1] + cb[1]; u++) {
                const ey = Ey[eIdx(lb, T, ca[1], cb[1], u)];
                if (ey === 0) continue;
                for (let w = 0; w <= ca[2] + cb[2]; w++) {
                  const ez = Ez[eIdx(lb, T, ca[2], cb[2], w)];
                  if (ez === 0) continue;
                  let r = 0;
                  for (let c = 0; c < atoms.length; c++) r -= atoms[c].Z * Rs[c][(t * D + u) * D + w];
                  v += ex * ey * ez * r;
                }
              }
            }
            V[mu * n + nu] += nrm * pr.c * 2 * Math.PI / pr.p * v;
          });
        });
      }
    }
  }
  for (const M of [S, Tm, V, Dx, Dy, Dz]) {
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) M[j * n + i] = M[i * n + j];
  }
  return { S, T: Tm, V, dipole: [Dx, Dy, Dz] };
}

// ---------------------------------------------------------------------------
// Integrali bielettronici (μν|λσ), con simmetria a 8 e screening di Schwarz
// ---------------------------------------------------------------------------

/** Blocco di integrali per un quartetto di gusci: out[((ia·nb + ib)·nc + ic)·nd + id] */
export function eriQuartet(bra, ket) {
  const { sa, sb } = bra;
  const { sa: sc, sb: sd } = ket;
  const na = sa.comps.length, nb = sb.comps.length, nc = sc.comps.length, nd = sd.comps.length;
  const out = new Float64Array(na * nb * nc * nd);
  const Lab = sa.l + sb.l;
  const Lcd = sc.l + sd.l;
  const L = Lab + Lcd;
  const D = L + 1;
  const nHab = (Lab + 1) ** 3;
  const G = new Float64Array(nc * nd * nHab);
  const lbB = bra.lb, TB = bra.T, lbK = ket.lb, TK = ket.T;
  const PI52 = 2 * Math.pow(Math.PI, 2.5);
  for (const pb of bra.prims) {
    if (pb.K < 1e-18) continue;
    for (const pk of ket.prims) {
      if (pb.K * pk.K < 1e-15) continue;
      const p = pb.p;
      const q = pk.p;
      const alpha = p * q / (p + q);
      const R = hermiteR(L, alpha, pb.P[0] - pk.P[0], pb.P[1] - pk.P[1], pb.P[2] - pk.P[2]);
      const pref = PI52 / (p * q * Math.sqrt(p + q)) * pb.c * pk.c;
      // G[ket comp][t,u,v] = Σ_τνφ (−1)^{τ+ν+φ} E^{cd} R_{t+τ,u+ν,v+φ}
      G.fill(0);
      let kc = 0;
      for (let ic = 0; ic < nc; ic++) {
        const cc = sc.comps[ic];
        for (let id = 0; id < nd; id++, kc++) {
          const cd = sd.comps[id];
          const base = kc * nHab;
          for (let tau = 0; tau <= cc[0] + cd[0]; tau++) {
            const ex = pk.Ex[eIdx(lbK, TK, cc[0], cd[0], tau)];
            if (ex === 0) continue;
            for (let nu = 0; nu <= cc[1] + cd[1]; nu++) {
              const ey = pk.Ey[eIdx(lbK, TK, cc[1], cd[1], nu)];
              if (ey === 0) continue;
              for (let phi = 0; phi <= cc[2] + cd[2]; phi++) {
                const ez = pk.Ez[eIdx(lbK, TK, cc[2], cd[2], phi)];
                if (ez === 0) continue;
                const e = ((tau + nu + phi) & 1 ? -1 : 1) * ex * ey * ez;
                for (let t = 0; t <= Lab; t++) {
                  for (let u = 0; u <= Lab - t; u++) {
                    const rowR = ((t + tau) * D + (u + nu)) * D + phi;
                    const rowG = base + (t * (Lab + 1) + u) * (Lab + 1);
                    for (let v = 0; v <= Lab - t - u; v++) G[rowG + v] += e * R[rowR + v];
                  }
                }
              }
            }
          }
        }
      }
      // contrazione con E^{ab}
      let o = 0;
      for (let ia = 0; ia < na; ia++) {
        const ca = sa.comps[ia];
        for (let ib = 0; ib < nb; ib++) {
          const cb = sb.comps[ib];
          const obase = (ia * nb + ib) * nc * nd;
          for (let t = 0; t <= ca[0] + cb[0]; t++) {
            const ex = pb.Ex[eIdx(lbB, TB, ca[0], cb[0], t)];
            if (ex === 0) continue;
            for (let u = 0; u <= ca[1] + cb[1]; u++) {
              const ey = pb.Ey[eIdx(lbB, TB, ca[1], cb[1], u)];
              if (ey === 0) continue;
              for (let v = 0; v <= ca[2] + cb[2]; v++) {
                const ez = pb.Ez[eIdx(lbB, TB, ca[2], cb[2], v)];
                if (ez === 0) continue;
                const e = pref * ex * ey * ez;
                const h = (t * (Lab + 1) + u) * (Lab + 1) + v;
                for (let k = 0; k < nc * nd; k++) out[obase + k] += e * G[k * nHab + h];
              }
            }
          }
          o++;
        }
      }
    }
  }
  // normalizzazione delle componenti cartesiane
  let idx = 0;
  for (let ia = 0; ia < na; ia++) for (let ib = 0; ib < nb; ib++) for (let ic = 0; ic < nc; ic++) for (let id = 0; id < nd; id++) {
    out[idx++] *= sa.compNorm[ia] * sb.compNorm[ib] * sc.compNorm[ic] * sd.compNorm[id];
  }
  return out;
}

/**
 * Tutti gli integrali bielettronici unici (i≥j, k≥l, ij≥kl) sopra la soglia.
 * Restituisce { idx: Uint16Array (4 per integrale), val: Float64Array, count }.
 */
export function twoElectron(basis, { threshold = 1e-12, onProgress } = {}) {
  const { shells } = basis;
  const ns = shells.length;
  const pairs = [];
  for (let a = 0; a < ns; a++) for (let b = 0; b <= a; b++) pairs.push(shellPair(shells[a], shells[b]));
  // screening di Schwarz: Q_ab = √max|(ab|ab)|
  const Q = pairs.map(pr => {
    const blk = eriQuartet(pr, pr);
    const na = pr.sa.comps.length, nb = pr.sb.comps.length;
    let m = 0;
    for (let i = 0; i < na; i++) for (let j = 0; j < nb; j++) {
      const k = ((i * nb + j) * na + i) * nb + j;
      m = Math.max(m, Math.abs(blk[k]));
    }
    return Math.sqrt(m);
  });
  let cap = 1 << 16;
  let idx = new Uint16Array(cap * 4);
  let val = new Float64Array(cap);
  let count = 0;
  const push = (i, j, k, l, v) => {
    if (count === cap) {
      cap *= 2;
      const ni = new Uint16Array(cap * 4); ni.set(idx); idx = ni;
      const nv = new Float64Array(cap); nv.set(val); val = nv;
    }
    idx[4 * count] = i; idx[4 * count + 1] = j; idx[4 * count + 2] = k; idx[4 * count + 3] = l;
    val[count++] = v;
  };
  for (let P = 0; P < pairs.length; P++) {
    const bra = pairs[P];
    for (let K = 0; K <= P; K++) {
      if (Q[P] * Q[K] < threshold) continue;
      const ket = pairs[K];
      const blk = eriQuartet(bra, ket);
      const { sa, sb } = bra;
      const { sa: sc, sb: sd } = ket;
      const na = sa.comps.length, nb = sb.comps.length, nc = sc.comps.length, nd = sd.comps.length;
      let o = 0;
      for (let ia = 0; ia < na; ia++) {
        const i = sa.offset + ia;
        for (let ib = 0; ib < nb; ib++) {
          const j = sb.offset + ib;
          for (let ic = 0; ic < nc; ic++) {
            const k = sc.offset + ic;
            for (let id = 0; id < nd; id++, o++) {
              const l = sd.offset + id;
              if (j > i || l > k) continue;
              const ij = i * (i + 1) / 2 + j;
              const kl = k * (k + 1) / 2 + l;
              const v = blk[o];
              if (kl > ij) {
                // nello stesso blocco diagonale il simmetrico (kl|ij) compare già; altrimenti va salvato qui
                if (P === K) continue;
                if (Math.abs(v) >= threshold) push(k, l, i, j, v);
                continue;
              }
              if (Math.abs(v) < threshold) continue;
              push(i, j, k, l, v);
            }
          }
        }
      }
    }
    if (onProgress && (P & 15) === 0) onProgress(P / pairs.length);
  }
  return { idx: idx.subarray(0, 4 * count), val: val.subarray(0, count), count, schwarz: Q };
}

/**
 * Costruisce le matrici di Coulomb J e di scambio K per una o due densità.
 * dens: array di matrici densità; restituisce J (dalla densità totale) e K per ciascuna densità.
 */
export function buildJK(eri, dens, n) {
  const J = new Float64Array(n * n);
  const Ks = dens.map(() => new Float64Array(n * n));
  const Pt = new Float64Array(n * n);
  for (const P of dens) for (let i = 0; i < n * n; i++) Pt[i] += P[i];
  const { idx, val, count } = eri;
  const nd = dens.length;
  for (let q = 0; q < count; q++) {
    const i = idx[4 * q], j = idx[4 * q + 1], k = idx[4 * q + 2], l = idx[4 * q + 3];
    let v = val[q];
    if (i === j) v *= 0.5;
    if (k === l) v *= 0.5;
    if (i === k && j === l) v *= 0.5;
    // Coulomb: J_ij += (ij|kl) P_kl sulle 8 permutazioni
    const jij = v * (Pt[k * n + l] + Pt[l * n + k]);
    const jkl = v * (Pt[i * n + j] + Pt[j * n + i]);
    J[i * n + j] += jij; J[j * n + i] += jij;
    J[k * n + l] += jkl; J[l * n + k] += jkl;
    // Scambio: K_ik += (ij|kl) P_jl
    for (let d = 0; d < nd; d++) {
      const P = dens[d];
      const K = Ks[d];
      K[i * n + k] += v * P[j * n + l];
      K[j * n + k] += v * P[i * n + l];
      K[i * n + l] += v * P[j * n + k];
      K[j * n + l] += v * P[i * n + k];
      K[k * n + i] += v * P[l * n + j];
      K[l * n + i] += v * P[k * n + j];
      K[k * n + j] += v * P[l * n + i];
      K[l * n + j] += v * P[k * n + i];
    }
  }
  return { J, K: Ks };
}

/** Valori delle funzioni di base in un punto (per disegnare orbitali e densità). */
export function basisValues(basis, x, y, z, out) {
  const { shells } = basis;
  let f = 0;
  for (const sh of shells) {
    const dx = x - sh.center[0];
    const dy = y - sh.center[1];
    const dz = z - sh.center[2];
    const r2 = dx * dx + dy * dy + dz * dz;
    let radial = 0;
    for (let k = 0; k < sh.exps.length; k++) {
      const e = sh.exps[k] * r2;
      if (e < 40) radial += sh.coefs[k] * Math.exp(-e);
    }
    for (let c = 0; c < sh.comps.length; c++) {
      const [i, j, k] = sh.comps[c];
      let ang = sh.compNorm[c];
      for (let a = 0; a < i; a++) ang *= dx;
      for (let a = 0; a < j; a++) ang *= dy;
      for (let a = 0; a < k; a++) ang *= dz;
      out[f++] = radial * ang;
    }
  }
  return out;
}
