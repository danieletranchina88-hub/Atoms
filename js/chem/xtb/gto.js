// Integrali su gaussiane contratte per GFN2-xTB: sovrapposizione, dipolo e quadrupolo (operatori centrati
// sull'atomo della funzione "ket") e le loro derivate rispetto alla posizione dell'atomo "bra".
// Le funzioni sono armoniche sferiche reali nell'ordine di tblite (m = −l … l): p = (y, z, x),
// d = (xy, yz, z², xz, x²−y²). Ogni funzione contratta è rinormalizzata sulla propria autosovrapposizione,
// come fanno tblite e dxtb.

const CART = [
  [[0, 0, 0]],
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
  [[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]],
];
const S3 = Math.sqrt(3);
// righe: funzioni sferiche; colonne: cartesiane di CART[l]
const TRAFO = [
  [[1]],
  [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
  [
    [0, 0, 0, S3, 0, 0],          // xy
    [0, 0, 0, 0, 0, S3],          // yz
    [-0.5, -0.5, 1, 0, 0, 0],     // z²
    [0, 0, 0, 0, S3, 0],          // xz
    [S3 / 2, -S3 / 2, 0, 0, 0, 0], // x² − y²
  ],
];
export const nSph = (l) => 2 * l + 1;

// coefficienti binomiali fino a 8
const BINOM = Array.from({ length: 9 }, (_, n) => {
  const r = [1];
  for (let k = 1; k <= n; k++) r.push(r[k - 1] * (n - k + 1) / k);
  return r;
});
// momenti gaussiani ∫ t^k e^{−p t²} dt / √(π/p) = (k−1)!!/(2p)^{k/2} per k pari
function gaussMoments(p, kmax, out) {
  out[0] = 1;
  const h = 1 / (2 * p);
  for (let k = 1; k <= kmax; k++) out[k] = k % 2 ? 0 : out[k - 2] * (k - 1) * h;
}

const MAXI = 4, MAXJ = 5; // i ≤ la+1 ≤ 3, j ≤ lb+2 ≤ 4
const G = new Float64Array(10);
const powA = new Float64Array(MAXI + 1), powB = new Float64Array(MAXJ + 1);
// s1d[d][i*(MAXJ+1)+j] = ∫ (x−A)^i (x−B)^j e^{−p(x−P)²} dx / √(π/p)
function table1d(PA, PB, imax, jmax, out) {
  powA[0] = powB[0] = 1;
  for (let i = 1; i <= imax; i++) powA[i] = powA[i - 1] * PA;
  for (let j = 1; j <= jmax; j++) powB[j] = powB[j - 1] * PB;
  for (let i = 0; i <= imax; i++) {
    for (let j = 0; j <= jmax; j++) {
      let s = 0;
      for (let m = 0; m <= i; m++) {
        const cm = BINOM[i][m] * powA[i - m];
        for (let n = m % 2; n <= j; n += 2) s += cm * BINOM[j][n] * powB[j - n] * G[m + n];
        // m+n deve essere pari: n parte con la stessa parità di m
      }
      out[i * (MAXJ + 1) + j] = s;
    }
  }
}

/**
 * Base di una molecola. elements: dati per elemento (shells con l, alpha, coeff) indicizzati per Z.
 * Restituisce { shells: [{atom, l, alpha, coeff, ao, norm}], nao, aoAtom, aoShell, atomAO }.
 */
export function buildBasis(Z, elements) {
  const shells = [], atomShells = [];
  let ao = 0;
  for (let A = 0; A < Z.length; A++) {
    const el = elements[Z[A] - 1];
    const list = [];
    for (const sh of el.shells) {
      list.push(shells.length);
      shells.push({ atom: A, l: sh.l, alpha: sh.alpha, coeff: sh.coeff, ao, norm: shellNorm(sh), ref: sh });
      ao += nSph(sh.l);
    }
    atomShells.push(list);
  }
  const aoAtom = new Int32Array(ao), aoShell = new Int32Array(ao);
  shells.forEach((s, k) => { for (let m = 0; m < nSph(s.l); m++) { aoAtom[s.ao + m] = s.atom; aoShell[s.ao + m] = k; } });
  return { shells, atomShells, nao: ao, aoAtom, aoShell };
}

const normCache = new WeakMap();
function shellNorm(sh) {
  if (normCache.has(sh)) return normCache.get(sh);
  const blk = new Float64Array(25);
  shellPair(sh.l, sh.alpha, sh.coeff, sh.l, sh.alpha, sh.coeff, [0, 0, 0], { S: blk });
  const n = 1 / Math.sqrt(blk[0]);
  normCache.set(sh, n);
  return n;
}

const tab = [new Float64Array((MAXI + 1) * (MAXJ + 1)), new Float64Array((MAXI + 1) * (MAXJ + 1)), new Float64Array((MAXI + 1) * (MAXJ + 1))];
// blocchi cartesiani temporanei (6×6 al massimo, 3 componenti di dipolo, 6 di quadrupolo, derivate in 3 direzioni)
const cS = new Float64Array(36), cD = new Float64Array(3 * 36), cQ = new Float64Array(6 * 36);
const cdS = new Float64Array(3 * 36), cdD = new Float64Array(9 * 36), cdQ = new Float64Array(18 * 36);
const QIDX = [[0, 0], [0, 1], [1, 1], [0, 2], [1, 2], [2, 2]]; // xx, xy, yy, xz, yz, zz

// somma i contributi di una coppia di primitive: S = X0Y0Z0, D_c con un fattore di ordine 1 nella direzione c,
// Q_uw con i fattori di ordine 1 o 2; g seleziona il blocco (0 per gli integrali, 0–2 per le derivate)
function accumulate(w, X0, X1, X2, Y0, Y1, Y2, Z0, Z1, Z2, oS, oD, oQ, g, nc, k, wantD, wantQ) {
  oS[g * nc + k] += w * X0 * Y0 * Z0;
  if (wantD) {
    const b = 3 * g * nc + k;
    oD[b] += w * X1 * Y0 * Z0; oD[b + nc] += w * X0 * Y1 * Z0; oD[b + 2 * nc] += w * X0 * Y0 * Z1;
  }
  if (wantQ) {
    const b = 6 * g * nc + k;
    oQ[b] += w * X2 * Y0 * Z0; oQ[b + nc] += w * X1 * Y1 * Z0; oQ[b + 2 * nc] += w * X0 * Y2 * Z0;
    oQ[b + 3 * nc] += w * X1 * Y0 * Z1; oQ[b + 4 * nc] += w * X0 * Y1 * Z1; oQ[b + 5 * nc] += w * X0 * Y0 * Z2;
  }
}

/**
 * Integrali fra due shell contratte; R = B − A (vettore dal centro bra al centro ket).
 * out.S [na×nb], out.D [3][na×nb] e out.Q [6][na×nb] (operatori centrati su B), e, se richieste,
 * out.dS [3][…], out.dD [3·3][…], out.dQ [3·6][…]: derivate rispetto ad A (indice: direzione*ncomp + comp).
 * Valori non normalizzati: la normalizzazione delle shell è applicata da chi chiama.
 */
export function shellPair(la, alphaA, coeffA, lb, alphaB, coeffB, R, out) {
  const ca = CART[la], cb = CART[lb], na = ca.length, nb = cb.length;
  const wantD = !!out.D, wantQ = !!out.Q, wantG = !!out.dS;
  const ncart = na * nb;
  cS.fill(0, 0, ncart);
  if (wantD) cD.fill(0, 0, 3 * ncart);
  if (wantQ) cQ.fill(0, 0, 6 * ncart);
  if (wantG) { cdS.fill(0, 0, 3 * ncart); if (wantD) cdD.fill(0, 0, 9 * ncart); if (wantQ) cdQ.fill(0, 0, 18 * ncart); }
  const r2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
  const imax = la + (wantG ? 1 : 0), jmax = lb + (wantQ ? 2 : wantD ? 1 : 0);
  const J1 = MAXJ + 1;
  for (let p1 = 0; p1 < alphaA.length; p1++) {
    const a = alphaA[p1];
    for (let p2 = 0; p2 < alphaB.length; p2++) {
      const b = alphaB[p2], p = a + b, mu = a * b / p;
      if (mu * r2 > 40) continue;
      const pref = coeffA[p1] * coeffB[p2] * Math.exp(-mu * r2) * (Math.PI / p) ** 1.5;
      gaussMoments(p, imax + jmax, G);
      // P − A = (b/p) R, P − B = −(a/p) R
      for (let d = 0; d < 3; d++) table1d(b / p * R[d], -a / p * R[d], imax, jmax, tab[d]);
      const tx = tab[0], ty = tab[1], tz = tab[2];
      for (let i = 0; i < na; i++) {
        const ax = ca[i][0], ay = ca[i][1], az = ca[i][2];
        for (let j = 0; j < nb; j++) {
          const bx = cb[j][0], by = cb[j][1], bz = cb[j][2], k = i * nb + j;
          const ix = ax * J1 + bx, iy = ay * J1 + by, iz = az * J1 + bz;
          const X0 = tx[ix], X1 = tx[ix + 1], X2 = tx[ix + 2];
          const Y0 = ty[iy], Y1 = ty[iy + 1], Y2 = ty[iy + 2];
          const Z0 = tz[iz], Z1 = tz[iz + 1], Z2 = tz[iz + 2];
          accumulate(pref, X0, X1, X2, Y0, Y1, Y2, Z0, Z1, Z2, cS, cD, cQ, 0, ncart, k, wantD, wantQ);
          if (wantG) {
            // ∂/∂A_g di (x−A)^e e^{−a(x−A)²} = 2a (x−A)^{e+1} − e (x−A)^{e−1}
            const a2 = 2 * a;
            const ux = ix + J1, lx = ix - J1, uy = iy + J1, ly = iy - J1, uz = iz + J1, lz = iz - J1;
            const dX0 = a2 * tx[ux] - (ax ? ax * tx[lx] : 0), dX1 = a2 * tx[ux + 1] - (ax ? ax * tx[lx + 1] : 0), dX2 = a2 * tx[ux + 2] - (ax ? ax * tx[lx + 2] : 0);
            const dY0 = a2 * ty[uy] - (ay ? ay * ty[ly] : 0), dY1 = a2 * ty[uy + 1] - (ay ? ay * ty[ly + 1] : 0), dY2 = a2 * ty[uy + 2] - (ay ? ay * ty[ly + 2] : 0);
            const dZ0 = a2 * tz[uz] - (az ? az * tz[lz] : 0), dZ1 = a2 * tz[uz + 1] - (az ? az * tz[lz + 1] : 0), dZ2 = a2 * tz[uz + 2] - (az ? az * tz[lz + 2] : 0);
            accumulate(pref, dX0, dX1, dX2, Y0, Y1, Y2, Z0, Z1, Z2, cdS, cdD, cdQ, 0, ncart, k, wantD, wantQ);
            accumulate(pref, X0, X1, X2, dY0, dY1, dY2, Z0, Z1, Z2, cdS, cdD, cdQ, 1, ncart, k, wantD, wantQ);
            accumulate(pref, X0, X1, X2, Y0, Y1, Y2, dZ0, dZ1, dZ2, cdS, cdD, cdQ, 2, ncart, k, wantD, wantQ);
          }
        }
      }
    }
  }
  // trasformazione cartesiane → sferiche
  const Ta = TRAFO[la], Tb = TRAFO[lb], sa = Ta.length, sb = Tb.length;
  const tr = (src, off, dst) => {
    for (let i = 0; i < sa; i++) for (let j = 0; j < sb; j++) {
      let v = 0;
      for (let p = 0; p < na; p++) { const t = Ta[i][p]; if (!t) continue; for (let q = 0; q < nb; q++) { const u = Tb[j][q]; if (u) v += t * u * src[off + p * nb + q]; } }
      dst[i * sb + j] = v;
    }
  };
  tr(cS, 0, out.S);
  if (wantD) for (let c = 0; c < 3; c++) tr(cD, c * ncart, out.D[c]);
  if (wantQ) for (let c = 0; c < 6; c++) tr(cQ, c * ncart, out.Q[c]);
  if (wantG) {
    for (let g = 0; g < 3; g++) tr(cdS, g * ncart, out.dS[g]);
    if (wantD) for (let c = 0; c < 9; c++) tr(cdD, c * ncart, out.dD[c]);
    if (wantQ) for (let c = 0; c < 18; c++) tr(cdQ, c * ncart, out.dQ[c]);
  }
  return out;
}

const blk = (n) => new Float64Array(n);
/** Contenitore riutilizzabile per i blocchi di una coppia di shell (fino a d–d). */
export function pairBuffers(withGrad) {
  const o = { S: blk(25), D: [blk(25), blk(25), blk(25)], Q: Array.from({ length: 6 }, () => blk(25)) };
  if (withGrad) {
    o.dS = [blk(25), blk(25), blk(25)];
    o.dD = Array.from({ length: 9 }, () => blk(25));
    o.dQ = Array.from({ length: 18 }, () => blk(25));
  }
  return o;
}

/**
 * Matrici complete: S (nao²), D[3] e Q[6] con D[c][μ·n+ν] = ⟨μ|(r−R_ν)_c|ν⟩ e Q non a traccia nulla.
 * Le coppie di atomi più lontane di cutoff (bohr) sono trascurate.
 */
export function moleculeIntegrals(basis, pos, { cutoff = 40, withGrad = false } = {}) {
  const n = basis.nao, sh = basis.shells;
  const S = new Float64Array(n * n);
  const D = [0, 1, 2].map(() => new Float64Array(n * n));
  const Q = Array.from({ length: 6 }, () => new Float64Array(n * n));
  const buf = pairBuffers(withGrad);
  // derivate (rispetto all'atomo bra) delle coppie di shell su atomi diversi, già normalizzate, per il gradiente
  const derivs = withGrad ? [] : null;
  const c2 = cutoff * cutoff;
  for (let I = 0; I < sh.length; I++) {
    const si = sh[I], A = si.atom;
    for (let J = 0; J <= I; J++) {
      const sj = sh[J], B = sj.atom;
      const R = [pos[3 * B] - pos[3 * A], pos[3 * B + 1] - pos[3 * A + 1], pos[3 * B + 2] - pos[3 * A + 2]];
      if (R[0] ** 2 + R[1] ** 2 + R[2] ** 2 > c2) continue;
      const g = withGrad && A !== B;
      shellPair(si.l, si.alpha, si.coeff, sj.l, sj.alpha, sj.coeff, R, g ? buf : { S: buf.S, D: buf.D, Q: buf.Q });
      const nb = nSph(sj.l), nrm = si.norm * sj.norm, m = nSph(si.l) * nb;
      if (g) {
        const copy = (arr) => arr.map(b => { const o = new Float64Array(m); for (let k = 0; k < m; k++) o[k] = b[k] * nrm; return o; });
        derivs.push({ I, J, R, dS: copy(buf.dS), dD: copy(buf.dD), dQ: copy(buf.dQ) });
      }
      for (let i = 0; i < nSph(si.l); i++) for (let j = 0; j < nb; j++) {
        const mu = si.ao + i, nu = sj.ao + j, k = i * nb + j;
        const s = buf.S[k] * nrm;
        S[mu * n + nu] = S[nu * n + mu] = s;
        // ⟨μ|(r−B)|ν⟩ e, spostando il centro, ⟨ν|(r−A)|μ⟩ = ⟨μ|(r−B)|ν⟩ + R s
        const d = [buf.D[0][k] * nrm, buf.D[1][k] * nrm, buf.D[2][k] * nrm];
        for (let c = 0; c < 3; c++) { D[c][mu * n + nu] = d[c]; D[c][nu * n + mu] = d[c] + R[c] * s; }
        for (let c = 0; c < 6; c++) {
          const [u, w] = QIDX[c], q = buf.Q[c][k] * nrm;
          Q[c][mu * n + nu] = q;
          Q[c][nu * n + mu] = q + R[u] * d[w] + d[u] * R[w] + R[u] * R[w] * s;
        }
      }
    }
  }
  return { S, D, Q, derivs };
}

export { QIDX };

/**
 * Valori delle funzioni di base della shell sh (normalizzate) nel punto (dx, dy, dz) relativo all'atomo, in bohr.
 * Scrive 2l+1 valori in out a partire da off.
 */
export function shellValues(sh, dx, dy, dz, out, off) {
  const r2 = dx * dx + dy * dy + dz * dz;
  let rad = 0;
  for (let p = 0; p < sh.alpha.length; p++) rad += sh.coeff[p] * Math.exp(-sh.alpha[p] * r2);
  rad *= sh.norm;
  if (sh.l === 0) { out[off] = rad; return; }
  if (sh.l === 1) { out[off] = rad * dy; out[off + 1] = rad * dz; out[off + 2] = rad * dx; return; }
  out[off] = rad * S3 * dx * dy;
  out[off + 1] = rad * S3 * dy * dz;
  out[off + 2] = rad * (dz * dz - 0.5 * (dx * dx + dy * dy));
  out[off + 3] = rad * S3 * dx * dz;
  out[off + 4] = rad * S3 / 2 * (dx * dx - dy * dy);
}
