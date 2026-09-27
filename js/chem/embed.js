// Geometria 3D iniziale da un grafo molecolare, con un campo di forze ispirato alla teoria VSEPR:
//  • legami: molle armoniche con lunghezza dai raggi covalenti di Pyykkö (dipende dall'ordine);
//  • domini elettronici (atomi legati e coppie solitarie) attorno a ogni atomo si respingono sulla
//    sfera unitaria, E = w Σ 1/|û_i − û_j|² (le coppie solitarie respingono di più: regola di Gillespie);
//  • legami doppi e tripli: i sostituenti restano complanari (termine in sin² del diedro);
//  • atomi non legati: repulsione a corto raggio.
// La geometria ottenuta è poi rifinita dall'ottimizzazione quantistica Hartree–Fock.

import { covalentRadius } from './elementData.js';
import { makeRandom } from '../physics/wavefunction.js';
import { eigh } from './linalg.js';

const ANG_TO_BOHR = 1 / 0.52917721090;

/**
 * @param mol { atoms: [{Z}], bonds: [{a,b,order}] }, lonePairs: array (coppie solitarie per atomo, anche frazionarie)
 * @returns coordinate in bohr [{Z, xyz}]
 */
export function embedMolecule(mol, lonePairs, { seeds = 6, iterations = 2500 } = {}) {
  const n = mol.atoms.length;
  if (n === 1) return [{ Z: mol.atoms[0].Z, xyz: [0, 0, 0] }];
  let best = null;
  for (let s = 0; s < seeds; s++) {
    const res = relax(mol, lonePairs, 1234 + 977 * s, iterations);
    if (!best || res.energy < best.energy) best = res;
  }
  // centra nell'origine
  const c = [0, 0, 0];
  for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) c[k] += best.x[3 * i + k] / n;
  return mol.atoms.map((a, i) => ({ Z: a.Z, xyz: [0, 1, 2].map(k => (best.x[3 * i + k] - c[k]) * ANG_TO_BOHR) }));
}

function relax(mol, lonePairs, seed, iterations) {
  const rand = makeRandom(seed);
  const { atoms, bonds } = mol;
  const n = atoms.length;
  // punti: atomi + coppie solitarie fittizie
  const lpOwner = [];
  atoms.forEach((a, i) => { for (let k = 0; k < Math.round(lonePairs[i] ?? 0); k++) lpOwner.push(i); });
  const np = n + lpOwner.length;
  const x = new Float64Array(3 * np);
  const init = distanceGeometry(atoms, bonds, lonePairs, rand);
  x.set(init);
  lpOwner.forEach((o, k) => { for (let c = 0; c < 3; c++) x[3 * (n + k) + c] = x[3 * o + c] + (rand() - 0.5) * 0.6; });

  const bondLength = (b) => {
    const r = (o) => covalentRadius(atoms[b.a].Z, o) + covalentRadius(atoms[b.b].Z, o);
    // legame aromatico: a metà tra singolo e doppio
    return (b.aromatic ? 0.5 * (r(1) + r(2)) : r(b.order)) / 100;
  };
  const springs = bonds.map(b => ({ i: b.a, j: b.b, r0: bondLength(b), k: 2 }));
  lpOwner.forEach((o, k) => springs.push({ i: o, j: n + k, r0: 0.5, k: 2 }));
  // domini attorno a ogni atomo
  const domains = Array.from({ length: n }, () => []);
  bonds.forEach(b => { domains[b.a].push({ p: b.b, lp: false }); domains[b.b].push({ p: b.a, lp: false }); });
  lpOwner.forEach((o, k) => domains[o].push({ p: n + k, lp: true }));
  const bonded = new Set();
  const key = (i, j) => (i < j ? i * 100000 + j : j * 100000 + i);
  bonds.forEach(b => bonded.add(key(b.a, b.b)));
  // coppie 1-3 (escluse dalla repulsione non legante)
  const near = new Set(bonded);
  domains.forEach(d => { for (const u of d) for (const v of d) if (u.p !== v.p && u.p < n && v.p < n) near.add(key(u.p, v.p)); });
  const vdw = (Z) => (Z === 1 ? 1.0 : 1.35);
  // diedri attorno ai legami multipli
  const planarity = [];
  bonds.forEach(b => {
    if (b.order < 2 && !b.aromatic) return;
    const na = domains[b.a].filter(d => !d.lp && d.p !== b.b).map(d => d.p);
    const nbb = domains[b.b].filter(d => !d.lp && d.p !== b.a).map(d => d.p);
    for (const i of na) for (const l of nbb) planarity.push([i, b.a, b.b, l]);
  });

  const grad = new Float64Array(3 * np);
  const energyAndGrad = () => {
    grad.fill(0);
    let E = 0;
    for (const s of springs) {
      const d = [0, 1, 2].map(c => x[3 * s.i + c] - x[3 * s.j + c]);
      const r = Math.hypot(...d) || 1e-6;
      const dr = r - s.r0;
      E += s.k * dr * dr;
      for (let c = 0; c < 3; c++) {
        const g = 2 * s.k * dr * d[c] / r;
        grad[3 * s.i + c] += g;
        grad[3 * s.j + c] -= g;
      }
    }
    // repulsione dei domini sulla sfera unitaria
    for (let c0 = 0; c0 < n; c0++) {
      const dom = domains[c0];
      if (dom.length < 2) continue;
      const u = dom.map(d => {
        const v = [0, 1, 2].map(c => x[3 * d.p + c] - x[3 * c0 + c]);
        const r = Math.hypot(...v) || 1e-6;
        return { v: v.map(q => q / r), r };
      });
      for (let a = 0; a < dom.length; a++) {
        for (let b = a + 1; b < dom.length; b++) {
          const w = 0.6 * (dom[a].lp ? 1.35 : 1) * (dom[b].lp ? 1.35 : 1);
          const diff = [0, 1, 2].map(c => u[a].v[c] - u[b].v[c]);
          const d2 = diff[0] ** 2 + diff[1] ** 2 + diff[2] ** 2 + 1e-4;
          E += w / d2;
          // dE/dû_a = −2w diff / d2²
          const gu = diff.map(q => -2 * w * q / (d2 * d2));
          for (const [idx, sign] of [[a, 1], [b, -1]]) {
            const U = u[idx];
            const g = gu.map(q => sign * q);
            const dotg = g[0] * U.v[0] + g[1] * U.v[1] + g[2] * U.v[2];
            for (let c = 0; c < 3; c++) {
              const gx = (g[c] - dotg * U.v[c]) / U.r;
              grad[3 * dom[idx].p + c] += gx;
              grad[3 * c0 + c] -= gx;
            }
          }
        }
      }
    }
    // repulsione tra atomi non legati
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (near.has(key(i, j))) continue;
        const d = [0, 1, 2].map(c => x[3 * i + c] - x[3 * j + c]);
        const r = Math.hypot(...d) || 1e-6;
        const dmin = vdw(atoms[i].Z) + vdw(atoms[j].Z);
        if (r >= dmin) continue;
        const dr = dmin - r;
        E += 0.5 * dr * dr;
        for (let c = 0; c < 3; c++) {
          const g = -dr * d[c] / r;
          grad[3 * i + c] += g;
          grad[3 * j + c] -= g;
        }
      }
    }
    return E;
  };
  // E = ½ sin²φ per ogni diedro attorno ai legami multipli; gradiente analitico (Blondel e Karplus, 1996)
  const addPlanarGrad = (weight) => {
    let E = 0;
    for (const [i, j, k, l] of planarity) {
      const P = (a) => [x[3 * a], x[3 * a + 1], x[3 * a + 2]];
      const ri = P(i), rj = P(j), rk = P(k), rl = P(l);
      const F = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
      const G = [rj[0] - rk[0], rj[1] - rk[1], rj[2] - rk[2]];
      const H = [rl[0] - rk[0], rl[1] - rk[1], rl[2] - rk[2]];
      const cr = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
      const dt = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
      const A = cr(F, G), B = cr(H, G);
      const A2 = dt(A, A), B2 = dt(B, B), Gn = Math.sqrt(dt(G, G));
      if (A2 < 1e-10 || B2 < 1e-10 || Gn < 1e-10) continue;
      const cosp = dt(A, B) / Math.sqrt(A2 * B2);
      const sinp = dt(cr(B, A), G) / (Math.sqrt(A2 * B2) * Gn);
      E += weight * 0.5 * sinp * sinp;
      const dEdphi = weight * sinp * cosp;
      const FG = dt(F, G), HG = dt(H, G);
      for (let c = 0; c < 3; c++) {
        const gi = -Gn / A2 * A[c];
        const gl = Gn / B2 * B[c];
        const gj = Gn / A2 * A[c] + FG / (A2 * Gn) * A[c] - HG / (B2 * Gn) * B[c];
        const gk = -Gn / B2 * B[c] - FG / (A2 * Gn) * A[c] + HG / (B2 * Gn) * B[c];
        grad[3 * i + c] += dEdphi * gi;
        grad[3 * j + c] += dEdphi * gj;
        grad[3 * k + c] += dEdphi * gk;
        grad[3 * l + c] += dEdphi * gl;
      }
    }
    return E;
  };

  // minimizzazione FIRE (Bitzek et al., 2006)
  const v = new Float64Array(3 * np);
  let dt = 0.02;
  let alpha = 0.1;
  let Npos = 0;
  let E = 0;
  for (let it = 0; it < iterations; it++) {
    // la planarità dei legami multipli entra gradualmente nella seconda metà della minimizzazione
    E = energyAndGrad() + addPlanarGrad(Math.max(0, 2 * it / iterations - 1));
    let P = 0, vn = 0, fn = 0;
    for (let i = 0; i < 3 * np; i++) { P += -grad[i] * v[i]; vn += v[i] * v[i]; fn += grad[i] * grad[i]; }
    vn = Math.sqrt(vn); fn = Math.sqrt(fn) || 1e-12;
    if (fn < 1e-6) break;
    for (let i = 0; i < 3 * np; i++) v[i] = (1 - alpha) * v[i] + alpha * (-grad[i]) / fn * vn;
    if (P >= 0) {
      Npos++;
      if (Npos > 5) { dt = Math.min(dt * 1.1, 0.1); alpha *= 0.99; }
    } else {
      Npos = 0;
      dt = Math.max(dt * 0.5, 2e-3);
      alpha = 0.1;
      v.fill(0);
    }
    // spostamento massimo per passo: 0,1 Å
    let maxStep = 0;
    for (let i = 0; i < 3 * np; i++) {
      v[i] += -grad[i] * dt;
      maxStep = Math.max(maxStep, Math.abs(v[i] * dt));
    }
    const scale = maxStep > 0.1 ? 0.1 / maxStep : 1;
    for (let i = 0; i < 3 * np; i++) x[i] += v[i] * dt * scale;
  }
  return { x: x.subarray(0, 3 * n), energy: E };
}

/**
 * Geometria di partenza con la geometria delle distanze (Crippen e Havel):
 * limiti inferiori/superiori per ogni coppia di atomi (1-2 dai raggi covalenti, 1-3 dagli angoli VSEPR,
 * altri dalla repulsione di van der Waals), lisciati con la disuguaglianza triangolare; si sceglie a caso
 * una distanza per coppia e si ricavano le coordinate dagli autovettori della matrice metrica.
 */
function distanceGeometry(atoms, bonds, lonePairs, rand) {
  const n = atoms.length;
  const L = Array.from({ length: n }, () => new Float64Array(n));
  const U = Array.from({ length: n }, () => new Float64Array(n).fill(100));
  const nb = Array.from({ length: n }, () => []);
  const len = (b) => {
    const r = (o) => covalentRadius(atoms[b.a].Z, o) + covalentRadius(atoms[b.b].Z, o);
    return (b.aromatic ? 0.5 * (r(1) + r(2)) : r(b.order)) / 100;
  };
  const d12 = Array.from({ length: n }, () => new Float64Array(n));
  bonds.forEach(b => { nb[b.a].push(b.b); nb[b.b].push(b.a); d12[b.a][b.b] = d12[b.b][b.a] = len(b); });
  for (let i = 0; i < n; i++) {
    U[i][i] = 0;
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      const vd = (Z) => (Z === 1 ? 1.0 : 1.35);
      L[i][j] = Math.min(2.2, 0.85 * (vd(atoms[i].Z) + vd(atoms[j].Z)));
    }
  }
  // 1-3: angolo ideale dal numero sterico
  for (let c = 0; c < n; c++) {
    const sn = nb[c].length + Math.round(lonePairs[c] ?? 0);
    const ideal = { 2: 180, 3: 120, 4: 109.47, 5: 105, 6: 90 }[sn] ?? 109.47;
    for (const a of nb[c]) {
      for (const b of nb[c]) {
        if (a >= b) continue;
        const da = d12[a][c], db = d12[b][c];
        const tri = (deg) => Math.sqrt(da * da + db * db - 2 * da * db * Math.cos(deg * Math.PI / 180));
        const lo = sn >= 5 ? tri(85) : tri(ideal - 8);
        const hi = sn >= 5 ? tri(180) : tri(Math.min(180, ideal + 8));
        L[a][b] = L[b][a] = Math.max(lo, 0.5);
        U[a][b] = U[b][a] = hi;
      }
    }
  }
  bonds.forEach(b => {
    const d = d12[b.a][b.b];
    L[b.a][b.b] = L[b.b][b.a] = d * 0.98;
    U[b.a][b.b] = U[b.b][b.a] = d * 1.02;
  });
  // lisciatura triangolare (Floyd–Warshall) dei limiti superiori, poi degli inferiori
  for (let k = 0; k < n; k++) for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    if (U[i][j] > U[i][k] + U[k][j]) U[i][j] = U[i][k] + U[k][j];
  }
  for (let k = 0; k < n; k++) for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    if (L[i][j] < L[i][k] - U[k][j]) L[i][j] = L[i][k] - U[k][j];
    if (L[i][j] < L[j][k] - U[k][i]) L[i][j] = L[j][k] - U[k][i];
  }
  // distanze casuali entro i limiti
  const D2 = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const lo = Math.min(L[i][j], U[i][j]);
    const d = lo + rand() * (U[i][j] - lo);
    D2[i][j] = D2[j][i] = d * d;
  }
  // matrice metrica rispetto al baricentro: G_ij = (d²_i0 + d²_j0 − d²_ij)/2
  const d0 = new Float64Array(n);
  let sumAll = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) sumAll += D2[i][j];
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += D2[i][j];
    d0[i] = s / n - sumAll / (2 * n * n);
  }
  const G = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) G[i * n + j] = 0.5 * (d0[i] + d0[j] - D2[i][j]);
  const { values, vectors } = eigh(G, n);
  const out = new Float64Array(3 * n);
  for (let c = 0; c < 3; c++) {
    const k = n - 1 - c;
    const lam = Math.sqrt(Math.max(values[k] ?? 0, 1e-4));
    for (let i = 0; i < n; i++) out[3 * i + c] = lam * (vectors[i * n + k] ?? 0) + (rand() - 0.5) * 0.05;
  }
  return out;
}
