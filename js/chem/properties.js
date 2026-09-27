// Proprietà molecolari dalla funzione d'onda Hartree–Fock.

import { matmul, symFunction } from './linalg.js';
import { AU_TO_DEBYE, HARTREE_EV } from './hf.js';

/** Popolazioni di Mulliken e Löwdin, cariche parziali, ordini di legame di Mayer, valenze. */
export function populationAnalysis(res) {
  const { basis, S, Pa, Pb, atoms } = res;
  const n = basis.nbf;
  const nat = atoms.length;
  const Pt = new Float64Array(n * n);
  const Ps = new Float64Array(n * n); // densità di spin
  for (let i = 0; i < n * n; i++) { Pt[i] = Pa[i] + Pb[i]; Ps[i] = Pa[i] - Pb[i]; }
  const PS = matmul(Pt, S, n);
  const PaS = matmul(Pa, S, n);
  const PbS = matmul(Pb, S, n);
  const atomOf = basis.functions.map(f => f.atom);

  // Mulliken: q_A = Z_A − Σ_{μ∈A} (PS)_μμ
  const mulliken = atoms.map(a => a.Z);
  const spinMulliken = new Array(nat).fill(0);
  const PsS = matmul(Ps, S, n);
  for (let m = 0; m < n; m++) {
    mulliken[atomOf[m]] -= PS[m * n + m];
    spinMulliken[atomOf[m]] += PsS[m * n + m];
  }

  // Löwdin: q_A = Z_A − Σ_{μ∈A} (S^½ P S^½)_μμ
  const Sh = symFunction(S, n, s => Math.sqrt(Math.max(s, 0)));
  const L = matmul(matmul(Sh, Pt, n), Sh, n);
  const lowdin = atoms.map(a => a.Z);
  for (let m = 0; m < n; m++) lowdin[atomOf[m]] -= L[m * n + m];

  // Ordine di legame di Mayer (1983): B_AB = 2 Σ_{μ∈A, ν∈B} [(P^αS)_μν (P^αS)_νμ + (P^βS)_μν (P^βS)_νμ]
  const bondOrder = Array.from({ length: nat }, () => new Float64Array(nat));
  for (let m = 0; m < n; m++) {
    for (let v = 0; v < n; v++) {
      const A = atomOf[m], B = atomOf[v];
      if (A === B) continue;
      bondOrder[A][B] += 2 * (PaS[m * n + v] * PaS[v * n + m] + PbS[m * n + v] * PbS[v * n + m]);
    }
  }
  const valence = bondOrder.map(row => row.reduce((s, x) => s + x, 0));
  return { mulliken, lowdin, spin: spinMulliken, bondOrder, valence };
}

/** Momento di dipolo elettrico μ = Σ_A Z_A R_A − Σ_μν P_μν ⟨μ|r|ν⟩ (in unità atomiche e debye). */
export function dipoleMoment(res) {
  const { atoms, one, Pa, Pb, n } = res;
  const mu = [0, 0, 0];
  for (const a of atoms) for (let k = 0; k < 3; k++) mu[k] += a.Z * a.xyz[k];
  for (let k = 0; k < 3; k++) {
    const D = one.dipole[k];
    let s = 0;
    for (let i = 0; i < n * n; i++) s += (Pa[i] + Pb[i]) * D[i];
    mu[k] -= s;
  }
  const au = Math.hypot(mu[0], mu[1], mu[2]);
  return { vector: mu, au, debye: au * AU_TO_DEBYE };
}

/**
 * Energia di correlazione MP2 (Møller–Plesset al secondo ordine).
 * Shell chiuse:  E⁽²⁾ = Σ_ijab (ia|jb) [2(ia|jb) − (ib|ja)] / (ε_i + ε_j − ε_a − ε_b)
 * Shell aperte (UMP2): somma dei contributi αα, ββ (con scambio) e αβ.
 * Gli integrali sono trasformati dalla base atomica a quella molecolare (costo O(N⁵)).
 */
export function mp2Energy(res, { frozenCore = true } = {}) {
  const { n, eri, atoms } = res;
  let ncore = 0;
  if (frozenCore) for (const a of atoms) ncore += a.Z > 18 ? 9 : a.Z > 10 ? 5 : a.Z > 2 ? 1 : 0;
  ncore = Math.min(ncore, res.nbeta > 0 ? res.nbeta - 1 : 0, res.nalpha - 1);
  ncore = Math.max(ncore, 0);
  if (res.nel === 0) return { energy: 0, ncore: 0 };
  const full = fullERI(eri, n);
  const range = (a, b) => { const r = []; for (let k = a; k < b; k++) r.push(k); return r; };
  const oa = range(ncore, res.nalpha), va = range(res.nalpha, n);
  if (!res.unrestricted) {
    const X = transform(full, n, res.Ca, oa, va, res.Ca, oa, va);
    const no = oa.length, nv = va.length;
    const eps = res.epsA;
    let e2 = 0;
    const at = (i, a, j, b) => X[((i * nv + a) * no + j) * nv + b];
    for (let i = 0; i < no; i++) for (let j = 0; j < no; j++) for (let a = 0; a < nv; a++) for (let b = 0; b < nv; b++) {
      const x = at(i, a, j, b);
      e2 += x * (2 * x - at(i, b, j, a)) / (eps[oa[i]] + eps[oa[j]] - eps[va[a]] - eps[va[b]]);
    }
    return { energy: e2, ncore };
  }
  const ob = range(ncore, res.nbeta), vb = range(res.nbeta, n);
  let e2 = 0;
  const same = (C, eps, o, v) => {
    if (o.length < 2) return 0;
    const X = transform(full, n, C, o, v, C, o, v);
    const no = o.length, nv = v.length;
    const at = (i, a, j, b) => X[((i * nv + a) * no + j) * nv + b];
    let e = 0;
    for (let i = 0; i < no; i++) for (let j = i + 1; j < no; j++) for (let a = 0; a < nv; a++) for (let b = a + 1; b < nv; b++) {
      const x = at(i, a, j, b) - at(i, b, j, a);
      e += x * x / (eps[o[i]] + eps[o[j]] - eps[v[a]] - eps[v[b]]);
    }
    return e;
  };
  e2 += same(res.Ca, res.epsA, oa, va);
  e2 += same(res.Cb, res.epsB, ob, vb);
  if (ob.length) {
    const X = transform(full, n, res.Ca, oa, va, res.Cb, ob, vb);
    const noa = oa.length, nva = va.length, nob = ob.length, nvb = vb.length;
    for (let i = 0; i < noa; i++) for (let a = 0; a < nva; a++) for (let j = 0; j < nob; j++) for (let b = 0; b < nvb; b++) {
      const x = X[((i * nva + a) * nob + j) * nvb + b];
      e2 += x * x / (res.epsA[oa[i]] + res.epsB[ob[j]] - res.epsA[va[a]] - res.epsB[vb[b]]);
    }
  }
  return { energy: e2, ncore };
}

function fullERI(eri, n) {
  const full = new Float64Array(n * n * n * n);
  const { idx, val, count } = eri;
  const put = (i, j, k, l, v) => { full[((i * n + j) * n + k) * n + l] = v; };
  for (let q = 0; q < count; q++) {
    const i = idx[4 * q], j = idx[4 * q + 1], k = idx[4 * q + 2], l = idx[4 * q + 3];
    const v = val[q];
    put(i, j, k, l, v); put(j, i, k, l, v); put(i, j, l, k, v); put(j, i, l, k, v);
    put(k, l, i, j, v); put(l, k, i, j, v); put(k, l, j, i, v); put(l, k, j, i, v);
  }
  return full;
}

/** (ia|jb) = Σ C_μi C_νa C_λj C_σb (μν|λσ), in quattro trasformazioni successive. */
function transform(full, n, C1, occ1, vir1, C2, occ2, vir2) {
  const no1 = occ1.length, nv1 = vir1.length, no2 = occ2.length, nv2 = vir2.length;
  const n2 = n * n, n3 = n2 * n;
  const t1 = new Float64Array(no1 * n3);
  for (let ii = 0; ii < no1; ii++) {
    for (let mu = 0; mu < n; mu++) {
      const c = C1[mu * n + occ1[ii]];
      if (c === 0) continue;
      const src = mu * n3, dst = ii * n3;
      for (let r = 0; r < n3; r++) t1[dst + r] += c * full[src + r];
    }
  }
  const t2 = new Float64Array(no1 * nv1 * n2);
  for (let ii = 0; ii < no1; ii++) for (let aa = 0; aa < nv1; aa++) {
    const dst = (ii * nv1 + aa) * n2;
    for (let nu = 0; nu < n; nu++) {
      const c = C1[nu * n + vir1[aa]];
      if (c === 0) continue;
      const src = (ii * n + nu) * n2;
      for (let r = 0; r < n2; r++) t2[dst + r] += c * t1[src + r];
    }
  }
  const t3 = new Float64Array(no1 * nv1 * no2 * n);
  for (let ia = 0; ia < no1 * nv1; ia++) for (let jj = 0; jj < no2; jj++) {
    const dst = (ia * no2 + jj) * n;
    for (let la = 0; la < n; la++) {
      const c = C2[la * n + occ2[jj]];
      if (c === 0) continue;
      const src = ia * n2 + la * n;
      for (let sg = 0; sg < n; sg++) t3[dst + sg] += c * t2[src + sg];
    }
  }
  const out = new Float64Array(no1 * nv1 * no2 * nv2);
  for (let q = 0; q < no1 * nv1 * no2; q++) for (let bb = 0; bb < nv2; bb++) {
    let s = 0;
    const b = vir2[bb];
    for (let sg = 0; sg < n; sg++) s += C2[sg * n + b] * t3[q * n + sg];
    out[q * nv2 + bb] = s;
  }
  return out;
}

/** Riepilogo degli orbitali molecolari: energie, occupazioni, HOMO/LUMO, IE di Koopmans. */
export function orbitalSummary(res) {
  const { epsA, epsB, nalpha, nbeta, unrestricted, n } = res;
  const list = [];
  if (!unrestricted) {
    for (let k = 0; k < n; k++) list.push({ index: k, spin: null, e: epsA[k], occ: k < nalpha ? 2 : 0 });
  } else {
    for (let k = 0; k < n; k++) list.push({ index: k, spin: 'α', e: epsA[k], occ: k < nalpha ? 1 : 0 });
    for (let k = 0; k < n; k++) list.push({ index: k, spin: 'β', e: epsB[k], occ: k < nbeta ? 1 : 0 });
  }
  const homoE = Math.max(epsA[nalpha - 1], nbeta > 0 ? epsB[nbeta - 1] : -Infinity);
  const lumoE = Math.min(epsA[nalpha] ?? Infinity, epsB[nbeta] ?? Infinity);
  return {
    list,
    homo: homoE,
    lumo: lumoE,
    gapEV: (lumoE - homoE) * HARTREE_EV,
    koopmansIE: -homoE * HARTREE_EV,
  };
}

/**
 * Carattere di un orbitale molecolare: peso di ciascun atomo e tipo (s, p, d) secondo Mulliken,
 * e simmetria rispetto a un asse di legame (σ se non ha piani nodali contenenti l'asse, π altrimenti).
 */
export function orbitalComposition(res, k, spin = 'α') {
  const { basis, S, n } = res;
  const C = spin === 'β' ? res.Cb : res.Ca;
  const SC = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += S[i * n + j] * C[j * n + k];
    SC[i] = s;
  }
  const byAtom = new Array(res.atoms.length).fill(0);
  const byL = [0, 0, 0, 0];
  for (let i = 0; i < n; i++) {
    const w = C[i * n + k] * SC[i];
    byAtom[basis.functions[i].atom] += w;
    byL[basis.functions[i].l] += w;
  }
  return { byAtom, byL };
}
