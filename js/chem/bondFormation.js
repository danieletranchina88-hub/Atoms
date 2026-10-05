// Come nasce un legame: due atomi che si avvicinano, risolti con la meccanica quantistica a ogni distanza.
//
// Per ogni distanza R si calcolano:
//  · Hartree–Fock ristretto (RHF, Roothaan 1951) e non ristretto (UHF, Pople–Nesbet 1954). La soluzione UHF parte
//    dagli atomi separati con gli spin degli elettroni spaiati opposti sui due atomi e viene seguita verso distanze
//    minori: dove coincide con RHF si trova il punto di Coulson–Fischer (Phil. Mag. 40, 386, 1949).
//  · per H₂, l'interazione di configurazioni completa (FCI): la soluzione esatta dell'equazione di Schrödinger nella
//    base, che dissocia correttamente in due atomi H (Szabo e Ostlund, "Modern Quantum Chemistry", cap. 4).
//  · GFN2-xTB (Bannwarth, Ehlert, Grimme, J. Chem. Theory Comput. 15, 1652, 2019).
// Dalla funzione d'onda: orbitali molecolari scomposti negli orbitali degli atomi liberi (1s, 2s, 2pz…), carattere
// legante o antilegante (popolazione di sovrapposizione di Mulliken), ordine di legame di Mayer, energia cinetica
// ⟨T⟩ e potenziale ⟨V⟩ = E − ⟨T⟩ (analisi di Ruedenberg, Rev. Mod. Phys. 34, 326, 1962), teorema del viriale
// molecolare 2⟨T⟩ + ⟨V⟩ = −R dE/dR (Slater, J. Chem. Phys. 1, 687, 1933).

import { runHF, ANGSTROM_TO_BOHR, HARTREE_EV } from './hf.js';
import { buildBasis, oneElectron, twoElectron, buildJK } from './integrals.js';
import { eigh, symFunction } from './linalg.js';
import { GFN2xTB } from './xtb/gfn2.js';
import { experimentalWell, reducedMass } from './diatomicData.js';

const BOHR_A = 1 / ANGSTROM_TO_BOHR;
const AMU_ME = 1822.888486209;            // unità di massa atomica in masse dell'elettrone
const HARTREE_CM = 219474.6313632;
const SYMBOL = { 1: 'H', 2: 'He', 3: 'Li', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 17: 'Cl' };

/** Atomi A e B sull'asse z, con il punto medio nell'origine (bohr). */
export const geometry = (Z, Rb) => [{ Z: Z[0], xyz: [0, 0, -Rb / 2] }, { Z: Z[1], xyz: [0, 0, Rb / 2] }];

/** Distanze della scansione (Å), dalla più grande alla più piccola: fitte vicino all'equilibrio, rade lontano. */
export function scanGrid(mol, { near = 22, far = 14 } = {}) {
  const re = mol.exp.re;
  const lo = mol.scanMin ?? 0.6 * re, mid = Math.max(1.6 * re, lo + 0.8), hi = Math.max(re + 3.4, 2.4 * re);
  const R = [];
  for (let k = 0; k < near; k++) R.push(lo + (mid - lo) * k / near);
  for (let k = 0; k <= far; k++) R.push(mid * Math.pow(hi / mid, k / far));
  return R.reverse();
}

/** Metodi disponibili per la molecola. */
export function methodsFor(mol) {
  const nel = mol.Z[0] + mol.Z[1] - mol.charge;
  if (nel === 1) return ['exact', 'gfn2'];
  const out = [];
  if (mol.mult === 1) out.push('rhf');
  if (mol.id !== 'He2') out.push('uhf');
  if (nel === 2 && mol.mult === 1) out.push('fci');
  out.push('gfn2');
  return out;
}

export const METHOD_INFO = {
  exact: { name: 'esatto nella base', short: 'esatto', ref: 'un solo elettrone: Hartree–Fock coincide con la soluzione esatta nella base' },
  rhf: { name: 'Hartree–Fock ristretto (RHF)', short: 'RHF', ref: 'C. C. J. Roothaan, Rev. Mod. Phys. 23, 69 (1951)' },
  uhf: { name: 'Hartree–Fock non ristretto (UHF)', short: 'UHF', ref: 'J. A. Pople, R. K. Nesbet, J. Chem. Phys. 22, 571 (1954)' },
  fci: { name: 'interazione di configurazioni completa (FCI)', short: 'FCI', ref: 'soluzione esatta nella base; Szabo e Ostlund, Modern Quantum Chemistry, cap. 4' },
  gfn2: { name: 'GFN2-xTB', short: 'GFN2-xTB', ref: 'C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 15, 1652 (2019)' },
};

// ---------------------------------------------------------------------------
// Atomi liberi: energia, densità media sferica e orbitali atomici nella stessa base
// ---------------------------------------------------------------------------

const parityClass = (c) => `${c[0] & 1}${c[1] & 1}${c[2] & 1}`;
const AXIS_OF = { 100: 'x', '010': 'y', '001': 'z' };

/** Media sferica (cubica) della densità di un atomo: p e d diventano isotropi, i blocchi misti s–p si annullano. */
export function sphericalAverage(P, basis) {
  const n = basis.nbf, out = new Float64Array(n * n), sh = basis.shells;
  for (const a of sh) for (const b of sh) {
    const na = a.comps.length, nb = b.comps.length;
    const get = (i, j) => P[(a.offset + i) * n + b.offset + j];
    const set = (i, j, v) => { out[(a.offset + i) * n + b.offset + j] = v; };
    if (a.l === 0 && b.l === 0) { set(0, 0, get(0, 0)); continue; }
    if (a.l === 1 && b.l === 1) { const v = (get(0, 0) + get(1, 1) + get(2, 2)) / 3; for (let k = 0; k < 3; k++) set(k, k, v); continue; }
    const diagIdx = (s) => s.comps.map((c, i) => (c[0] === 2 || c[1] === 2 || c[2] === 2 ? i : -1)).filter(i => i >= 0);
    const offIdx = (s) => s.comps.map((c, i) => (c[0] < 2 && c[1] < 2 && c[2] < 2 ? i : -1)).filter(i => i >= 0);
    if (a.l === 2 && b.l === 2) {
      const da = diagIdx(a), db = diagIdx(b), oa = offIdx(a), ob = offIdx(b);
      let dd = 0, dx = 0, oo = 0;
      for (let k = 0; k < 3; k++) { dd += get(da[k], db[k]); oo += get(oa[k], ob[k]); for (let m = 0; m < 3; m++) if (m !== k) dx += get(da[k], db[m]); }
      for (let k = 0; k < 3; k++) { set(oa[k], ob[k], oo / 3); for (let m = 0; m < 3; m++) set(da[k], db[m], k === m ? dd / 3 : dx / 6); }
      continue;
    }
    if (a.l + b.l === 2 && (a.l === 0 || b.l === 0)) {
      // s con le componenti xx, yy, zz
      if (a.l === 0) { const db = diagIdx(b); const v = db.reduce((s, k) => s + get(0, k), 0) / 3; for (const k of db) set(0, k, v); }
      else { const da = diagIdx(a); const v = da.reduce((s, k) => s + get(k, 0), 0) / 3; for (const k of da) set(k, 0, v); }
    }
  }
  return out;
}

/** Autovettori del problema generalizzato F v = ε S v ristretto agli indici idx. */
function subEig(F, S, n, idx) {
  const m = idx.length, Fs = new Float64Array(m * m), Ss = new Float64Array(m * m);
  for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) { Fs[i * m + j] = F[idx[i] * n + idx[j]]; Ss[i * m + j] = S[idx[i] * n + idx[j]]; }
  const X = symFunction(Ss, m, (s) => (s > 1e-10 ? 1 / Math.sqrt(s) : 0));
  const XF = new Float64Array(m * m), Fp = new Float64Array(m * m);
  for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) { let v = 0; for (let k = 0; k < m; k++) v += X[i * m + k] * Fs[k * m + j]; XF[i * m + j] = v; }
  for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) { let v = 0; for (let k = 0; k < m; k++) v += XF[i * m + k] * X[k * m + j]; Fp[i * m + j] = v; }
  const { values, vectors } = eigh(Fp, m);
  const out = [];
  for (let k = 0; k < m; k++) {
    const v = new Float64Array(m);
    for (let i = 0; i < m; i++) { let s = 0; for (let j = 0; j < m; j++) s += X[i * m + j] * vectors[j * m + k]; v[i] = s; }
    out.push({ e: values[k], v });
  }
  return out;
}

/**
 * Atomo libero (UHF ad alto spin) nella base indicata. Restituisce l'energia, ⟨T⟩, la densità media sferica e gli
 * orbitali atomici (dalla Fock α costruita con le densità di spin sferiche, diagonalizzata per classi di simmetria) con nome
 * (1s, 2s, 2pz…), energia e occupazione.
 */
export function freeAtom(Z, { basis = '6-31G**', mult = 1, charge = 0 } = {}) {
  const atoms = [{ Z, xyz: [0, 0, 0] }];
  const B = buildBasis(atoms, basis);
  const n = B.nbf;
  const one = oneElectron(B, atoms);
  const eri = twoElectron(B);
  const nel = Z - charge;
  let E = 0, T = 0, Pa = new Float64Array(n * n), Pb = new Float64Array(n * n);
  if (nel > 0) {
    const r = runHF(atoms, { basis, basisObj: B, one, eri, multiplicity: mult, charge, unrestricted: true, conv: 1e-11 });
    E = r.energy; Pa = r.Pa; Pb = r.Pb;
    for (let i = 0; i < n * n; i++) T += (Pa[i] + Pb[i]) * one.T[i];
  }
  const Ptot = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) Ptot[i] = Pa[i] + Pb[i];
  const Psph = sphericalAverage(Ptot, B);
  // Fock degli elettroni α (spin maggioritario) con le densità di spin mediate sfericamente: F = H + J[Pα + Pβ] − K[Pα].
  // Per l'orbitale occupato più alto ε ≈ −energia di ionizzazione (Koopmans); per l'idrogeno ε(1s) = −½ Eh.
  const { J, K } = buildJK(eri, [sphericalAverage(Pa, B), sphericalAverage(Pb, B)], n);
  const F = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) F[i] = one.T[i] + one.V[i] + J[i] - K[0][i];
  // classi di parità (x, y, z): la Fock sferica non le mescola
  const classes = {};
  B.functions.forEach((f, mu) => { const k = parityClass(f.comp); (classes[k] ??= []).push(mu); });
  const orbitals = [];
  const count = {};
  const S = one.S;
  for (const [cls, idx] of Object.entries(classes)) {
    for (const { e, v } of subEig(F, S, n, idx)) {
      const full = new Float64Array(n);
      idx.forEach((mu, i) => { full[mu] = v[i]; });
      // carattere s/p/d (popolazioni di Mulliken) e occupazione nella densità sferica
      const Sv = new Float64Array(n);
      for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += S[i * n + j] * full[j]; Sv[i] = s; }
      const wl = [0, 0, 0];
      for (let mu = 0; mu < n; mu++) wl[Math.min(B.functions[mu].l, 2)] += full[mu] * Sv[mu];
      const l = wl.indexOf(Math.max(...wl));
      let occ = 0;
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) occ += Sv[i] * Psph[i * n + j] * Sv[j];
      orbitals.push({ cls, l, e, v: full, occ, axis: l === 1 ? AXIS_OF[cls] : '' });
    }
  }
  orbitals.sort((a, b) => a.e - b.e);
  for (const o of orbitals) {
    const key = `${o.l}${o.axis}`;
    count[key] = (count[key] ?? 0) + 1;
    const nq = o.l + count[key];
    o.n = nq;
    o.label = o.l === 2 ? 'd' : `${nq}${'spd'[o.l]}${o.axis}`;
    // gli orbitali oltre la valenza, nelle basi estese, servono solo a dare flessibilità
    o.extra = o.l === 2 || (o.occ < 0.02 && nq > valenceN(Z));
    o.core = !o.extra && Z > 2 && nq < valenceN(Z);
  }
  return { Z, symbol: SYMBOL[Z] ?? `Z${Z}`, mult, charge, basis: B, n, E, T, Pa, Pb, Psph, orbitals, nel };
}

function valenceN(Z) { return Z <= 2 ? 1 : Z <= 10 ? 2 : Z <= 18 ? 3 : 4; }

// ---------------------------------------------------------------------------
// Analisi di una funzione d'onda molecolare
// ---------------------------------------------------------------------------

/** Contesto condiviso dalla scansione: base, atomi liberi, orbitali atomici disposti nella base molecolare. */
export function makeContext(mol, basisName = '6-31G**') {
  const charges = mol.atomCharges ?? [0, 0];
  const A = freeAtom(mol.Z[0], { basis: basisName, mult: mol.atoms[0], charge: charges[0] });
  const B = mol.Z[0] === mol.Z[1] && mol.atoms[0] === mol.atoms[1] && charges[0] === charges[1] ? A : freeAtom(mol.Z[1], { basis: basisName, mult: mol.atoms[1], charge: charges[1] });
  const n = A.n + B.n;
  // orbitali atomici nella base molecolare (A occupa le prime A.n funzioni)
  const aos = [];
  [A, B].forEach((at, k) => {
    const off = k ? A.n : 0;
    for (const o of at.orbitals) {
      const v = new Float64Array(n);
      v.set(o.v, off);
      aos.push({ atom: k, label: o.label, l: o.l, axis: o.axis, e: o.e, occ: o.occ, extra: o.extra, v, off, nloc: at.n, vloc: o.v });
    }
  });
  const Eatoms = A.E + B.E, Tatoms = A.T + B.T;
  // densità dei due atomi liberi e sferici (promolecola), nella base molecolare
  const Ppro = new Float64Array(n * n);
  [A, B].forEach((at, k) => {
    const off = k ? A.n : 0;
    for (let i = 0; i < at.n; i++) for (let j = 0; j < at.n; j++) Ppro[(off + i) * n + off + j] = at.Psph[i * at.n + j];
  });
  return { mol, basis: basisName, A, B, n, aos, Eatoms, Tatoms, Ppro, homonuclear: mol.Z[0] === mol.Z[1] };
}

const matVec = (M, v, n) => { const o = new Float64Array(n); for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += M[i * n + j] * v[j]; o[i] = s; } return o; };
const column = (C, n, k) => { const v = new Float64Array(n); for (let i = 0; i < n; i++) v[i] = C[i * n + k]; return v; };

/**
 * Descrive gli orbitali molecolari: simmetria (σ/π/δ, g/u), carattere legante (popolazione di sovrapposizione
 * fra A e B), composizione negli orbitali atomici liberi. C: coefficienti per colonne; occ: occupazioni.
 */
export function describeOrbitals(ctx, basis, S, C, eps, occ, { maxVirtual = 4 } = {}) {
  const n = ctx.n, nA = ctx.A.n;
  const lastOcc = occ.reduce((m, o, k) => (o > 1e-6 ? k : m), -1);
  const upto = Math.min(n, lastOcc + 1 + maxVirtual);
  const out = [];
  const counters = {};
  for (let k = 0; k < upto; k++) {
    const c = column(C, n, k), Sc = matVec(S, c, n);
    // simmetria rispetto all'asse z
    const w = { s: 0, p: 0, d: 0 };
    let xw = 0, yw = 0;
    basis.functions.forEach((f, mu) => {
      const q = c[mu] * Sc[mu];
      const ox = f.comp[0] & 1, oy = f.comp[1] & 1;
      if (ox && oy) w.d += q; else if (ox || oy) { w.p += q; if (ox) xw += q; else yw += q; } else w.s += q;
    });
    const type = w.p >= w.s && w.p >= w.d ? 'π' : w.d > w.s ? 'δ' : 'σ';
    let gu = '';
    if (ctx.homonuclear) {
      // ⟨φ|î φ⟩ (î = inversione): ±1 per un orbitale g/u; nelle soluzioni UHF a simmetria rotta gli orbitali si
      // localizzano su un atomo e la parità perde significato
      let par = 0;
      for (let mu = 0; mu < nA; mu++) {
        const f = basis.functions[mu];
        par += 2 * c[mu] * Sc[mu + nA] * ((f.l & 1) ? -1 : 1);
      }
      gu = par > 0.8 ? 'g' : par < -0.8 ? 'u' : '';
    }
    const sym = type + gu;
    const key = sym + (type === 'π' ? (xw >= yw ? 'x' : 'y') : '');
    counters[key] = (counters[key] ?? 0) + 1;
    // popolazione di sovrapposizione 2 Σ_{μ∈A, ν∈B} c_μ c_ν S_μν
    let op = 0;
    for (let mu = 0; mu < nA; mu++) for (let nu = nA; nu < n; nu++) op += 2 * c[mu] * c[nu] * S[mu * n + nu];
    // composizione negli orbitali atomici: d_a = u_aᵀ S_AA c_A, peso d_a (Uᵀ S c)_a
    const comp = [];
    for (const ao of ctx.aos) {
      let da = 0, ga = 0;
      for (let i = 0; i < ao.nloc; i++) {
        const mu = ao.off + i;
        if (ao.vloc[i] === 0) continue;
        let sAA = 0;
        for (let j = 0; j < ao.nloc; j++) sAA += S[mu * n + ao.off + j] * c[ao.off + j];
        da += ao.vloc[i] * sAA;
        ga += ao.vloc[i] * Sc[mu];
      }
      const wgt = da * ga;
      if (Math.abs(wgt) > 0.005) comp.push({ atom: ao.atom, label: ao.label, extra: ao.extra, coef: da, w: wgt });
    }
    comp.sort((a, b) => Math.abs(b.w) - Math.abs(a.w));
    const character = op > 0.03 ? 'legante' : op < -0.03 ? 'antilegante' : 'non legante';
    out.push({ index: k, e: eps[k], occ: occ[k], sym, type, gu, label: `${counters[key]}${sym}`, star: op < -0.03, op, character, comp, axis: type === 'π' ? (xw >= yw ? 'x' : 'y') : '' });
  }
  return out;
}

/** Ordine di legame di Mayer (Chem. Phys. Lett. 97, 270, 1983) per densità α e β. */
export function mayerOrder(Pa, Pb, S, n, nA) {
  const PS = (P) => { const o = new Float64Array(n * n); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { let s = 0; for (let k = 0; k < n; k++) s += P[i * n + k] * S[k * n + j]; o[i * n + j] = s; } return o; };
  const a = PS(Pa), b = PS(Pb);
  let B = 0;
  for (let mu = 0; mu < nA; mu++) for (let nu = nA; nu < n; nu++) B += 2 * (a[mu * n + nu] * a[nu * n + mu] + b[mu * n + nu] * b[nu * n + mu]);
  return B;
}

/** Cariche di Mulliken degli atomi. */
function mullikenCharges(P, S, n, nA, Z) {
  let qA = 0;
  for (let mu = 0; mu < nA; mu++) for (let nu = 0; nu < n; nu++) qA += P[mu * n + nu] * S[nu * n + mu];
  let tot = 0;
  for (let i = 0; i < n * n; i++) tot += P[i] * S[i];
  return [Z[0] - qA, Z[1] - (tot - qA)];
}

// ---------------------------------------------------------------------------
// FCI per due elettroni (singoletto)
// ---------------------------------------------------------------------------

/**
 * Interazione di configurazioni completa per due elettroni in singoletto, negli orbitali RHF C (colonne).
 * Ψ = Σ_pq c_pq φ_p(1) φ_q(2) con c simmetrica; (Hc)_pq = Σ_r h_pr c_rq + Σ_s c_ps h_sq + Σ_rs (pr|qs) c_rs.
 * Restituisce energia, densità AO (spin sommato), occupazioni naturali e pesi delle configurazioni principali.
 */
export function fci2(hf) {
  const { n, C, H, eri, Enuc } = hf;
  // integrali sugli orbitali molecolari
  const h = new Float64Array(n * n);
  const HC = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let q = 0; q < n; q++) { let s = 0; for (let j = 0; j < n; j++) s += H[i * n + j] * C[j * n + q]; HC[i * n + q] = s; }
  for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) { let s = 0; for (let i = 0; i < n; i++) s += C[i * n + p] * HC[i * n + q]; h[p * n + q] = s; }
  const g = new Float64Array(n * n * n * n);   // (pq|rs)
  const D = new Float64Array(n * n);
  for (let r = 0; r < n; r++) for (let s = 0; s <= r; s++) {
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) D[i * n + j] = 0.5 * (C[i * n + r] * C[j * n + s] + C[i * n + s] * C[j * n + r]);
    const { J } = buildJK(eri, [D], n);
    const JC = new Float64Array(n * n);
    for (let i = 0; i < n; i++) for (let q = 0; q < n; q++) { let v = 0; for (let j = 0; j < n; j++) v += J[i * n + j] * C[j * n + q]; JC[i * n + q] = v; }
    for (let p = 0; p < n; p++) for (let q = 0; q <= p; q++) {
      let v = 0;
      for (let i = 0; i < n; i++) v += C[i * n + p] * JC[i * n + q];
      for (const [a, b, c2, d] of [[p, q, r, s], [q, p, r, s], [p, q, s, r], [q, p, s, r], [r, s, p, q], [s, r, p, q], [r, s, q, p], [s, r, q, p]]) g[((a * n + b) * n + c2) * n + d] = v;
    }
  }
  // base simmetrica: coppie p ≥ q
  const pairs = [];
  for (let p = 0; p < n; p++) for (let q = 0; q <= p; q++) pairs.push([p, q]);
  const m = pairs.length;
  const Hm = new Float64Array(m * m);
  const vec = (p, q) => { const B = new Float64Array(n * n); if (p === q) B[p * n + p] = 1; else { B[p * n + q] = B[q * n + p] = Math.SQRT1_2; } return B; };
  const apply = (c) => {
    const o = new Float64Array(n * n);
    for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) {
      let v = 0;
      for (let r = 0; r < n; r++) v += h[p * n + r] * c[r * n + q] + c[p * n + r] * h[r * n + q];
      for (let r = 0; r < n; r++) for (let s = 0; s < n; s++) { const crs = c[r * n + s]; if (crs) v += g[((p * n + r) * n + q) * n + s] * crs; }
      o[p * n + q] = v;
    }
    return o;
  };
  for (let b = 0; b < m; b++) {
    const Hb = apply(vec(...pairs[b]));
    for (let a = 0; a <= b; a++) {
      const [p, q] = pairs[a];
      const v = p === q ? Hb[p * n + p] : Math.SQRT1_2 * (Hb[p * n + q] + Hb[q * n + p]);
      Hm[a * m + b] = Hm[b * m + a] = v;
    }
  }
  const { values, vectors } = eigh(Hm, m);
  let k0 = 0;
  for (let k = 1; k < m; k++) if (values[k] < values[k0]) k0 = k;
  const c = new Float64Array(n * n);
  for (let a = 0; a < m; a++) {
    const [p, q] = pairs[a], x = vectors[a * m + k0];
    if (p === q) c[p * n + p] = x; else { c[p * n + q] += x * Math.SQRT1_2; c[q * n + p] += x * Math.SQRT1_2; }
  }
  // matrice densità ridotta a un corpo (spin sommato) γ = 2 c c
  const gam = new Float64Array(n * n);
  for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) { let s = 0; for (let r = 0; r < n; r++) s += c[p * n + r] * c[r * n + q]; gam[p * n + q] = 2 * s; }
  const nat = eigh(gam, n);
  const natOcc = Array.from(nat.values).sort((a, b) => b - a);
  // densità AO P = C γ Cᵀ
  const P = new Float64Array(n * n);
  const Cg = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let q = 0; q < n; q++) { let s = 0; for (let p = 0; p < n; p++) s += C[i * n + p] * gam[p * n + q]; Cg[i * n + q] = s; }
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { let s = 0; for (let q = 0; q < n; q++) s += Cg[i * n + q] * C[j * n + q]; P[i * n + j] = s; }
  // configurazioni principali (coppie di orbitali RHF), peso = |c_pq|² + |c_qp|²
  const configs = pairs.map(([p, q]) => ({ p, q, w: p === q ? c[p * n + p] ** 2 : 2 * c[p * n + q] ** 2, sign: Math.sign(p === q ? c[p * n + p] : c[p * n + q]) }))
    .sort((a, b) => b.w - a.w).slice(0, 4);
  // occupazioni degli orbitali RHF (diagonale di γ)
  const moOcc = Array.from({ length: n }, (_, p) => gam[p * n + p]);
  return { energy: values[k0] + Enuc, P, natOcc, configs, moOcc, dim: m };
}

// ---------------------------------------------------------------------------
// Scansione
// ---------------------------------------------------------------------------

/** Densità di stima per UHF a simmetria rotta: atomi liberi con gli spin dell'atomo B scambiati. */
function brokenSymmetryGuess(ctx) {
  const { A, B, n } = ctx;
  const Pa = new Float64Array(n * n), Pb = new Float64Array(n * n);
  const put = (dst, src, off, m) => { for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) dst[(off + i) * n + off + j] = src[i * m + j]; };
  put(Pa, A.Pa, 0, A.n); put(Pb, A.Pb, 0, A.n);
  put(Pa, B.Pb, A.n, B.n); put(Pb, B.Pa, A.n, B.n);
  return { Pa, Pb };
}

function pointFromHF(ctx, mol, Rang, r, kind) {
  const n = r.n, nA = ctx.A.n, S = r.S;
  let T = 0;
  const Pt = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) { Pt[i] = r.Pa[i] + r.Pb[i]; T += Pt[i] * r.one.T[i]; }
  const occA = Array.from({ length: n }, (_, k) => (k < r.nalpha ? 1 : 0));
  const occB = Array.from({ length: n }, (_, k) => (k < r.nbeta ? 1 : 0));
  const restricted = !r.unrestricted;
  const mos = restricted
    ? describeOrbitals(ctx, r.basis, S, r.Ca, r.epsA, occA.map((o, k) => o + occB[k]))
    : describeOrbitals(ctx, r.basis, S, r.Ca, r.epsA, occA);
  const mosB = restricted ? null : describeOrbitals(ctx, r.basis, S, r.Cb, r.epsB, occB);
  let spin = 0;
  if (!restricted) {
    // popolazione di spin di Mulliken sull'atomo A
    for (let mu = 0; mu < nA; mu++) for (let nu = 0; nu < n; nu++) spin += (r.Pa[mu * n + nu] - r.Pb[mu * n + nu]) * S[nu * n + mu];
  }
  return {
    method: kind, R: Rang, E: r.energy, T, V: r.energy - T, converged: r.converged, iterations: r.iterations,
    S2: r.S2, spinA: spin, bondOrder: mayerOrder(r.Pa, r.Pb, S, n, nA), charges: mullikenCharges(Pt, S, n, nA, mol.Z),
    mos, mosB, Pa: r.Pa, Pb: r.Pb, Ca: r.Ca, Cb: restricted ? null : r.Cb, S,
    overlap: valenceOverlaps(ctx, S),
  };
}

/** Integrali di sovrapposizione fra gli orbitali atomici di valenza di A e B (quanto "si toccano"). */
function valenceOverlaps(ctx, S) {
  const n = ctx.n;
  const val = ctx.aos.filter(a => !a.extra);
  const pick = (atom) => val.filter(a => a.atom === atom && (a.axis === '' || a.axis === 'z' || a.axis === 'x'));
  const out = [];
  for (const a of pick(0)) for (const b of pick(1)) {
    if ((a.axis === 'x') !== (b.axis === 'x')) continue;
    let s = 0;
    const Sb = matVec(S, b.v, n);
    for (let i = 0; i < n; i++) s += a.v[i] * Sb[i];
    // il segno dipende dalla fase arbitraria: si riporta |S|
    // nelle molecole omonucleari S(2s_A, 2pz_B) = S(2pz_A, 2s_B): basta una volta
    const key = ctx.homonuclear ? [a.label, b.label].sort().join('|') : `${a.label}|${b.label}`;
    if (Math.abs(s) > 1e-4 && !out.some(o => o.key === key)) out.push({ a: a.label, b: b.label, S: Math.abs(s), key });
  }
  return out.sort((x, y) => y.S - x.S).slice(0, 6);
}

function hfAt(ctx, mol, Rang, { unrestricted, guess, bias }) {
  const atoms = geometry(mol.Z, Rang * ANGSTROM_TO_BOHR);
  const basisObj = buildBasis(atoms, ctx.basis);
  const one = oneElectron(basisObj, atoms);
  const eri = twoElectron(basisObj);
  const opts = { basis: ctx.basis, basisObj, one, eri, charge: mol.charge, multiplicity: mol.mult, unrestricted, conv: 1e-10, maxIter: 200 };
  if (bias) {
    // prima una SCF con la polarizzazione di spin imposta, poi quella vera a partire da lì
    const pre = runHF(atoms, { ...opts, guess, bias, maxIter: 60, conv: 1e-6 });
    guess = { Pa: pre.Pa, Pb: pre.Pb };
  }
  return runHF(atoms, { ...opts, guess });
}

/** Stima che rompe la simmetria: l'HOMO α ruotato verso il LUMO e l'HOMO β in verso opposto (α su un atomo, β sull'altro). */
function perturbedGuess(r, theta = Math.PI / 7) {
  const n = r.n;
  const dens = (C, nocc, sgn) => {
    const P = new Float64Array(n * n);
    const h = nocc - 1, c = Math.cos(theta), s = sgn * Math.sin(theta);
    for (let k = 0; k < nocc; k++) {
      for (let i = 0; i < n; i++) {
        const ci = k === h && nocc < n ? c * C[i * n + h] + s * C[i * n + nocc] : C[i * n + k];
        if (!ci) continue;
        for (let j = 0; j < n; j++) {
          const cj = k === h && nocc < n ? c * C[j * n + h] + s * C[j * n + nocc] : C[j * n + k];
          P[i * n + j] += ci * cj;
        }
      }
    }
    return P;
  };
  return { Pa: dens(r.Ca, r.nalpha, 1), Pb: dens(r.Cb, r.nbeta, -1) };
}

/** Stima per O₂ e simili (molteplicità > 1) agli atomi separati: A con il suo spin, B con la densità media α = β. */
function highSpinGuess(ctx) {
  const { A, B, n } = ctx;
  const Pa = new Float64Array(n * n), Pb = new Float64Array(n * n);
  for (let i = 0; i < A.n; i++) for (let j = 0; j < A.n; j++) { Pa[i * n + j] = A.Pa[i * A.n + j]; Pb[i * n + j] = A.Pb[i * A.n + j]; }
  for (let i = 0; i < B.n; i++) for (let j = 0; j < B.n; j++) {
    const v = 0.5 * B.Psph[i * B.n + j];
    Pa[(A.n + i) * n + A.n + j] = v; Pb[(A.n + i) * n + A.n + j] = v;
  }
  return { Pa, Pb };
}

function makePoint(ctx, kind, Rang, r) {
  const mol = ctx.mol;
  if (kind !== 'fci') return pointFromHF(ctx, mol, Rang, r, kind);
  const f = fci2({ n: r.n, C: r.Ca, H: r.H, eri: r.eri, Enuc: r.Enuc });
  let T = 0;
  for (let k = 0; k < r.n * r.n; k++) T += f.P[k] * r.one.T[k];
  const half = f.P.map(v => v / 2);
  return {
    method: 'fci', R: Rang, E: f.energy, T, V: f.energy - T, converged: r.converged, natOcc: f.natOcc, configs: f.configs,
    bondOrder: mayerOrder(half, half, r.S, r.n, ctx.A.n), charges: mullikenCharges(f.P, r.S, r.n, ctx.A.n, mol.Z),
    mos: describeOrbitals(ctx, r.basis, r.S, r.Ca, r.epsA, f.moOcc), mosB: null, Pa: half, Pb: half, Ca: r.Ca, Cb: null, S: r.S,
    overlap: valenceOverlaps(ctx, r.S), Ehf: r.energy, dim: f.dim,
  };
}

/**
 * Scansione con Hartree–Fock (kind = 'rhf' | 'uhf' | 'exact') o FCI ('fci'). onPoint(point, index) riceve i punti
 * nell'ordine in cui sono calcolati (per UHF anche due candidati per la stessa distanza: vale quello di energia
 * minore). Restituisce i punti ordinati dalla distanza maggiore alla minore.
 *
 * Ogni distanza parte dalla densità della distanza vicina già risolta, così si segue con continuità lo stesso stato.
 * UHF si cerca in due passate e si tiene la soluzione più bassa: dagli atomi separati verso l'interno (spin opposti
 * sui due atomi) e dall'equilibrio verso l'esterno, rompendo ogni volta la simmetria fra HOMO e LUMO.
 */
export function scanHF(ctx, kind, { grid, onPoint } = {}) {
  const mol = ctx.mol;
  const R = grid ?? scanGrid(mol);
  const unrestricted = kind === 'uhf' || kind === 'exact' || mol.mult !== 1;
  const i0 = R.reduce((b, r, i) => (Math.abs(r - mol.exp.re) < Math.abs(R[b] - mol.exp.re) ? i : b), 0);
  const fromEq = [];
  for (let i = i0; i >= 0; i--) fromEq.push(i);
  for (let i = i0 + 1; i < R.length; i++) fromEq.push(i);
  const best = new Array(R.length);

  const pass = (order, { first, perturb }) => {
    let prev = null, prevI = -1, start = null;
    for (const i of order) {
      let guess = null, bias = null;
      const base = prev && Math.abs(prevI - i) === 1 ? prev : start;
      if (base) guess = perturb ? perturbedGuess(base) : { Pa: base.Pa, Pb: base.Pb };
      else if (first) ({ guess, bias } = first(R[i]));
      const r = hfAt(ctx, mol, R[i], { unrestricted, guess, bias });
      const point = makePoint(ctx, kind, R[i], r);
      if (!best[i] || (point.converged && (!best[i].converged || point.E < best[i].E - 1e-9))) {
        best[i] = point;
        onPoint?.(point, i);
      }
      prev = r; prevI = i;
      start ??= r;
    }
  };

  if (kind === 'uhf') {
    const inward = R.map((_, i) => i);
    if (mol.mult === 1) pass(inward, { first: (Rang) => ({ guess: brokenSymmetryGuess(ctx), bias: spinBias(ctx, Rang) }) });
    else pass(inward, { first: () => ({ guess: highSpinGuess(ctx) }) });
    pass(fromEq, { perturb: true });
  } else {
    pass(fromEq, {});
  }
  return best;
}

/** Fock aggiuntiva che spinge gli elettroni α verso A e i β verso B (solo per la prima stima a simmetria rotta). */
function spinBias(ctx, Rang) {
  const mol = ctx.mol;
  const atoms = geometry(mol.Z, Rang * ANGSTROM_TO_BOHR);
  const B = buildBasis(atoms, ctx.basis);
  const { S } = oneElectron(B, atoms);
  const n = B.nbf, nA = ctx.A.n, lam = 0.15;
  const a = new Float64Array(n * n), b = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const w = 0.5 * ((i < nA ? -1 : 1) + (j < nA ? -1 : 1)) * lam * S[i * n + j];
    a[i * n + j] = w; b[i * n + j] = -w;
  }
  return { a, b };
}

/** Scansione GFN2-xTB: energie (Eh) e riferimento degli atomi liberi. */
export function scanGFN2(mol, { grid, onPoint } = {}) {
  const R = grid ?? scanGrid(mol);
  const calc = new GFN2xTB();
  const charges = mol.atomCharges ?? [0, 0];
  let Eatoms = 0;
  for (let k = 0; k < 2; k++) {
    Eatoms += calc.compute([mol.Z[k]], [0, 0, 0], { charge: charges[k], uhf: mol.atoms[k] - 1, gradient: false }).energy;
  }
  const pts = [];
  R.forEach((Rang, i) => {
    const Rb = Rang * ANGSTROM_TO_BOHR;
    const r = calc.compute(mol.Z, [0, 0, -Rb / 2, 0, 0, Rb / 2], { charge: mol.charge, uhf: mol.mult - 1, gradient: false });
    const p = { method: 'gfn2', R: Rang, E: r.energy, converged: r.converged, charges: r.charges ? Array.from(r.charges) : null };
    pts.push(p);
    onPoint?.(p, i);
  });
  return { points: pts, Eatoms };
}

// ---------------------------------------------------------------------------
// Analisi delle curve
// ---------------------------------------------------------------------------

/**
 * Minimo della curva E(R) (Å, Eh): interpolazione cubica attorno al punto più basso. Restituisce re (Å), Emin (Eh),
 * la curvatura k (Eh/bohr²) e ωe (cm⁻¹) per la massa ridotta μ (u). null se il minimo è al bordo.
 */
export function curveMinimum(R, E, muAmu) {
  const idx = R.map((_, i) => i).sort((a, b) => R[a] - R[b]);
  const x = idx.map(i => R[i]), y = idx.map(i => E[i]);
  let k = 0;
  for (let i = 1; i < y.length; i++) if (y[i] < y[k]) k = i;
  if (k === 0 || k === y.length - 1) return null;
  // cubica per i quattro punti attorno al minimo (il punto più basso, i due vicini e il successivo più basso), in bohr
  const extra = k + 2 < y.length && (k - 2 < 0 || y[k + 2] < y[k - 2]) ? k + 2 : k - 2;
  const sel = [k - 1, k, k + 1, extra].filter(i => i >= 0 && i < y.length).sort((a, b) => a - b);
  const x0 = x[k] * ANGSTROM_TO_BOHR;
  const A = [], b = [];
  for (const i of sel) { const t = x[i] * ANGSTROM_TO_BOHR - x0; A.push([1, t, t * t, t * t * t]); b.push(y[i]); }
  const c = leastSquares(A, b);
  // minimo: c1 + 2 c2 t + 3 c3 t² = 0
  let t;
  if (Math.abs(c[3]) < 1e-14) t = -c[1] / (2 * c[2]);
  else {
    const disc = 4 * c[2] * c[2] - 12 * c[3] * c[1];
    const r1 = (-2 * c[2] + Math.sqrt(Math.max(disc, 0))) / (6 * c[3]), r2 = (-2 * c[2] - Math.sqrt(Math.max(disc, 0))) / (6 * c[3]);
    t = [r1, r2].filter(v => 2 * c[2] + 6 * c[3] * v > 0).sort((p, q) => Math.abs(p) - Math.abs(q))[0] ?? 0;
  }
  const kk = 2 * c[2] + 6 * c[3] * t;
  const Emin = c[0] + c[1] * t + c[2] * t * t + c[3] * t ** 3;
  const we = kk > 0 && muAmu ? Math.sqrt(kk / (muAmu * AMU_ME)) * HARTREE_CM : null;
  return { re: (x0 + t) * BOHR_A, Emin, k: kk, we };
}

function leastSquares(A, b) {
  const m = A[0].length, N = new Float64Array(m * m), r = new Float64Array(m);
  for (let i = 0; i < A.length; i++) for (let p = 0; p < m; p++) { r[p] += A[i][p] * b[i]; for (let q = 0; q < m; q++) N[p * m + q] += A[i][p] * A[i][q]; }
  // eliminazione di Gauss con pivot
  const M = Array.from({ length: m }, (_, i) => [...N.slice(i * m, i * m + m), r[i]]);
  for (let c = 0; c < m; c++) {
    let p = c;
    for (let i = c + 1; i < m; i++) if (Math.abs(M[i][c]) > Math.abs(M[p][c])) p = i;
    [M[c], M[p]] = [M[p], M[c]];
    for (let i = c + 1; i < m; i++) { const f = M[i][c] / M[c][c]; for (let j = c; j <= m; j++) M[i][j] -= f * M[c][j]; }
  }
  const x = new Array(m).fill(0);
  for (let i = m - 1; i >= 0; i--) { let s = M[i][m]; for (let j = i + 1; j < m; j++) s -= M[i][j] * x[j]; x[i] = s / M[i][i]; }
  return x;
}

/**
 * Derivata dE/dR (Eh/bohr) su una griglia non uniforme (differenze finite a tre punti di Lagrange).
 * R in Å, E in Eh; restituisce un array allineato a R.
 */
export function derivative(R, E) {
  const n = R.length, d = new Array(n).fill(NaN);
  const x = R.map(r => r * ANGSTROM_TO_BOHR);
  for (let i = 0; i < n; i++) {
    const [a, b, c] = i === 0 ? [0, 1, 2] : i === n - 1 ? [n - 3, n - 2, n - 1] : [i - 1, i, i + 1];
    const xa = x[a], xb = x[b], xc = x[c], t = x[i];
    d[i] = E[a] * (2 * t - xb - xc) / ((xa - xb) * (xa - xc)) + E[b] * (2 * t - xa - xc) / ((xb - xa) * (xb - xc)) + E[c] * (2 * t - xa - xb) / ((xc - xa) * (xc - xb));
  }
  return d;
}

/**
 * Teorema del viriale molecolare (Slater 1933), esatto per autofunzioni dell'hamiltoniano elettronico:
 * 2T + V = −R dE/dR, E = T + V  ⇒  T = −E − R dE/dR,  V = 2E + R dE/dR.
 * Da una curva di energia (anche sperimentale) dà l'energia cinetica e potenziale degli elettroni a ogni R.
 */
export function virialSplit(R, E, dEdR) {
  return R.map((r, i) => {
    const Rb = r * ANGSTROM_TO_BOHR;
    return { T: -E[i] - Rb * dEdR[i], V: 2 * E[i] + Rb * dEdR[i] };
  });
}

/** Riassunto di una curva calcolata contro l'esperimento. */
export function curveSummary(mol, R, E, Eatoms) {
  const mu = reducedMass(mol.Z);
  const min = curveMinimum(R, E, mu);
  const exp = experimentalWell(mol);
  if (!min) return { bound: false, exp };
  const De = (Eatoms - min.Emin) * HARTREE_EV;
  // un minimo più debole di un decimo di quello sperimentale (o sopra gli atomi) non è un legame: errore di
  // sovrapposizione della base o residuo di dispersione
  const bound = De > 0.1 * exp.De;
  return { bound, re: min.re, De, we: min.we, k: min.k, exp };
}

export { HARTREE_EV, ANGSTROM_TO_BOHR, BOHR_A };
