// Metodo di Hartree–Fock (Roothaan–Hall per shell chiuse, Pople–Nesbet per shell aperte).
//
//   F C = S C ε,     F_μν = H_μν + Σ_λσ P_λσ [ (μν|λσ) − ½ (μλ|νσ) ]            (RHF)
//   F^α = H + J[P^α + P^β] − K[P^α],   F^β = H + J[P^α + P^β] − K[P^β]            (UHF)
//   E = ½ Σ_μν [ P^α (H + F^α) + P^β (H + F^β) ]_μν + V_NN
//
// Convergenza accelerata con DIIS (Pulay, 1980): si estrapola la matrice di Fock che minimizza
// il commutatore FPS − SPF.

import { buildBasis, oneElectron, twoElectron, buildJK } from './integrals.js';
import { eigh, matmul, symFunction, solve } from './linalg.js';

export const ANGSTROM_TO_BOHR = 1 / 0.52917721090;
export const HARTREE_EV = 27.211386245988;
export const HARTREE_KJMOL = 2625.4996394799;
export const AU_TO_DEBYE = 2.541746473;

export function nuclearRepulsion(atoms) {
  let e = 0;
  for (let i = 0; i < atoms.length; i++) {
    for (let j = 0; j < i; j++) {
      const a = atoms[i].xyz, b = atoms[j].xyz;
      e += atoms[i].Z * atoms[j].Z / Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
    }
  }
  return e;
}

/** Diagonalizza F nella base ortonormale X = S^{−1/2}: restituisce ε e C = X C'. */
function diagonalize(F, X, n) {
  const Fp = matmul(matmul(X, F, n), X, n); // X simmetrica
  const { values, vectors } = eigh(Fp, n);
  const C = matmul(X, vectors, n);
  return { eps: values, C };
}

function density(C, n, nocc) {
  const P = new Float64Array(n * n);
  for (let k = 0; k < nocc; k++) {
    for (let i = 0; i < n; i++) {
      const ci = C[i * n + k];
      if (ci === 0) continue;
      for (let j = 0; j < n; j++) P[i * n + j] += ci * C[j * n + k];
    }
  }
  return P;
}

/** Errore DIIS: e = F P S − S P F, espresso nella base ortonormale. */
function diisError(F, P, S, X, n) {
  const FPS = matmul(matmul(F, P, n), S, n);
  const e = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) e[i * n + j] = FPS[i * n + j] - FPS[j * n + i];
  return matmul(matmul(X, e, n), X, n);
}

function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

/**
 * Calcolo Hartree–Fock.
 * @param atoms [{Z, xyz (bohr)}]
 * @param opts { basis, charge, multiplicity, maxIter, conv, eri (riuso), onProgress, bias: {a, b} aggiunte alle Fock α e β }
 */
export function runHF(atoms, opts = {}) {
  const basisName = opts.basis ?? 'STO-3G';
  const charge = opts.charge ?? 0;
  const nel = atoms.reduce((s, a) => s + a.Z, 0) - charge;
  const mult = opts.multiplicity ?? (nel % 2 === 0 ? 1 : 2);
  if ((nel + mult - 1) % 2 !== 0 || mult < 1 || mult - 1 > nel) throw new Error(`Molteplicità ${mult} incompatibile con ${nel} elettroni`);
  const nalpha = (nel + mult - 1) / 2;
  const nbeta = nel - nalpha;
  const unrestricted = opts.unrestricted ?? (mult !== 1);

  const basis = opts.basisObj ?? buildBasis(atoms, basisName);
  const n = basis.nbf;
  const t0 = Date.now();
  const one = opts.one ?? oneElectron(basis, atoms);
  const eri = opts.eri ?? twoElectron(basis, { onProgress: opts.onProgress });
  const tInts = Date.now() - t0;
  const { S } = one;
  const H = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) H[i] = one.T[i] + one.V[i];
  // ortogonalizzazione simmetrica (Löwdin), scartando le combinazioni quasi dipendenti
  const X = symFunction(S, n, (s) => (s > 1e-8 ? 1 / Math.sqrt(s) : 0));
  const Enuc = nuclearRepulsion(atoms);

  // Stima iniziale: Hückel generalizzato (Wolfsberg–Helmholz), F_μν = 1,75 S_μν (H_μμ + H_νν)/2
  let Pa, Pb;
  if (opts.guess) {
    Pa = Float64Array.from(opts.guess.Pa);
    Pb = Float64Array.from(opts.guess.Pb);
  } else {
    const F0 = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        F0[i * n + j] = i === j ? H[i * n + i] : 0.875 * S[i * n + j] * (H[i * n + i] + H[j * n + j]);
      }
    }
    const { C } = diagonalize(F0, X, n);
    Pa = density(C, n, nalpha);
    Pb = density(C, n, nbeta);
    if (unrestricted && mult === 1 && opts.breakSymmetry) {
      // miscela HOMO/LUMO per rompere la simmetria di spin (UHF a shell aperta)
      const Cm = Float64Array.from(C);
      const h = nalpha - 1, l = nalpha;
      const th = Math.PI / 4;
      for (let i = 0; i < n; i++) {
        const a = C[i * n + h], b = C[i * n + l];
        Cm[i * n + h] = Math.cos(th) * a + Math.sin(th) * b;
      }
      Pa = density(Cm, n, nalpha);
    }
  }

  const maxIter = opts.maxIter ?? 128;
  const conv = opts.conv ?? 1e-10;
  const diisF = [];
  const diisE = [];
  const DIIS_MAX = 8;
  let E = 0;
  let Eold = 0;
  let converged = false;
  let iter = 0;
  let Fa, Fb, resA, resB, errNorm = 1;
  const history = [];

  for (iter = 1; iter <= maxIter; iter++) {
    const { J, K } = buildJK(eri, unrestricted ? [Pa, Pb] : [Pa, Pb], n);
    Fa = new Float64Array(n * n);
    Fb = new Float64Array(n * n);
    for (let i = 0; i < n * n; i++) {
      Fa[i] = H[i] + J[i] - K[0][i];
      Fb[i] = H[i] + J[i] - K[1][i];
    }
    // energia elettronica
    let Eel = 0;
    for (let i = 0; i < n * n; i++) Eel += 0.5 * (Pa[i] * (H[i] + Fa[i]) + Pb[i] * (H[i] + Fb[i]));
    E = Eel + Enuc;
    // polarizzazione di spin imposta (solo per preparare una stima iniziale a simmetria rotta): entra nella Fock da
    // diagonalizzare, non nell'energia
    if (opts.bias) for (let i = 0; i < n * n; i++) { Fa[i] += opts.bias.a[i]; Fb[i] += opts.bias.b[i]; }

    const ea = diisError(Fa, Pa, S, X, n);
    const eb = unrestricted ? diisError(Fb, Pb, S, X, n) : null;
    errNorm = Math.sqrt(dot(ea, ea) + (eb ? dot(eb, eb) : 0));
    history.push({ iter, E, err: errNorm });
    if (iter > 1 && Math.abs(E - Eold) < conv && errNorm < Math.sqrt(conv) * 10) {
      converged = true;
      break;
    }
    Eold = E;

    // DIIS
    let Fna = Fa, Fnb = Fb;
    diisF.push(unrestricted ? [Fa, Fb] : [Fa]);
    diisE.push(unrestricted ? [ea, eb] : [ea]);
    if (diisF.length > DIIS_MAX) { diisF.shift(); diisE.shift(); }
    if (diisF.length >= 2) {
      const m = diisF.length;
      const B = new Float64Array((m + 1) * (m + 1));
      for (let i = 0; i < m; i++) {
        for (let j = 0; j <= i; j++) {
          let v = 0;
          for (let s = 0; s < diisE[i].length; s++) v += dot(diisE[i][s], diisE[j][s]);
          B[i * (m + 1) + j] = B[j * (m + 1) + i] = v;
        }
        B[i * (m + 1) + m] = B[m * (m + 1) + i] = -1;
      }
      const rhs = new Float64Array(m + 1);
      rhs[m] = -1;
      const c = solve(B, rhs, m + 1);
      if (c.every(Number.isFinite)) {
        Fna = new Float64Array(n * n);
        Fnb = new Float64Array(n * n);
        for (let k = 0; k < m; k++) {
          const fk = diisF[k];
          for (let i = 0; i < n * n; i++) {
            Fna[i] += c[k] * fk[0][i];
            if (unrestricted) Fnb[i] += c[k] * fk[1][i];
          }
        }
        if (!unrestricted) Fnb = Fna;
      }
    } else if (!unrestricted) {
      Fnb = Fna;
    }
    resA = diagonalize(Fna, X, n);
    resB = unrestricted ? diagonalize(Fnb, X, n) : resA;
    const Pa_new = density(resA.C, n, nalpha);
    const Pb_new = density(resB.C, n, nbeta);
    // leggero smorzamento nelle prime iterazioni
    const damp = iter < 4 && !opts.guess ? 0.3 : 0;
    for (let i = 0; i < n * n; i++) {
      Pa[i] = damp * Pa[i] + (1 - damp) * Pa_new[i];
      Pb[i] = damp * Pb[i] + (1 - damp) * Pb_new[i];
    }
  }

  // Orbitali finali dalla Fock convergente
  const fa = diagonalize(Fa, X, n);
  const fb = unrestricted ? diagonalize(Fb, X, n) : fa;
  Pa = density(fa.C, n, nalpha);
  Pb = density(fb.C, n, nbeta);

  // ⟨S²⟩ per UHF: S(S+1) + N_β − Σ_ij |⟨i_α|j_β⟩|²
  let S2 = null;
  if (unrestricted) {
    const Sz = (nalpha - nbeta) / 2;
    const SCb = matmul(S, fb.C, n);
    let ov = 0;
    for (let i = 0; i < nalpha; i++) {
      for (let j = 0; j < nbeta; j++) {
        let s = 0;
        for (let m = 0; m < n; m++) s += fa.C[m * n + i] * SCb[m * n + j];
        ov += s * s;
      }
    }
    S2 = Sz * (Sz + 1) + nbeta - ov;
  }

  return {
    atoms, basis, n, nel, nalpha, nbeta, charge, multiplicity: mult, unrestricted,
    energy: E, Enuc, converged, iterations: iter, history, errNorm,
    S, H, one, eri, X,
    Pa, Pb, Ca: fa.C, Cb: fb.C, epsA: fa.eps, epsB: fb.eps, Fa, Fb,
    S2,
    timeIntegrals: tInts,
    time: Date.now() - t0,
  };
}
