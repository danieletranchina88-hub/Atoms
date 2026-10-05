// GFN2-xTB come motore delle forze della sandbox: converte unità (Å, eV) e carica/spin/campo della simulazione.

import { GFN2xTB, gfn2Supports, GFN2_MAX_Z } from './gfn2.js';
import { GFN2 as PAR } from './gfn2Data.js';

const BOHR_ANG = 0.52917721090;
const HARTREE_EV = 27.211386245988;
const VA_AU = 1 / 51.422067476; // 1 V/Å in hartree/(e·bohr)
export const GFN2_MAX_ATOMS = 120;

/** Elettroni di valenza GFN2 (occupazioni di riferimento delle shell). */
export function gfn2Electrons(Z, charge = 0) {
  return Z.reduce((s, z) => s + gfn2Valence(z), 0) - charge;
}
const gfn2Valence = (z) => PAR.elements[z - 1].shells.reduce((s, sh) => s + sh.refocc, 0);

export function gfn2Feasible(Z) {
  if (!Z.length) return { ok: false, reason: 'La scatola è vuota.' };
  if (!gfn2Supports(Z)) return { ok: false, reason: `GFN2-xTB copre gli elementi da H (1) a Rn (${GFN2_MAX_Z}).` };
  if (Z.length > GFN2_MAX_ATOMS) return { ok: false, reason: `GFN2-xTB è limitato a ${GFN2_MAX_ATOMS} atomi nella sandbox (ce ne sarebbero ${Z.length}).` };
  return { ok: true };
}

/**
 * multiplicity: 2S+1 fissato, oppure null (spin più basso compatibile con il numero di elettroni);
 * field: campo elettrico uniforme [Ex, Ey, Ez] in V/Å; solvent: 'water' (solvente implicito ALPB) oppure null.
 */
export function makeGFN2Provider({ Tel = 300, multiplicity = null, field = null, solvent = null } = {}) {
  // tolleranza SCF per la dinamica: con 1e-5 la deriva NVE è la stessa che con 1e-6 (domina l'integratore), con il 30 % di iterazioni in meno
  const calc = new GFN2xTB({ Tel, etol: 1e-9, ptol: 1e-5, solvent });
  let charge = 0, uhfUsed = null, fieldVA = null;
  const setField = (f) => {
    fieldVA = f && f.some(v => v) ? f.slice() : null;
    calc.field = fieldVA ? fieldVA.map(v => v * VA_AU) : null;
  };
  setField(field);
  const provider = {
    kind: 'gfn2',
    info: {},
    multiplicity,
    reset() { calc.reset(); },
    setField,
    get field() { return fieldVA; },
    setSolvent(sv) { calc.solvent = sv || null; calc.reset(); },
    _setTol(ptol, etol = Math.min(1e-9, ptol * ptol * 10)) { calc.ptol = ptol; calc.etol = etol; },
    get solvent() { return calc.solvent; },
    compute(Z, pos, F, ctx = {}) {
      if ((ctx.charge ?? 0) !== charge) { charge = ctx.charge ?? 0; calc.reset(); }
      const nel = gfn2Electrons(Z, charge);
      let uhf;
      if (provider.multiplicity) {
        uhf = provider.multiplicity - 1;
        if (nel < 0 || (nel - uhf) % 2 !== 0 || uhf > nel) throw new Error(`Molteplicità ${provider.multiplicity} impossibile con ${nel} elettroni di valenza: scegli ${nel % 2 ? 'doppietto, quartetto…' : 'singoletto, tripletto…'}`);
      } else uhf = ((nel % 2) + 2) % 2;
      if (uhf !== uhfUsed) { uhfUsed = uhf; calc.reset(); }
      const xb = Float64Array.from(pos, v => v / BOHR_ANG);
      let r = calc.compute(Z, xb, { charge, uhf });
      if (!r.converged) { calc.reset(); r = calc.compute(Z, xb, { charge, uhf }); }
      if (!r.converged || !Number.isFinite(r.energy)) throw new Error('GFN2-xTB: SCF non convergente, dinamica sospesa.');
      const f = HARTREE_EV / BOHR_ANG;
      for (let k = 0; k < F.length; k++) F[k] = -r.gradient[k] * f;
      const bonds = mayerBonds(calc.last, Z.length);
      const { spin } = spinDensity(calc.last, Z.length);
      const mu = r.dipole.map(v => v * BOHR_ANG / 0.20819434); // e·bohr → debye
      provider.info = {
        method: 'GFN2-xTB', nbf: r.nao, iterations: r.iterations, converged: r.converged, charge, Tel,
        multiplicity: uhf + 1, homo: r.homoLumo.homo * HARTREE_EV, lumo: r.homoLumo.lumo * HARTREE_EV, gap: r.homoLumo.gap * HARTREE_EV,
        dipole: r.dipole.map(v => v * BOHR_ANG), field: fieldVA, solvent: calc.solvent,
        solvation: calc.solvent ? r.parts.solvation * HARTREE_EV : null, Sz: uhf / 2, dipoleDebye: Math.hypot(...mu),
      };
      const parts = { gfn2: r.energy * HARTREE_EV, repulsion: r.parts.repulsion * HARTREE_EV, dispersion3: r.parts.atm * HARTREE_EV };
      return { E: r.energy * HARTREE_EV, parts, q: Float64Array.from(r.charges), bonds, hbonds: [], spin };
    },
    /** Funzione d'onda per il disegno: densità di valenza, HOMO e LUMO nella base STO-nG di GFN2. */
    /**
     * Funzione d'onda per il disegno. Con { occupied: true } aggiunge gli orbitali doppiamente occupati (per la
     * localizzazione); levels contiene sempre tutte le energie orbitali e le occupazioni (diagramma dei livelli).
     */
    wavefunction({ occupied = false } = {}) {
      const w = calc.last;
      if (!w) return null;
      const n = w.n;
      let ho = -1, lu = -1;
      for (let k = 0; k < n; k++) { if (w.f[k] > 0.5) ho = k; else if (lu < 0) lu = k; }
      const vec = (k) => (k < 0 || k >= w.nC ? null : { e: w.e[k], spin: w.f[k] > 1.5 || w.f[k] < 0.5 ? 'α+β' : 'α', c: w.Ct.slice(k * n, (k + 1) * n) });
      const { Ps } = spinDensity(w, w.Z.length);
      const orbitals = [];
      for (let k = Math.max(0, ho - 5); k < Math.min(w.nC, ho + 7); k++) orbitals.push({ ...vec(k), id: `MO:${k}`, index: k, occ: w.f[k] });
      const levels = { e: Float64Array.from(w.e), f: Float64Array.from(w.f), homo: ho };
      let occ = null;
      if (occupied) {
        let m = 0;
        while (m < w.nC && w.f[m] > 1.5) m++;
        occ = { n, m, Ct: w.Ct.slice(0, m * n) };
      }
      return { kind: 'xtb', Ps, orbitals, levels, occ, Z: w.Z.slice(), pos: Float64Array.from(w.pos, v => v * BOHR_ANG), n, P: w.P, homo: vec(ho), lumo: vec(lu), eUnit: HARTREE_EV };
    },
  };
  return provider;
}

/**
 * Ordini di legame di Mayer B_AB = Σ_{μ∈A, ν∈B} (PS)_μν (PS)_νμ, solo per le coppie entro 4,5 Å
 * (oltre l'ordine di legame è trascurabile): costa O(coppie vicine · n) invece del prodotto completo P·S.
 */
function mayerBonds(w, N) {
  if (!w?.P || !w.S) return [];
  const { P, S, n, basis, pos } = w;
  const first = new Int32Array(N + 1);
  for (let mu = 0; mu < n; mu++) first[basis.aoAtom[mu] + 1] = mu + 1;
  for (let A = 1; A <= N; A++) first[A] = Math.max(first[A], first[A - 1]);
  const cut2 = (4.5 / BOHR_ANG) ** 2;
  const ps = (mu, nu) => { let v = 0; const r = mu * n; for (let k = 0; k < n; k++) v += P[r + k] * S[k * n + nu]; return v; };
  const bonds = [];
  for (let A = 0; A < N; A++) for (let C = A + 1; C < N; C++) {
    const dx = pos[3 * A] - pos[3 * C], dy = pos[3 * A + 1] - pos[3 * C + 1], dz = pos[3 * A + 2] - pos[3 * C + 2];
    if (dx * dx + dy * dy + dz * dz > cut2) continue;
    let b = 0;
    for (let mu = first[A]; mu < first[A + 1]; mu++) for (let nu = first[C]; nu < first[C + 1]; nu++) b += ps(mu, nu) * ps(nu, mu);
    if (b > 0.1) bonds.push({ i: A, j: C, n: b, w: 1 });
  }
  return bonds;
}

/** Restricted spatial orbitals with two Fermi occupation channels, as in GFN2 (not UHF). */
function spinDensity(w, N) {
  const { n, nC, Ct, fa, fb, S, basis } = w;
  const Ps = new Float64Array(n * n), spin = new Float64Array(N);
  for (let k = 0; k < nC; k++) {
    const f = fa[k] - fb[k];
    if (Math.abs(f) < 1e-14) continue;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Ps[i * n + j] += f * Ct[k * n + i] * Ct[k * n + j];
  }
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) spin[basis.aoAtom[i]] += Ps[i * n + j] * S[j * n + i];
  return { Ps, spin };
}
