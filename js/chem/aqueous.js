// Chimica in soluzione acquosa: equilibrio chimico completo di un becher.
//
//  • Bilanci di massa per ogni componente e legge di azione di massa per ogni specie:
//      c_s = β_s Π {c}^{a_sc} / γ_s,   T_c = Σ_s a_sc c_s + Σ_p a_pc m_p
//    risolti con il metodo di Newton–Raphson nei logaritmi delle concentrazioni libere.
//  • Coefficienti di attività dall'equazione di Davies: log γ = −A z² [√I/(1+√I) − 0,3 I].
//  • Precipitati: un solido entra nel sistema quando il suo indice di saturazione log(Q/K_sp) supera 0
//    ed esce se la sua quantità diventa negativa (regola delle fasi di Gibbs).
//  • Metalli: reazioni redox spinte finché la f.e.m. (equazione di Nernst) si annulla.
//  • Gas: CO₂ si libera oltre la solubilità di Henry, H₂ lascia subito la soluzione.
//  • Calore: legge di Hess con le entalpie di formazione; ΔT = q / (m c_p).
//  • Colore: legge di Lambert–Beer, A(λ) = Σ ε(λ) c l, e funzioni colorimetriche CIE 1931.

import { COMPONENTS, SPECIES, SOLIDS, METALS, REAGENTS, H2_OVERPOTENTIAL } from './aqueousData.js';

const LN10 = Math.LN10;
const A_DAVIES = 0.509;          // a 25 °C
const NERNST = 0.05916;          // RT ln10 / F a 25 °C (V)
const HENRY_CO2 = 0.0339;        // mol/(L·atm) a 25 °C
const DFH_WATER = -285.83;       // kJ/mol
const DFH_CO2 = -393.51;
const CP_WATER = 4.184;          // J/(g K)

const COMP_IDS = Object.keys(COMPONENTS);

// specie e solidi in forma numerica
const SP = SPECIES.map(([id, label, z, comp, logb, w, dfh, bands]) => ({ id, label, z, comp, lnb: logb * LN10, logb, w, dfh, bands: bands ?? null }));
const SPI = new Map(SP.map((s, i) => [s.id, i]));
const SO = SOLIDS.map(([id, label, comp, logK, w, dfh, M, color]) => ({ id, label, comp, lnK: logK * LN10, logK, w, dfh, M, color }));
const SOI = new Map(SO.map((s, i) => [s.id, i]));
// entalpie mancanti: si assume che la formazione dalle componenti non scambi calore (stima dichiarata)
const freeDfh = {};
for (const s of SP) {
  const keys = Object.keys(s.comp);
  if (keys.length === 1 && s.comp[keys[0]] === 1 && s.dfh !== null) freeDfh[keys[0]] = s.dfh;
}
freeDfh.H = 0;
const estimateDfh = (comp, w) => Object.entries(comp).reduce((h, [c, a]) => h + a * (freeDfh[c] ?? 0), 0) + w * DFH_WATER;
for (const s of SP) { s.dfhEstimated = s.dfh === null; if (s.dfh === null) s.dfh = estimateDfh(s.comp, s.w); }
for (const s of SO) { s.dfhEstimated = s.dfh === null; if (s.dfh === null) s.dfh = estimateDfh(s.comp, s.w); }

export const speciesInfo = (id) => SP[SPI.get(id)];
export const solidInfo = (id) => SO[SOI.get(id)];

function lnGamma(z, sqrtI, I) {
  if (!z) return 0;
  return -A_DAVIES * z * z * (sqrtI / (1 + sqrtI) - 0.3 * I) * LN10;
}

/** Risolve A x = b (eliminazione di Gauss con pivot parziale). */
function linsolve(A, b) {
  const n = b.length;
  const M = A.map((r, i) => [...r, b[i]]);
  for (let k = 0; k < n; k++) {
    let p = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
    if (Math.abs(M[p][k]) < 1e-300) return null;
    [M[k], M[p]] = [M[p], M[k]];
    for (let i = k + 1; i < n; i++) {
      const f = M[i][k] / M[k][k];
      if (f === 0) continue;
      for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

/**
 * Equilibrio chimico: totali delle componenti (mol) nel volume V (L).
 * Restituisce concentrazioni delle specie, quantità di solidi, pH, forza ionica.
 */
export function equilibrate(totals, V, guess = null) {
  const comps = COMP_IDS.filter(c => c === 'H' || Math.abs(totals[c] ?? 0) > 1e-15);
  const ci = new Map(comps.map((c, i) => [c, i]));
  const nc = comps.length;
  const T = comps.map(c => (totals[c] ?? 0) / V);
  const species = SP.map((s, k) => ({ k, s, a: Object.entries(s.comp).map(([c, v]) => [ci.get(c), v]) }))
    .filter(e => e.a.every(([i]) => i !== undefined));
  const solids = SO.map((s, k) => ({ k, s, a: Object.entries(s.comp).map(([c, v]) => [ci.get(c), v]) }))
    .filter(e => e.a.every(([i]) => i !== undefined));
  // stima iniziale
  const x = comps.map((c, i) => {
    if (guess?.x?.[c] !== undefined) return guess.x[c];
    if (c === 'H') return Math.log(1e-7);
    return Math.log(Math.max(T[i], 1e-14));
  });
  let lnG = new Array(nc).fill(0);
  let I = 0;
  let active = Object.keys(guess?.solids ?? {}).map(id => solids.find(e => e.s.id === id)).filter(Boolean);
  let m = new Array(active.length).fill(0);
  const conc = new Float64Array(species.length);
  const lnGs = new Float64Array(species.length);

  const speciesConc = () => {
    for (let j = 0; j < species.length; j++) {
      const e = species[j];
      let l = e.s.lnb - lnGs[j];
      for (const [i, v] of e.a) l += v * (x[i] + lnG[i]);
      conc[j] = Math.exp(Math.max(-700, Math.min(l, 50)));
    }
  };

  const newton = () => {
    for (let it = 0; it < 200; it++) {
      speciesConc();
      const n = nc + active.length;
      const R = new Array(n).fill(0);
      const J = Array.from({ length: n }, () => new Array(n).fill(0));
      const scale = new Array(nc).fill(0);
      for (let j = 0; j < species.length; j++) {
        const e = species[j], c = conc[j];
        for (const [i, v] of e.a) {
          R[i] += v * c;
          scale[i] += Math.abs(v) * c;
          for (const [i2, v2] of e.a) J[i][i2] += v * v2 * c;
        }
      }
      active.forEach((e, p) => {
        for (const [i, v] of e.a) {
          R[i] += v * m[p];
          scale[i] += Math.abs(v * m[p]);
          J[i][nc + p] += v;
          J[nc + p][i] += v;
          R[nc + p] += v * (x[i] + lnG[i]);
        }
        R[nc + p] -= e.s.lnK;
      });
      let done = true;
      for (let i = 0; i < nc; i++) {
        R[i] -= T[i];
        if (Math.abs(R[i]) > 1e-13 + 1e-10 * (scale[i] + Math.abs(T[i]))) done = false;
      }
      for (let p = 0; p < active.length; p++) if (Math.abs(R[nc + p]) > 1e-9) done = false;
      if (done) return true;
      const d = linsolve(J, R.map(v => -v));
      if (!d) return false;
      let f = 1;
      for (let i = 0; i < nc; i++) f = Math.min(f, 2.3 / Math.max(Math.abs(d[i]), 1e-300));
      for (let i = 0; i < nc; i++) x[i] += f * d[i];
      for (let p = 0; p < active.length; p++) m[p] += f * d[nc + p];
    }
    return false;
  };

  const solveWithActivity = () => {
    for (let outer = 0; outer < 40; outer++) {
      newton();
      let Inew = 0;
      for (let j = 0; j < species.length; j++) Inew += 0.5 * conc[j] * species[j].s.z * species[j].s.z;
      const Iold = I;
      I = 0.5 * (I + Math.min(Inew, 5));
      const sq = Math.sqrt(I);
      lnG = comps.map(c => lnGamma(COMPONENTS[c].z, sq, I));
      for (let j = 0; j < species.length; j++) lnGs[j] = lnGamma(species[j].s.z, sq, I);
      if (Math.abs(I - Iold) < 1e-9 + 1e-7 * I) break;
    }
    newton();
  };

  // insieme delle fasi solide
  for (let phase = 0; phase < 40; phase++) {
    solveWithActivity();
    const neg = m.map((v, p) => [v, p]).filter(([v]) => v < 0).sort((a, b) => a[0] - b[0]);
    if (neg.length) {
      const p = neg[0][1];
      active.splice(p, 1); m.splice(p, 1);
      continue;
    }
    let best = null;
    for (const e of solids) {
      if (active.includes(e)) continue;
      let si = -e.s.lnK;
      for (const [i, v] of e.a) si += v * (x[i] + lnG[i]);
      if (si > 1e-7 && (!best || si > best.si)) best = { e, si };
    }
    if (!best || active.length >= nc) break;
    active.push(best.e);
    m.push(0);
  }

  const out = { V, I, comps, species: {}, solids: {}, x: {} };
  comps.forEach((c, i) => { out.x[c] = x[i]; });
  species.forEach((e, j) => { out.species[e.s.id] = conc[j]; });
  active.forEach((e, p) => { if (m[p] > 0) out.solids[e.s.id] = m[p] * V; });
  const iH = ci.get('H');
  out.pH = -(x[iH] + lnG[iH]) / LN10;
  out.gammaH = Math.exp(lnG[iH]);
  out.activity = (id) => {
    const j = species.findIndex(e => e.s.id === id);
    return j < 0 ? 0 : conc[j] * Math.exp(lnGs[j]);
  };
  // indici di saturazione (per mostrare quanto si è vicini a precipitare)
  out.SI = {};
  for (const e of solids) {
    let si = -e.s.lnK;
    for (const [i, v] of e.a) si += v * (x[i] + lnG[i]);
    out.SI[e.s.id] = si / LN10;
  }
  return out;
}

// ---------------------------------------------------------------------------
// Colore: CIE 1931 (approssimazione analitica di Wyman, Sloan e Shirley, 2013)
// ---------------------------------------------------------------------------

const g2 = (l, mu, s1, s2) => { const t = (l - mu) / (l < mu ? s1 : s2); return Math.exp(-0.5 * t * t); };
export const cieX = (l) => 1.056 * g2(l, 599.8, 37.9, 31.0) + 0.362 * g2(l, 442.0, 16.0, 26.7) - 0.065 * g2(l, 501.1, 20.4, 26.2);
export const cieY = (l) => 0.821 * g2(l, 568.8, 46.9, 40.5) + 0.286 * g2(l, 530.9, 16.3, 31.1);
export const cieZ = (l) => 1.217 * g2(l, 437.0, 11.8, 36.0) + 0.681 * g2(l, 459.0, 26.0, 13.8);

/** Assorbanza A(λ) della soluzione (cammino ottico in cm). */
export function absorbance(speciesConc, lambda, path = 1) {
  let A = 0;
  const nu = 1e7 / lambda;
  for (const [id, c] of Object.entries(speciesConc)) {
    const s = SP[SPI.get(id)];
    if (!s?.bands || c <= 0) continue;
    for (const [lmax, eps, fwhm] of s.bands) {
      const d = nu - 1e7 / lmax;
      A += eps * c * path * Math.exp(-4 * Math.LN2 * d * d / (fwhm * fwhm));
    }
  }
  return A;
}

function xyzToSrgb(X, Y, Z) {
  const r = 3.2406 * X - 1.5372 * Y - 0.4986 * Z;
  const g = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
  const b = 0.0557 * X - 0.2040 * Y + 1.0570 * Z;
  return [r, g, b];
}
const gammaEnc = (v) => { const c = Math.max(0, Math.min(1, v)); return c <= 0.0031308 ? 12.92 * c : 1.055 * Math.pow(c, 1 / 2.4) - 0.055; };

let WHITE = null;
/** Colore percepito della soluzione (luce bianca trasmessa): [r, g, b] 0–255 e trasmittanza media. */
export function solutionColor(speciesConc, path = 1.5) {
  const integrate = (Tfun) => {
    let X = 0, Y = 0, Z = 0;
    for (let l = 380; l <= 780; l += 5) {
      const t = Tfun(l);
      X += t * cieX(l); Y += t * cieY(l); Z += t * cieZ(l);
    }
    return xyzToSrgb(X, Y, Z);
  };
  if (!WHITE) WHITE = integrate(() => 1);
  const lin = integrate((l) => Math.pow(10, -absorbance(speciesConc, l, path)));
  const rgb = lin.map((v, k) => Math.round(255 * gammaEnc(v / WHITE[k])));
  return rgb;
}

// ---------------------------------------------------------------------------
// Il becher
// ---------------------------------------------------------------------------

const EMPTY_STATE = () => ({
  V: 0, totals: {}, metals: {}, T: 298.15, gasH2: 0, gasCO2: 0, heat: 0, eq: null, nW: 0,
});

export class Beaker {
  constructor() { this.reset(); }

  reset() {
    this.st = EMPTY_STATE();
    this.log = [];
    this.history = [];
    this.lastEquation = null;
  }

  get volumeML() { return this.st.V * 1000; }

  /** Entalpia del contenuto (kJ), secondo le specie presenti all'equilibrio. */
  enthalpy(eq, nW, extraNominal = null) {
    let H = nW * DFH_WATER;
    if (eq) {
      for (const [id, c] of Object.entries(eq.species)) H += c * eq.V * SP[SPI.get(id)].dfh;
      for (const [id, n] of Object.entries(eq.solids)) H += n * SO[SOI.get(id)].dfh;
    }
    if (extraNominal) {
      for (const [id, n] of Object.entries(extraNominal.species)) H += n * SP[SPI.get(id)].dfh;
      for (const [id, n] of Object.entries(extraNominal.solids)) H += n * SO[SOI.get(id)].dfh;
    }
    return H;
  }

  /**
   * Aggiunge un reagente: `amount` in mL per le soluzioni, in grammi per i solidi.
   */
  add(reagentId, amount) {
    const r = REAGENTS.find(x => x.id === reagentId);
    if (!r || !(amount > 0)) return null;
    const st = this.st;
    const before = {
      eq: st.eq, metals: { ...st.metals }, nW: st.nW, gasH2: st.gasH2, gasCO2: st.gasCO2,
    };
    const nominal = { species: {}, solids: {}, metals: {} };
    let massAdded = 0;
    let dV = 0;
    if (r.kind === 'sol' || r.kind === 'ind') {
      dV = amount / 1000;
      const n = r.c * dV;
      for (const [id, k] of r.species) nominal.species[id] = (nominal.species[id] ?? 0) + k * n;
      massAdded = amount; // g (densità ≈ 1)
    } else if (r.kind === 'metal') {
      nominal.metals[r.metal] = amount / METALS[r.metal].M;
    } else if (r.kind === 'solid') {
      nominal.solids[r.solid] = amount / SO[SOI.get(r.solid)].M;
    }
    if (st.V + dV <= 0) return null;
    // totali delle componenti
    for (const [id, n] of Object.entries(nominal.species)) {
      for (const [c, a] of Object.entries(SP[SPI.get(id)].comp)) st.totals[c] = (st.totals[c] ?? 0) + a * n;
    }
    for (const [id, n] of Object.entries(nominal.solids)) {
      for (const [c, a] of Object.entries(SO[SOI.get(id)].comp)) st.totals[c] = (st.totals[c] ?? 0) + a * n;
    }
    for (const [M, n] of Object.entries(nominal.metals)) st.metals[M] = (st.metals[M] ?? 0) + n;
    const H0 = this.enthalpy(before.eq, before.nW, nominal);
    // l'acqua della soluzione aggiunta, in moli (per il bilancio del calore basta la variazione)
    const massOld = st.V * 1000;
    st.V += dV;
    // temperatura di miscela (i reagenti sono a 25 °C)
    const massNew = st.V * 1000;
    if (massNew > 0) st.T = (massOld * st.T + massAdded * 298.15) / Math.max(massNew, 1e-9);
    if (st.V <= 0) return null;

    let eq = equilibrate(st.totals, st.V, st.eq);
    // bilancio dell'acqua: formare una specie con w > 0 consuma acqua
    const waterChange = (eqA, eqB, nomB) => {
      let dw = 0;
      const nA = (id) => (eqA?.species[id] ?? 0) * (eqA?.V ?? 0);
      const sids = new Set([...Object.keys(eqA?.species ?? {}), ...Object.keys(eqB.species), ...Object.keys(nomB?.species ?? {})]);
      for (const id of sids) dw -= SP[SPI.get(id)].w * (nA(id) - (eqB.species[id] ?? 0) * eqB.V - (nomB?.species[id] ?? 0));
      const pids = new Set([...Object.keys(eqA?.solids ?? {}), ...Object.keys(eqB.solids), ...Object.keys(nomB?.solids ?? {})]);
      for (const id of pids) dw -= SO[SOI.get(id)].w * ((eqA?.solids[id] ?? 0) - (eqB.solids[id] ?? 0) - (nomB?.solids[id] ?? 0));
      return dw;
    };
    // redox con i metalli, poi fuga dei gas
    const redoxEvents = this.redox(eq);
    eq = this.lastEq ?? eq;
    const co2 = this.degas();
    eq = this.lastEq ?? eq;
    // prima dell'aggiunta (specie precedenti + nominali) → dopo
    const prevEqLike = { species: {}, solids: {}, V: 1 };
    for (const [id, c] of Object.entries(before.eq?.species ?? {})) prevEqLike.species[id] = c * before.eq.V;
    for (const [id, n] of Object.entries(before.eq?.solids ?? {})) prevEqLike.solids[id] = n;
    st.nW = before.nW + waterChange(eq, prevEqLike, nominal) + co2;
    st.eq = eq;
    // calore: legge di Hess (i gas escono portando la loro entalpia)
    const H1 = this.enthalpy(eq, st.nW) + co2 * DFH_CO2;
    const q = H0 - H1; // kJ ceduti alla soluzione
    st.heat += q;
    const mass = st.V * 1000;
    if (mass > 0) st.T += q * 1000 / (mass * CP_WATER);
    // equazione netta della trasformazione
    this.lastEquation = this.netEquation(prevEqLike, nominal, eq, before, co2);
    const entry = {
      reagent: r, amount, V: this.volumeML, pH: eq.pH, T: st.T, q, equation: this.lastEquation, redox: redoxEvents,
      co2, h2: st.gasH2 - before.gasH2,
    };
    this.log.push(entry);
    this.history.push({ V: this.volumeML, pH: eq.pH, T: st.T - 273.15, reagent: r.id });
    return entry;
  }

  /** Reazioni redox tra i metalli solidi e gli ossidanti in soluzione, fino a f.e.m. nulla. */
  redox(eq) {
    const st = this.st;
    this.lastEq = eq;
    const events = [];
    const lg = (v) => Math.log10(Math.max(v, 1e-14));
    const ionAct = (e, comp) => {
      const sp = SP.find(s => Object.keys(s.comp).length === 1 && s.comp[comp] === 1);
      return sp ? e.activity(sp.id) : 0;
    };
    const oxidants = () => {
      const list = [{ id: 'H', label: 'H⁺', n: 2, E: (e) => -NERNST * e.pH - H2_OVERPOTENTIAL, avail: () => true }];
      for (const [M, info] of Object.entries(METALS)) {
        list.push({ id: M, label: SP.find(s => s.comp[info.ion] === 1 && Object.keys(s.comp).length === 1)?.label ?? M, n: info.n, metal: M,
          E: (e) => info.E0 + NERNST / info.n * lg(ionAct(e, info.ion)), avail: () => (st.totals[info.ion] ?? 0) > 1e-12 });
      }
      list.push({ id: 'Fe3', label: 'Fe³⁺', n: 1, E: (e) => 0.771 + NERNST * (lg(ionAct(e, 'Fe3')) - lg(Math.max(ionAct(e, 'Fe2'), 1e-14))), avail: () => (st.totals.Fe3 ?? 0) > 1e-12 });
      return list;
    };
    for (let round = 0; round < 12; round++) {
      let best = null;
      const e0 = this.lastEq;
      for (const [M, amt] of Object.entries(st.metals)) {
        if (amt <= 1e-12) continue;
        const info = METALS[M];
        const Em = info.E0 + NERNST / info.n * lg(ionAct(e0, info.ion));
        for (const ox of oxidants()) {
          if (ox.metal === M || !ox.avail()) continue;
          const dE = ox.E(e0) - Em;
          if (dE > 0.002 && (!best || dE > best.dE)) best = { M, ox, dE };
        }
      }
      if (!best) break;
      const { M, ox } = best;
      const info = METALS[M];
      const oxInfo = ox.metal ? METALS[ox.metal] : null;
      const snapshot = JSON.parse(JSON.stringify({ totals: st.totals, metals: st.metals, gasH2: st.gasH2 }));
      const apply = (xi) => {
        st.totals = { ...snapshot.totals }; st.metals = { ...snapshot.metals }; st.gasH2 = snapshot.gasH2;
        st.metals[M] -= xi / info.n;
        st.totals[info.ion] = (st.totals[info.ion] ?? 0) + xi / info.n;
        if (ox.id === 'H') { st.totals.H -= xi; st.gasH2 += xi / 2; }
        else if (ox.id === 'Fe3') { st.totals.Fe3 -= xi; st.totals.Fe2 = (st.totals.Fe2 ?? 0) + xi; }
        else { st.totals[oxInfo.ion] -= xi / oxInfo.n; st.metals[ox.metal] = (st.metals[ox.metal] ?? 0) + xi / oxInfo.n; }
        const e = equilibrate(st.totals, st.V, this.lastEq);
        const Em = info.E0 + NERNST / info.n * lg(ionAct(e, info.ion));
        return { e, dE: ox.E(e) - Em };
      };
      let hi = snapshot.metals[M] * info.n;
      if (ox.id === 'Fe3') hi = Math.min(hi, snapshot.totals.Fe3);
      else if (ox.id !== 'H') hi = Math.min(hi, snapshot.totals[oxInfo.ion] * oxInfo.n);
      hi *= 0.999999;
      let res = apply(hi);
      let xi = hi;
      if (res.dE < 0) {
        let lo = 0;
        for (let k = 0; k < 40; k++) {
          const mid = 0.5 * (lo + hi);
          const r2 = apply(mid);
          if (r2.dE > 0) lo = mid; else hi = mid;
        }
        xi = lo;
        res = apply(xi);
      }
      this.lastEq = res.e;
      if (xi <= 1e-12) break;
      events.push({ metal: M, oxidant: ox.id, xi, E: best.dE });
    }
    return events;
  }

  /** Fuga di CO₂ oltre la solubilità di Henry (P_CO₂ = 1 atm). */
  degas() {
    const st = this.st;
    let escaped = 0;
    for (let k = 0; k < 20; k++) {
      const c = this.lastEq.species.H2CO3 ?? 0;
      if (c <= HENRY_CO2 * 1.0001) break;
      const d = (c - HENRY_CO2) * st.V;
      st.totals.CO3 -= d;
      st.totals.H -= 2 * d;
      escaped += d;
      this.lastEq = equilibrate(st.totals, st.V, this.lastEq);
    }
    st.gasCO2 += escaped;
    return escaped;
  }

  /** Equazione ionica netta dalla variazione delle quantità (prima + aggiunta → dopo). */
  netEquation(prev, nominal, eq, before, co2) {
    const d = new Map();
    const addD = (key, label, v) => { const x = d.get(key) ?? { label, v: 0 }; x.v += v; d.set(key, x); };
    const ids = new Set([...Object.keys(prev.species), ...Object.keys(nominal.species), ...Object.keys(eq.species)]);
    for (const id of ids) {
      const s = SP[SPI.get(id)];
      addD(id, s.label, (eq.species[id] ?? 0) * eq.V - (prev.species[id] ?? 0) - (nominal.species[id] ?? 0));
    }
    const pids = new Set([...Object.keys(prev.solids), ...Object.keys(nominal.solids), ...Object.keys(eq.solids)]);
    for (const id of pids) addD(id, `${SO[SOI.get(id)].label}(s)`, (eq.solids[id] ?? 0) - (prev.solids[id] ?? 0) - (nominal.solids[id] ?? 0));
    const mids = new Set([...Object.keys(before.metals), ...Object.keys(nominal.metals), ...Object.keys(this.st.metals)]);
    for (const M of mids) addD(`metal:${M}`, `${METALS[M].label}(s)`, (this.st.metals[M] ?? 0) - (before.metals[M] ?? 0) - (nominal.metals[M] ?? 0));
    const h2 = this.st.gasH2 - before.gasH2;
    if (h2 > 0) addD('H2', 'H₂(g)↑', h2);
    if (co2 > 0) addD('CO2', 'CO₂(g)↑', co2);
    // acqua: dal bilancio dell'idrogeno delle specie coinvolte (w)
    let dw = 0;
    for (const id of ids) dw -= SP[SPI.get(id)].w * (d.get(id)?.v ?? 0);
    for (const id of pids) dw -= SO[SOI.get(id)].w * (d.get(id)?.v ?? 0);
    dw += co2;
    if (Math.abs(dw) > 0) addD('H2O', 'H₂O', dw);
    const items = [...d.values()].filter(x => Math.abs(x.v) > 1e-12);
    if (!items.length) return null;
    const vmax = Math.max(...items.map(x => Math.abs(x.v)));
    // trasformazione trascurabile rispetto a quanto aggiunto (per esempio il solo spostamento di un equilibrio per diluizione)
    const added = [...Object.values(nominal.species), ...Object.values(nominal.solids), ...Object.values(nominal.metals)].reduce((a, v) => a + Math.abs(v), 0);
    if (vmax < 0.02 * added + 2e-6) return null;
    const kept = items.filter(x => Math.abs(x.v) > 0.08 * vmax);
    const reac = kept.filter(x => x.v < 0), prod = kept.filter(x => x.v > 0);
    if (!reac.length || !prod.length) return null;
    const unit = Math.min(...kept.map(x => Math.abs(x.v)));
    const coef = (v) => {
      const r = Math.abs(v) / unit;
      const k = Math.round(r * 2) / 2;
      if (Math.abs(r - Math.round(r)) < 0.12) return Math.round(r) === 1 ? '' : `${Math.round(r)} `;
      return `${k.toLocaleString('it-IT')} `;
    };
    const side = (arr) => arr.sort((a, b) => Math.abs(b.v) - Math.abs(a.v)).map(x => `${coef(x.v)}${x.label}`).join(' + ');
    return { text: `${side(reac)} → ${side(prod)}`, extent: unit };
  }

  /** Stato per l'interfaccia. */
  summary() {
    const st = this.st;
    const eq = st.eq;
    return {
      V: this.volumeML, T: st.T, heat: st.heat, pH: eq?.pH ?? 7, I: eq?.I ?? 0, gasH2: st.gasH2, gasCO2: st.gasCO2,
      species: eq ? Object.entries(eq.species).filter(([id]) => id !== 'H2O').map(([id, c]) => ({ id, label: SP[SPI.get(id)].label, c, z: SP[SPI.get(id)].z, estimated: SP[SPI.get(id)].dfhEstimated })).sort((a, b) => b.c - a.c) : [],
      solids: eq ? Object.entries(eq.solids).map(([id, n]) => ({ id, label: SO[SOI.get(id)].label, n, mass: n * SO[SOI.get(id)].M, color: SO[SOI.get(id)].color, Ksp: SO[SOI.get(id)].logK })) : [],
      metals: Object.entries(st.metals).filter(([, n]) => n > 1e-9).map(([M, n]) => ({ M, n, mass: n * METALS[M].M, color: METALS[M].color, label: METALS[M].label })),
      SI: eq?.SI ?? {},
      color: eq ? solutionColor(eq.species) : [235, 240, 245],
    };
  }
}
