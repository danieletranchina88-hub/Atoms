// Chimica in soluzione acquosa: equilibrio chimico completo del contenuto di un becher.
//
// Dati termodinamici: database MINTEQ v4 (U.S. EPA), lo stesso usato da PHREEQC (USGS), più poche specie
// documentate in aqueousExtra.js. Nessuna costante è stimata dal programma.
//
//  • Ogni specie si forma dalle componenti (ioni principali, H⁺, elettroni e⁻):
//      c_s = K_s(T) Π_i (γ_i c_i)^{a_si} / γ_s
//    con K_s(T) dalla legge di van 't Hoff: log K(T) = log K(298) − ΔrH/(R ln10) (1/T − 1/298,15).
//  • Bilanci di massa per ogni componente, compresi gli elettroni: le coppie redox (Fe²⁺/Fe³⁺, MnO₄⁻/Mn²⁺,
//    Cr₂O₇²⁻/Cr³⁺, Cu²⁺/Cu⁺, Co³⁺/Co²⁺, metalli/ioni) raggiungono l'equilibrio insieme. Il potenziale redox
//    della soluzione è pe = −log a(e⁻), Eh = 0,05916 pe (a 25 °C).
//  • Attività: Debye–Hückel esteso (parametri di WATEQ dal database) o Davies; A e B dalla costante
//    dielettrica e dalla densità dell'acqua a temperatura T. Specie neutre: log γ = 0,1 I.
//  • Fasi solide: precipitano quando l'indice di saturazione log(Q/K_sp) supera 0, si sciolgono se la loro
//    quantità diventa negativa (regola delle fasi di Gibbs).
//  • Gas: H₂ si libera quando pe < −pH − η/0,05916 (η = sovratensione dell'idrogeno); CO₂ quando la sua
//    pressione di equilibrio supera 1 atm (legge di Henry).
//  • Calore: legge di Hess con le entalpie di reazione del database; ΔT = q/(m c_p).
//  • Colore: legge di Lambert–Beer A(λ) = Σ ε(λ) c ℓ e funzioni colorimetriche CIE 1931.

import { DB_COMPONENTS, DB_SPECIES, DB_PHASES } from './aqueousDB.js';
import { EXTRA_COMPONENTS, EXTRA_SPECIES, EXTRA_PHASES, SPECTRA, PHASE_INFO, ORGANIC_ABBR } from './aqueousExtra.js';

const LN10 = Math.LN10;
const R_KJ = 8.314462618e-3;          // kJ/(mol K)
const T0 = 298.15;
const H2_OVERPOTENTIAL = 0; // Equilibrio termodinamico: nessuna sovratensione empirica o previsione cinetica.
const CP_WATER = 4.184;               // J/(g K)

export const COMPONENTS = [...DB_COMPONENTS, ...EXTRA_COMPONENTS];
const speciesMap = new Map();
for (const s of [...DB_SPECIES, ...EXTRA_SPECIES]) speciesMap.set(s.n, s);
export const SPECIES = [...speciesMap.values()];
export const PHASES = [...DB_PHASES, ...EXTRA_PHASES];
const SPI = new Map(SPECIES.map((s, i) => [s.n, i]));
const PHI = new Map(PHASES.map((p, i) => [p.n, i]));
const COMP_CHARGE = Object.fromEntries(COMPONENTS.map(c => [c, c === 'e-' ? -1 : (speciesMap.get(c)?.z ?? 0)]));
const REDOX_BASIS = new Set(['Fe+3', 'Mn+3', 'CrO4-2', 'Cu+2', 'Co+3']);

export const speciesByName = (n) => SPECIES[SPI.get(n)];
export const phaseByName = (n) => PHASES[PHI.get(n)];

// ---------------------------------------------------------------------------
// Nomi leggibili
// ---------------------------------------------------------------------------

const SUB = '₀₁₂₃₄₅₆₇₈₉';
const SUPD = '⁰¹²³⁴⁵⁶⁷⁸⁹';
export function prettyName(name) {
  if (name === 'H2CO3' || name === 'CO2') return 'CO₂(aq)';
  let s = name;
  for (const [k, v] of Object.entries(ORGANIC_ABBR)) s = s.split(k).join(v);
  s = s.replace(/^([A-Z][a-z]?\d*)\(([A-Za-z]+)\)(?=[+-]|$)/, '$1$2'); // H(Ac) → HAc, Na(Ac) → NaAc
  const m = s.match(/([+-])(\d*)$/);
  let charge = '';
  if (m) {
    s = s.slice(0, m.index);
    const n = m[2] ? parseInt(m[2], 10) : 1;
    charge = (n > 1 ? String(n).split('').map(d => SUPD[+d]).join('') : '') + (m[1] === '+' ? '⁺' : '⁻');
  }
  // pedici: cifre dopo una lettera o una parentesi chiusa
  s = s.replace(/([A-Za-z)])(\d+)/g, (_, a, d) => a + d.split('').map(x => SUB[+x]).join(''));
  return s + charge;
}
export const phaseLabel = (n) => PHASE_INFO[n]?.[0] ?? prettyName(PHASES[PHI.get(n)]?.f ?? n);
export const phaseName = (n) => PHASE_INFO[n]?.[1] ?? n;
export const phaseColor = (n) => PHASE_INFO[n]?.[2] ?? '#f2f2f2';

// ---------------------------------------------------------------------------
// Proprietà dell'acqua e costanti in funzione della temperatura
// ---------------------------------------------------------------------------

/** Costante dielettrica dell'acqua (Malmberg e Maryott, J. Res. NBS 1956), t in °C. */
export function waterDielectric(T) {
  const t = T - 273.15;
  return 87.740 - 0.40008 * t + 9.398e-4 * t * t - 1.410e-6 * t * t * t;
}
/** Densità dell'acqua (formula di Tilton e Taylor, 1937), g/cm³. */
export function waterDensity(T) {
  const t = T - 273.15;
  return 1 - (t - 3.9863) ** 2 * (t + 288.9414) / (508929.2 * (t + 68.12963));
}
/** Parametri di Debye–Hückel A (kg^½ mol^−½) e B (Å⁻¹ kg^½ mol^−½). */
export function debyeHuckel(T) {
  const eps = waterDielectric(T), rho = waterDensity(T);
  return { A: 1.82483e6 * Math.sqrt(rho) / Math.pow(eps * T, 1.5), B: 50.2916 * Math.sqrt(rho) / Math.sqrt(eps * T) };
}
/** log K alla temperatura T (van 't Hoff con ΔrH costante). */
export function logKAt(lk, dh, T) {
  if (dh === null || dh === undefined || T === T0) return lk;
  return lk - dh / (R_KJ * LN10) * (1 / T - 1 / T0);
}

function lnGamma(z, g, sqI, I, dh) {
  if (!z) return 0.1 * I * LN10;
  if (g && g[0] > 0) return (-dh.A * z * z * sqI / (1 + dh.B * g[0] * sqI) + g[1] * I) * LN10;
  return -dh.A * z * z * (sqI / (1 + sqI) - 0.3 * I) * LN10;
}

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

// ---------------------------------------------------------------------------
// Equilibrio
// ---------------------------------------------------------------------------

const PE_LIMIT = 60;

/**
 * @param totals quantità totali delle componenti (mol), comprese H⁺ ed e⁻ (possono essere negative)
 * @param V volume (L), T temperatura (K)
 * @param opts { guess, exclude: Set di fasi escluse, include: Set di fasi (se dato, solo queste) }
 */
export function equilibrate(totals, V, T = T0, opts = {}) {
  if (!(V > 0) || !Number.isFinite(V) || !Number.isFinite(T) || T < 273.15 || T > 373.15 || Object.values(totals).some(v=>!Number.isFinite(v))) throw new Error('Volume, temperatura o composizione fuori dal dominio del modello acquoso.');
  const present = COMPONENTS.filter(c => c !== 'H+' && c !== 'e-' && Math.abs(totals[c] ?? 0) > 1e-15);
  // gli elettroni sono una componente se ci sono coppie redox o metalli
  const hasRedox = totals['e-'] !== undefined;
  const comps = ['H+', ...(hasRedox ? ['e-'] : []), ...present];
  const ci = new Map(comps.map((c, i) => [c, i]));
  const nc = comps.length;
  const Tot = comps.map(c => (totals[c] ?? 0) / V);
  const dhp = debyeHuckel(T);
  const species = [];
  for (const s of SPECIES) {
    const a = [];
    let ok = true;
    for (const [c, v] of Object.entries(s.c)) {
      const i = ci.get(c);
      if (i === undefined) { ok = false; break; }
      a.push([i, v]);
    }
    if (ok) species.push({ s, a, lnK: logKAt(s.lk, s.dh, T) * LN10 });
  }
  const solids = [];
  for (const p of PHASES) {
    if (p.n === 'CO2(g)') continue;
    if (opts.include && !opts.include.has(p.n)) continue;
    if (opts.exclude?.has(p.n)) continue;
    const a = [];
    let ok = true;
    for (const [c, v] of Object.entries(p.c)) {
      const i = ci.get(c);
      if (i === undefined) { ok = false; break; }
      a.push([i, v]);
    }
    if (ok) solids.push({ p, a, lnK: logKAt(p.lk, p.dh, T) * LN10 });
  }
  const g = opts.guess;
  const x = comps.map((c, i) => {
    if (g?.x?.[c] !== undefined && Number.isFinite(g.x[c])) return g.x[c];
    if (c === 'H+') return Math.log(1e-7);
    if (c === 'e-') return 0;
    return Math.log(Math.max(Tot[i], 1e-14));
  });
  let lnG = new Array(nc).fill(0);
  const lnGs = new Float64Array(species.length);
  let I = g?.I ?? 0;
  let active = Object.keys(g?.solids ?? {}).map(n => solids.find(e => e.p.n === n)).filter(Boolean);
  let m = active.map(e => (g.solids[e.p.n] ?? 0) / V);
  // se gli elettroni non stanno in nessuna specie disciolta (per esempio zinco metallico in HCl), all'inizio
  // si trovano nel metallo: la fase metallica parte già attiva, altrimenti il metodo di Newton non converge
  if (hasRedox) {
    const ie = ci.get('e-');
    // capacità delle specie disciolte: quanti elettroni può trattenere al massimo ogni componente (per esempio il
    // rame in soluzione arriva al più a Cu⁺, un elettrone per atomo). Se gli elettroni totali la superano, una
    // parte deve restare nel metallo e la fase metallica parte già attiva
    let capacity = 0;
    for (const [c, i] of ci) {
      if (c === 'H+' || c === 'e-') continue;
      let best = 0;
      for (const e of species) {
        const others = e.a.filter(([k]) => k !== ie && ci.get('H+') !== k);
        if (others.length !== 1 || others[0][0] !== i) continue;
        const ve = e.a.find(([k]) => k === ie)?.[1] ?? 0;
        best = Math.max(best, ve / others[0][1]);
      }
      capacity += best * Math.max(Tot[i], 0);
    }
    if (Tot[ie] > capacity * (1 + 1e-9) + 1e-15) {
      const metals = solids.filter(e => !active.includes(e) && e.a.some(([i, v]) => i === ie && v > 0));
      for (const e of metals) {
        const nu = e.a.find(([i]) => i === ie)[1];
        active.push(e);
        m.push((Tot[ie] - capacity) / nu / metals.length);
      }
    }
  }
  const conc = new Float64Array(species.length);

  const setGammas = () => {
    const sq = Math.sqrt(I);
    lnG = comps.map(c => (c === 'e-' ? 0 : lnGamma(COMP_CHARGE[c], speciesMap.get(c)?.g, sq, I, dhp)));
    for (let j = 0; j < species.length; j++) lnGs[j] = lnGamma(species[j].s.z, species[j].s.g, sq, I, dhp);
  };
  const speciesConc = () => {
    for (let j = 0; j < species.length; j++) {
      const e = species[j];
      let l = e.lnK - lnGs[j];
      for (const [i, v] of e.a) l += v * (x[i] + lnG[i]);
      conc[j] = Math.exp(Math.max(-700, Math.min(l, 60)));
    }
  };
  let converged = false;
  const newton = () => {
    for (let it = 0; it < 300; it++) {
      speciesConc();
      const n = nc + active.length;
      const Rv = new Array(n).fill(0);
      const J = Array.from({ length: n }, () => new Array(n).fill(0));
      const scale = new Array(nc).fill(0);
      for (let j = 0; j < species.length; j++) {
        const e = species[j], c = conc[j];
        for (const [i, v] of e.a) {
          Rv[i] += v * c;
          scale[i] += Math.abs(v) * c;
          for (const [i2, v2] of e.a) J[i][i2] += v * v2 * c;
        }
      }
      active.forEach((e, p) => {
        for (const [i, v] of e.a) {
          Rv[i] += v * m[p];
          scale[i] += Math.abs(v * m[p]);
          J[i][nc + p] += v;
          J[nc + p][i] += v;
          Rv[nc + p] += v * (x[i] + lnG[i]);
        }
        Rv[nc + p] -= e.lnK;
      });
      let done = true;
      for (let i = 0; i < nc; i++) {
        Rv[i] -= Tot[i];
        if (Math.abs(Rv[i]) > 1e-15 + 1e-10 * (scale[i] + Math.abs(Tot[i]))) done = false;
      }
      for (let p = 0; p < active.length; p++) if (Math.abs(Rv[nc + p]) > 1e-9) done = false;
      if (done) return true;
      // regolarizzazione minima: una componente senza specie (per esempio e⁻ senza forma ridotta) non blocca il sistema
      for (let i = 0; i < nc; i++) J[i][i] += 1e-30;
      let d = linsolve(J, Rv.map(v => -v));
      // jacobiana singolare (per esempio solo Mn²⁺, senza Mn³⁺ né MnO₄⁻: la direzione Mn³⁺ ↔ e⁻ è indeterminata):
      // si riprova con uno smorzamento di Levenberg–Marquardt relativo alla diagonale
      if (!d || d.some(v => !Number.isFinite(v))) {
        for (let i = 0; i < nc; i++) J[i][i] *= 1 + 1e-10;
        d = linsolve(J, Rv.map(v => -v));
      }
      if (!d || d.some(v => !Number.isFinite(v))) return false;
      // deriva senza fine lungo una direzione indeterminata: si lascia decidere alla regola delle fasi
      if (x.some(v => Math.abs(v) > 460)) return false;
      let f = 1;
      for (let i = 0; i < nc; i++) f = Math.min(f, 2.3 / Math.max(Math.abs(d[i]), 1e-300));
      for (let i = 0; i < nc; i++) x[i] += f * d[i];
      for (let p = 0; p < active.length; p++) m[p] += f * d[nc + p];
      // senza coppia redox (solo Fe²⁺, solo Mn²⁺…) il pe non ha un valore determinato e scivolerebbe all'infinito:
      // lo si ferma a |pe| = 60, ben oltre il campo di qualunque coppia del database (Mg²⁺/Mg: pe ≈ −40)
      if (hasRedox) { const ie = ci.get('e-'); x[ie] = Math.max(-PE_LIMIT * LN10, Math.min(PE_LIMIT * LN10, x[ie])); }
    }
    return false;
  };
  const solveWithActivity = () => {
    for (let outer = 0; outer < 60; outer++) {
      setGammas();
      converged = newton();
      let Inew = 0;
      for (let j = 0; j < species.length; j++) Inew += 0.5 * conc[j] * species[j].s.z * species[j].s.z;
      const Iold = I;
      I = 0.5 * (I + Math.min(Inew, 6));
      if (Math.abs(I - Iold) < 1e-10 + 1e-8 * I) break;
    }
    setGammas();
    converged = newton();
  };
  for (let phase = 0; phase < 60; phase++) {
    solveWithActivity();
    const neg = m.map((v, p) => [v, p]).filter(([v]) => v < 0).sort((a, b) => a[0] - b[0]);
    if (neg.length) { active.splice(neg[0][1], 1); m.splice(neg[0][1], 1); continue; }
    let best = null;
    for (const e of solids) {
      if (active.includes(e)) continue;
      let si = -e.lnK;
      for (const [i, v] of e.a) si += v * (x[i] + lnG[i]);
      if (si > 1e-7 && (!best || si > best.si)) best = { e, si };
    }
    if (!best || active.length >= nc) break;
    active.push(best.e);
    m.push(0);
  }

  const out = { V, T, I, converged, comps, species: {}, solids: {}, x: {}, SI: {}, gamma: {} };
  comps.forEach((c, i) => { out.x[c] = x[i]; });
  species.forEach((e, j) => { out.species[e.s.n] = conc[j]; out.gamma[e.s.n] = Math.exp(lnGs[j]); });
  active.forEach((e, p) => { if (m[p] > 0) out.solids[e.p.n] = m[p] * V; });
  const iH = ci.get('H+');
  out.pH = -(x[iH] + lnG[iH]) / LN10;
  out.pe = hasRedox ? -x[ci.get('e-')] / LN10 : null;
  out.activity = (n) => (out.species[n] ?? 0) * (out.gamma[n] ?? 1);
  // Il pe ha senso fisico (un elettrodo di platino misura un potenziale stabile) solo se una coppia redox è
  // "tamponata": la stessa componente è presente in almeno due stati di ossidazione, ognuno ≥ 1 µmol/L
  // (o come metallo). Senza coppia il pe calcolato è solo un numero formale e non viene mostrato.
  out.poised = false;
  if (hasRedox) {
    const ie = ci.get('e-'), iH2 = ci.get('H+');
    const levels = new Map();
    const addLevel = (a, c) => {
      const red = a.filter(([i]) => i !== ie && i !== iH2);
      if (red.length !== 1) return;
      const [ic, vc] = red[0];
      const lev = Math.round(100 * (a.find(([i]) => i === ie)?.[1] ?? 0) / vc);
      if (!levels.has(ic)) levels.set(ic, new Map());
      const L = levels.get(ic);
      L.set(lev, (L.get(lev) ?? 0) + vc * c);
    };
    for (let j = 0; j < species.length; j++) addLevel(species[j].a, conc[j]);
    for (const c of present) if (!species.some(e => e.s.n === c)) addLevel([[ci.get(c), 1]], Math.exp(x[ci.get(c)]));
    active.forEach((e, p) => addLevel(e.a, Math.max(m[p], 0)));
    for (const L of levels.values()) if ([...L.values()].filter(v => v >= 1e-6).length >= 2) out.poised = true;
  }
  for (const e of [...solids, ...PHASES.filter(p => p.n === 'CO2(g)').map(p => ({ p, a: Object.entries(p.c).map(([c, v]) => [ci.get(c), v]), lnK: logKAt(p.lk, p.dh, T) * LN10 }))]) {
    if (e.a.some(([i]) => i === undefined)) continue;
    let si = -e.lnK;
    for (const [i, v] of e.a) si += v * (x[i] + lnG[i]);
    out.SI[e.p.n] = si / LN10;
  }
  out.residuals = {};
  let balancesOK = true, maxRelative = 0;
  for (const c of comps) {
    let found = 0, scale = Math.abs(totals[c] ?? 0);
    for (const e of species) { const n = (e.s.c[c] ?? 0)*out.species[e.s.n]*V; found+=n; scale+=Math.abs(n); }
    for (const [id,n] of Object.entries(out.solids)) { const v=(phaseByName(id).c[c]??0)*n;found+=v;scale+=Math.abs(v); }
    const r=found-(totals[c]??0);out.residuals[c]=r;
    if(Math.abs(r)>1e-12+1e-8*scale)balancesOK=false;
    maxRelative=Math.max(maxRelative,Math.abs(r)/(1e-12+scale));
  }
  const physicalI=species.reduce((v,e)=>v+.5*e.s.z**2*out.species[e.s.n],0);
  out.chargeResidual=species.reduce((v,e)=>v+e.s.z*out.species[e.s.n]*V,0)-comps.reduce((v,c)=>v+COMP_CHARGE[c]*(totals[c]??0),0);
  out.maxRelativeResidual=maxRelative;
  out.converged=balancesOK && Math.abs(physicalI-I)<1e-7*(1+physicalI)
    && Object.entries(out.SI).every(([id,si])=>id==='CO2(g)' || (out.solids[id] ? Math.abs(si)<1e-5 : si<1e-5));
  out.I=physicalI;
  if(!out.converged && opts.retry!==false && solids.length<=20){
    for(const candidate of solids){
      const alternative=equilibrate(totals,V,T,{...opts,retry:false,guess:{solids:{[candidate.p.n]:1e-12}}});
      if(alternative.converged)return alternative;
    }
  }
  out.warnings=[];
  if(!out.converged)out.warnings.push('Equilibrio non convergente: risultato non utilizzabile.');
  if(physicalI>.5)out.warnings.push('I > 0,5 mol/L: attività fuori dal dominio diluito. Diluire.');
  return out;
}

// ---------------------------------------------------------------------------
// Colore: CIE 1931 (approssimazione analitica di Wyman, Sloan e Shirley, JCGT 2013)
// ---------------------------------------------------------------------------

const g2 = (l, mu, s1, s2) => { const t = (l - mu) / (l < mu ? s1 : s2); return Math.exp(-0.5 * t * t); };
const cieX = (l) => 1.056 * g2(l, 599.8, 37.9, 31.0) + 0.362 * g2(l, 442.0, 16.0, 26.7) - 0.065 * g2(l, 501.1, 20.4, 26.2);
const cieY = (l) => 0.821 * g2(l, 568.8, 46.9, 40.5) + 0.286 * g2(l, 530.9, 16.3, 31.1);
const cieZ = (l) => 1.217 * g2(l, 437.0, 11.8, 36.0) + 0.681 * g2(l, 459.0, 26.0, 13.8);

export function absorbance(speciesConc, lambda, path = 1) {
  let A = 0;
  const nu = 1e7 / lambda;
  for (const [n, c] of Object.entries(speciesConc)) {
    const bands = SPECTRA[n];
    if (!bands || c <= 0) continue;
    for (const [lmax, eps, fwhm] of bands) {
      const d = nu - 1e7 / lmax;
      A += eps * c * path * Math.exp(-4 * Math.LN2 * d * d / (fwhm * fwhm));
    }
  }
  return A;
}

function xyzToSrgb(X, Y, Z) {
  return [3.2406 * X - 1.5372 * Y - 0.4986 * Z, -0.9689 * X + 1.8758 * Y + 0.0415 * Z, 0.0557 * X - 0.2040 * Y + 1.0570 * Z];
}
const gammaEnc = (v) => { const c = Math.max(0, Math.min(1, v)); return c <= 0.0031308 ? 12.92 * c : 1.055 * Math.pow(c, 1 / 2.4) - 0.055; };
let WHITE = null;
export function solutionColor(speciesConc, path = 1.5) {
  const integrate = (Tf) => {
    let X = 0, Y = 0, Z = 0;
    for (let l = 380; l <= 780; l += 5) { const t = Tf(l); X += t * cieX(l); Y += t * cieY(l); Z += t * cieZ(l); }
    return xyzToSrgb(X, Y, Z);
  };
  if (!WHITE) WHITE = integrate(() => 1);
  const lin = integrate((l) => Math.pow(10, -absorbance(speciesConc, l, path)));
  return lin.map((v, k) => Math.round(255 * gammaEnc(v / WHITE[k])));
}

// ---------------------------------------------------------------------------
// Il becher
// ---------------------------------------------------------------------------

/** Entalpia (kJ) di formazione dalle componenti di tutto ciò che è nel becher (specie + solidi). */
function enthalpy(eq) {
  let H = 0;
  if (!eq) return H;
  for (const [n, c] of Object.entries(eq.species)) { const s = speciesMap.get(n); if (s?.dh) H += c * eq.V * s.dh; }
  for (const [n, amt] of Object.entries(eq.solids)) { const p = PHASES[PHI.get(n)]; if (p?.dh) H -= amt * p.dh; }
  return H;
}

/**
 * Una "ricetta" descrive che cosa si versa:
 *  { label, kind: 'solution' | 'solid' | 'metal' | 'phase' | 'water', species: [[nome, coeff], ...], phase, c (mol/L), M (g/mol), dsol (kJ/mol) }
 */
/** Composizione nelle componenti del database (più l'acqua) di una voce del registro delle variazioni. */
function compositionOf(key) {
  if (key === 'H2O') return { c: {}, w: 1 };
  if (key === 'H2') return { c: { 'H+': 2, 'e-': 2 }, w: 0 };
  if (key === 'CO2g') { const p = PHASES[PHI.get('CO2(g)')]; return { c: p.c, w: p.w ?? 0 }; }
  if (key.startsWith('p:')) { const p = PHASES[PHI.get(key.slice(2))]; return p ? { c: p.c, w: p.w ?? 0 } : null; }
  const sp = speciesMap.get(key);
  return sp ? { c: sp.c, w: speciesWater(key) } : null;
}

// H₂CO₃* del database è mostrato come CO₂(aq) (la forma prevalente): CO₂(aq) = H₂CO₃ − H₂O
function speciesWater(name) { return name === 'H2CO3' ? -1 : speciesMap.get(name)?.w ?? 0; }

/**
 * Cerca coefficienti interi (≤ 12) proporzionali alle variazioni in moli, con scarto ≤ 4 %, e accetta l'equazione
 * solo se bilancia esattamente ogni componente (quindi atomi e cariche) e l'acqua.
 */
function balancedEquation(entries) {
  const rows = entries.map(([key, x]) => ({ ...x, comp: compositionOf(key) }));
  if (rows.length < 2 || rows.some(r => !r.comp) || !rows.some(r => r.v < 0) || !rows.some(r => r.v > 0)) return null;
  const unit = Math.min(...rows.map(r => Math.abs(r.v)));
  for (let k = 1; k <= 12; k++) {
    const coef = rows.map(r => Math.round(Math.abs(r.v) / unit * k));
    if (coef.some(c => c < 1 || c > 12) || rows.some((r, i) => Math.abs(Math.abs(r.v) / unit * k - coef[i]) > 0.04 * coef[i])) continue;
    const bal = {};
    rows.forEach((r, i) => {
      const n = Math.sign(r.v) * coef[i];
      for (const [comp, v] of Object.entries(r.comp.c)) bal[comp] = (bal[comp] ?? 0) + n * v;
      bal.H2O = (bal.H2O ?? 0) + n * r.comp.w;
    });
    if (Object.values(bal).some(v => Math.abs(v) > 1e-9)) continue;
    const g = coef.reduce((a, b) => { while (b) [a, b] = [b, a % b]; return a; });
    const side = (sign) => rows.map((r, i) => ({ r, c: coef[i] / g })).filter(x => Math.sign(x.r.v) === sign)
      .map(x => `${x.c === 1 ? '' : `${x.c} `}${x.r.label}`).join(' + ');
    return `${side(-1)} → ${side(1)}`;
  }
  return null;
}

export class Beaker {
  constructor() { this.reset(); }

  reset() {
    this.st = { V: 0, totals: {}, T: T0, heat: 0, gasH2: 0, gasCO2: 0, eq: null, thermostat: false, Tset: T0 };
    this.log = [];
    this.history = [];
    this.undoStack = [];
  }

  get volumeML() { return this.st.V * 1000; }

  setThermostat(on, T) {
    if(T!==undefined && (!Number.isFinite(T) || T<273.15 || T>368.15))throw new Error('Temperatura del bagno fuori intervallo 0–95 °C.');
    const saved={...this.st,totals:{...this.st.totals}},lastGas=this.lastGas;
    try {
      this.st.thermostat = on;
      if (T !== undefined) this.st.Tset = T;
      if (on) { this.st.T = this.st.Tset; this.reequilibrate(); }
      if(this.st.eq && !this.st.eq.converged)throw new Error('Equilibrio non convergente: temperatura ripristinata.');
    }catch(error){this.st=saved;this.lastGas=lastGas;throw error;}
  }

  reequilibrate() {
    if (this.st.V <= 0) return;
    this.st.eq = this.solve(this.st.eq);
  }

  solve(guess) {
    const st = this.st;
    let eq = equilibrate(st.totals, st.V, st.T, { guess });
    const h2 = this.evolveH2(eq); eq = h2.eq;
    const co2 = this.degasCO2(eq); eq = co2.eq;
    this.lastGas = { h2: h2.n, co2: co2.n };
    return eq;
  }

  /** Sviluppo di idrogeno: finché pe < −pH − η/0,05916, gli elettroni in eccesso riducono H⁺ a H₂ che se ne va. */
  evolveH2(eq) {
    const st = this.st;
    if (eq.pe === null) return { eq, n: 0 };
    // solo un metallo (Zn, Fe, Mg…) ha abbastanza potere riducente per sviluppare H₂ in tempi di laboratorio;
    // senza questo controllo una soluzione di sola Fe²⁺ (pe formalmente −∞ perché manca Fe³⁺) "svilupperebbe" H₂
    if (!Object.keys(eq.solids).some(n => /metal$/.test(n))) return { eq, n: 0 };
    const thr = (e) => e.pe + e.pH + H2_OVERPOTENTIAL / (0.05916 * st.T / T0);
    if (thr(eq) >= 0) return { eq, n: 0 };
    const e0 = st.totals['e-'] ?? 0;
    const h0 = st.totals['H+'] ?? 0;
    const apply = (n) => { st.totals['e-'] = e0 - 2 * n; st.totals['H+'] = h0 - 2 * n; return equilibrate(st.totals, st.V, st.T, { guess: eq }); };
    let lo = 0, hi = Math.max(0, e0 / 2);
    let best = apply(hi);
    if (thr(best) < 0) { st.gasH2 += hi; return { eq: best, n: hi }; }
    for (let k = 0; k < 50; k++) {
      const mid = 0.5 * (lo + hi);
      const e = apply(mid);
      if (thr(e) < 0) lo = mid; else { hi = mid; best = e; }
    }
    best = apply(hi);
    st.gasH2 += hi;
    return { eq: best, n: hi };
  }

  /** CO₂ oltre 1 atm di pressione di equilibrio lascia la soluzione (bolle). */
  degasCO2(eq) {
    const st = this.st;
    const p = PHASES[PHI.get('CO2(g)')];
    let n = 0;
    for (let k = 0; k < 30; k++) {
      const si = eq.SI['CO2(g)'];
      if (si === undefined || si <= 1e-4) break;
      // quantità da togliere: si stima dalla CO₂ disciolta (specie con CO₃ e 2 H⁺)
      const dissolved = (eq.species.CO2 ?? eq.species.H2CO3 ?? 0) * st.V;
      const d = dissolved * (1 - Math.pow(10, -si));
      if (d <= 0) break;
      for (const [c, v] of Object.entries(p.c)) st.totals[c] = (st.totals[c] ?? 0) - v * d;
      n += d;
      eq = equilibrate(st.totals, st.V, st.T, { guess: eq });
    }
    st.gasCO2 += n;
    return { eq, n };
  }

  /** Versa: `amount` in mL per le soluzioni, in grammi per i solidi. Restituisce la voce di registro. */
  add(recipe, amount) {
    if (!Number.isFinite(amount) || amount<=0) return null;
    const saved={...this.st,totals:{...this.st.totals}};
    const logLength=this.log.length,historyLength=this.history.length,lastGas=this.lastGas;
    try {
      const result=this._add(recipe,amount);
      if(result && !this.st.eq.converged)throw new Error('Equilibrio non convergente: aggiunta annullata. Prova una miscela più diluita.');
      if(result)this.undoStack.push({st:saved,logLength,historyLength,lastGas});
      return result;
    } catch(error){this.st=saved;this.log.length=logLength;this.history.length=historyLength;this.lastGas=lastGas;throw error;}
  }

  undo() {
    const previous=this.undoStack.pop();if(!previous)return false;
    this.st=previous.st;this.log.length=previous.logLength;this.history.length=previous.historyLength;this.lastGas=previous.lastGas;return true;
  }

  _add(recipe, amount) {
    const st = this.st;
    const prevEq = st.eq;
    const nominal = { species: {}, solids: {} };
    let dV = 0, massAdded = 0, dsolHeat = 0;
    if (recipe.kind === 'solution' || recipe.kind === 'water') {
      dV = amount / 1000;
      massAdded = amount;
      const n = (recipe.c ?? 0) * dV;
      for (const [s, k] of recipe.species ?? []) nominal.species[s] = (nominal.species[s] ?? 0) + k * n;
    } else if (recipe.kind === 'solid') {
      const n = amount / recipe.M;
      for (const [s, k] of recipe.species) nominal.species[s] = (nominal.species[s] ?? 0) + k * n;
      // il solido ha entalpia più bassa degli ioni sciolti di ΔsolH: H(solido) = H(ioni) − ΔsolH
      if (recipe.dsol !== null && recipe.dsol !== undefined) dsolHeat = -recipe.dsol * n;
    } else if (recipe.kind === 'phase' || recipe.kind === 'metal') {
      nominal.solids[recipe.phase] = amount / recipe.M;
    }
    if (st.V + dV <= 0) return null;
    for (const [s, n] of Object.entries(nominal.species)) {
      const sp = speciesMap.get(s);
      for (const [c, v] of Object.entries(sp.c)) st.totals[c] = (st.totals[c] ?? 0) + v * n;
    }
    for (const [pn, n] of Object.entries(nominal.solids)) {
      const p = PHASES[PHI.get(pn)];
      for (const [c, v] of Object.entries(p.c)) st.totals[c] = (st.totals[c] ?? 0) + v * n;
    }
    const unsupported = message => { const e=new Error(message);e.code='UNSUPPORTED_CHEMISTRY';throw e; };
    if((st.totals['Cu+2']??0)>1e-10 && (st.totals['I-']??0)>1e-10)
      unsupported('Cu/ioduro richiede anche I₂ e la sua redox, assente dal database: combinazione non calcolabile.');
    const hasMetal=[...Object.keys(nominal.solids),...Object.keys(prevEq?.solids??{})].some(n=>/metal/.test(n));
    if(hasMetal && (st.totals['NO3-']??0)>1e-10 && (st.totals['H+']??0)>1e-7)
      unsupported('Metallo in nitrato acido: riduzione a NOₓ non modellata. Combinazione non calcolabile.');
    if (REDOX_BASIS && Object.keys(st.totals).some(c => REDOX_BASIS.has(c))) st.totals['e-'] = st.totals['e-'] ?? 0;
    // entalpia prima della reazione: contenuto precedente + ciò che si aggiunge (a 25 °C)
    let Hnom = dsolHeat;
    for (const [s, n] of Object.entries(nominal.species)) { const sp = speciesMap.get(s); if (sp.dh) Hnom += n * sp.dh; }
    for (const [pn, n] of Object.entries(nominal.solids)) { const p = PHASES[PHI.get(pn)]; if (p.dh) Hnom -= n * p.dh; }
    const H0 = enthalpy(prevEq) + Hnom;
    const massOld = st.V * 1000;
    st.V += dV;
    const mass = st.V * 1000;
    if (!st.thermostat && mass > 0) st.T = (massOld * st.T + massAdded * T0) / mass;
    const gasBefore = { h2: st.gasH2, co2: st.gasCO2 };
    // equilibrio e bilancio termico (si itera perché le costanti dipendono dalla temperatura)
    let eq = null, q = 0;
    const Tstart = st.T;
    const baseTotals = { ...st.totals };
    for (let it = 0; it < (st.thermostat ? 1 : 5); it++) {
      // i gas che se ne vanno cambiano i totali: a ogni iterazione si riparte dal miscuglio appena versato
      st.totals = { ...baseTotals }; st.gasH2 = gasBefore.h2; st.gasCO2 = gasBefore.co2;
      eq = this.solve(eq ?? prevEq);
      const pCO2 = PHASES[PHI.get('CO2(g)')];
      const H1 = enthalpy(eq) - (st.gasCO2 - gasBefore.co2) * (pCO2.dh ?? 0);
      q = H0 - H1;
      if (st.thermostat) break;
      const Tnew = Tstart + q * 1000 / (mass * CP_WATER);
      const changed = Math.abs(Tnew - st.T) > 1e-3;
      st.T = Tnew;
      if (!changed) break;
    }
    st.eq = eq;
    st.heat += q;
    const entry = {
      recipe, amount, V: this.volumeML, pH: eq.pH, pe: eq.poised ? eq.pe : null, T: st.T, q,
      h2: st.gasH2 - gasBefore.h2, co2: st.gasCO2 - gasBefore.co2,
      equation: this.netEquation(prevEq, nominal, eq, st.gasH2 - gasBefore.h2, st.gasCO2 - gasBefore.co2),
      heatKnown: this.heatKnown(eq, nominal) && (recipe.kind !== 'solid' || recipe.dsol !== null && recipe.dsol !== undefined),
    };
    this.log.push(entry);
    this.history.push({ V: this.volumeML, pH: eq.pH, T: st.T - 273.15, label: recipe.label });
    return entry;
  }

  /** Vero se tutte le specie e le fasi principali hanno l'entalpia nel database. */
  heatKnown(eq, nominal) {
    const names = [...Object.keys(nominal.species), ...Object.entries(eq.species).filter(([, c]) => c * eq.V > 1e-6).map(([n]) => n)];
    return names.every(n => speciesMap.get(n)?.dh !== null) && Object.keys(eq.solids).every(n => PHASES[PHI.get(n)].dh !== null);
  }

  netEquation(prevEq, nominal, eq, h2, co2) {
    const d = new Map();
    const add = (key, label, v) => { const x = d.get(key) ?? { label, v: 0, w: 0 }; x.v += v; d.set(key, x); };
    const names = new Set([...Object.keys(prevEq?.species ?? {}), ...Object.keys(nominal.species), ...Object.keys(eq.species)]);
    for (const n of names) {
      const before = (prevEq?.species[n] ?? 0) * (prevEq?.V ?? 0) + (nominal.species[n] ?? 0);
      add(n, prettyName(n), (eq.species[n] ?? 0) * eq.V - before);
    }
    const pnames = new Set([...Object.keys(prevEq?.solids ?? {}), ...Object.keys(nominal.solids), ...Object.keys(eq.solids)]);
    for (const n of pnames) add(`p:${n}`, `${phaseLabel(n)}(s)`, (eq.solids[n] ?? 0) - (prevEq?.solids[n] ?? 0) - (nominal.solids[n] ?? 0));
    if (h2 > 0) add('H2', 'H₂(g)↑', h2);
    if (co2 > 0) add('CO2g', 'CO₂(g)↑', co2);
    // acqua prodotta: −Σ w Δn (formare una specie con w > 0 consuma acqua)
    let dw = 0;
    for (const n of names) dw -= speciesWater(n) * d.get(n).v;
    for (const n of pnames) dw -= (PHASES[PHI.get(n)]?.w ?? 0) * d.get(`p:${n}`).v;
    if (co2 > 0) dw -= (PHASES[PHI.get('CO2(g)')].w ?? 0) * co2;
    if (Math.abs(dw) > 0) add('H2O', 'H₂O', dw);
    const items = [...d.values()].filter(x => Math.abs(x.v) > 1e-13);
    if (!items.length) return null;
    const vmax = Math.max(...items.map(x => Math.abs(x.v)));
    const added = [...Object.values(nominal.species), ...Object.values(nominal.solids)].reduce((a, v) => a + Math.abs(v), 0);
    if (vmax < 0.02 * added + 2e-7) return null;
    // una diluizione sposta un poco tutti gli equilibri: si mostra solo se si scioglie o si forma un solido
    if (added === 0 && ![...pnames].some(n => Math.abs(d.get(`p:${n}`).v) > 1e-6)) return null;
    const kept = items.filter(x => Math.abs(x.v) > 0.08 * vmax);
    const reac = kept.filter(x => x.v < 0), prod = kept.filter(x => x.v > 0);
    if (!reac.length || !prod.length) return null;
    const detail = [...kept].sort((a, b) => Math.abs(b.v) - Math.abs(a.v))
      .map(x => `${x.label}: ${x.v > 0 ? '+' : ''}${(x.v * 1000).toLocaleString('it-IT', { maximumSignificantDigits: 4 })} mmol`).join('; ');
    const changes = items.map(x => ({ label: x.label, mol: x.v }));
    // equazione netta solo se coefficienti interi piccoli bilanciano esattamente componenti, cariche e acqua;
    // altrimenti (molti complessi, reazioni parallele) si mostrano le variazioni in mmol
    for (const frac of [0.08, 0.03]) {
      const eqn = balancedEquation([...d.entries()].filter(([, x]) => Math.abs(x.v) > frac * vmax));
      if (eqn) return { text: eqn, balanced: true, detail, changes };
    }
    return { text: detail, balanced: false, detail, changes };
  }

  summary() {
    const st = this.st, eq = st.eq;
    const species = eq ? Object.entries(eq.species).map(([n, c]) => ({ n, label: prettyName(n), c, z: speciesMap.get(n)?.z ?? 0, a: c * (eq.gamma[n] ?? 1), src: speciesMap.get(n)?.src })).sort((a, b) => b.c - a.c) : [];
    const solids = eq ? Object.entries(eq.solids).map(([n, amt]) => {
      const p = PHASES[PHI.get(n)];
      return { n, label: phaseLabel(n), name: phaseName(n), amt, mass: amt * (molarMassOfPhase(p) ?? 0), color: phaseColor(n), logKsp: logKAt(p.lk, p.dh, st.T), src: p.src, metal: /metal$/.test(n) };
    }) : [];
    return {
      V: this.volumeML, T: st.T, heat: st.heat, pH: eq?.pH ?? 7, pe: eq?.poised ? eq.pe : null, I: eq?.I ?? 0,
      Eh: eq?.poised ? 0.05916 * (st.T / T0) * eq.pe : null,
      gasH2: st.gasH2, gasCO2: st.gasCO2, species, solids, SI: eq?.SI ?? {},
      color: eq ? solutionColor(eq.species) : [235, 240, 245], converged: eq?.converged ?? true,
      thermostat: st.thermostat, diagnostics: eq ? {converged:eq.converged,residual:eq.maxRelativeResidual,charge:eq.chargeResidual,warnings:eq.warnings}:null,
    };
  }
}

// masse atomiche (IUPAC 2021, valori abbreviati) per le masse molari dei solidi
const ATOMIC = { H: 1.008, C: 12.011, N: 14.007, O: 15.999, F: 18.998, Na: 22.990, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974, S: 32.06, Cl: 35.45, K: 39.098, Ca: 40.078, Cr: 51.996, Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38, Br: 79.904, Sr: 87.62, Ag: 107.868, Cd: 112.414, I: 126.904, Ba: 137.327, Pb: 207.2, Li: 6.94 };

/** Massa molare da una formula come "Ca3(PO4)2" o "CaSO4:2H2O". */
export function molarMass(formula) {
  const parts = formula.split(/[:·]/);
  let total = 0;
  for (const part of parts) {
    const m = part.match(/^(\d+)(.*)$/);
    const mult = m ? parseInt(m[1], 10) : 1;
    const f = m ? m[2] : part;
    const stack = [0];
    let i = 0;
    while (i < f.length) {
      const ch = f[i];
      if (ch === '(') { stack.push(0); i++; continue; }
      if (ch === ')') {
        i++;
        let num = '';
        while (i < f.length && /\d/.test(f[i])) num += f[i++];
        const v = stack.pop() * (num ? parseInt(num, 10) : 1);
        stack[stack.length - 1] += v;
        continue;
      }
      const em = f.slice(i).match(/^([A-Z][a-z]?)(\d*\.?\d*)/);
      if (!em) return null;
      const w = ATOMIC[em[1]];
      if (w === undefined) return null;
      stack[stack.length - 1] += w * (em[2] ? parseFloat(em[2]) : 1);
      i += em[0].length;
    }
    total += mult * stack[0];
  }
  return total;
}
function molarMassOfPhase(p) { return molarMass(p.f); }
