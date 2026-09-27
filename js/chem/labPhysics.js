// Leggi della chimica fisica per il laboratorio: gas reali, equilibri acido-base, elettrochimica, cinetica.

import { R_LBAR, R_GAS, FARADAY } from './labData.js';

// ---------------------------------------------------------------------------
// Gas reali: equazione di van der Waals  (P + a/V²)(V − b) = RT     (V molare)
// ---------------------------------------------------------------------------

export function vdwPressure(V, T, a, b) {
  return R_LBAR * T / (V - b) - a / (V * V);
}

/** Punto critico: dP/dV = d²P/dV² = 0  →  Vc = 3b, Tc = 8a/27Rb, Pc = a/27b² */
export function vdwCritical(a, b) {
  return { Vc: 3 * b, Tc: 8 * a / (27 * R_LBAR * b), Pc: a / (27 * b * b), Zc: 3 / 8 };
}

/** Radici reali dell'equazione cubica in V: P V³ − (Pb + RT) V² + a V − a b = 0 */
export function vdwVolumes(P, T, a, b) {
  const A = P, B = -(P * b + R_LBAR * T), C = a, D = -a * b;
  // metodo di Cardano sulla cubica monica
  const p = B / A, q = C / A, r = D / A;
  const Q = (3 * q - p * p) / 9;
  const Rr = (9 * p * q - 27 * r - 2 * p * p * p) / 54;
  const disc = Q * Q * Q + Rr * Rr;
  const roots = [];
  if (disc > 0) {
    const s = Math.cbrt(Rr + Math.sqrt(disc));
    const t = Math.cbrt(Rr - Math.sqrt(disc));
    roots.push(s + t - p / 3);
  } else {
    const th = Math.acos(Rr / Math.sqrt(-Q * Q * Q));
    for (let k = 0; k < 3; k++) roots.push(2 * Math.sqrt(-Q) * Math.cos((th + 2 * Math.PI * k) / 3) - p / 3);
  }
  return roots.filter(v => v > b).sort((x, y) => x - y);
}

/**
 * Costruzione di Maxwell delle aree uguali (T < Tc): pressione di vapore P_sat tale che
 * ∫ P dV tra V_liquido e V_gas = P_sat (V_gas − V_liquido), cioè
 * RT ln[(V_g − b)/(V_l − b)] + a (1/V_g − 1/V_l) − P_sat (V_g − V_l) = 0
 */
export function maxwellConstruction(T, a, b) {
  const { Tc } = vdwCritical(a, b);
  if (T >= Tc) return null;
  // estremi locali (spinodali) di P(V): dP/dV = 0  ⇔  RT V³ = 2a (V − b)²
  const f = (V) => R_LBAR * T * V * V * V - 2 * a * (V - b) * (V - b);
  const Vc = 3 * b;
  const bisect = (lo, hi) => {
    for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); if (Math.sign(f(m)) === Math.sign(f(lo))) lo = m; else hi = m; }
    return 0.5 * (lo + hi);
  };
  const V1 = bisect(b * 1.0001, Vc);        // minimo locale (lato liquido)
  const V2 = bisect(Vc, 1e4 * b);           // massimo locale (lato gas)
  let Plo = Math.max(vdwPressure(V1, T, a, b), 1e-9);
  let Phi = vdwPressure(V2, T, a, b);
  const area = (P) => {
    const V = vdwVolumes(P, T, a, b);
    if (V.length < 3) return null;
    const Vl = V[0], Vg = V[2];
    return { g: R_LBAR * T * Math.log((Vg - b) / (Vl - b)) + a * (1 / Vg - 1 / Vl) - P * (Vg - Vl), Vl, Vg };
  };
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (Plo + Phi);
    const r = area(m);
    if (!r) { Plo = m; continue; }
    if (r.g > 0) Plo = m; else Phi = m;
  }
  const P = 0.5 * (Plo + Phi);
  const r = area(P);
  return r ? { P, Vl: r.Vl, Vg: r.Vg, spinodal: [V1, V2] } : null;
}

// ---------------------------------------------------------------------------
// Equilibri acido-base in acqua (25 °C): pH esatto dal bilancio di carica
// ---------------------------------------------------------------------------

export const KW = 1e-14;

/** Frazioni delle specie H_nA, H_(n−1)A⁻, … , A^(n−) a una data [H⁺]. */
export function fractions(H, pKa) {
  const Ka = pKa.map(p => Math.pow(10, -p));
  const n = Ka.length;
  const terms = [];
  let prod = 1;
  for (let i = 0; i <= n; i++) {
    if (i > 0) prod *= Ka[i - 1];
    terms.push(Math.pow(H, n - i) * prod);
  }
  const D = terms.reduce((s, x) => s + x, 0);
  return terms.map(t => t / D);
}

/** Risolve il bilancio di carica f(pH) = 0 per bisezione su pH ∈ [−1, 15]. */
function solvePH(balance) {
  let lo = -1, hi = 15;
  const g = (pH) => balance(Math.pow(10, -pH));
  for (let i = 0; i < 100; i++) {
    const m = 0.5 * (lo + hi);
    if (g(m) > 0) lo = m; else hi = m;
  }
  return 0.5 * (lo + hi);
}

/**
 * Titolazione di un acido (Ca, Va, pKa) con base forte (Cb), oppure di una base debole
 * (Cb0, Vb0, pKa del coniugato) con acido forte. Volumi in mL, concentrazioni in mol/L.
 * Bilancio di carica:  [H⁺] + [cationi] = [OH⁻] + [anioni]
 */
export function titrationPH({ analyte, titrant, Canalyte, Vanalyte, Ctitrant, Vadded }) {
  const Vt = Vanalyte + Vadded;
  if (analyte.kind === 'acid') {
    const CA = Canalyte * Vanalyte / Vt;
    const Na = Ctitrant * Vadded / Vt;
    return solvePH((H) => {
      const a = fractions(H, analyte.pKa);
      const charge = a.reduce((s, x, i) => s + i * x, 0) * CA;
      return H + Na - KW / H - charge;
    });
  }
  // base debole B (acido coniugato BH⁺ con pKa) titolata con acido forte HCl
  const CB = Canalyte * Vanalyte / Vt;
  const Cl = Ctitrant * Vadded / Vt;
  const Ka = Math.pow(10, -analyte.pKaConj);
  return solvePH((H) => H + CB * H / (H + Ka) - KW / H - Cl);
}

// ---------------------------------------------------------------------------
// Elettrochimica: pila galvanica
// ---------------------------------------------------------------------------

function gcd(a, b) { return b ? gcd(b, a % b) : a; }

/**
 * Catodo: coppia con potenziale di riduzione maggiore. n = m.c.m. degli elettroni scambiati.
 * E = E° − (RT/nF) ln Q  (equazione di Nernst),  ΔG° = −nFE°,  K = exp(nFE°/RT)
 */
export function galvanicCell(h1, h2, { cOxCathode = 1, cOxAnode = 1, T = 298.15 } = {}) {
  const [cat, an] = h1.E >= h2.E ? [h1, h2] : [h2, h1];
  const n = cat.n * an.n / gcd(cat.n, an.n);
  const E0 = cat.E - an.E;
  const dG0 = -n * FARADAY * E0 / 1000;                  // kJ/mol
  const lnK = n * FARADAY * E0 / (R_GAS * T);
  // Q per coppie metallo/ione: l'anodo produce ioni, il catodo li consuma
  const Q = Math.pow(cOxAnode, n / an.n) / Math.pow(cOxCathode, n / cat.n);
  const E = E0 - R_GAS * T / (n * FARADAY) * Math.log(Q);
  return { cathode: cat, anode: an, n, E0, E, dG0, dG: -n * FARADAY * E / 1000, lnK, Q };
}

// ---------------------------------------------------------------------------
// Cinetica chimica: leggi integrate
// ---------------------------------------------------------------------------

export function integratedRate(order, A0, k, t) {
  if (order === 0) return Math.max(0, A0 - k * t);
  if (order === 1) return A0 * Math.exp(-k * t);
  return 1 / (1 / A0 + k * t);
}

export function halfLife(order, A0, k) {
  if (order === 0) return A0 / (2 * k);
  if (order === 1) return Math.LN2 / k;
  return 1 / (k * A0);
}

/** Reazioni consecutive A → B → C (soluzione analitica di Bateman). */
export function consecutive(A0, k1, k2, t) {
  const A = A0 * Math.exp(-k1 * t);
  const B = Math.abs(k2 - k1) < 1e-12
    ? A0 * k1 * t * Math.exp(-k1 * t)
    : A0 * k1 / (k2 - k1) * (Math.exp(-k1 * t) - Math.exp(-k2 * t));
  return { A, B, C: A0 - A - B };
}
