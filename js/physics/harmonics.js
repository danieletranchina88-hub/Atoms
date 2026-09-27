// Armoniche sferiche reali Y_lm(θ, φ), la parte angolare degli orbitali atomici.
//   m > 0:  √2 · K_lm · P_l^m(cos θ) · cos(mφ)
//   m = 0:        K_l0 · P_l^0(cos θ)
//   m < 0:  √2 · K_l|m| · P_l^|m|(cos θ) · sin(|m|φ)
// con K_lm = √[(2l+1)/(4π) · (l−m)!/(l+m)!]. Senza la fase di Condon–Shortley, così che
// il lobo positivo di p_x punti verso +x (convenzione dei testi di chimica).

import { L_LETTERS } from './configuration.js';

function factorial(n) {
  let f = 1;
  for (let i = 2; i <= n; i++) f *= i;
  return f;
}

const normCache = new Map();
function norm(l, m) {
  const key = l * 100 + m;
  let k = normCache.get(key);
  if (k === undefined) {
    k = Math.sqrt((2 * l + 1) / (4 * Math.PI) * factorial(l - m) / factorial(l + m));
    if (m > 0) k *= Math.SQRT2;
    normCache.set(key, k);
  }
  return k;
}

/**
 * Polinomio associato di Legendre P_l^m(x)/(1−x²)^{m/2}, cioè senza il fattore sin^m θ
 * (che viene reintrodotto in coordinate cartesiane tramite (x ± iy)^m).
 */
function legendreReduced(l, m, x) {
  // P_m^m / sin^m = (2m−1)!!
  let pmm = 1;
  for (let i = 1; i <= m; i++) pmm *= 2 * i - 1;
  if (l === m) return pmm;
  let pm1 = x * (2 * m + 1) * pmm;
  if (l === m + 1) return pm1;
  let pll = 0;
  for (let ll = m + 2; ll <= l; ll++) {
    pll = ((2 * ll - 1) * x * pm1 - (ll + m - 1) * pmm) / (ll - m);
    pmm = pm1;
    pm1 = pll;
  }
  return pll;
}

/**
 * Valuta Y_lm nel punto (x, y, z) (non serve che sia normalizzato: si usa la direzione).
 * Restituisce 0 all'origine.
 */
export function realYlm(l, m, x, y, z) {
  const r = Math.sqrt(x * x + y * y + z * z);
  if (r === 0) return l === 0 ? norm(0, 0) : 0;
  const ct = z / r;
  const am = Math.abs(m);
  const P = legendreReduced(l, am, ct);
  if (am === 0) return norm(l, 0) * P;
  // sin^m θ · cos(mφ) = Re[(x + iy)^m] / r^m ; sin^m θ · sin(mφ) = Im[(x + iy)^m] / r^m
  let re = 1;
  let im = 0;
  const xr = x / r;
  const yr = y / r;
  for (let i = 0; i < am; i++) {
    const t = re * xr - im * yr;
    im = re * yr + im * xr;
    re = t;
  }
  return norm(l, am) * P * (m > 0 ? re : im);
}

/** Crea una funzione veloce (x, y, z) → Y_lm per l ed m fissati. */
export function makeYlm(l, m) {
  return (x, y, z) => realYlm(l, m, x, y, z);
}

// Nomi cartesiani degli orbitali reali.
const NAMES = {
  0: { 0: '' },
  1: { '-1': 'y', 0: 'z', 1: 'x' },
  2: { '-2': 'xy', '-1': 'yz', 0: 'z²', 1: 'xz', 2: 'x²−y²' },
  3: { '-3': 'y(3x²−y²)', '-2': 'xyz', '-1': 'yz²', 0: 'z³', 1: 'xz²', 2: 'z(x²−y²)', 3: 'x(x²−3y²)' },
};

/** Etichetta dell'orbitale reale, es. "3d" + pedice "xy". */
export function orbitalName(n, l, m) {
  const sub = NAMES[l]?.[m];
  return { base: `${n}${L_LETTERS[l]}`, sub: sub ?? `m=${m}` };
}

export function orbitalNameText(n, l, m) {
  const { base, sub } = orbitalName(n, l, m);
  return sub ? `${base}${sub.startsWith('m=') ? ` (${sub})` : sub}` : base;
}

/** Ordine degli m usato per mostrare gli orbitali (convenzione dei libri: p_x, p_y, p_z; d_xy … ). */
export function mOrder(l) {
  if (l === 1) return [1, -1, 0];
  if (l === 2) return [-2, -1, 1, 2, 0];
  const list = [];
  for (let m = -l; m <= l; m++) list.push(m);
  return list;
}

/** Superfici nodali angolari, descritte a parole. */
export function angularNodesDescription(l, m) {
  if (l === 0) return 'nessuna';
  const key = `${l},${m}`;
  const map = {
    '1,1': 'piano yz (x = 0)',
    '1,-1': 'piano xz (y = 0)',
    '1,0': 'piano xy (z = 0)',
    '2,-2': 'piani xz e yz',
    '2,-1': 'piani xy e xz',
    '2,1': 'piani xy e yz',
    '2,2': 'piani y = x e y = −x',
    '2,0': 'due coni a θ = 54,7° e 125,3°',
    '3,0': 'piano xy e due coni a θ = 39,2° e 140,8°',
    '3,1': 'piano yz e due coni a θ = 63,4° e 116,6°',
    '3,-1': 'piano xz e due coni a θ = 63,4° e 116,6°',
    '3,2': 'piano xy e piani y = ±x',
    '3,-2': 'piani xy, xz e yz',
    '3,3': 'tre piani verticali a φ = 30°, 90°, 150°',
    '3,-3': 'tre piani verticali a φ = 0°, 60°, 120°',
  };
  return map[key] ?? `${l} superfici nodali`;
}
