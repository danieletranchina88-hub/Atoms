// Funzioni d'onda tridimensionali: ψ(x,y,z) = R_nl(r) · Y_lm(θ,φ), combinazioni lineari (LCAO, ibridi)
// e campionamento Monte Carlo esatto della densità di probabilità |ψ|².

import { realYlm } from './harmonics.js';

/** Funzione radiale R(r) = u(r)/r interpolata sulla griglia logaritmica del calcolo. */
export class RadialFunction {
  constructor(grid, u, l) {
    this.l = l;
    this.xmin = grid.xmin;
    this.h = grid.h;
    this.Z = grid.Z;
    this.N = u.length;
    const r0 = Math.exp(grid.xmin) / grid.Z;
    this.r0 = r0;
    this.R = new Float64Array(this.N);
    this.r = new Float64Array(this.N);
    for (let i = 0; i < this.N; i++) {
      const r = Math.exp(grid.xmin + i * grid.h) / grid.Z;
      this.r[i] = r;
      this.R[i] = u[i] / r;
    }
    this.u = u;
    // Funzione di ripartizione della probabilità radiale ∫ u² dr, per il campionamento.
    this.cdf = new Float64Array(this.N);
    let acc = 0;
    for (let i = 1; i < this.N; i++) {
      const dr = this.r[i] - this.r[i - 1];
      acc += 0.5 * (u[i] * u[i] + u[i - 1] * u[i - 1]) * dr;
      this.cdf[i] = acc;
    }
    for (let i = 0; i < this.N; i++) this.cdf[i] /= acc;
  }

  at(r) {
    if (r <= this.r0) return this.R[0] * Math.pow(r / this.r0, this.l);
    const t = (Math.log(r * this.Z) - this.xmin) / this.h;
    const i = Math.floor(t);
    if (i >= this.N - 1) return 0;
    const f = t - i;
    return this.R[i] * (1 - f) + this.R[i + 1] * f;
  }

  /** Raggio campionato dalla distribuzione P(r) = u(r)². */
  sampleR(rand) {
    const q = rand();
    let lo = 0;
    let hi = this.N - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (this.cdf[mid] < q) lo = mid; else hi = mid;
    }
    const span = this.cdf[hi] - this.cdf[lo];
    const f = span > 0 ? (q - this.cdf[lo]) / span : 0;
    return this.r[lo] + f * (this.r[hi] - this.r[lo]);
  }

  /** Raggio che racchiude la frazione indicata della probabilità. */
  radiusEnclosing(fraction) {
    for (let i = 0; i < this.N; i++) if (this.cdf[i] >= fraction) return this.r[i];
    return this.r[this.N - 1];
  }
}

/**
 * Un termine di una combinazione lineare: coeff · R_nl(|r − c|) · Y_lm(r − c).
 * `radial` è una RadialFunction; `center` in raggi di Bohr.
 */
export function term(radial, l, m, coeff = 1, center = [0, 0, 0]) {
  return { radial, l, m, coeff, center };
}

/** Valutazione di ψ = Σ c_k φ_k nel punto (x, y, z). */
export function evaluate(terms, x, y, z) {
  let psi = 0;
  for (let k = 0; k < terms.length; k++) {
    const t = terms[k];
    const dx = x - t.center[0];
    const dy = y - t.center[1];
    const dz = z - t.center[2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
    psi += t.coeff * t.radial.at(r) * realYlm(t.l, t.m, dx, dy, dz);
  }
  return psi;
}

// Generatore pseudo-casuale riproducibile (mulberry32).
export function makeRandom(seed = 12345) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Massimo di |Y_lm|² sulla sfera (stimato su una griglia fine, con margine).
const ymaxCache = new Map();
function maxY2(l, m) {
  const key = `${l},${m}`;
  if (ymaxCache.has(key)) return ymaxCache.get(key);
  let max = 0;
  const nt = 90;
  const np = 180;
  for (let i = 0; i <= nt; i++) {
    const th = Math.PI * i / nt;
    for (let j = 0; j < np; j++) {
      const ph = 2 * Math.PI * j / np;
      const y = realYlm(l, m, Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th));
      if (y * y > max) max = y * y;
    }
  }
  max *= 1.05;
  ymaxCache.set(key, max);
  return max;
}

function randomDirection(rand) {
  const z = 2 * rand() - 1;
  const phi = 2 * Math.PI * rand();
  const s = Math.sqrt(1 - z * z);
  return [s * Math.cos(phi), s * Math.sin(phi), z];
}

/** Punto campionato esattamente da |R_nl Y_lm|² (radiale per inversione, angolare per rigetto). */
function sampleSingle(t, rand) {
  const ymax = maxY2(t.l, t.m);
  let d;
  for (;;) {
    d = randomDirection(rand);
    const y = realYlm(t.l, t.m, d[0], d[1], d[2]);
    if (rand() * ymax <= y * y) break;
  }
  const r = t.radial.sampleR(rand);
  return [t.center[0] + r * d[0], t.center[1] + r * d[1], t.center[2] + r * d[2]];
}

/**
 * Campiona `count` punti dalla densità |Σ c_k φ_k|².
 * Metodo del rigetto con proposta q = (1/K) Σ |φ_k|²: per la disuguaglianza di Cauchy–Schwarz
 * |Σ c_k φ_k|² ≤ (Σ c_k²)(Σ |φ_k|²) = (Σ c_k²)·K·q, quindi il campionamento è esatto.
 * Restituisce posizioni e segno di ψ in ciascun punto.
 */
export function samplePoints(terms, count, seed = 1) {
  const rand = makeRandom(seed);
  const positions = new Float32Array(count * 3);
  const signs = new Int8Array(count);
  const K = terms.length;
  const c2 = terms.reduce((s, t) => s + t.coeff * t.coeff, 0);
  let n = 0;
  let guard = 0;
  while (n < count && guard < count * 400) {
    guard++;
    const t = terms[Math.floor(rand() * K)];
    const p = sampleSingle(t, rand);
    let psi = 0;
    if (K > 1) {
      let q = 0;
      for (const s of terms) {
        const dx = p[0] - s.center[0];
        const dy = p[1] - s.center[1];
        const dz = p[2] - s.center[2];
        const phi = s.radial.at(Math.sqrt(dx * dx + dy * dy + dz * dz)) * realYlm(s.l, s.m, dx, dy, dz);
        psi += s.coeff * phi;
        q += phi * phi;
      }
      if (rand() * c2 * q > psi * psi) continue;
    } else {
      psi = evaluate(terms, p[0], p[1], p[2]);
    }
    positions[3 * n] = p[0];
    positions[3 * n + 1] = p[1];
    positions[3 * n + 2] = p[2];
    signs[n] = psi >= 0 ? 1 : -1;
    n++;
  }
  return { positions: positions.subarray(0, 3 * n), signs: signs.subarray(0, n), count: n };
}

/**
 * Valuta ψ su una griglia cubica di lato 2·half centrata nell'origine.
 * Restituisce i valori e la soglia |ψ| la cui isosuperficie racchiude la frazione `enclosed` della probabilità.
 */
export function sampleGrid(terms, res, half, enclosed = 0.9) {
  const values = new Float32Array(res * res * res);
  const step = (2 * half) / (res - 1);
  let idx = 0;
  for (let k = 0; k < res; k++) {
    const z = -half + k * step;
    for (let j = 0; j < res; j++) {
      const y = -half + j * step;
      for (let i = 0; i < res; i++) {
        values[idx++] = evaluate(terms, -half + i * step, y, z);
      }
    }
  }
  return { values, res, half, step, iso: isoForFraction(values, enclosed) };
}

/** Soglia su |ψ| tale che la regione |ψ| ≥ soglia contenga la frazione richiesta di ∫|ψ|². */
export function isoForFraction(values, fraction) {
  const n = values.length;
  const sq = new Float32Array(n);
  let total = 0;
  for (let i = 0; i < n; i++) {
    sq[i] = values[i] * values[i];
    total += sq[i];
  }
  sq.sort();
  let acc = 0;
  for (let i = n - 1; i >= 0; i--) {
    acc += sq[i];
    if (acc >= fraction * total) return Math.sqrt(sq[i]);
  }
  return 0;
}

/**
 * Nuvola elettronica dell'intero atomo. `parts` = [{ radial, l, m | null, weight, group }]:
 * con m = null la direzione è uniforme (media sferica del sottolivello), altrimenti si segue |Y_lm|².
 * I punti sono ripartiti in proporzione al numero di elettroni (weight).
 */
export function sampleAtom(parts, count, seed = 7) {
  const rand = makeRandom(seed);
  const total = parts.reduce((s, p) => s + p.weight, 0);
  const positions = new Float32Array(count * 3);
  const groups = new Uint8Array(count);
  let n = 0;
  parts.forEach((p, idx) => {
    const k = idx === parts.length - 1 ? count - n : Math.round(count * p.weight / total);
    const t = { radial: p.radial, l: p.l, m: p.m ?? 0, center: [0, 0, 0] };
    for (let i = 0; i < k && n < count; i++) {
      let pos;
      if (p.m === null || p.m === undefined) {
        const d = randomDirection(rand);
        const r = p.radial.sampleR(rand);
        pos = [r * d[0], r * d[1], r * d[2]];
      } else {
        pos = sampleSingle(t, rand);
      }
      positions[3 * n] = pos[0];
      positions[3 * n + 1] = pos[1];
      positions[3 * n + 2] = pos[2];
      groups[n] = p.group;
      n++;
    }
  });
  return { positions: positions.subarray(0, 3 * n), groups: groups.subarray(0, n), count: n };
}
