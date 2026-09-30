// Strumenti di chimica fisica per la dinamica molecolare della sandbox:
//  • funzione di distribuzione radiale g(r), normalizzata sulla geometria vera della scatola (cubo con pareti,
//    senza condizioni periodiche): g(r) = istogramma / (coppie × densità delle distanze di punti uniformi nel cubo);
//  • spostamento quadratico medio e coefficiente di diffusione (relazione di Einstein, MSD = 6 D t);
//  • capacità termica dalle fluttuazioni dell'energia nell'insieme canonico, C_V = (⟨E²⟩ − ⟨E⟩²)/(k_B T²).

import { KB_EV } from './reactiveData.js';

// direzioni quasi uniformi sulla sfera (spirale di Fibonacci) per la media sulle orientazioni
const DIRS = (() => {
  const n = 1200, out = [];
  const ga = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < n; i++) {
    const z = 1 - (2 * i + 1) / n, r = Math.sqrt(1 - z * z), phi = i * ga;
    out.push([Math.abs(r * Math.cos(phi)), Math.abs(r * Math.sin(phi)), Math.abs(z)]);
  }
  return out;
})();

/**
 * Densità di probabilità della distanza fra due punti uniformi in un cubo di lato L:
 * p(r) = 4π r² γ̄(r) / L⁶, dove γ̄(r) = ⟨(L − r|u_x|)(L − r|u_y|)(L − r|u_z|)⟩ è il covariogramma medio.
 */
export function cubePairDensity(r, L) {
  let g = 0;
  for (const [a, b, c] of DIRS) {
    const x = L - r * a, y = L - r * b, z = L - r * c;
    if (x > 0 && y > 0 && z > 0) g += x * y * z;
  }
  g /= DIRS.length;
  return 4 * Math.PI * r * r * g / Math.pow(L, 6);
}

export class RDF {
  constructor(rmax = 10, nb = 100) { this.rmax = rmax; this.nb = nb; this.reset(); }
  reset(pair = null) { this.hist = new Float64Array(this.nb); this.samples = 0; this.pairs = 0; this.Lsum = 0; this.pair = pair; }
  /** Aggiunge una configurazione. pair = [ZA, ZB]; si escludono le coppie nella stessa molecola. */
  accumulate(sim) {
    if (!this.pair) return;
    const [za, zb] = this.pair;
    const A = [], B = [];
    for (let i = 0; i < sim.N; i++) { if (sim.Z[i] === za) A.push(i); if (sim.Z[i] === zb) B.push(i); }
    const same = za === zb;
    const frag = sim.fragOf;
    const dr = this.rmax / this.nb;
    let np = 0;
    for (let ia = 0; ia < A.length; ia++) {
      const i = A[ia];
      const list = same ? A.slice(ia + 1) : B;
      for (const j of list) {
        if (i === j) continue;
        if (frag && frag[i] && frag[i] === frag[j]) continue;
        np++;
        const dx = sim.pos[3 * i] - sim.pos[3 * j], dy = sim.pos[3 * i + 1] - sim.pos[3 * j + 1], dz = sim.pos[3 * i + 2] - sim.pos[3 * j + 2];
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r < this.rmax) this.hist[Math.floor(r / dr)]++;
      }
    }
    this.pairs += np;
    this.samples++;
    this.Lsum += sim.box;
  }
  result() {
    if (!this.samples || !this.pairs) return null;
    const L = this.Lsum / this.samples;
    const dr = this.rmax / this.nb;
    const r = [], g = [];
    for (let k = 0; k < this.nb; k++) {
      const rc = (k + 0.5) * dr;
      // densità attesa integrata sul bin (regola di Simpson)
      const p = (cubePairDensity(k * dr, L) + 4 * cubePairDensity(rc, L) + cubePairDensity((k + 1) * dr, L)) / 6 * dr;
      r.push(rc);
      g.push(p > 0 ? this.hist[k] / (this.pairs * p) : 0);
    }
    return { r, g, samples: this.samples, pair: this.pair, L };
  }
}

/** Spostamento quadratico medio degli atomi dall'istante di riferimento. */
export class MSD {
  constructor() { this.reset(null); }
  reset(sim) {
    this.ref = sim ? Float64Array.from(sim.pos) : null;
    this.t0 = sim ? sim.time : 0;
    this.N = sim ? sim.N : 0;
    this.points = [];
  }
  sample(sim) {
    if (!this.ref || sim.N !== this.N) { this.reset(sim); return; }
    let s = 0;
    for (let k = 0; k < 3 * sim.N; k++) { const d = sim.pos[k] - this.ref[k]; s += d * d; }
    this.points.push([sim.time - this.t0, sim.N ? s / sim.N : 0]);
    if (this.points.length > 1500) this.points = this.points.filter((_, i) => i % 2 === 0);
  }
  /**
   * D dalla pendenza della parte lineare (Einstein, MSD = 6 D t). Si usa l'intervallo fra il 30 % e il 70 % del
   * tempo finché MSD resta sotto il 30 % del valore di saturazione L²/2 imposto dalle pareti. L'esponente
   * α = d ln MSD / d ln t distingue il moto balistico (α = 2, fra un urto e l'altro) da quello diffusivo (α = 1):
   * D si riporta solo se 0,75 < α < 1,25.
   */
  result(box) {
    const pts = this.points;
    if (pts.length < 8) return { points: pts, D: null, alpha: null };
    const tmax = pts[pts.length - 1][0];
    const sat = box * box / 2;
    const sel = pts.filter(([t, m]) => t >= 0.3 * tmax && t <= 0.7 * tmax && m < 0.3 * sat && t > 0 && m > 0);
    let D = null, fit = null, alpha = null;
    if (sel.length >= 5) {
      const lin = (xs, ys) => {
        let sx = 0, sy = 0, sxx = 0, sxy = 0;
        const n = xs.length;
        for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
        const den = n * sxx - sx * sx;
        if (!(den > 0)) return null;
        const b = (n * sxy - sx * sy) / den;
        return { a: (sy - b * sx) / n, b };
      };
      const lg = lin(sel.map(p => Math.log(p[0])), sel.map(p => Math.log(p[1])));
      alpha = lg?.b ?? null;
      const l = lin(sel.map(p => p[0]), sel.map(p => p[1]));
      if (l) {
        fit = { a: l.a, b: l.b, t0: sel[0][0], t1: sel[sel.length - 1][0] };
        if (alpha > 0.75 && alpha < 1.25) D = l.b / 6 * 0.1; // Å²/fs → cm²/s
      }
    }
    return { points: pts, D, fit, alpha, saturation: sat };
  }
}

/** Capacità termica a volume costante dalle fluttuazioni dell'energia totale (insieme canonico). */
export class HeatCapacity {
  constructor() { this.reset(); }
  reset(key = null) { this.E = []; this.key = key; }
  sample(sim, key) {
    if (key !== this.key) this.reset(key);
    this.E.push(sim.totalEnergy());
    if (this.E.length > 20000) this.E.splice(0, this.E.length - 20000);
  }
  /** Restituisce C_V in eV/K con l'errore dalla media a blocchi (10 blocchi). */
  result(T) {
    const E = this.E, n = E.length;
    if (n < 40 || !(T > 0)) return null;
    const cv = (arr) => {
      let m = 0; for (const e of arr) m += e; m /= arr.length;
      let v = 0; for (const e of arr) v += (e - m) ** 2; v /= arr.length - 1;
      return v / (KB_EV * T * T);
    };
    const all = cv(E);
    const nb = 10, bl = Math.floor(n / nb);
    const blocks = [];
    for (let b = 0; b < nb; b++) blocks.push(cv(E.slice(b * bl, (b + 1) * bl)));
    const mb = blocks.reduce((a, x) => a + x, 0) / nb;
    const err = Math.sqrt(blocks.reduce((a, x) => a + (x - mb) ** 2, 0) / (nb * (nb - 1)));
    // stazionarietà: l'energia media della prima e della seconda metà non deve differire più di quanto
    // consentono le fluttuazioni (reazioni in corso o riscaldamento rendono C_V privo di significato)
    const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;
    const h = Math.floor(n / 2);
    const sd = Math.sqrt(all * KB_EV * T * T);
    const drift = Math.abs(mean(E.slice(0, h)) - mean(E.slice(h)));
    const stationary = drift < 1.5 * sd;
    return { Cv: all, err, samples: n, stationary };
  }
}
