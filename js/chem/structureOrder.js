// Ordine strutturale locale: parametri di Steinhardt q₆ (Steinhardt, Nelson, Ronchetti, Phys. Rev. B 28, 784, 1983)
// e criterio dei "legami solidi" di ten Wolde, Ruiz-Montero e Frenkel (J. Chem. Phys. 104, 9932, 1996).
// Ogni atomo riceve il vettore q₆ₘ medio sui vicini; due vicini sono "connessi da un legame solido" se il prodotto
// scalare normalizzato dei loro q₆ supera 0,7. Un atomo con almeno 7 legami solidi è di tipo cristallino.
// Il criterio non dipende dall'orientazione del reticolo né dalla posizione iniziale degli atomi, quindi riconosce
// anche un cristallo che si riforma dopo la fusione. Non distingue FCC, HCP e BCC fra loro.

const L = 6;
// fattori di normalizzazione √((l−m)!/(l+m)!) per m = 0…6 (il fattore comune √((2l+1)/4π) si semplifica)
const NORM = Array.from({ length: L + 1 }, (_, m) => {
  let r = 1;
  for (let k = L - m + 1; k <= L + m; k++) r /= k;
  return Math.sqrt(r);
});

// Polinomi associati di Legendre P₆ᵐ(x), m = 0…6, per ricorrenza.
function legendre6(x, out) {
  const s = Math.sqrt(Math.max(0, 1 - x * x));
  for (let m = 0; m <= L; m++) {
    let pmm = 1;
    for (let k = 1; k <= m; k++) pmm *= -(2 * k - 1) * s;
    if (m === L) { out[m] = pmm; continue; }
    let pm1 = x * (2 * m + 1) * pmm;
    if (m + 1 === L) { out[m] = pm1; continue; }
    let pl = 0;
    for (let l = m + 2; l <= L; l++) {
      pl = ((2 * l - 1) * x * pm1 - (l + m - 1) * pmm) / (l - m);
      pmm = pm1; pm1 = pl;
    }
    out[m] = pl;
  }
}

/**
 * pos: coordinate (Float64Array 3N); box: lato per condizioni periodiche (minima immagine) oppure null;
 * cutoff: raggio della prima sfera di coordinazione (per il Lennard–Jones ≈ 1,5σ, primo minimo di g(r)).
 * include: indici degli atomi da analizzare (predefinito: tutti).
 */
export function localOrder(pos, { box = null, cutoff, include = null, dotThreshold = 0.7, minSolidBonds = 7 } = {}) {
  const idx = include ?? Array.from({ length: pos.length / 3 }, (_, i) => i);
  const n = idx.length;
  const re = new Float64Array(n * (L + 1)), im = new Float64Array(n * (L + 1));
  const nb = Array.from({ length: n }, () => []);
  const c2 = cutoff * cutoff, P = new Float64Array(L + 1);
  for (let a = 0; a < n; a++) {
    const i = idx[a];
    for (let b = a + 1; b < n; b++) {
      const j = idx[b];
      let dx = pos[3 * j] - pos[3 * i], dy = pos[3 * j + 1] - pos[3 * i + 1], dz = pos[3 * j + 2] - pos[3 * i + 2];
      if (box) { dx -= box * Math.round(dx / box); dy -= box * Math.round(dy / box); dz -= box * Math.round(dz / box); }
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 >= c2 || r2 === 0) continue;
      nb[a].push(b); nb[b].push(a);
      const r = Math.sqrt(r2);
      legendre6(dz / r, P);
      const phi = Math.atan2(dy, dx);
      for (let m = 0; m <= L; m++) {
        const y = NORM[m] * P[m], cr = y * Math.cos(m * phi), ci = y * Math.sin(m * phi);
        // Y₆ₘ(−r̂) = (−1)^l Y₆ₘ(r̂) = Y₆ₘ(r̂) per l pari: stesso contributo ai due atomi
        re[a * (L + 1) + m] += cr; im[a * (L + 1) + m] += ci;
        re[b * (L + 1) + m] += cr; im[b * (L + 1) + m] += ci;
      }
    }
  }
  // norme: Σ_m |q₆ₘ|² = |q₆₀|² + 2 Σ_{m>0} |q₆ₘ|²
  const norm = new Float64Array(n);
  for (let a = 0; a < n; a++) {
    let s = 0;
    for (let m = 0; m <= L; m++) { const k = a * (L + 1) + m; s += (m ? 2 : 1) * (re[k] ** 2 + im[k] ** 2); }
    norm[a] = Math.sqrt(s);
  }
  const solidBonds = new Int16Array(n);
  for (let a = 0; a < n; a++) {
    if (!norm[a]) continue;
    for (const b of nb[a]) {
      if (b < a || !norm[b]) continue;
      let d = 0;
      for (let m = 0; m <= L; m++) {
        const ka = a * (L + 1) + m, kb = b * (L + 1) + m;
        d += (m ? 2 : 1) * (re[ka] * re[kb] + im[ka] * im[kb]);
      }
      if (d / (norm[a] * norm[b]) > dotThreshold) { solidBonds[a]++; solidBonds[b]++; }
    }
  }
  // Q₆ globale (media di tutti i q₆ₘ di legame): 0,575 per FCC perfetto, ≈ 0 per un liquido grande
  let bonds = 0;
  const sumRe = new Float64Array(L + 1), sumIm = new Float64Array(L + 1);
  for (let a = 0; a < n; a++) {
    bonds += nb[a].length;
    for (let m = 0; m <= L; m++) { sumRe[m] += re[a * (L + 1) + m]; sumIm[m] += im[a * (L + 1) + m]; }
  }
  let Q = 0;
  if (bonds) for (let m = 0; m <= L; m++) Q += (m ? 2 : 1) * (sumRe[m] ** 2 + sumIm[m] ** 2);
  // con la normalizzazione completa Y = √(13/4π)·NORM·P·e^{imφ}: Q₆ = √(4π/13 Σ|Q₆ₘ|²) = √(Σ|NORM·…|²)/legami
  const Q6 = bonds ? Math.sqrt(Q) / bonds : 0;
  const kind = new Array(n);
  const counts = { crystal: 0, liquid: 0, vapor: 0 };
  for (let a = 0; a < n; a++) {
    kind[a] = solidBonds[a] >= minSolidBonds ? 'crystal' : nb[a].length >= 3 ? 'liquid' : 'vapor';
    counts[kind[a]]++;
  }
  return {
    kind, neighbors: nb.map(x => x.length), solidBonds, Q6,
    fraction: { crystal: counts.crystal / (n || 1), liquid: counts.liquid / (n || 1), vapor: counts.vapor / (n || 1) },
  };
}

/**
 * Lettura qualitativa delle frazioni di atomi cristallini, condensati e isolati.
 * mobile: true se le gabbie di primi vicini si rinnovano (vedi Mobility), false se gli atomi vibrano soltanto,
 * null se non è ancora misurato. Distingue il liquido da un solido disordinato (vetro, aggregato non cristallino).
 */
export function phaseVerdict(f, mobile = null) {
  const cond = f.crystal + f.liquid;
  if (f.crystal >= 0.6) return { id: 'solido', name: 'Cristallo', why: `${pct(f.crystal)} degli atomi ha un intorno ordinato (q₆)` };
  if (f.crystal >= 0.15 && cond >= 0.5) return { id: 'misto', name: 'Solido + liquido', why: `${pct(f.crystal)} cristallino, ${pct(f.liquid)} disordinato ma condensato` };
  if (f.vapor >= 0.85) return { id: 'gas', name: 'Gas', why: `${pct(f.vapor)} degli atomi ha meno di 3 vicini` };
  if (cond >= 0.85 && mobile === false) return { id: 'solido', name: 'Solido disordinato', why: `${pct(cond)} condensato, senza ordine cristallino e senza diffusione: vetro o aggregato non cristallino` };
  if (f.vapor < 0.03) return { id: 'liquido', name: 'Liquido', why: `${pct(cond)} condensato senza ordine cristallino${mobile ? ', atomi che diffondono' : ''}` };
  return { id: 'misto', name: 'Liquido + vapore', why: `${pct(cond)} condensato, ${pct(f.vapor)} isolato` };
}

/**
 * Mobilità indipendente da traslazioni e rotazioni d'insieme (una goccia che deriva nella scatola non è un liquido):
 * fra due istantanee distanti almeno 3τ si contano le coppie di primi vicini (r < 1,35σ) che si sono separate
 * oltre 1,65σ. In un liquido le gabbie di vicini si rinnovano in circa τ; in un solido, anche disordinato, no.
 * Misurato sull'argon LJ in sandbox: 0,5–0,6 nel liquido a 90 K, 0,05–0,12 nell'aggregato a 40 K (solo la superficie
 * si muove), 0,29 nel liquido LJ di bulk vicino al punto triplo (T* = 0,75); soglia 0,2.
 */
export class Mobility {
  constructor() { this.snaps = []; }
  /** t e tau nella stessa unità; sigma nell'unità delle posizioni; box per la minima immagine (o null). */
  sample(pos, t, tau, sigma, box = null) {
    if (this.snaps.length && (this.snaps.at(-1).pos.length !== pos.length || t < this.snaps.at(-1).t)) this.snaps = [];
    this.snaps.push({ t, pos: Float64Array.from(pos) });
    while (this.snaps.length > 2 && t - this.snaps[1].t >= 3 * tau) this.snaps.shift();
    const old = this.snaps[0];
    if (t - old.t < 3 * tau) return null;
    const n = pos.length / 3, near = (1.35 * sigma) ** 2, far = (1.65 * sigma) ** 2;
    const d2 = (p, i, j) => {
      let s = 0;
      for (let a = 0; a < 3; a++) { let v = p[3 * i + a] - p[3 * j + a]; if (box) v -= box * Math.round(v / box); s += v * v; }
      return s;
    };
    let pairs = 0, lost = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      if (d2(old.pos, i, j) >= near) continue;
      pairs++;
      if (d2(pos, i, j) > far) lost++;
    }
    this.lostFraction = pairs ? lost / pairs : null;
    return pairs ? lost / pairs > 0.2 : null;
  }
}

const pct = (x) => `${Math.round(100 * x)}%`;
