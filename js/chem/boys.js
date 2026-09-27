// Funzione di Boys F_n(T) = ∫₀¹ t²ⁿ e^(−T t²) dt, cuore degli integrali coulombiani tra gaussiane.
// Tabulata su una griglia fine e interpolata con uno sviluppo di Taylor a 6 termini:
//   F_n(T) = Σ_k F_{n+k}(T₀) (−ΔT)^k / k!

const NMAX = 24;          // ordine massimo richiesto (fino a funzioni f e derivate)
const TAYLOR = 7;
const STEP = 0.05;
const TMAX = 40;
const NPTS = Math.round(TMAX / STEP) + 1;
const TABLE = new Float64Array(NPTS * (NMAX + TAYLOR + 1));
const W = NMAX + TAYLOR + 1;

// F_n(T) accurata: serie per la massima n, poi ricorrenza discendente (stabile).
function boysExact(T, nmax, out) {
  const e = Math.exp(-T);
  // F_nmax(T) = e^{−T} Σ_k (2T)^k / [(2n+1)(2n+3)…(2n+2k+1)]
  let term = 1 / (2 * nmax + 1);
  let sum = term;
  for (let k = 1; k < 400; k++) {
    term *= 2 * T / (2 * nmax + 2 * k + 1);
    sum += term;
    if (term < sum * 1e-17) break;
  }
  out[nmax] = e * sum;
  for (let n = nmax - 1; n >= 0; n--) out[n] = (2 * T * out[n + 1] + e) / (2 * n + 1);
}

{
  const tmp = new Float64Array(W);
  for (let i = 0; i < NPTS; i++) {
    boysExact(i * STEP, W - 1, tmp);
    TABLE.set(tmp, i * W);
  }
}

const INV_FACT = [1, 1, 1 / 2, 1 / 6, 1 / 24, 1 / 120, 1 / 720, 1 / 5040];

/** Riempie out[0..nmax] con F_n(T). */
export function boys(T, nmax, out) {
  if (T < TMAX - 1) {
    const i = Math.round(T / STEP);
    const dT = i * STEP - T;
    const base = i * W;
    for (let n = 0; n <= nmax; n++) {
      let s = 0;
      let p = 1;
      for (let k = 0; k < TAYLOR; k++) {
        s += TABLE[base + n + k] * p * INV_FACT[k];
        p *= dT;
      }
      out[n] = s;
    }
    return out;
  }
  // T grande: F_0 = ½√(π/T) (erf ≈ 1), ricorrenza ascendente (stabile per T grande)
  const e = Math.exp(-T);
  out[0] = 0.5 * Math.sqrt(Math.PI / T);
  for (let n = 0; n < nmax; n++) out[n + 1] = ((2 * n + 1) * out[n] - e) / (2 * T);
  return out;
}

export { boysExact };
