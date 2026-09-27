// Calcolo autoconsistente di campo centrale per atomi (Kohn–Sham, approssimazione LDA).
//
// Per ogni sottolivello (n, l) si risolve l'equazione di Schrödinger radiale
//     −½ u''(r) + [ l(l+1)/(2r²) + V(r) ] u(r) = ε u(r),   ψ = u(r)/r · Y_lm(θ, φ)
// con V(r) = −Z/r + V_Hartree(r) + v_xc(r). La densità elettronica ottenuta genera un nuovo
// potenziale e si itera fino all'autoconsistenza (metodo di Hartree–Fock–Slater / Herman–Skillman,
// qui con il funzionale di scambio-correlazione LDA di Dirac + Vosko–Wilk–Nusair).
// Unità atomiche di Hartree: lunghezze in a₀ (raggi di Bohr), energie in Hartree.

import { groundStateConfiguration } from './configuration.js';

// ---------------------------------------------------------------------------
// Griglia radiale logaritmica: r_i = exp(x_min + i·h) / Z
// ---------------------------------------------------------------------------

export function makeGrid(Z, { xmin = -8, h = 0.008, rmax = 250 } = {}) {
  const N = Math.ceil((Math.log(rmax * Z) - xmin) / h) + 1;
  const r = new Float64Array(N);
  const r2 = new Float64Array(N);
  const sqr = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    r[i] = Math.exp(xmin + i * h) / Z;
    r2[i] = r[i] * r[i];
    sqr[i] = Math.sqrt(r[i]);
  }
  return { Z, N, h, xmin, r, r2, sqr };
}

/** ∫ f(r) dr su griglia logaritmica (dr = r·h), regola di Simpson. */
export function integrate(grid, f) {
  const { N, h, r } = grid;
  let s = 0;
  for (let i = 0; i < N; i++) {
    const w = (i === 0 || i === N - 1) ? 1 : (i % 2 === 1 ? 4 : 2);
    s += w * f[i] * r[i];
  }
  return s * h / 3;
}

// ---------------------------------------------------------------------------
// Risolutore radiale di Numerov con ricerca dell'autovalore per shooting.
// Con u(r) = √r · y(x), x = ln r:   y'' = g(x) y,   g = (l+½)² + 2r²(V − ε)
// ---------------------------------------------------------------------------

export function solveRadial(grid, V, n, l, eGuess, Zn) {
  const { N, h, r, r2, sqr } = grid;
  const h2 = h * h / 12;
  const lhalf2 = (l + 0.5) * (l + 0.5);
  const nodesWanted = n - l - 1;

  // Limiti iniziali: minimo del potenziale efficace e valore al bordo.
  let elw = Infinity;
  for (let i = 0; i < N; i++) {
    const veff = V[i] + l * (l + 1) / (2 * r2[i]);
    if (veff < elw) elw = veff;
  }
  let eup = V[N - 1] + l * (l + 1) / (2 * r2[N - 1]);
  if (eup > 0) eup = 0;
  let e = (eGuess !== undefined && eGuess > elw && eGuess < eup) ? eGuess : 0.5 * (elw + eup);

  const f = new Float64Array(N);
  const y = new Float64Array(N);
  let converged = false;
  let icl = -1;
  let imax = N - 1;

  for (let iter = 0; iter < 200; iter++) {
    // Funzione di Numerov e punto di inversione classico più esterno.
    icl = -1;
    for (let i = 0; i < N; i++) {
      f[i] = h2 * (lhalf2 + 2 * r2[i] * (V[i] - e));
    }
    for (let i = N - 2; i >= 1; i--) {
      if (f[i] < 0 && f[i + 1] >= 0) { icl = i + 1; break; }
    }
    if (icl < 0 || icl > N - 10) {
      // nessuna regione classicamente permessa (e troppo bassa) o stato non legato (e troppo alta)
      if (icl < 0) elw = e; else eup = e;
      e = 0.5 * (elw + eup);
      if (eup - elw < 1e-12) break;
      continue;
    }
    // Punto di partenza dell'integrazione verso l'interno: dove la funzione d'onda è trascurabile.
    imax = N - 1;
    let decay = 0;
    for (let i = icl; i < N - 1; i++) {
      const g = f[i] / h2;
      decay += Math.sqrt(Math.max(g, 0)) * h;
      if (decay > 45) { imax = i; break; }
    }
    for (let i = 0; i < N; i++) f[i] = 1 - f[i];

    // Integrazione verso l'esterno (0 → icl), comportamento r^(l+1) vicino al nucleo.
    y[0] = Math.pow(r[0], l + 1) * (1 - Zn * r[0] / (l + 1)) / sqr[0];
    y[1] = Math.pow(r[1], l + 1) * (1 - Zn * r[1] / (l + 1)) / sqr[1];
    let nodes = 0;
    for (let i = 1; i < icl; i++) {
      y[i + 1] = ((12 - 10 * f[i]) * y[i] - f[i - 1] * y[i - 1]) / f[i + 1];
      if (y[i] * y[i + 1] < 0) nodes++;
      if (Math.abs(y[i + 1]) > 1e150) {
        for (let k = 0; k <= i + 1; k++) y[k] *= 1e-150;
      }
    }
    const ycl = y[icl];

    if (nodes !== nodesWanted) {
      if (nodes > nodesWanted) eup = e; else elw = e;
      e = 0.5 * (eup + elw);
      if (eup - elw < 1e-12) break;
      continue;
    }

    // Integrazione verso l'interno (imax → icl).
    for (let i = imax + 1; i < N; i++) y[i] = 0;
    y[imax] = h;
    y[imax - 1] = (12 - 10 * f[imax]) * y[imax] / f[imax - 1];
    for (let i = imax - 1; i > icl; i--) {
      y[i - 1] = ((12 - 10 * f[i]) * y[i] - f[i + 1] * y[i + 1]) / f[i - 1];
      if (Math.abs(y[i - 1]) > 1e150) {
        for (let k = i - 1; k <= imax; k++) y[k] *= 1e-150;
      }
    }
    // Raccordo continuo in icl.
    const scale = ycl / y[icl];
    for (let i = icl; i <= imax; i++) y[i] *= scale;

    // Normalizzazione: ∫u² dr = ∫ y² r² dx
    let norm = 0;
    for (let i = 0; i <= imax; i++) {
      const w = (i === 0 || i === imax) ? 1 : (i % 2 === 1 ? 4 : 2);
      norm += w * y[i] * y[i] * r2[i];
    }
    norm = Math.sqrt(norm * h / 3);
    for (let i = 0; i <= imax; i++) y[i] /= norm;

    // Correzione perturbativa dell'energia dalla discontinuità della derivata in icl.
    const ycusp = (y[icl - 1] * f[icl - 1] + y[icl + 1] * f[icl + 1] + 10 * f[icl] * y[icl]) / 12;
    const dfcusp = f[icl] * (y[icl] / ycusp - 1);
    const de = 0.5 * dfcusp / h2 * ycusp * ycusp * h;
    if (de > 0) elw = e;
    if (de < 0) eup = e;
    e += de;
    if (e > eup || e < elw) e = 0.5 * (eup + elw);
    if (Math.abs(de) < 1e-11) { converged = true; break; }
  }

  // u(r) = √r · y, con segno positivo vicino al nucleo (convenzione dei polinomi di Laguerre).
  const u = new Float64Array(N);
  let sign = 1;
  for (let i = 0; i < N; i++) {
    if (Math.abs(y[i]) > 1e-12) { sign = Math.sign(y[i]); break; }
  }
  for (let i = 0; i < N; i++) u[i] = sign * sqr[i] * y[i];
  return { e, u, converged, icl };
}

// ---------------------------------------------------------------------------
// Scambio e correlazione LDA: scambio di Dirac-Slater, correlazione di Vosko–Wilk–Nusair (VWN5).
// ---------------------------------------------------------------------------

const VWN = { A: 0.0310907, x0: -0.10498, b: 3.72744, c: 12.9352 };

export function lda(rho) {
  if (rho < 1e-30) return { exc: 0, vxc: 0 };
  const kf = Math.cbrt(3 * rho / Math.PI);
  const ex = -0.75 * kf;
  const vx = -kf;
  const rs = Math.cbrt(3 / (4 * Math.PI * rho));
  const { A, x0, b, c } = VWN;
  const x = Math.sqrt(rs);
  const X = x * x + b * x + c;
  const X0 = x0 * x0 + b * x0 + c;
  const Q = Math.sqrt(4 * c - b * b);
  const atanTerm = Math.atan(Q / (2 * x + b));
  const ec = A * (Math.log(x * x / X) + 2 * b / Q * atanTerm
    - b * x0 / X0 * (Math.log((x - x0) * (x - x0) / X) + 2 * (b + 2 * x0) / Q * atanTerm));
  const den = (2 * x + b) * (2 * x + b) + Q * Q;
  const dec = A * (2 / x - (2 * x + b) / X - 4 * b / den
    - b * x0 / X0 * (2 / (x - x0) - (2 * x + b) / X - 4 * (b + 2 * x0) / den));
  const vc = ec - x / 6 * dec;
  return { exc: ex + ec, vxc: vx + vc };
}

// ---------------------------------------------------------------------------
// Potenziale di Hartree: V_H(r) = (1/r)∫₀ʳ 4πr'²ρ dr' + ∫ᵣ^∞ 4πr'ρ dr'
// ---------------------------------------------------------------------------

function hartree(grid, rho) {
  const { N, h, r } = grid;
  const inner = new Float64Array(N);
  const outer = new Float64Array(N);
  // Integrazione cumulativa con la regola dei trapezi corretta (termine di Simpson a 3 punti).
  const a = new Float64Array(N); // integrando in x di Q(r): 4π r³ ρ
  const b = new Float64Array(N); // integrando in x del termine esterno: 4π r² ρ
  for (let i = 0; i < N; i++) {
    a[i] = 4 * Math.PI * r[i] * r[i] * r[i] * rho[i];
    b[i] = 4 * Math.PI * r[i] * r[i] * rho[i];
  }
  // Q(r0): per r → 0 la densità è ~costante, ∫ 4πr²ρ dr ≈ 4π r³ ρ /3
  inner[0] = a[0] / 3;
  for (let i = 1; i < N; i++) {
    // formula a 3 punti di ordine superiore dove possibile
    if (i >= 2) {
      inner[i] = inner[i - 1] + h / 12 * (5 * a[i] + 8 * a[i - 1] - a[i - 2]);
    } else {
      inner[i] = inner[i - 1] + h / 2 * (a[i] + a[i - 1]);
    }
  }
  outer[N - 1] = 0;
  for (let i = N - 2; i >= 0; i--) {
    if (i <= N - 3) {
      outer[i] = outer[i + 1] + h / 12 * (5 * b[i] + 8 * b[i + 1] - b[i + 2]);
    } else {
      outer[i] = outer[i + 1] + h / 2 * (b[i] + b[i + 1]);
    }
  }
  const VH = new Float64Array(N);
  for (let i = 0; i < N; i++) VH[i] = inner[i] / r[i] + outer[i];
  return { VH, charge: inner };
}

// ---------------------------------------------------------------------------
// Potenziale iniziale: modello di Thomas–Fermi (funzione di schermo di Tietz)
// ---------------------------------------------------------------------------

function thomasFermiPotential(grid, Z, Nel) {
  const { N, r } = grid;
  const V = new Float64Array(N);
  const b = 0.8853 * Math.pow(Z, -1 / 3);
  const ion = Z - Nel;
  for (let i = 0; i < N; i++) {
    const x = r[i] / b;
    const phi = 1 / ((1 + 0.53625 * x) * (1 + 0.53625 * x));
    V[i] = -Math.max(Z * phi, ion + 1) / r[i];
  }
  return V;
}

// ---------------------------------------------------------------------------
// Ciclo autoconsistente
// ---------------------------------------------------------------------------

/**
 * @param {number} Z numero atomico
 * @param {Array<{n,l,occ}>} [config] configurazione (default: stato fondamentale)
 * @returns risultato con orbitali, energie, densità e potenziale
 */
export function runSCF(Z, config = groundStateConfiguration(Z), opts = {}) {
  const grid = makeGrid(Z, opts.grid);
  const { N, r } = grid;
  const Nel = config.reduce((s, c) => s + c.occ, 0);
  const tol = opts.tol ?? 1e-8;
  const maxIter = opts.maxIter ?? 300;

  let Vin = thomasFermiPotential(grid, Z, Nel);
  const orbitals = config.map(c => ({ ...c, e: undefined, u: null }));
  const Vext = new Float64Array(N);
  for (let i = 0; i < N; i++) Vext[i] = -Z / r[i];

  // Mescolamento di Anderson (Broyden a un passo) sui potenziali.
  let prevIn = null;
  let prevRes = null;
  let alpha = opts.alpha ?? 0.4;
  let rho = new Float64Array(N);
  let VH = null;
  let dv = Infinity;
  let iter = 0;

  for (iter = 1; iter <= maxIter; iter++) {
    rho = new Float64Array(N);
    for (const orb of orbitals) {
      const res = solveRadial(grid, Vin, orb.n, orb.l, orb.e, Z);
      orb.e = res.e;
      orb.u = res.u;
      for (let i = 0; i < N; i++) rho[i] += orb.occ * res.u[i] * res.u[i] / (4 * Math.PI * r[i] * r[i]);
    }
    const hres = hartree(grid, rho);
    VH = hres.VH;
    const Vout = new Float64Array(N);
    for (let i = 0; i < N; i++) Vout[i] = Vext[i] + VH[i] + lda(rho[i]).vxc;

    // Residuo pesato con r (il potenziale diverge come 1/r vicino al nucleo)
    const res = new Float64Array(N);
    dv = 0;
    for (let i = 0; i < N; i++) {
      res[i] = Vout[i] - Vin[i];
      const w = Math.abs(res[i]) * r[i];
      if (w > dv) dv = w;
    }
    if (dv < tol) break;

    const next = new Float64Array(N);
    if (prevIn && iter > 2) {
      // β di Anderson: minimizza ‖(1−β)R_n + βR_{n−1}‖
      let num = 0;
      let den = 0;
      for (let i = 0; i < N; i++) {
        const d = res[i] - prevRes[i];
        const w = r[i] * r[i] * r[i];
        num += res[i] * d * w;
        den += d * d * w;
      }
      let beta = den > 0 ? num / den : 0;
      if (!Number.isFinite(beta)) beta = 0;
      beta = Math.max(-0.5, Math.min(0.9, beta));
      for (let i = 0; i < N; i++) {
        const vbar = (1 - beta) * Vin[i] + beta * prevIn[i];
        const rbar = (1 - beta) * res[i] + beta * prevRes[i];
        next[i] = vbar + alpha * rbar;
      }
    } else {
      for (let i = 0; i < N; i++) next[i] = Vin[i] + alpha * res[i];
    }
    prevIn = Vin;
    prevRes = res;
    Vin = next;
  }

  // Energia totale di Kohn–Sham
  const eigenSum = orbitals.reduce((s, o) => s + o.occ * o.e, 0);
  const tmp = new Float64Array(N);
  const r2rho = (fn) => { for (let i = 0; i < N; i++) tmp[i] = 4 * Math.PI * r[i] * r[i] * rho[i] * fn(i); return integrate(grid, tmp); };
  const { VH: VHout } = hartree(grid, rho);
  const excArr = new Float64Array(N);
  for (let i = 0; i < N; i++) excArr[i] = lda(rho[i]).exc;
  const kinetic = eigenSum - r2rho(i => Vin[i]);
  const eExt = r2rho(i => Vext[i]);
  const eH = 0.5 * r2rho(i => VHout[i]);
  const eXC = r2rho(i => excArr[i]);
  const total = kinetic + eExt + eH + eXC;

  return {
    Z,
    config,
    grid,
    orbitals: orbitals.map(o => ({ n: o.n, l: o.l, occ: o.occ, e: o.e, u: o.u })),
    V: Vin,
    rho,
    energy: { total, kinetic, external: eExt, hartree: eH, xc: eXC },
    iterations: iter,
    residual: dv,
    converged: dv < tol * 10,
  };
}

/**
 * Potenziale per stati non occupati (orbitali "virtuali"): potenziale autoconsistente con la coda
 * asintotica −(Z−N+1)/r (correzione di Latter), che un elettrone lontano vede realmente.
 */
export function latterPotential(scf) {
  const { grid, V, config, Z } = scf;
  const Nel = config.reduce((s, c) => s + c.occ, 0);
  const zion = Z - Nel + 1;
  const out = new Float64Array(grid.N);
  for (let i = 0; i < grid.N; i++) out[i] = Math.min(V[i], -zion / grid.r[i]);
  return out;
}

/** Orbitale radiale per un (n, l) qualsiasi nel potenziale dell'atomo. */
export function orbitalFor(scf, n, l) {
  const occ = scf.orbitals.find(o => o.n === n && o.l === l);
  if (occ) return { n, l, e: occ.e, u: occ.u, occupied: true, occ: occ.occ };
  if (!scf._latter) scf._latter = latterPotential(scf);
  const res = solveRadial(scf.grid, scf._latter, n, l, undefined, scf.Z);
  return { n, l, e: res.e, u: res.u, occupied: false, occ: 0, converged: res.converged };
}

/** Grandezze radiali: ⟨r⟩, r di massima probabilità, nodi radiali, raggio che racchiude una frazione della probabilità. */
export function radialStats(grid, u) {
  const { N, r, h } = grid;
  let rAvg = 0;
  let pmax = 0;
  let rmax = 0;
  let nodes = [];
  for (let i = 0; i < N; i++) {
    const p = u[i] * u[i];
    rAvg += p * r[i] * r[i] * h;
    if (p > pmax) { pmax = p; rmax = r[i]; }
  }
  // nodi: cambi di segno dove la funzione è significativa
  let umax = 0;
  for (let i = 0; i < N; i++) umax = Math.max(umax, Math.abs(u[i]));
  for (let i = 1; i < N; i++) {
    if (u[i - 1] * u[i] < 0 && Math.max(Math.abs(u[i - 1]), Math.abs(u[i])) > 1e-7 * umax) {
      const t = u[i - 1] / (u[i - 1] - u[i]);
      nodes.push(r[i - 1] + t * (r[i] - r[i - 1]));
    }
  }
  return { rAvg, rMaxProb: rmax, nodes };
}

export function enclosingRadius(grid, u, fraction) {
  const { N, r, h } = grid;
  let total = 0;
  for (let i = 0; i < N; i++) total += u[i] * u[i] * r[i] * h;
  let acc = 0;
  for (let i = 0; i < N; i++) {
    acc += u[i] * u[i] * r[i] * h;
    if (acc >= fraction * total) return r[i];
  }
  return r[N - 1];
}
