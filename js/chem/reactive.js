// Campo di forze reattivo per la dinamica molecolare della sandbox.
//
// I legami si formano e si rompono in modo continuo. Nessuna connettività è fissata in anticipo:
// l'energia dipende solo dalle posizioni degli atomi ed è derivabile ovunque, così le forze sono
// esatte e l'energia totale si conserva.
//
//  1. Legame covalente (Abell–Tersoff–Brenner):
//       E = f_c(r) [ V_R(r) − b(n) V_A(r) ],  V_R = D/(S−1) e^{−λ₁(r−r₁)},  V_A = DS/(S−1) e^{−λ₂(r−r₁)}
//     Con b = n^{p(S−1)/S}, un legame di ordine n ha energia D·nᵖ e lunghezza r₁ − c·ln n: esattamente
//     la relazione di Pauling tra ordine e lunghezza di legame, con D, p e c dalle tabelle sperimentali.
//  2. Ordine di legame dalla valenza: n_ij è la valenza lasciata libera dagli altri legami di i e di j.
//     I legami π in eccesso (sistemi coniugati, benzene) sono ridotti in proporzione: n = 1,5.
//  3. Angoli dal teorema di Coulson: cos θ₀ = −1/λ, con λ = (numero sterico − 1) come negli ibridi
//     spᵏ; ogni doppietto solitario stringe l'angolo (regola VSEPR: 109,5° → 107° → 104,5°).
//  4. Cariche parziali per equalizzazione dell'elettronegatività (split-charge equilibration,
//     Nistor et al. 2006): E_Q = Σ (χq + ½Jq²) + Σ q_i q_j γ_ij(r) + Σ_legami ½κ x²/w.
//     Le cariche si spostano solo lungo i legami, quindi i frammenti separati restano neutri.
//  5. van der Waals (UFF): E = D[(x/r)¹² − 2(x/r)⁶], escluso tra atomi legati e tra atomi 1-3.
//
// Unità interne: Å, eV, amu, fs.

import { ATOM_PARAMS, bondData, KJ_PER_EV, KE, S_BRENNER } from './reactiveData.js';
import { PAULING } from './elementData.js';

const RC = 8;                  // raggio di taglio di Coulomb e van der Waals (Å)
const SWITCH_IN = 0.35;        // inizio dello spegnimento del legame: r₁ + 0,35 Å
const SWITCH_OUT = 0.75;       // fine: r₁ + 0,75 Å
const SHIELD = 2;               // raggio di schermatura di Coulomb: a = 2 k_e/√(J_i J_j)
const KAPPA = 1.5;              // durezza di legame dello split-charge (eV)
const HB_D = 0.19;              // legame a idrogeno (DREIDING): profondità (eV) …
const HB_R = 2.75;              // … e distanza donatore–accettore di equilibrio (Å)
const HB_ELEMENTS = new Set([7, 8, 9]);            // durezza di legame dello split-charge (eV)
const ANGLE_K = 3.0;           // costante di piegamento (eV)
const LP_SQUEEZE = 0.042;      // Δcos θ₀ per doppietto solitario (109,5° → 104,5° con due doppietti)
const KCAL_EV = 0.0433641;
const G_MID = 0.25, G_W = 0.07;  // peso di conteggio: funzione di Fermi centrata a r₁ + 0,25 Å

const countWeight = (y) => 1 / (1 + Math.exp((y - G_MID) / G_W));

// funzioni lisce: max(0,x), min(a,b), max(a,b)
const sp = (x, d) => 0.5 * (x + Math.sqrt(x * x + d * d));
const spd = (x, d) => 0.5 * (1 + x / Math.sqrt(x * x + d * d));

// spegnimento dolce tra R_ON e RC (polinomio di 5° grado: continuo con derivate prima e seconda)
const R_ON = 6;
function taper(r) {
  if (r <= R_ON) return [1, 0];
  if (r >= RC) return [0, 0];
  const t = (r - R_ON) / (RC - R_ON), t2 = t * t, t3 = t2 * t;
  return [1 - 10 * t3 + 15 * t3 * t - 6 * t3 * t2, (-30 * t2 + 60 * t3 - 30 * t3 * t) / (RC - R_ON)];
}

export class ReactiveFF {
  constructor({ kappa = KAPPA } = {}) {
    this.kappa = kappa;
    this.cache = new Map();
    this.splitPrev = new Map();
    this.tolerance = 1e-7;
  }

  /** Tipi atomici e tabelle di coppia (van der Waals UFF, schermatura di Coulomb). */
  types(Z) {
    if (!this.typeOf) {
      const els = Object.keys(ATOM_PARAMS).map(Number);
      this.ntypes = els.length;
      this.typeOf = new Map(els.map((z, i) => [z, i]));
      const nt = els.length;
      this.ljX = new Float64Array(nt * nt); this.ljD = new Float64Array(nt * nt); this.a3t = new Float64Array(nt * nt);
      for (let a = 0; a < nt; a++) {
        for (let b = 0; b < nt; b++) {
          const pa = ATOM_PARAMS[els[a]], pb = ATOM_PARAMS[els[b]];
          this.ljX[a * nt + b] = Math.sqrt(pa[5] * pb[5]);
          // H···N/O/F: il contatto è un legame a idrogeno (termine dedicato); il raggio UFF di H,
          // pensato per C–H, lo renderebbe fortemente repulsivo a 1,9 Å
          const hPolar = (els[a] === 1 && HB_ELEMENTS.has(els[b])) || (els[b] === 1 && HB_ELEMENTS.has(els[a]));
          if (hPolar) this.ljX[a * nt + b] = 2.4;
          this.ljD[a * nt + b] = Math.sqrt(pa[6] * pb[6]) * KCAL_EV;
          this.a3t[a * nt + b] = Math.pow(SHIELD * KE / Math.sqrt(pa[4] * pb[4]), 3);
        }
      }
    }
    if (!this.tcache || this.tcache.length !== Z.length || this.tcacheZ !== Z) {
      this.tcache = Int32Array.from(Z, z => this.typeOf.get(z));
      this.tcacheZ = Z;
    }
    return this.tcache;
  }

  /** Parametri di coppia (in eV e Å), calibrati sull'energia di legame misurata. */
  pairParams(Za, Zb) {
    const key = Za <= Zb ? Za * 200 + Zb : Zb * 200 + Za;
    let pp = this.cache.get(key);
    if (pp) return pp;
    const pa = ATOM_PARAMS[Za], pb = ATOM_PARAMS[Zb];
    const bd = bondData(Za, Zb, PAULING);
    const S = S_BRENNER;
    const r1 = bd.r1;
    // parte elettrostatica del legame nella molecola biatomica isolata (split-charge a due siti):
    // la si sottrae dall'energia misurata, perché nel modello la fornisce il termine di Coulomb
    const g1 = countWeight(0);
    const a3 = Math.pow(SHIELD * KE / Math.sqrt(pa[4] * pb[4]), 3);
    const gam = KE / Math.cbrt(r1 * r1 * r1 + a3) * taper(r1)[0];
    const dchi = pa[3] - pb[3];
    const esStab = dchi * dchi / (2 * (pa[4] + pb[4] - 2 * gam + this.kappa / g1));
    const Dexp = bd.D1 / KJ_PER_EV;
    const D = Math.max(0.3 * Dexp, Dexp - esStab);
    pp = {
      D, Dexp, r1, p: bd.p, c: bd.c, source: bd.source,
      A: D / (S - 1), B: D * S / (S - 1),
      l1: bd.p / bd.c, l2: bd.p / (S * bd.c),
      gb: bd.p * (S - 1) / S,
      // legami parziali (n < 1): esponente almeno 1,05, altrimenti tanti legami frazionari "renderebbero"
      // più di pochi legami interi e il carbonio formerebbe grumi sovracoordinati invece di catene e anelli
      gs: Math.max(bd.p, 1.05) * (S - 1) / S,
      rin: r1 + SWITCH_IN, rout: r1 + SWITCH_OUT,
    };
    this.cache.set(key, pp);
    return pp;
  }

  /**
   * Energia e forze.
   * @param Z Int32Array/array di numeri atomici; pos Float64Array 3N (Å); F Float64Array 3N (eV/Å, uscita)
   */
  compute(Z, pos, F) {
    const N = Z.length;
    F.fill(0);
    const V = new Float64Array(N), EV = new Float64Array(N), CHI = new Float64Array(N), JJ = new Float64Array(N);
    const XV = new Float64Array(N), DV = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      const pa = ATOM_PARAMS[Z[i]];
      V[i] = pa[0]; EV[i] = pa[1]; CHI[i] = pa[3]; JJ[i] = pa[4]; XV[i] = pa[5]; DV[i] = pa[6] * KCAL_EV;
    }

    // --- coppie entro il raggio di taglio: van der Waals e interazione di Coulomb schermata ------
    // (tutto in array tipizzati riutilizzati tra un passo e l'altro: è il ciclo più costoso)
    const cap = (N * (N - 1)) >> 1;
    if (!this.buf || this.buf.cap < cap) {
      this.buf = {
        cap, I: new Int32Array(cap), J: new Int32Array(cap), R: new Float64Array(cap),
        LJ: new Float64Array(cap), dLJ: new Float64Array(cap), G: new Float64Array(cap), dG: new Float64Array(cap),
        W: new Float64Array(cap),
      };
    }
    const { I: allI, J: allJ, R: allR, LJ, dLJ, G: gam, dG: dgam, W: Wv } = this.buf;
    const T = this.types(Z);
    const nt = this.ntypes;
    const { ljX, ljD, a3t } = this;
    let nAll = 0;
    const bI = [], bJ = [], bR = [], bP = [];
    const RC2 = RC * RC;
    const S6C = 5.618655692729767; // (1/0,75)⁶
    for (let i = 0; i < N; i++) {
      const xi = pos[3 * i], yi = pos[3 * i + 1], zi = pos[3 * i + 2];
      const ti = T[i] * nt;
      for (let j = i + 1; j < N; j++) {
        const dx = pos[3 * j] - xi, dy = pos[3 * j + 1] - yi, dz = pos[3 * j + 2] - zi;
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= RC2) continue;
        const r = Math.sqrt(r2);
        const k = nAll++;
        allI[k] = i; allJ[k] = j; allR[k] = r; Wv[k] = 1;
        // spegnimento dolce tra R_ON e RC
        let Tp = 1, dTp = 0;
        if (r > R_ON) {
          const t = (r - R_ON) / (RC - R_ON), t2 = t * t, t3 = t2 * t;
          Tp = 1 - 10 * t3 + 15 * t3 * t - 6 * t3 * t2;
          dTp = (-30 * t2 + 60 * t3 - 30 * t3 * t) / (RC - R_ON);
        }
        const tp = ti + T[j];
        // van der Waals UFF (con continuazione lineare sotto 0,75 x, contro forze infinite negli urti)
        const X = ljX[tp], D = ljD[tp], rc = 0.75 * X;
        let E, dE;
        if (r >= rc) {
          const q2 = X * X / r2, s6 = q2 * q2 * q2;
          E = D * (s6 * s6 - 2 * s6); dE = -12 * D * (s6 * s6 - s6) / r;
        } else {
          const dEc = -12 * D * (S6C * S6C - S6C) / rc;
          E = D * (S6C * S6C - 2 * S6C) + dEc * (r - rc); dE = dEc;
        }
        LJ[k] = E * Tp; dLJ[k] = dE * Tp + E * dTp;
        // Coulomb schermato: γ = k_e / (r³ + a³)^{1/3}, con a = k_e/√(J_i J_j) (γ → √(J_i J_j) per r → 0)
        const u = r2 * r + a3t[tp], cb = Math.cbrt(u);
        const g0 = KE / cb;
        gam[k] = g0 * Tp; dgam[k] = -g0 * r2 / u * Tp + g0 * dTp;
        if (V[i] > 0 && V[j] > 0) {
          const pp = this.pairParams(Z[i], Z[j]);
          if (r < pp.rout) { bI.push(i); bJ.push(j); bR.push(r); bP.push(pp); }
        }
      }
    }
    const M = bI.length;
    const unit = (i, j, r) => [(pos[3 * j] - pos[3 * i]) / r, (pos[3 * j + 1] - pos[3 * i + 1]) / r, (pos[3 * j + 2] - pos[3 * i + 2]) / r];
    const addRadial = (i, j, r, dEdr) => {
      // F_i = −∂E/∂x_i = +dE/dr · u,  F_j = −dE/dr · u   (u da i verso j)
      const f = dEdr / r;
      const dx = pos[3 * j] - pos[3 * i], dy = pos[3 * j + 1] - pos[3 * i + 1], dz = pos[3 * j + 2] - pos[3 * i + 2];
      F[3 * i] += f * dx; F[3 * i + 1] += f * dy; F[3 * i + 2] += f * dz;
      F[3 * j] -= f * dx; F[3 * j + 1] -= f * dy; F[3 * j + 2] -= f * dz;
    };

    // --- 1. termini di coppia dei legami ---------------------------------------------------
    // Due pesi per ogni coppia:
    //  • gs = ordine σ "geometrico" di Pauling, e^{−(r−r₁)/c} limitato a 1: conserva la valenza lungo un
    //    trasferimento di atomo (A + B–C → A–B + C, come nel modello BEBO di Johnston e Parr);
    //  • w = occupazione di un vicino, ≈ 1 fino a r₁ + 0,25 Å (funzione di Fermi): decide quanta valenza
    //    resta per i legami π, così un legame σ appena allungato non "libera" legami multipli altrove.
    const fc = new Float64Array(M), dfc = new Float64Array(M), w = new Float64Array(M), dw = new Float64Array(M);
    const gs = new Float64Array(M), dgs = new Float64Array(M);
    const VR = new Float64Array(M), dVR = new Float64Array(M), VA = new Float64Array(M), dVA = new Float64Array(M);
    const adj = Array.from({ length: N }, () => []);
    const Ncnt = new Float64Array(N), Nsig = new Float64Array(N);
    for (let p = 0; p < M; p++) {
      const pp = bP[p], r = bR[p];
      if (r <= pp.rin) { fc[p] = 1; dfc[p] = 0; } else {
        const L = pp.rout - pp.rin, t = Math.PI * (r - pp.rin) / L;
        fc[p] = 0.5 * (1 + Math.cos(t)); dfc[p] = -0.5 * Math.PI / L * Math.sin(t);
      }
      const y = r - pp.r1;
      const g = countWeight(y), dg = -g * (1 - g) / G_W;
      w[p] = fc[p] * g; dw[p] = dfc[p] * g + fc[p] * dg;
      const u = y / pp.c;
      const e = Math.exp(-sp(u, 0.05)), de = -e * spd(u, 0.05) / pp.c;
      gs[p] = fc[p] * e; dgs[p] = dfc[p] * e + fc[p] * de;
      VR[p] = pp.A * Math.exp(-pp.l1 * y); dVR[p] = -pp.l1 * VR[p];
      VA[p] = pp.B * Math.exp(-pp.l2 * y); dVA[p] = -pp.l2 * VA[p];
      adj[bI[p]].push(p); adj[bJ[p]].push(p);
      Ncnt[bI[p]] += w[p]; Ncnt[bJ[p]] += w[p];
      Nsig[bI[p]] += gs[p]; Nsig[bJ[p]] += gs[p];
    }

    // --- 2. ordini di legame dalla valenza ---------------------------------------------------
    // σ = min(valenza libera di i, di j) ∈ [0, 1];  π = min(valenza libera − 1) ∈ [0, 2]
    const D1 = 0.02;
    const smin = (a1, a2, e2) => {
      const dd = a1 - a2, sq = Math.sqrt(dd * dd + e2);
      return [0.5 * (a1 + a2 - sq), 0.5 * (1 - dd / sq), 0.5 * (1 + dd / sq)];
    };
    const sig = new Float64Array(M), dsI = new Float64Array(M), dsJ = new Float64Array(M);
    const pi0 = new Float64Array(M), dpI = new Float64Array(M), dpJ = new Float64Array(M);
    for (let p = 0; p < M; p++) {
      const i = bI[p], j = bJ[p];
      const [ms, mI, mJ] = smin(V[i] - Nsig[i] + gs[p], V[j] - Nsig[j] + gs[p], 0.0016);
      sig[p] = sp(ms, D1) - sp(ms - 1, D1);
      const ds = spd(ms, D1) - spd(ms - 1, D1);
      dsI[p] = ds * mI; dsJ[p] = ds * mJ;
      const [mp, qI, qJ] = smin(V[i] - Ncnt[i] + w[p] - 1, V[j] - Ncnt[j] + w[p] - 1, 0.0016);
      pi0[p] = sp(mp, D1) - sp(mp - 2, D1);
      const dp = spd(mp, D1) - spd(mp - 2, D1);
      dpI[p] = dp * qI; dpJ[p] = dp * qJ;
    }
    // conservazione della valenza σ: se Σ σ supera la valenza (tre centri, come H₃), si riscalano
    const sig0 = Float64Array.from(sig);
    const Tsig = new Float64Array(N), usc = new Float64Array(N), duT = new Float64Array(N);
    for (let p = 0; p < M; p++) { Tsig[bI[p]] += fc[p] * sig0[p]; Tsig[bJ[p]] += fc[p] * sig0[p]; }
    for (let i = 0; i < N; i++) {
      const dd = Tsig[i] - V[i], sq = Math.sqrt(dd * dd + 0.0004);
      const den = 0.5 * (Tsig[i] + V[i] + sq);
      usc[i] = V[i] / den;
      duT[i] = -V[i] / (den * den) * 0.5 * (1 + dd / sq);
    }
    const uu = new Float64Array(M), duI = new Float64Array(M), duJ = new Float64Array(M);
    for (let p = 0; p < M; p++) {
      [uu[p], duI[p], duJ[p]] = smin(usc[bI[p]], usc[bJ[p]], 0.0004);
      sig[p] = sig0[p] * uu[p];
    }
    // legami π in conflitto (coniugazione): ridotti in proporzione al bilancio di valenza dell'atomo
    const Bud = new Float64Array(N), dBud = new Float64Array(N), Pi = new Float64Array(N);
    for (let p = 0; p < M; p++) { Pi[bI[p]] += w[p] * pi0[p]; Pi[bJ[p]] += w[p] * pi0[p]; }
    const s = new Float64Array(N), dsPi = new Float64Array(N), dsB = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      Bud[i] = sp(V[i] - Ncnt[i], D1); dBud[i] = spd(V[i] - Ncnt[i], D1);
      const dd = Pi[i] - Bud[i], sq = Math.sqrt(dd * dd + 0.0025);
      const den = 0.5 * (Pi[i] + Bud[i] + sq);
      const dden_dPi = 0.5 * (1 + dd / sq), dden_dB = 0.5 * (1 - dd / sq);
      s[i] = Bud[i] / den;
      dsPi[i] = -Bud[i] / (den * den) * dden_dPi;
      dsB[i] = 1 / den - Bud[i] / (den * den) * dden_dB;
    }
    const tt = new Float64Array(M), dtI = new Float64Array(M), dtJ = new Float64Array(M);
    const n = new Float64Array(M), b = new Float64Array(M), db = new Float64Array(M);
    let Ebond = 0;
    const Sn = new Float64Array(N);
    for (let p = 0; p < M; p++) {
      const i = bI[p], j = bJ[p], pp = bP[p];
      [tt[p], dtI[p], dtJ[p]] = smin(s[i], s[j], 0.0004);
      n[p] = sig[p] + pi0[p] * tt[p];
      // b = (n² + ε²)^{G/2} normalizzato a b(1) = 1; G passa con continuità dal regime σ (n < 1) a quello π
      const q2 = (n[p] * n[p] + 0.0025) / 1.0025, lq = Math.log(q2);
      const sw = 1 / (1 + Math.exp((n[p] - 1) / 0.08)), dsw = -sw * (1 - sw) / 0.08;
      const G = pp.gb + (pp.gs - pp.gb) * sw, dG = (pp.gs - pp.gb) * dsw;
      b[p] = Math.exp(0.5 * G * lq);
      db[p] = b[p] * (G * n[p] / (n[p] * n[p] + 0.0025) + 0.5 * lq * dG);
      Ebond += fc[p] * (VR[p] - b[p] * VA[p]);
      Sn[i] += w[p] * n[p]; Sn[j] += w[p] * n[p];
    }

    // --- 3. angoli (Coulson + VSEPR) -------------------------------------------------------
    const LP = new Float64Array(N), dLP = new Float64Array(N), lam = new Float64Array(N), dlam = new Float64Array(N);
    const c0 = new Float64Array(N), dc0lp = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      const arg = 0.5 * (EV[i] - Sn[i]);
      LP[i] = sp(arg, 0.02); dLP[i] = spd(arg, 0.02);
      const SN = Ncnt[i] + LP[i];
      lam[i] = 1 + sp(SN - 2, 0.02) - sp(SN - 4, 0.02);
      dlam[i] = spd(SN - 2, 0.02) - spd(SN - 4, 0.02);
      const dd = LP[i] - 2, sq = Math.sqrt(dd * dd + 0.01);
      const lpc = 0.5 * (LP[i] + 2 - sq);
      dc0lp[i] = LP_SQUEEZE * 0.5 * (1 - dd / sq);
      c0[i] = -1 / lam[i] + LP_SQUEEZE * lpc;
    }
    const a_w = new Float64Array(M), a_fc = new Float64Array(M), a_r = new Float64Array(M), a_n = new Float64Array(M);
    const a_c0 = new Float64Array(N);
    let Eangle = 0;
    for (let i = 0; i < N; i++) {
      const L = adj[i];
      if (L.length < 2) continue;
      for (let x = 0; x < L.length; x++) {
        const p = L[x];
        const j = bI[p] === i ? bJ[p] : bI[p];
        const uj = unit(i, j, bR[p]);
        for (let y = x + 1; y < L.length; y++) {
          const q = L[y];
          const k = bI[q] === i ? bJ[q] : bI[q];
          const uk = unit(i, k, bR[q]);
          const cs = uj[0] * uk[0] + uj[1] * uk[1] + uj[2] * uk[2];
          const diff = cs - c0[i];
          const ww = w[p] * w[q];
          Eangle += ANGLE_K * ww * diff * diff;
          a_w[p] += ANGLE_K * w[q] * diff * diff;
          a_w[q] += ANGLE_K * w[p] * diff * diff;
          a_c0[i] -= 2 * ANGLE_K * ww * diff;
          const G = 2 * ANGLE_K * ww * diff;
          for (let c = 0; c < 3; c++) {
            const gj = (uk[c] - cs * uj[c]) / bR[p];
            const gk = (uj[c] - cs * uk[c]) / bR[q];
            F[3 * j + c] -= G * gj;
            F[3 * k + c] -= G * gk;
            F[3 * i + c] += G * (gj + gk);
          }
        }
      }
    }

    // --- 3b. legame a idrogeno D–H···A (D, A = N, O, F), forma di DREIDING (Mayo, Olafson, Goddard 1990):
    //     E = D_hb [5 (R₀/R)¹² − 6 (R₀/R)¹⁰] cos⁴θ,  R = distanza D···A,  θ = angolo D–H···A
    const bondIndex = new Map();
    for (let p = 0; p < M; p++) bondIndex.set(bI[p] * N + bJ[p], p);
    const acceptors = [];
    for (let i = 0; i < N; i++) if (HB_ELEMENTS.has(Z[i])) acceptors.push(i);
    let Ehb = 0;
    const hbonds = [];
    if (acceptors.length > 1) {
      for (let p = 0; p < M; p++) {
        let h, d;
        if (Z[bI[p]] === 1 && HB_ELEMENTS.has(Z[bJ[p]])) { h = bI[p]; d = bJ[p]; }
        else if (Z[bJ[p]] === 1 && HB_ELEMENTS.has(Z[bI[p]])) { h = bJ[p]; d = bI[p]; }
        else continue;
        if (w[p] < 1e-6) continue;
        const rhd = bR[p];
        const uhd = [(pos[3 * d] - pos[3 * h]) / rhd, (pos[3 * d + 1] - pos[3 * h + 1]) / rhd, (pos[3 * d + 2] - pos[3 * h + 2]) / rhd];
        for (const a of acceptors) {
          if (a === d) continue;
          const dx = pos[3 * a] - pos[3 * d], dy = pos[3 * a + 1] - pos[3 * d + 1], dz = pos[3 * a + 2] - pos[3 * d + 2];
          const R = Math.sqrt(dx * dx + dy * dy + dz * dz);
          if (R > 4.5) continue;
          const hx = pos[3 * a] - pos[3 * h], hy = pos[3 * a + 1] - pos[3 * h + 1], hz = pos[3 * a + 2] - pos[3 * h + 2];
          const rha = Math.sqrt(hx * hx + hy * hy + hz * hz);
          const uha = [hx / rha, hy / rha, hz / rha];
          const cs = uhd[0] * uha[0] + uhd[1] * uha[1] + uhd[2] * uha[2];
          if (cs >= 0) continue;
          // H–A non devono essere legati (altrimenti è un legame covalente, non un legame a idrogeno)
          const qa = bondIndex.get(Math.min(h, a) * N + Math.max(h, a));
          const fha = qa === undefined ? 0 : fc[qa];
          // distanza: 12-10 con spegnimento tra 4,0 e 4,5 Å
          const s2 = (HB_R * HB_R) / (R * R), s10 = s2 * s2 * s2 * s2 * s2;
          let V = HB_D * (5 * s10 * s2 - 6 * s10), dV = HB_D * (-60 * s10 * s2 + 60 * s10) / R;
          if (R > 4) {
            const t = (R - 4) / 0.5, t3 = t * t * t;
            const S = 1 - 10 * t3 + 15 * t3 * t - 6 * t3 * t * t, dS = (-30 * t * t + 60 * t3 - 30 * t3 * t) / 0.5;
            dV = dV * S + V * dS; V *= S;
          }
          const c4 = cs * cs * cs * cs, dc4 = 4 * cs * cs * cs;
          const W = w[p] * (1 - fha);
          const E = V * c4 * W;
          Ehb += E;
          if (E < -0.04) hbonds.push(h, a);
          a_w[p] += V * c4 * (1 - fha);
          if (qa !== undefined) a_fc[qa] += -V * c4 * w[p];
          // radiale D···A
          const fR = dV * c4 * W / R;
          F[3 * d] += fR * dx; F[3 * d + 1] += fR * dy; F[3 * d + 2] += fR * dz;
          F[3 * a] -= fR * dx; F[3 * a + 1] -= fR * dy; F[3 * a + 2] -= fR * dz;
          // angolare: cos θ = u_hd · u_ha
          const G = V * dc4 * W;
          for (let c = 0; c < 3; c++) {
            const gd = (uha[c] - cs * uhd[c]) / rhd;
            const ga = (uhd[c] - cs * uha[c]) / rha;
            F[3 * d + c] -= G * gd;
            F[3 * a + c] -= G * ga;
            F[3 * h + c] += G * (gd + ga);
          }
        }
      }
    }

    // --- 4. coppie "speciali": legate (1-2) o con un vicino comune (1-3) -------------------------
    // van der Waals è escluso tra atomi legati e 1-3 (fattori lisci). Coulomb agisce tra tutte le
    // coppie: togliere dei termini renderebbe la matrice delle cariche non definita positiva.
    const keyOf = (i, j) => (i < j ? i * N + j : j * N + i);
    const special = new Map();
    for (let p = 0; p < M; p++) special.set(keyOf(bI[p], bJ[p]), { i: bI[p], j: bJ[p], p, c13: 0, tri: [], k: -1 });
    for (let k = 0; k < N; k++) {
      const L = adj[k];
      for (let x = 0; x < L.length; x++) {
        const p = L[x], i = bI[p] === k ? bJ[p] : bI[p];
        for (let y = x + 1; y < L.length; y++) {
          const q = L[y], j = bI[q] === k ? bJ[q] : bI[q];
          const key = keyOf(i, j);
          let e = special.get(key);
          if (!e) { e = { i: Math.min(i, j), j: Math.max(i, j), p: -1, c13: 0, tri: [], k: -1 }; special.set(key, e); }
          e.c13 += w[p] * w[q];
          e.tri.push(p, q);
        }
      }
    }
    // indice di ogni coppia speciale nell'elenco di tutte le coppie (le righe sono ordinate per i)
    const byRow = Array.from({ length: N }, () => []);
    for (const e of special.values()) byRow[e.i].push(e);
    const stamp = new Int32Array(N).fill(-1);
    for (let k = 0, i = -1; k <= nAll; k++) {
      if (k === nAll || allI[k] !== i) {
        if (i >= 0) {
          for (const e of byRow[i]) e.k = stamp[e.j];
          for (let t = k - 1; t >= 0 && allI[t] === i; t--) stamp[allJ[t]] = -1;
        }
        if (k === nAll) break;
        i = allI[k];
      }
      stamp[allJ[k]] = k;
    }
    const spm1 = sp(-1, 0.05);
    const specList = [];
    for (const e of special.values()) {
      if (e.k < 0) continue;
      e.f = e.p >= 0 ? fc[e.p] : 0;
      e.C = e.tri.length ? e.c13 - sp(e.c13 - 1, 0.05) + spm1 : 0;
      e.dC = e.tri.length ? 1 - spd(e.c13 - 1, 0.05) : 0;
      Wv[e.k] = (1 - e.f) * (1 - e.C);
      specList.push(e);
    }

    let Evdw = 0;
    for (let k = 0; k < nAll; k++) Evdw += Wv[k] * LJ[k];

    // --- 6. cariche: split-charge equilibration ------------------------------------------------
    // la carica si sposta solo lungo i legami veri: peso ω = w·n (un contatto tra molecole sature ha n ≈ 0)
    const act = [], omega = new Float64Array(M);
    for (let p = 0; p < M; p++) { omega[p] = w[p] * n[p]; if (omega[p] > 1e-6) act.push(p); }
    const nb = act.length;
    const kap = new Float64Array(nb), x = new Float64Array(nb), diag = new Float64Array(nb);
    for (let a = 0; a < nb; a++) {
      const p = act[a];
      const e = special.get(keyOf(bI[p], bJ[p]));
      kap[a] = this.kappa / omega[p];
      x[a] = this.splitPrev.get(keyOf(bI[p], bJ[p])) ?? 0;
      diag[a] = JJ[bI[p]] + JJ[bJ[p]] - 2 * (e.k >= 0 ? gam[e.k] : 0) + kap[a];
    }
    const q = new Float64Array(N), mu = new Float64Array(N);
    const chargesFrom = (xv) => {
      q.fill(0);
      for (let a = 0; a < nb; a++) { const p = act[a]; q[bI[p]] -= xv[a]; q[bJ[p]] += xv[a]; }
      return q;
    };
    const potential = (qv, withChi) => {
      for (let i = 0; i < N; i++) mu[i] = (withChi ? CHI[i] : 0) + JJ[i] * qv[i];
      for (let k = 0; k < nAll; k++) { mu[allI[k]] += gam[k] * qv[allJ[k]]; mu[allJ[k]] += gam[k] * qv[allI[k]]; }
      return mu;
    };
    // ∂E/∂x_a = μ_j − μ_i + κ x_a = 0: gradiente coniugato precondizionato (partendo dal passo precedente)
    const residual = (xv, withChi, out) => {
      potential(chargesFrom(xv), withChi);
      for (let a = 0; a < nb; a++) { const p = act[a]; out[a] = mu[bJ[p]] - mu[bI[p]] + kap[a] * xv[a]; }
      return out;
    };
    if (nb > 0) {
      const g = residual(x, true, new Float64Array(nb));
      const z = new Float64Array(nb), d = new Float64Array(nb), Ad = new Float64Array(nb);
      for (let a = 0; a < nb; a++) { z[a] = g[a] / diag[a]; d[a] = -z[a]; }
      let rz = 0;
      for (let a = 0; a < nb; a++) rz += g[a] * z[a];
      for (let it = 0; it < 80; it++) {
        let gn = 0;
        for (let a = 0; a < nb; a++) gn += g[a] * g[a];
        if (Math.sqrt(gn / nb) < this.tolerance) break;
        residual(d, false, Ad);
        let dAd = 0;
        for (let a = 0; a < nb; a++) dAd += d[a] * Ad[a];
        if (dAd <= 0) break;
        const alpha = rz / dAd;
        let rz2 = 0;
        for (let a = 0; a < nb; a++) { x[a] += alpha * d[a]; g[a] += alpha * Ad[a]; z[a] = g[a] / diag[a]; rz2 += g[a] * z[a]; }
        const beta = rz2 / rz;
        rz = rz2;
        for (let a = 0; a < nb; a++) d[a] = -z[a] + beta * d[a];
      }
    }
    this.splitPrev.clear();
    for (let a = 0; a < nb; a++) { const p = act[a]; this.splitPrev.set(keyOf(bI[p], bJ[p]), x[a]); }
    chargesFrom(x);
    const charges = Float64Array.from(q);
    let Ees = 0;
    for (let i = 0; i < N; i++) Ees += CHI[i] * q[i] + 0.5 * JJ[i] * q[i] * q[i];
    for (let k = 0; k < nAll; k++) {
      const i = allI[k], j = allJ[k];
      const qq = q[i] * q[j];
      Ees += qq * gam[k];
      const f = (qq * dgam[k] + Wv[k] * dLJ[k]) / allR[k];
      if (f === 0) continue;
      const dx = pos[3 * j] - pos[3 * i], dy = pos[3 * j + 1] - pos[3 * i + 1], dz = pos[3 * j + 2] - pos[3 * i + 2];
      F[3 * i] += f * dx; F[3 * i + 1] += f * dy; F[3 * i + 2] += f * dz;
      F[3 * j] -= f * dx; F[3 * j + 1] -= f * dy; F[3 * j + 2] -= f * dz;
    }
    for (let a = 0; a < nb; a++) {
      const p = act[a];
      Ees += 0.5 * kap[a] * x[a] * x[a];
      const g = -0.5 * this.kappa * x[a] * x[a] / (omega[p] * omega[p]);
      a_w[p] += g * n[p];
      a_n[p] += g * w[p];
    }
    // derivate dei fattori di esclusione rispetto ai pesi dei legami
    for (const e of specList) {
      const k = e.k;
      if (e.p >= 0) a_fc[e.p] += -(1 - e.C) * LJ[k];
      if (e.tri.length) {
        const a_C = -(1 - e.f) * LJ[k];
        const a_c13 = a_C * e.dC;
        for (let t = 0; t < e.tri.length; t += 2) {
          const p = e.tri[t], qd = e.tri[t + 1];
          a_w[p] += a_c13 * w[qd];
          a_w[qd] += a_c13 * w[p];
        }
      }
    }

    // --- 7. propagazione all'indietro delle derivate (regola della catena) -------------------------
    const a_b = new Float64Array(M);
    for (let p = 0; p < M; p++) {
      a_fc[p] += VR[p] - b[p] * VA[p];
      a_r[p] += fc[p] * (dVR[p] - b[p] * dVA[p]);
      a_b[p] = -fc[p] * VA[p];
    }
    const a_N = new Float64Array(N), a_LP = new Float64Array(N), a_Sn = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      const a_lam = a_c0[i] / (lam[i] * lam[i]);
      a_LP[i] += a_c0[i] * dc0lp[i];
      const a_SN = a_lam * dlam[i];
      a_N[i] += a_SN; a_LP[i] += a_SN;
      a_Sn[i] = -0.5 * a_LP[i] * dLP[i];
    }
    const a_s = new Float64Array(N), a_Nsig = new Float64Array(N), a_gs = new Float64Array(M);
    const a_pi0 = new Float64Array(M), a_sigF = new Float64Array(M), a_usc = new Float64Array(N), a_T = new Float64Array(N);
    for (let p = 0; p < M; p++) {
      const i = bI[p], j = bJ[p];
      a_w[p] += (a_Sn[i] + a_Sn[j]) * n[p];
      a_n[p] += (a_Sn[i] + a_Sn[j]) * w[p] + a_b[p] * db[p];
      // n = σ·u + π·t
      a_sigF[p] = a_n[p];
      a_pi0[p] += a_n[p] * tt[p];
      const a_t = a_n[p] * pi0[p];
      a_s[i] += a_t * dtI[p]; a_s[j] += a_t * dtJ[p];
      const a_u = a_n[p] * sig0[p];
      a_usc[i] += a_u * duI[p]; a_usc[j] += a_u * duJ[p];
    }
    for (let i = 0; i < N; i++) a_T[i] = a_usc[i] * duT[i];
    for (let p = 0; p < M; p++) {
      const i = bI[p], j = bJ[p];
      // σ₀ entra in n (×u) e nella somma T (×f_c)
      const a_sig = a_sigF[p] * uu[p] + (a_T[i] + a_T[j]) * fc[p];
      a_fc[p] += (a_T[i] + a_T[j]) * sig0[p];
      a_Nsig[i] -= a_sig * dsI[p]; a_Nsig[j] -= a_sig * dsJ[p];
      a_gs[p] += a_sig * (dsI[p] + dsJ[p]);
    }
    const a_Pi = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      a_Pi[i] = a_s[i] * dsPi[i];
      a_N[i] += -a_s[i] * dsB[i] * dBud[i];
    }
    for (let p = 0; p < M; p++) {
      const i = bI[p], j = bJ[p];
      const aP = a_Pi[i] + a_Pi[j];
      a_w[p] += aP * pi0[p];
      a_pi0[p] += aP * w[p];
      const aI = a_pi0[p] * dpI[p], aJ = a_pi0[p] * dpJ[p];
      a_N[i] -= aI; a_N[j] -= aJ;
      a_w[p] += aI + aJ;
    }
    for (let p = 0; p < M; p++) {
      a_w[p] += a_N[bI[p]] + a_N[bJ[p]];
      a_gs[p] += a_Nsig[bI[p]] + a_Nsig[bJ[p]];
      a_r[p] += a_w[p] * dw[p] + a_gs[p] * dgs[p] + a_fc[p] * dfc[p];
      addRadial(bI[p], bJ[p], bR[p], a_r[p]);
    }

    const bonds = [];
    for (let p = 0; p < M; p++) {
      if (w[p] * n[p] > 0.02) bonds.push({ i: bI[p], j: bJ[p], n: n[p], w: w[p] });
    }
    return {
      E: Ebond + Eangle + Evdw + Ees + Ehb,
      parts: { bond: Ebond, angle: Eangle, vdw: Evdw, es: Ees, hbond: Ehb },
      q: charges, bonds, hbonds,
    };
  }
}
