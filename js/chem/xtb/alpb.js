// Solvente implicito ALPB (analytical linearized Poisson–Boltzmann) per GFN2-xTB, con i parametri per l'acqua.
// Riferimenti: Ehlert, Stahn, Spicher, Grimme, J. Chem. Theory Comput. 17, 4250 (2021); Sigalov, Fenley, Onufriev,
// J. Chem. Phys. 124, 124902 (2006). Porting dell'implementazione di tblite 0.7.0 (solvation/alpb, born, surface, cds):
//  • polarizzazione: Born generalizzato con kernel P16, raggi di Born GBOBC e correzione ALPB di forma;
//  • termine non polare: area accessibile al solvente (griglia di Lebedev a 230 punti, bordo smussato) per tensioni
//    superficiali atomiche, più una correzione per i legami a idrogeno con il solvente proporzionale a q²;
//  • spostamento costante dell'energia libera di solvatazione.
// Unità: bohr e hartree. Le cariche sono quelle atomiche di Mulliken di GFN2-xTB.

import { ALPB_WATER as P } from './alpbData.js';

const ZETA = 1.028, ZETA16 = ZETA / 16, ALPHA_ALPB = 0.571412;
const OBC = [1.0, 0.8, 4.85];
const AATOAU = 1 / 0.529177210903;
const SMOOTH = 0.3 * AATOAU, SURF_OFFSET = 2.0 * AATOAU, TOL_SESP = 1e-6;

export class ALPBWater {
  constructor() {
    this.eps = P.epsilon;
    this.alpbet = ALPHA_ALPB / this.eps;
    this.keps = (1 / this.eps - 1) / (1 + this.alpbet);
  }

  /** Grandezze che dipendono solo dalla geometria: raggi di Born, matrice J, aree accessibili e derivate. */
  setup(Z, pos) {
    const N = Z.length;
    this.N = N; this.Z = Z; this.pos = pos;
    const rvdw = Float64Array.from(Z, z => P.vdwD3[z - 1]);
    this.rvdw = rvdw;
    this.bornRadii(rvdw);
    // matrice di interazione di Born generalizzata (kernel P16)
    const J = new Float64Array(N * N), b = this.brad;
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < i; j++) {
        const r = dist(pos, i, j), ab = Math.sqrt(b[i] * b[j]);
        let a = ab / (ab + ZETA16 * r); a *= a; a *= a; a *= a; a *= a;
        J[i * N + j] = J[j * N + i] = this.keps / (r + ab * a);
      }
      J[i * N + i] = this.keps / b[i];
    }
    if (this.alpbet > 0) {
      this.adet = aDet(N, pos, rvdw);
      const c = this.keps * this.alpbet / this.adet;
      for (let k = 0; k < N * N; k++) J[k] += c;
    }
    this.J = J;
    this.surfaceArea();
    this.tension = Float64Array.from(Z, z => P.tension[z - 1]);
    this.hb = Float64Array.from(Z, (z, i) => P.hbond[z - 1] / (4 * Math.PI * (rvdw[i] + P.probe) ** 2));
    // energia indipendente dalle cariche: tensione superficiale e spostamento costante
    let e = P.gshift;
    for (let i = 0; i < N; i++) e += this.surf[i] * this.tension[i];
    this.eClassical = e;
  }

  /** Potenziale atomico dE/dq. */
  potential(q, out) {
    const N = this.N, J = this.J;
    for (let i = 0; i < N; i++) {
      let s = 0;
      for (let j = 0; j < N; j++) s += J[i * N + j] * q[j];
      out[i] += s + 2 * this.hb[i] * this.surf[i] * q[i];
    }
  }

  /** Energia che dipende dalle cariche: ½ qᵀJq + Σ h_i S_i q_i². */
  energy(q) {
    const N = this.N, J = this.J;
    let e = 0;
    for (let i = 0; i < N; i++) {
      let s = 0;
      for (let j = 0; j < N; j++) s += J[i * N + j] * q[j];
      e += 0.5 * q[i] * s + this.hb[i] * this.surf[i] * q[i] * q[i];
    }
    return e;
  }

  /** Gradiente a cariche fissate (somma in g, hartree/bohr). */
  gradient(q, g) {
    const N = this.N, pos = this.pos, b = this.brad, keps = this.keps;
    const dEdb = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < i; j++) {
        const vx = pos[3 * i] - pos[3 * j], vy = pos[3 * i + 1] - pos[3 * j + 1], vz = pos[3 * i + 2] - pos[3 * j + 2];
        const r = Math.hypot(vx, vy, vz), qq = q[i] * q[j], ab = Math.sqrt(b[i] * b[j]);
        const a1 = ab / (ab + ZETA16 * r);
        let a16 = a1 * a1; a16 *= a16; a16 *= a16; a16 *= a16;
        const fgb = r + ab * a16, d2 = 1 / (fgb * fgb);
        const ap = (1 - ZETA * a1 * a16) * d2 * keps / r * qq;
        g[3 * i] -= ap * vx; g[3 * i + 1] -= ap * vy; g[3 * i + 2] -= ap * vz;
        g[3 * j] += ap * vx; g[3 * j + 1] += ap * vy; g[3 * j + 2] += ap * vz;
        const bp = -0.5 * (r * ZETA / ab * a1 + 1) / ab * a16 * d2;
        dEdb[i] += b[j] * bp * keps * qq;
        dEdb[j] += b[i] * bp * keps * qq;
      }
      dEdb[i] += -0.5 * keps * q[i] * q[i] / (b[i] * b[i]);
    }
    const br = this.brdr;
    for (let i = 0; i < N; i++) {
      const d = dEdb[i];
      if (!d) continue;
      for (let j = 0; j < N; j++) for (let c = 0; c < 3; c++) g[3 * j + c] += d * br[(i * N + j) * 3 + c];
    }
    if (this.alpbet > 0) aDetGradient(N, pos, this.rvdw, keps * this.alpbet, q, g);
    // area accessibile: tensione superficiale + legami a idrogeno
    const ds = this.dsdr;
    for (let i = 0; i < N; i++) {
      const w = this.tension[i] + this.hb[i] * q[i] * q[i];
      if (!w) continue;
      for (let j = 0; j < N; j++) for (let c = 0; c < 3; c++) g[3 * j + c] += w * ds[(i * N + j) * 3 + c];
    }
  }

  /** Raggi di Born con l'integrazione GBOBC (Onufriev, Bashford, Case) e le loro derivate brdr[(i·N+j)·3+c] = ∂b_i/∂R_jc. */
  bornRadii(rvdw) {
    const N = this.N, pos = this.pos;
    const rho = Float64Array.from(this.Z, (z, i) => rvdw[i] * P.descreening[z - 1]);
    const svdw = Float64Array.from(rvdw, r => r - P.bornOffset);
    const psi = new Float64Array(N), dp = new Float64Array(N * N * 3), dtr = new Float64Array(N * 3);
    const add = (i, j, sgnTr, sgnCross, dg, v) => {
      // contributo a ψ_i dalla coppia: ∂ψ_i/∂R_i += sgn·dg·v, ∂ψ_i/∂R_j −= sgn·dg·v
      for (let c = 0; c < 3; c++) { dtr[3 * i + c] += sgnTr * dg * v[c]; dp[(i * N + j) * 3 + c] += sgnCross * dg * v[c]; }
    };
    for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) {
      const v = [pos[3 * i] - pos[3 * j], pos[3 * i + 1] - pos[3 * j + 1], pos[3 * i + 2] - pos[3 * j + 2]];
      const r = Math.hypot(v[0], v[1], v[2]);
      if (r > 35 * AATOAU) continue;
      const ijov = r < rvdw[i] + rho[j], jiov = r < rho[i] + rvdw[j];
      // ψ_i riceve il contributo di j (raggio di descreening ρ_j), ψ_j quello di i; v = R_i − R_j
      const termI = () => {
        if (!ijov) { const [gi, dgi] = nonOverlap(r, rho[j]); psi[i] += gi; add(i, j, 1, -1, dgi, v); }
        else if (r + rho[j] > rvdw[i]) { const [gi, dgi] = overlap(r, rho[j], rvdw[i]); psi[i] += gi; add(i, j, 1, -1, dgi, v); }
      };
      const termJ = () => {
        if (!jiov) { const [gj, dgj] = nonOverlap(r, rho[i]); psi[j] += gj; add(j, i, -1, 1, dgj, v); }
        else if (r + rho[i] > rvdw[j]) { const [gj, dgj] = overlap(r, rho[i], rvdw[j]); psi[j] += gj; add(j, i, -1, 1, dgj, v); }
      };
      termI(); termJ();
    }
    for (let i = 0; i < N; i++) for (let c = 0; c < 3; c++) dp[(i * N + i) * 3 + c] = dtr[3 * i + c];
    const brad = new Float64Array(N), c1 = P.bornScale;
    for (let i = 0; i < N; i++) {
      const s1 = 1 / svdw[i], v1 = 1 / rvdw[i], s2 = 0.5 * svdw[i];
      const br = psi[i] * s2;
      const a2 = br * (OBC[2] * br - OBC[1]);
      const arg = br * (OBC[0] + a2);
      const darg = 2 * a2 + OBC[0] + OBC[2] * br * br;
      const th = Math.tanh(arg), ch = Math.cosh(arg);
      brad[i] = c1 / (s1 - v1 * th);
      const den = ch * (s1 - v1 * th);
      const dpsi = c1 * s2 * v1 * darg / (den * den);
      for (let j = 0; j < N; j++) for (let c = 0; c < 3; c++) dp[(i * N + j) * 3 + c] *= dpsi;
    }
    this.brad = brad; this.brdr = dp;
  }

  /** Area accessibile al solvente per atomo (bohr²) e derivate dsdr[(i·N+j)·3+c] = ∂S_i/∂R_jc. */
  surfaceArea() {
    const N = this.N, pos = this.pos, w = SMOOTH, w3 = w * w * w;
    const ah0 = 0.5, ah1 = 3 / (4 * w), ah3 = -1 / (4 * w3);
    const vdwsa = Float64Array.from(this.rvdw, r => r + P.probe);
    const t1 = Float64Array.from(vdwsa, v => (v - w) ** 2), t2 = Float64Array.from(vdwsa, v => (v + w) ** 2);
    const wrp = Float64Array.from(vdwsa, v => {
      const f = (r) => (0.25 / w + 3 * ah3 * (0.2 * r * r - 0.5 * r * v + v * v / 3)) * r * r * r;
      return f(v + w) - f(v - w);
    });
    let maxsa = 0;
    for (const v of vdwsa) maxsa = Math.max(maxsa, v);
    const cut2 = (2 * (w + maxsa) + SURF_OFFSET) ** 2;
    const grid = P.lebedev230;
    const surf = new Float64Array(N), ds = new Float64Array(N * N * 3);
    const nb = new Int32Array(N), gx = new Float64Array(3 * N), gi = new Int32Array(N);
    for (let i = 0; i < N; i++) {
      let nn = 0;
      for (let j = 0; j < N; j++) if (j !== i && dist2(pos, i, j) < cut2) nb[nn++] = j;
      const rs = vdwsa[i], xi = pos[3 * i], yi = pos[3 * i + 1], zi = pos[3 * i + 2];
      let si = 0;
      const base = i * N * 3;
      for (let p = 0; p < grid.length; p++) {
        const [ux, uy, uz, wp] = grid[p];
        const px = xi + rs * ux, py = yi + rs * uy, pz = zi + rs * uz;
        let sp = 1, ni = 0;
        for (let k = 0; k < nn; k++) {
          const j = nb[k];
          const tx = px - pos[3 * j], ty = py - pos[3 * j + 1], tz = pz - pos[3 * j + 2];
          const tj2 = tx * tx + ty * ty + tz * tz;
          if (tj2 >= t2[j]) continue;
          if (tj2 <= t1[j]) { sp = 0; break; }
          const sq = Math.sqrt(tj2), u = sq - vdwsa[j], h3 = ah3 * u * u;
          const sij = ah0 + (ah1 + h3) * u, dsij = (ah1 + 3 * h3) / (sij * sq);
          sp *= sij;
          gi[ni] = j; gx[3 * ni] = dsij * tx; gx[3 * ni + 1] = dsij * ty; gx[3 * ni + 2] = dsij * tz; ni++;
        }
        if (sp <= TOL_SESP) continue;
        const wsa = wp * 4 * Math.PI * wrp[i] * sp;
        si += wsa;
        for (let k = 0; k < ni; k++) {
          const j = gi[k];
          for (let c = 0; c < 3; c++) { const d = wsa * gx[3 * k + c]; ds[base + 3 * i + c] += d; ds[base + 3 * j + c] -= d; }
        }
      }
      surf[i] = si;
    }
    this.surf = surf; this.dsdr = ds;
  }
}

function dist2(pos, i, j) { const x = pos[3 * i] - pos[3 * j], y = pos[3 * i + 1] - pos[3 * j + 1], z = pos[3 * i + 2] - pos[3 * j + 2]; return x * x + y * y + z * z; }
function dist(pos, i, j) { return Math.sqrt(dist2(pos, i, j)); }

/** Contributo GBOBC di una sfera di raggio ρ che non si sovrappone: g e (dg/dr)/r. */
function nonOverlap(r, rho) {
  const r1 = 1 / r, ap = r + rho, am = r - rho, ab = ap * am, rhab = rho / ab;
  const lnab = 0.5 * Math.log(am / ap) * r1;
  return [rhab + lnab, -2 * rhab / ab + (rhab - lnab) * r1 * r1];
}
/** Contributo di una sfera di raggio ρ che si sovrappone alla sfera di van der Waals di raggio rv. */
function overlap(r, rho, rv) {
  const r1 = 1 / r, r12 = 0.5 * r1, r24 = r12 * r12, ap = r + rho, am = r - rho;
  const rh1 = 1 / rv, rhr1 = 1 / ap, aprh1 = ap * rh1, lnab = Math.log(aprh1);
  const g = rh1 - rhr1 + r12 * (0.5 * am * (rhr1 - rh1 * aprh1) - lnab);
  const dg = (rhr1 * rhr1 * (1 - 0.25 * am * r1 * (1 + aprh1 * aprh1)) + rho * r24 * (rhr1 - rh1 * aprh1) + r12 * (r1 * lnab - rhr1)) * r1;
  return [g, dg];
}

/** Tensore d'inerzia delle sfere di van der Waals e dimensione elettrostatica A_det della correzione ALPB. */
function inertia(N, pos, rad) {
  let vol = 0;
  const cen = [0, 0, 0];
  for (let i = 0; i < N; i++) { const r3 = rad[i] ** 3; vol += r3; for (let c = 0; c < 3; c++) cen[c] += pos[3 * i + c] * r3; }
  for (let c = 0; c < 3; c++) cen[c] /= vol;
  const I = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < N; i++) {
    const r2 = rad[i] ** 2, r3 = r2 * rad[i];
    const v = [pos[3 * i] - cen[0], pos[3 * i + 1] - cen[1], pos[3 * i + 2] - cen[2]];
    const d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) I[a][b] += r3 * ((a === b ? d2 + 0.4 * r2 : 0) - v[a] * v[b]);
  }
  const det = I[0][0] * (I[1][1] * I[2][2] - I[1][2] * I[2][1]) - I[0][1] * (I[1][0] * I[2][2] - I[1][2] * I[2][0]) + I[0][2] * (I[1][0] * I[2][1] - I[1][1] * I[2][0]);
  return { I, vol, cen, adet: Math.sqrt(Math.cbrt(det) / (0.4 * vol)) };
}
function aDet(N, pos, rad) { return inertia(N, pos, rad).adet; }
function aDetGradient(N, pos, rad, kEps, q, g) {
  const { I, vol, cen, adet } = inertia(N, pos, rad);
  let qt = 0;
  for (let i = 0; i < N; i++) qt += q[i];
  const f = (250 / (48 * vol ** 3 * adet ** 5)) * (-0.5 * kEps * qt * qt / (adet * adet));
  const D = [
    [I[0][0] * (I[1][1] + I[2][2]) - I[0][1] ** 2 - I[0][2] ** 2, I[0][1] * I[2][2] - I[0][2] * I[1][2], I[0][2] * I[1][1] - I[0][1] * I[2][1]],
    [I[0][1] * I[2][2] - I[0][2] * I[1][2], I[1][1] * (I[0][0] + I[2][2]) - I[0][1] ** 2 - I[1][2] ** 2, I[0][0] * I[1][2] - I[0][1] * I[0][2]],
    [I[0][2] * I[1][1] - I[0][1] * I[2][1], I[0][0] * I[1][2] - I[0][1] * I[0][2], I[2][2] * (I[0][0] + I[1][1]) - I[0][2] ** 2 - I[1][2] ** 2],
  ];
  for (let i = 0; i < N; i++) {
    const r3 = rad[i] ** 3, v = [pos[3 * i] - cen[0], pos[3 * i + 1] - cen[1], pos[3 * i + 2] - cen[2]];
    for (let a = 0; a < 3; a++) g[3 * i + a] += r3 * f * (D[a][0] * v[0] + D[a][1] * v[1] + D[a][2] * v[2]);
  }
}
