// Lennard–Jones per gas nobili, con parametri pubblicati e mixing di Lorentz–Berthelot.
// Taglio a 2,5σ con correzione lineare della forza (come NIST SRS per il fluido LJ): energia e forza
// sono nulle al taglio. Non è un potenziale reattivo e non descrive la chimica.
//
// σ, ε/kB: He, Ne, Ar da Hirschfelder / valori usati per l'argon nel confronto NIST;
// Kr e Xe dai valori classici di gas nobile (arrotondati, non metrologici).

import { KB_EV } from './reactiveData.js';

export const NOBLE_LJ = {
  2: { symbol: 'He', sigma: 2.556, epsilonK: 10.22 },
  10: { symbol: 'Ne', sigma: 2.789, epsilonK: 35.7 },
  18: { symbol: 'Ar', sigma: 3.405, epsilonK: 119.8 },
  36: { symbol: 'Kr', sigma: 3.60, epsilonK: 171 },
  54: { symbol: 'Xe', sigma: 4.055, epsilonK: 229 },
};

const RC_OVER_SIGMA = 2.5;

export function nobleOnly(Z) {
  return Z.length > 0 && Z.every(z => NOBLE_LJ[z]);
}

export function pairLJ(za, zb) {
  const a = NOBLE_LJ[za], b = NOBLE_LJ[zb];
  if (!a || !b) return null;
  const sigma = 0.5 * (a.sigma + b.sigma);
  const epsilon = Math.sqrt(a.epsilonK * b.epsilonK) * KB_EV; // eV
  return { sigma, epsilon };
}

function shifted(r, sigma, epsilon) {
  const rc = RC_OVER_SIGMA * sigma;
  if (!(r > 0) || r >= rc) return { u: 0, f: 0 };
  const sr = sigma / r, sr6 = sr ** 6, sr12 = sr6 * sr6;
  const src = sigma / rc, src6 = src ** 6, src12 = src6 * src6;
  const u = 4 * epsilon * (sr12 - sr6);
  const uc = 4 * epsilon * (src12 - src6);
  const f = 24 * epsilon * (2 * sr12 - sr6) / r;
  const fc = 24 * epsilon * (2 * src12 - src6) / rc;
  return { u: u - uc + (r - rc) * fc, f: f - fc };
}

export function shiftedEnergy(r, sigma, epsilon) {
  return shifted(r, sigma, epsilon).u;
}

export function makeNobleLJProvider() {
  const provider = {
    kind: 'lj',
    info: {},
    reset() {},
    compute(Z, pos, F) {
      F.fill(0);
      if (!nobleOnly(Z)) throw new Error('Il Lennard–Jones pubblicato copre solo He, Ne, Ar, Kr, Xe.');
      let E = 0;
      const n = Z.length;
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const dx = pos[3 * i] - pos[3 * j], dy = pos[3 * i + 1] - pos[3 * j + 1], dz = pos[3 * i + 2] - pos[3 * j + 2];
          const r = Math.hypot(dx, dy, dz);
          const p = pairLJ(Z[i], Z[j]);
          const s = shifted(r, p.sigma, p.epsilon);
          E += s.u;
          if (!s.f) continue;
          const fx = s.f * dx / r, fy = s.f * dy / r, fz = s.f * dz / r;
          F[3 * i] += fx; F[3 * i + 1] += fy; F[3 * i + 2] += fz;
          F[3 * j] -= fx; F[3 * j + 1] -= fy; F[3 * j + 2] -= fz;
        }
      }
      const eps = Z.reduce((s, z) => s + NOBLE_LJ[z].epsilonK, 0) / n;
      const sig = Z.reduce((s, z) => s + NOBLE_LJ[z].sigma, 0) / n;
      provider.info = { epsilonK: eps, sigma: sig, cutoff: RC_OVER_SIGMA };
      return { E, parts: { vdw: E, bond: 0, angle: 0, es: 0, hbond: 0 }, q: new Float64Array(n), bonds: [], hbonds: [] };
    },
  };
  return provider;
}
