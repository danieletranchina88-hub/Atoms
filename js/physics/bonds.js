// Legami: orbitali molecolari LCAO di molecole biatomiche omonucleari e orbitali ibridi.
// L'asse di legame è z; i nuclei A e B stanno in z = −d/2 e z = +d/2.

import { term } from './wavefunction.js';

// Distanze di legame sperimentali (pm) delle molecole biatomiche omonucleari (NIST / CRC Handbook).
export const BOND_LENGTH_PM = {
  1: 74.14, 3: 267.3, 4: 245, 5: 159.0, 6: 124.2, 7: 109.8, 8: 120.7, 9: 141.2,
  11: 307.9, 13: 270.1, 14: 224.6, 15: 189.3, 16: 188.9, 17: 198.8, 19: 390.5,
  29: 222, 33: 210.3, 34: 216.6, 35: 228.1, 37: 417, 47: 253, 51: 234.2, 52: 255.7,
  53: 266.6, 55: 465, 79: 247.2, 83: 266,
};

/** Tipi di orbitale molecolare formati da due orbitali atomici di valenza. */
export const MO_TYPES = [
  { id: 'sigma-s', label: 'σ', sub: 's', ao: 's', antibonding: false,
    text: 'σ legante: i due orbitali s si sommano, la densità si accumula tra i nuclei.' },
  { id: 'sigma-s*', label: 'σ*', sub: 's', ao: 's', antibonding: true,
    text: 'σ* antilegante: i due s si sottraggono; un piano nodale perpendicolare all\'asse passa a metà tra i nuclei.' },
  { id: 'sigma-p', label: 'σ', sub: 'p', ao: 'pz', antibonding: false,
    text: 'σ legante da p<sub>z</sub>: sovrapposizione frontale lungo l\'asse di legame, simmetria cilindrica.' },
  { id: 'sigma-p*', label: 'σ*', sub: 'p', ao: 'pz', antibonding: true,
    text: 'σ* antilegante da p<sub>z</sub>: piano nodale tra i nuclei, lobi esterni più grandi.' },
  { id: 'pi-p', label: 'π', sub: 'p', ao: 'px', antibonding: false,
    text: 'π legante da p<sub>x</sub>: sovrapposizione laterale; il piano nodale contiene l\'asse di legame.' },
  { id: 'pi-p*', label: 'π*', sub: 'p', ao: 'px', antibonding: true,
    text: 'π* antilegante: un piano nodale contiene l\'asse e un secondo lo taglia a metà.' },
];

/**
 * Termini LCAO dell'orbitale molecolare. Coefficienti non normalizzati (la forma non dipende
 * dalla normalizzazione; la normalizzazione è 1/√(2 ± 2S)).
 */
export function moTerms(type, sRadial, pRadial, dBohr) {
  const A = [0, 0, -dBohr / 2];
  const B = [0, 0, dBohr / 2];
  const sign = type.antibonding ? -1 : 1;
  if (type.ao === 's') return [term(sRadial, 0, 0, 1, A), term(sRadial, 0, 0, sign, B)];
  if (type.ao === 'pz') {
    // Il lobo positivo di p_z(A) punta verso B, quello di p_z(B) si allontana da A:
    // la combinazione legante è p_z(A) − p_z(B).
    return [term(pRadial, 1, 0, 1, A), term(pRadial, 1, 0, -sign, B)];
  }
  return [term(pRadial, 1, 1, 1, A), term(pRadial, 1, 1, sign, B)];
}

/** Integrale di sovrapposizione S = ⟨φ_A|φ_B⟩ calcolato numericamente su griglia cilindrica. */
export function overlapIntegral(terms, extent) {
  const [a, b] = terms;
  // S = ∫ φ_A φ_B dV con coefficienti unitari; simmetria cilindrica per σ, cos²φ per π (media = ½·2π).
  const nz = 240;
  const nr = 160;
  const zmax = Math.abs(a.center[2]) + extent;
  const rmax = extent;
  const dz = 2 * zmax / nz;
  const dr = rmax / nr;
  let s = 0;
  const isPi = a.m === 1 && a.l === 1;
  for (let i = 0; i < nz; i++) {
    const z = -zmax + (i + 0.5) * dz;
    for (let j = 0; j < nr; j++) {
      const rho = (j + 0.5) * dr;
      const fa = phiAt(a, rho, 0, z);
      const fb = phiAt(b, rho, 0, z);
      // integrazione in φ: per σ fattore 2π; per π (∝ cos φ · cos φ) fattore π
      s += fa * fb * rho * (isPi ? Math.PI : 2 * Math.PI);
    }
  }
  return s * dz * dr;
}

function phiAt(t, x, y, z) {
  const dx = x - t.center[0];
  const dy = y - t.center[1];
  const dz = z - t.center[2];
  const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
  const R = t.radial.at(r);
  // Y normalizzate: s = 1/√(4π); p_z = √(3/4π)·z/r; p_x = √(3/4π)·x/r
  if (t.l === 0) return R * 0.28209479177387814;
  if (r === 0) return 0;
  const c = 0.4886025119029199;
  return R * c * (t.m === 0 ? dz : dx) / r;
}

/**
 * Ordine qualitativo degli OM di valenza per le biatomiche del 2° periodo.
 * Fino all'azoto l'interazione s–p porta σ(2p) sopra π(2p).
 */
export function valenceMODiagram(Z, valenceElectrons) {
  const sp = Z <= 7;
  const levels = [
    { id: 'sigma-s', name: 'σs', cap: 2, bonding: true },
    { id: 'sigma-s*', name: 'σ*s', cap: 2, bonding: false },
    ...(sp
      ? [{ id: 'pi-p', name: 'πp', cap: 4, bonding: true }, { id: 'sigma-p', name: 'σp', cap: 2, bonding: true }]
      : [{ id: 'sigma-p', name: 'σp', cap: 2, bonding: true }, { id: 'pi-p', name: 'πp', cap: 4, bonding: true }]),
    { id: 'pi-p*', name: 'π*p', cap: 4, bonding: false },
    { id: 'sigma-p*', name: 'σ*p', cap: 2, bonding: false },
  ];
  let left = valenceElectrons;
  let bonding = 0;
  let anti = 0;
  let unpaired = 0;
  for (const lv of levels) {
    lv.electrons = Math.min(lv.cap, left);
    left -= lv.electrons;
    if (lv.bonding) bonding += lv.electrons; else anti += lv.electrons;
    const orbitals = lv.cap / 2;
    // Hund: prima un elettrone per orbitale degenere
    lv.unpaired = lv.electrons <= orbitals ? lv.electrons : lv.cap - lv.electrons;
    unpaired += lv.unpaired;
  }
  return { levels, bondOrder: (bonding - anti) / 2, unpaired, sMixing: sp };
}

// ---------------------------------------------------------------------------
// Orbitali ibridi. Ogni ibrido è una combinazione di s, p, d centrati sullo stesso atomo.
// Componenti: [l, m, coeff] con armoniche reali (m: p_x = 1, p_y = −1, p_z = 0; d_z² = 0, d_x²−y² = 2).
// ---------------------------------------------------------------------------

const r2 = Math.SQRT1_2;
const r3 = 1 / Math.sqrt(3);
const r6 = 1 / Math.sqrt(6);
const r12 = 1 / Math.sqrt(12);

export const HYBRIDS = {
  sp: {
    label: 'sp', geometry: 'lineare', angle: '180°', needsD: false,
    example: 'BeH₂, C in CO₂ e C₂H₂',
    set: [
      [[0, 0, r2], [1, 0, r2]],
      [[0, 0, r2], [1, 0, -r2]],
    ],
  },
  sp2: {
    label: 'sp²', geometry: 'trigonale planare', angle: '120°', needsD: false,
    example: 'BF₃, C in C₂H₄ e nel benzene',
    set: [
      [[0, 0, r3], [1, 1, Math.sqrt(2 / 3)]],
      [[0, 0, r3], [1, 1, -r6], [1, -1, r2]],
      [[0, 0, r3], [1, 1, -r6], [1, -1, -r2]],
    ],
  },
  sp3: {
    label: 'sp³', geometry: 'tetraedrica', angle: '109,47°', needsD: false,
    example: 'CH₄, NH₃, H₂O',
    set: [
      [[0, 0, 0.5], [1, 1, 0.5], [1, -1, 0.5], [1, 0, 0.5]],
      [[0, 0, 0.5], [1, 1, 0.5], [1, -1, -0.5], [1, 0, -0.5]],
      [[0, 0, 0.5], [1, 1, -0.5], [1, -1, 0.5], [1, 0, -0.5]],
      [[0, 0, 0.5], [1, 1, -0.5], [1, -1, -0.5], [1, 0, 0.5]],
    ],
  },
  sp3d: {
    label: 'sp³d', geometry: 'bipiramidale trigonale', angle: '120° e 90°', needsD: true,
    example: 'PCl₅',
    set: [
      [[0, 0, r3], [1, 1, Math.sqrt(2 / 3)]],
      [[0, 0, r3], [1, 1, -r6], [1, -1, r2]],
      [[0, 0, r3], [1, 1, -r6], [1, -1, -r2]],
      [[1, 0, r2], [2, 0, r2]],
      [[1, 0, -r2], [2, 0, r2]],
    ],
  },
  sp3d2: {
    label: 'sp³d²', geometry: 'ottaedrica', angle: '90°', needsD: true,
    example: 'SF₆',
    set: [
      [[0, 0, r6], [1, 1, r2], [2, 0, -r12], [2, 2, 0.5]],
      [[0, 0, r6], [1, 1, -r2], [2, 0, -r12], [2, 2, 0.5]],
      [[0, 0, r6], [1, -1, r2], [2, 0, -r12], [2, 2, -0.5]],
      [[0, 0, r6], [1, -1, -r2], [2, 0, -r12], [2, 2, -0.5]],
      [[0, 0, r6], [1, 0, r2], [2, 0, r3]],
      [[0, 0, r6], [1, 0, -r2], [2, 0, r3]],
    ],
  },
};

/** Termini di un ibrido dato l'insieme delle funzioni radiali {s, p, d}. */
export function hybridTerms(components, radials) {
  return components.map(([l, m, c]) => term(radials[l], l, m, c));
}

/** Energia di un ibrido: ⟨h|H|h⟩ = Σ c_k² ε_k (gli orbitali atomici sono autostati di H). */
export function hybridEnergy(components, energies) {
  return components.reduce((s, [l, , c]) => s + c * c * energies[l], 0);
}
