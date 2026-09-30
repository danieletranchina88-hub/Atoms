// IUPAC Nernst equation: https://goldbook.iupac.org/terms/view/09068
// Reduction stoichiometry, products positive. Electrons are omitted from Q.
import { R_GAS, FARADAY } from './labData.js';

const sp = (id, nu, phase = 'aq') => ({ id, nu, phase });
export function halfStoichiometry(h) {
  if (h.metal) return [sp(h.ox, -1), sp(h.red, 1, 's')];
  switch (h.id) {
    case 'H2': return [sp('H⁺', -2), sp('H₂', 1, 'g')];
    case 'O2': return [sp('O₂', -1, 'g'), sp('H⁺', -4), sp('H₂O', 2, 'l')];
    case 'I2': return [sp('I₂', -1, 's'), sp('I⁻', 2)];
    case 'Br2': return [sp('Br₂', -1, 'l'), sp('Br⁻', 2)];
    case 'Cl2': return [sp('Cl₂', -1, 'g'), sp('Cl⁻', 2)];
    case 'F2': return [sp('F₂', -1, 'g'), sp('F⁻', 2)];
    case 'Fe3': return [sp('Fe³⁺', -1), sp('Fe²⁺', 1)];
    case 'MnO4': return [sp('MnO₄⁻', -1), sp('H⁺', -8), sp('Mn²⁺', 1), sp('H₂O', 4, 'l')];
    default: throw new Error('Semireazione priva di stechiometria verificata.');
  }
}

export function electrodePotential(h, activities = {}, T = 298.15) {
  if (!Number.isFinite(T) || T <= 0) throw new Error('Temperatura assoluta non valida.');
  let lnQ = 0;
  for (const s of halfStoichiometry(h)) {
    if (s.phase === 's' || s.phase === 'l') continue; // pure standard phases, a = 1
    const a = activities[s.id] ?? 1;
    if (!Number.isFinite(a) || a <= 0) throw new Error(`Attività di ${s.id} non valida.`);
    lnQ += s.nu * Math.log(a);
  }
  return { E: h.E - R_GAS * T / (h.n * FARADAY) * lnQ, lnQ };
}
const gcd = (a, b) => b ? gcd(b, a % b) : a;

/** Fixed written direction: left oxidizes, right reduces. E may be negative.
 * No current/rate is inferred from equilibrium thermodynamics.
 * E° values are tabulated at 298.15 K; varying T without E°(T) is only a frozen-E° model.
 */
export function redoxCell(left, right, { leftActivities = {}, rightActivities = {}, T = 298.15 } = {}) {
  const n = left.n * right.n / gcd(left.n, right.n);
  const L = electrodePotential(left, leftActivities, T);
  const R = electrodePotential(right, rightActivities, T);
  const lnQ = n / right.n * R.lnQ - n / left.n * L.lnQ;
  const E0 = right.E - left.E, E = R.E - L.E;
  return { anode: left, cathode: right, left: L, right: R, n, E0, E, lnQ, Q: Math.exp(lnQ),
    dG0: -n * FARADAY * E0 / 1000, dG: -n * FARADAY * E / 1000,
    lnK: n * FARADAY * E0 / (R_GAS * T), reverse: E < -1e-10,
    atEquilibrium: Math.abs(E) <= 1e-10, leftFactor: n / left.n, rightFactor: n / right.n };
}

export function cellEquation(left, right) {
  const n = left.n * right.n / gcd(left.n, right.n);
  const terms = [];
  for (const [h, factor, side] of [[left, -n / left.n, 'sx'], [right, n / right.n, 'dx']]) {
    for (const s of halfStoichiometry(h)) terms.push({ ...s, nu: s.nu * factor, side });
  }
  const format = t => `${Math.abs(t.nu) === 1 ? '' : Math.abs(t.nu)}${t.id}(${t.phase}; ${t.side})`;
  return `${terms.filter(t => t.nu < 0).map(format).join(' + ')} → ${terms.filter(t => t.nu > 0).map(format).join(' + ')}`;
}
