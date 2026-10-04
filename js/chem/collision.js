import { ELEMENTS } from '../physics/elements.js';
import { MV2 } from './md.js';

/** Two molecular COMs with prescribed relative translational energy and zero total momentum.
 * Molecule coordinates in Å. Energy in eV, impact parameter in Å. No reaction is imposed.
 */
export function prepareCollision(left, right, { energy = 0.6, distance = 7, impact = 0, rotation = 0 } = {}) {
  if (![energy, distance, impact, rotation].every(Number.isFinite) || energy < 0 || distance < 2 || Math.abs(impact) > distance)
    throw new Error('Urto: energia ≥0, separazione ≥2 Å e parametro d’urto entro la separazione.');
  const center = m => {
    if (!m.Z.length || m.pos.length !== 3 * m.Z.length || !Array.from(m.pos).every(Number.isFinite)) throw new Error('Geometria dell’urto non valida.');
    const mass = m.Z.map(z => Number(ELEMENTS[z - 1].mass)), M = mass.reduce((s, v) => s + v, 0);
    const c = [0, 1, 2].map(k => mass.reduce((s, v, i) => s + v * m.pos[3 * i + k], 0) / M);
    return { mass, M, c };
  };
  const a = center(left), b = center(right), total = a.M + b.M, reduced = a.M * b.M / total;
  const speed = Math.sqrt(2 * energy / (reduced * MV2));
  const longitudinal = Math.sqrt(Math.max(0, distance * distance - impact * impact));
  const atoms = [];
  for (const [m, com, other, side] of [[left, a, b, -1], [right, b, a, 1]]) {
    for (let i = 0; i < m.Z.length; i++) {
      let [x, y, z] = [0, 1, 2].map(k => m.pos[3 * i + k] - com.c[k]);
      if (side === 1) { const th = rotation * Math.PI / 180; [x, y] = [x * Math.cos(th) - y * Math.sin(th), x * Math.sin(th) + y * Math.cos(th)]; }
      atoms.push({ Z: m.Z[i], formal: m.formal?.[i] ?? 0,
        pos: [x + side * other.M / total * longitudinal, y + side * other.M / total * impact, z],
        vel: [-side * other.M / total * speed, 0, 0] });
    }
  }
  for (let i = 0; i < left.Z.length; i++) for (let j = left.Z.length; j < atoms.length; j++) {
    const r = Math.hypot(...atoms[i].pos.map((v, c) => v - atoms[j].pos[c]));
    if (r < 1.2) throw new Error('Le molecole si sovrappongono: aumenta la separazione iniziale.');
  }
  const extent = Math.max(...atoms.flatMap(a => a.pos.map(Math.abs)));
  return { atoms, energy, speed, reducedMass: reduced, massLeft: a.M, massRight: b.M, distance, impact, box: Math.max(16, 2 * extent + 10) };
}
