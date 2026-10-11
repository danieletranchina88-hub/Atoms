// Exact nonrelativistic hydrogen, infinite nuclear mass, atomic units.
// A coherent superposition is different from an excited-state occupation.
import { makeRandom } from './wavefunction.js';

export const ATOMIC_TIME_AS = 24.188843265857;
export const HYDROGEN_GAP_HA = 3 / 8;
export const HYDROGEN_PERIOD_AS = 2 * Math.PI / HYDROGEN_GAP_HA * ATOMIC_TIME_AS;

export function hydrogenAmplitudes(x, y, z) {
  const r = Math.hypot(x, y, z);
  return [Math.exp(-r) / Math.sqrt(Math.PI), z * Math.exp(-r / 2) / (4 * Math.sqrt(2 * Math.PI))];
}

export function hydrogenDensity(x, y, z, timeAU) {
  const [a, b] = hydrogenAmplitudes(x, y, z);
  return (a * a + b * b) / 2 + a * b * Math.cos(HYDROGEN_GAP_HA * timeAU);
}

/** Samples q=(|1s|²+|2pz|²)/2 analytically. GPU weights ρ/q preserve interference.
 * No positions are interpolated; a point is an integration sample, not an electron.
 * Gamma(k, θ) radii follow from r²|R|²: k=3, θ=1/2 (1s), k=5, θ=1 (2p).
 */
export function sampleHydrogenCoherence(count, seed = 901) {
  const rand = makeRandom(seed);
  const positions = new Float32Array(count * 3);
  const visibility = new Float32Array(count);
  for (let i = 0; i < count; i++) {
    const p = rand() >= 0.5;
    let r = 0;
    for (let k = 0; k < (p ? 5 : 3); k++) r -= Math.log(Math.max(rand(), 1e-15)) * (p ? 1 : 0.5);
    // pz angular distribution has marginal p(cos θ)=3 cos²θ/2.
    const cosTheta = p ? Math.cbrt(2 * rand() - 1) : 2 * rand() - 1;
    const phi = 2 * Math.PI * rand();
    const sinTheta = Math.sqrt(1 - cosTheta * cosTheta);
    const x = r * sinTheta * Math.cos(phi), y = r * sinTheta * Math.sin(phi), z = r * cosTheta;
    positions.set([x, y, z], 3 * i);
    const [a, b] = hydrogenAmplitudes(x, y, z);
    visibility[i] = 2 * a * b / (a * a + b * b || 1);
  }
  return { positions, visibility, count };
}
