import assert from 'node:assert/strict';
import { hydrogenDensity, hydrogenAmplitudes, sampleHydrogenCoherence, HYDROGEN_GAP_HA, HYDROGEN_PERIOD_AS } from '../js/physics/quantumDynamics.js';
import { makeGrid, solveRadial } from '../js/physics/scf.js';
import { prepareOrbitalScene } from '../js/render/orbitalScene.js';

const near = (got, want, tol, message) => assert.ok(Math.abs(got - want) <= tol, `${message}: ${got} vs ${want}`);
near(HYDROGEN_PERIOD_AS, 405.3, 0.1, 'Hydrogen beat period (as)');
const T = 2 * Math.PI / HYDROGEN_GAP_HA;
for (const t of [0, T / 4, T / 2, T]) {
  // Spatial quadrature: dV = 2π r² dr d(cos θ). Normalization and dipole are independent checks.
  let norm = 0, dipole = 0;
  const dr = 0.035, dc = 2 / 100;
  for (let r = dr / 2; r < 30; r += dr) for (let j = 0; j < 100; j++) {
    const c = -1 + (j + 0.5) * dc;
    const rho = hydrogenDensity(r * Math.sqrt(1 - c * c), 0, r * c, t);
    assert.ok(rho >= -1e-16);
    const mass = rho * 2 * Math.PI * r * r * dr * dc;
    norm += mass;
    dipole += r * c * mass;
  }
  near(norm, 1, 0.0001, 'Normalized at every phase');
  near(dipole, 256 / (243 * Math.sqrt(2)) * Math.cos(HYDROGEN_GAP_HA * t), 0.0002, 'Analytic <z>(t)');
}
for (const p of [[0, 0, 1], [1, 2, -3], [0, 0, 0]]) {
  const [a, b] = hydrogenAmplitudes(...p);
  near(hydrogenDensity(...p, T / 4), (a * a + b * b) / 2, 1e-15, 'Quarter cycle has no interference');
  near(hydrogenDensity(...p, 0), hydrogenDensity(...p, T), 1e-15, 'Periodicity');
  near(hydrogenDensity(...p, 0), hydrogenDensity(p[0], p[1], -p[2], T / 2), 1e-15, 'Reflection at half cycle');
}
const samples = sampleHydrogenCoherence(100000);
let radius = 0, weight = 0, dipole = 0;
for (let i = 0; i < samples.count; i++) {
  const p = samples.positions.subarray(3 * i, 3 * i + 3);
  radius += Math.hypot(...p);
  assert.ok(Math.abs(samples.visibility[i]) <= 1);
  weight += 1 + samples.visibility[i];
  dipole += p[2] * (1 + samples.visibility[i]);
}
near(radius / samples.count, 3.25, 0.025, 'Mixture radial mean');
near(weight / samples.count, 1, 0.006, 'GPU importance weights conserve probability');
near(dipole / samples.count, 256 / (243 * Math.sqrt(2)), 0.04, 'GPU weights recover coherent dipole');
assert.deepEqual(sampleHydrogenCoherence(10), sampleHydrogenCoherence(10), 'Deterministic samples');

const grid = makeGrid(1);
const V = Float64Array.from(grid.r, r => -1 / r);
const orbitals = [[1, 0, 0, 2], [2, 1, 0, 1]].map(([n, l, m, occ]) => ({
  n, l, m, occ, subshell: `${n},${l}`, grid, u: solveRadial(grid, V, n, l, undefined, 1).u,
}));
const scene = await prepareOrbitalScene({ orbitals, points: 6000, quality: 40, enclosed: 0.9, render: 'both' });
assert.equal(scene.cloud.count, 6000);
assert.equal(scene.cloud.groups.filter(g => g === 0).length, 4000, 'Occupancy-weighted cloud');
assert.ok(scene.surfaces.some(s => s.group === 0 && s.sign === 1), 'Compact 1s survives its own grid');
assert.ok(scene.surfaces.some(s => s.group === 1 && s.sign === -1), 'Negative pz lobe');
for (const mesh of scene.surfaces) {
  assert.ok(mesh.positions.every(Number.isFinite) && mesh.normals.every(Number.isFinite));
  if (mesh.group === 1) for (let i = 2; i < mesh.positions.length; i += 3) assert.ok(mesh.sign * mesh.positions[i] >= -0.01, 'pz nodal plane');
}
const empty = await prepareOrbitalScene({ orbitals: [], points: 6000, quality: 40, enclosed: 0.9, render: 'both' });
assert.equal(empty.cloud, null);
assert.equal(empty.surfaces.length, 0);
const cloud = await prepareOrbitalScene({ orbitals, points: 6000, quality: 40, enclosed: 0.9, render: 'cloud' });
assert.equal(cloud.surfaces.length, 0, 'Cloud mode does not calculate surfaces');
await assert.rejects(prepareOrbitalScene({ orbitals, points: 6000, quality: 40, enclosed: 0.9, render: 'surface' }, () => true), /sostituita/);
console.log('OK quantum dynamics: normalization, dipole, Born weights, orbital grids and cancellation');
