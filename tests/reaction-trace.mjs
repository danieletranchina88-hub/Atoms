import assert from 'node:assert/strict';
import { Worker } from 'node:worker_threads';
import { geometry, populationChange, TrajectoryRecorder } from '../js/chem/reactionTrace.js';
import { prepareCollision } from '../js/chem/collision.js';
import { MV2 } from '../js/chem/md.js';
import { ELEMENTS } from '../js/physics/elements.js';
import { makeHFProvider } from '../js/chem/aimd.js';
import { makeMindo3Provider } from '../js/chem/mindo3.js';
import { makeGFN2Provider } from '../js/chem/xtb/provider.js';
import { computeGrid, gridSlice } from '../js/chem/densityWorker.js';
import { PRESETS } from '../js/chem/sandboxData.js';

const near = (a, b, e = 1e-9) => assert.ok(Math.abs(a - b) < e, `${a} ≠ ${b}`);
const mass = z => Number(ELEMENTS[z - 1].mass);
const left = { Z: [8, 1, 1], pos: [0, 0, 0, 0.76, 0.59, 0, -0.76, 0.59, 0] };
const right = { Z: [17], pos: [0.3, -0.8, 1.2], formal: [-1] };
for (const energy of [0, 0.01, 0.6, 5]) for (const impact of [-2, 0, 2]) {
  const c = prepareCollision(left, right, { energy, distance: 8, impact, rotation: 73 });
  const momentum = [0, 1, 2].map(k => c.atoms.reduce((s, a) => s + mass(a.Z) * a.vel[k], 0));
  momentum.forEach(v => near(v, 0));
  const com = group => [0, 1, 2].map(k => group.reduce((s, a) => s + mass(a.Z) * a.pos[k], 0) / group.reduce((s, a) => s + mass(a.Z), 0));
  const ca = com(c.atoms.slice(0, 3)), cb = com(c.atoms.slice(3));
  near(Math.hypot(...ca.map((v, k) => cb[k] - v)), 8);
  near(cb[1] - ca[1], impact);
  com(c.atoms).forEach(v => near(v, 0));
  near(c.atoms.reduce((s, a) => s + 0.5 * mass(a.Z) * a.vel.reduce((t, v) => t + v * v, 0) * MV2, 0), energy);
  const distances = a => a.slice(1, 3).map(h => Math.hypot(...h.pos.map((v, k) => v - a[0].pos[k])));
  distances(c.atoms).forEach((v, i) => near(v, distances(left.Z.map((Z, i) => ({ Z, pos: left.pos.slice(3 * i, 3 * i + 3) })))[i]));
}
assert.throws(() => prepareCollision(left, right, { energy: -1 }));
assert.throws(() => prepareCollision(left, right, { impact: 8, distance: 3 }));
near(geometry([1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1], [0, 1, 2, 3]).angle, 90);
near(Math.abs(geometry([1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1], [0, 1, 2, 3]).dihedral), 90);
assert.equal(geometry([0, 0, 0, 0, 0, 0], [0, 1, 0]).angle, null);
const change = populationChange([0.6, -0.4, -0.2], [0.1, -0.1, 0]);
near(change.total, 0); near(change.gained, 0.5);

// Spin density must have the expected electron count in the actual orbital basis.
for (const make of [makeHFProvider, makeMindo3Provider, makeGFN2Provider]) {
  const provider = make({ multiplicity: 2, fixedMultiplicity: 2 });
  const F = new Float64Array(3); const r = provider.compute([1], new Float64Array(3), F);
  near(r.spin.reduce((s, v) => s + v, 0), 1, 1e-6);
  const w = provider.wavefunction();
  assert.ok(w.Ps && w.orbitals.length);
  const request = w.kind === 'gauss' ? { mode: 'gauss', atoms: w.atoms, basisName: w.basisName, Z: [1], pos: [0, 0, 0] }
    : { mode: w.kind, Z: w.Z, pos: w.pos, first: w.first };
  const g = computeGrid({ ...request, P: w.P, Ps: w.Ps, what: 'spin', box: 5, res: 49 });
  const integral = g.values.reduce((s, v) => s + v, 0) * (g.step / 0.52917721090) ** 3;
  near(integral, 1, 0.035);
  const d = computeGrid({ ...request, P: w.P, what: 'difference', box: 5, res: 29, reference: { ...request, P: w.P } });
  assert.ok(d.values.every(v => v === 0));
  const slice = gridSlice(g, { axis: 'x', offset: 0 });
  assert.equal(slice.values.length, 49 ** 2); near(slice.offset, 0);
}
const singlet = makeGFN2Provider({ multiplicity: 1 });
const sr = singlet.compute([1, 1], Float64Array.from([-0.37, 0, 0, 0.37, 0, 0]), new Float64Array(6));
near(sr.spin.reduce((s, v) => s + Math.abs(v), 0), 0);
assert.ok(singlet.wavefunction().Ps.every(v => v === 0));

// Worker test: record a real UHF exchange, resolve both changes by atomic identity,
// replay the original SCF matrices, and verify inspection cannot change the live state.
const url = new URL('../js/chem/sandboxWorker.js', import.meta.url).href;
const worker = new Worker(`const {parentPort}=require('node:worker_threads'); globalThis.onmessage=null;
globalThis.postMessage=m=>parentPort.postMessage(m);
import(${JSON.stringify(url)}).then(()=>parentPort.on('message',data=>globalThis.onmessage({data})));`, { eval: true });
function wait(predicate) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => { cleanup(); reject(new Error('trace worker timeout')); }, 30000);
    const handler = m => { if (predicate(m)) { cleanup(); resolve(m); } else if (m.type === 'error') { cleanup(); reject(new Error(m.text)); } };
    const cleanup = () => { clearTimeout(timer); worker.off('message', handler); };
    worker.on('message', handler);
  });
}
try {
  await wait(m => m.type === 'frame');
  let pending = wait(m => m.type === 'frame' && m.N === 3 && m.stats.t === 0);
  worker.postMessage({ type: 'preset', preset: PRESETS.find(p => p.id === 'aimd-h3') });
  const initial = await pending; const epoch = initial.trace.epoch;
  assert.equal(initial.trace.count, 1);
  pending = wait(m => m.type === 'frame' && m.stats.t > 49);
  worker.postMessage({ type: 'step', n: 250 });
  const current = await pending;
  const events = current.trace.events;
  assert.ok(events.some(e => e.kind === 'formed' && e.i === 0 && e.j === 1), JSON.stringify(events));
  assert.ok(events.some(e => e.kind === 'broken' && e.i === 1 && e.j === 2), JSON.stringify(events));
  assert.equal(current.trace.count, 251);
  assert.ok(current.spin && current.forces.length === 9 && current.vel.length === 9);
  pending = wait(m => m.type === 'trace-frame');
  worker.postMessage({ type: 'trace-inspect', epoch, id: 0, reqId: 1 });
  const replay = (await pending).frame;
  assert.equal(replay.stats.t, 0); assert.deepEqual(replay.pos, initial.pos); assert.ok(replay.wave.P.length);
  pending = wait(m => m.type === 'frame' && m.stats.t === current.stats.t);
  assert.deepEqual((await pending).pos, current.pos);
  pending = wait(m => m.type === 'wave');
  worker.postMessage({ type: 'wave', traceId: 0, reqId: 3 });
  assert.deepEqual((await pending).wave.P, replay.wave.P);
  pending = wait(m => m.type === 'trace-export'); worker.postMessage({ type: 'trace-export' });
  const exported = (await pending).data;
  assert.equal(exported.frames.length, 251); assert.equal(exported.units.forces, 'eV/angstrom');
  assert.ok(exported.frames.every(f => !('wave' in f)));
  const tiny = new TrajectoryRecorder({ maxFrames: 3, maxBytes: 10_000 });
  for (const f of exported.frames.slice(0, 8)) tiny.add(f);
  assert.equal(tiny.frames.length, 3); assert.equal(tiny.evicted, 5);
  pending = wait(m => m.type === 'frame' && m.trace.epoch !== epoch);
  worker.postMessage({ type: 'collision', left: '[H-]', right: '[H+]', energy: 0.6, distance: 6, impact: 0.4 });
  const collision = await pending;
  assert.equal(collision.N, 2); assert.equal(collision.stats.forceField, 'gfn2'); assert.equal(collision.stats.thermostat, false);
  near(collision.stats.Ekin, 0.6, 1e-8); assert.equal(collision.stats.charge, 0);
  pending = wait(m => m.type === 'frame' && m.N === 0);
  worker.postMessage({ type: 'clear' }); await pending;
  worker.postMessage({ type: 'set', paused: true });
  pending = wait(m => m.type === 'frame' && m.N === 40);
  worker.postMessage({ type: 'add', symbol: 'Ar', count: 40 });
  assert.equal((await pending).stats.forceField, 'lj');
  pending = wait(m => m.type === 'frame' && m.N === 42);
  worker.postMessage({ type: 'add', smiles: '[H][H]', count: 1 });
  assert.equal((await pending).stats.forceField, 'gfn2');
  pending = wait(m => m.type === 'error');
  worker.postMessage({ type: 'add', smiles: '[H][H]', count: 40 });
  assert.match((await pending).text, /120/);
  pending = wait(m => m.type === 'frame');
  assert.equal((await pending).N, 42);
  console.log('Reaction observatory: real UHF exchange, immutable replay, spin integrals, energy/momentum of collisions, geometry, density differences and bounded recording passed.');
} finally { await worker.terminate(); }
