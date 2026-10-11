// Pure numerical preparation, shared by the module worker and its fallback.
import { RadialFunction, term, sampleAtom, sampleGrid } from '../physics/wavefunction.js';
import { marchingCubes } from './marching.js';

export async function prepareOrbitalScene({ orbitals, points, quality, enclosed, render }, cancelled = () => false) {
  const radialCache = new Map();
  const list = orbitals.map(o => {
    if (!radialCache.has(o.subshell)) radialCache.set(o.subshell, new RadialFunction(o.grid, o.u, o.l));
    return { ...o, radial: radialCache.get(o.subshell) };
  });
  const cloud = render === 'surface' || !list.length ? null : sampleAtom(
    list.map((o, group) => ({ ...o, weight: o.occ, group })), points, 113,
  );
  const surfaces = [];
  // Keep memory and startup bounded on heavy atoms; each grid resolves its OWN orbital.
  const resolution = Math.max(32, Math.min(quality, Math.floor(Math.cbrt(6000000 / Math.max(list.length, 1)))));
  if (render !== 'cloud') {
    for (let i = 0; i < list.length; i++) {
      if (cancelled()) throw new Error('Preparazione sostituita da una nuova vista.');
      const o = list[i];
      const half = o.radial.radiusEnclosing(0.997) * 1.02;
      const field = sampleGrid([term(o.radial, o.l, o.m)], resolution, half, enclosed);
      for (const sign of [1, -1]) {
        const mesh = marchingCubes(field, field.iso, sign);
        if (mesh.positions.length) surfaces.push({ ...mesh, group: i, sign });
      }
      // Yield between grids so fallback rendering doesn't lock the UI for an entire atom.
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }
  return { cloud, surfaces, resolution };
}
