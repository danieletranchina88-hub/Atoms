// Worker della modalità "Come nasce un legame": esegue le scansioni e invia i punti man mano che sono pronti, così le
// curve crescono mentre gli atomi si avvicinano.

import { DIATOMICS } from './diatomicData.js';
import { makeContext, scanHF, scanGFN2, scanGrid, methodsFor } from './bondFormation.js';

self.onmessage = (ev) => {
  const m = ev.data;
  if (m.type !== 'run') return;
  const mol = DIATOMICS.find(d => d.id === m.id);
  try {
    const grid = scanGrid(mol, m.grid ?? {});
    const ctx = makeContext(mol, m.basis);
    const atomInfo = (at) => ({
      Z: at.Z, symbol: at.symbol, E: at.E, T: at.T, n: at.n,
      orbitals: at.orbitals.map(o => ({ label: o.label, l: o.l, axis: o.axis, e: o.e, occ: o.occ, extra: o.extra, core: o.core })),
    });
    self.postMessage({ type: 'context', token: m.token, grid, n: ctx.n, Eatoms: ctx.Eatoms, Tatoms: ctx.Tatoms, Ppro: ctx.Ppro, A: atomInfo(ctx.A), B: atomInfo(ctx.B) });
    const methods = m.methods ?? methodsFor(mol);
    for (const method of methods) {
      const t0 = Date.now();
      if (method === 'gfn2') {
        const g = scanGFN2(mol, { grid, onPoint: (p, i) => self.postMessage({ type: 'point', token: m.token, method, i, point: p }) });
        self.postMessage({ type: 'done', token: m.token, method, Eatoms: g.Eatoms, ms: Date.now() - t0 });
      } else {
        scanHF(ctx, method, { grid, onPoint: (p, i) => self.postMessage({ type: 'point', token: m.token, method, i, point: p }) });
        self.postMessage({ type: 'done', token: m.token, method, Eatoms: ctx.Eatoms, ms: Date.now() - t0 });
      }
    }
    self.postMessage({ type: 'finished', token: m.token });
  } catch (err) {
    self.postMessage({ type: 'error', token: m.token, message: err.message });
  }
};
