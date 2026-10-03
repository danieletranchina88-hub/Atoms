// GFN2-xTB in JavaScript contro tblite (tests/data/gfn2-reference.json, generato da tools/xtb/reference_gfn2.py):
// energia totale, cariche di Mulliken, dipolo e gradiente analitico per molecole neutre, ioni, un tripletto,
// metalli di transizione ed elementi pesanti fino al piombo; più un controllo del gradiente alle differenze finite.

import fs from 'fs';
import { GFN2xTB } from '../js/chem/xtb/gfn2.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};
const maxDiff = (a, b) => a.reduce((m, x, i) => Math.max(m, Math.abs(x - b[i])), 0);

const ref = JSON.parse(fs.readFileSync(new URL('./data/gfn2-reference.json', import.meta.url)));
for (const [name, m] of Object.entries(ref.molecules)) {
  const calc = new GFN2xTB();
  const t = performance.now();
  const r = calc.compute(m.Z, m.pos.flat(), { charge: m.charge, uhf: m.uhf });
  const ms = performance.now() - t;
  const dE = Math.abs(r.energy - m.energy), dq = maxDiff(m.charges, r.charges);
  const dg = maxDiff(m.gradient.flat(), r.gradient), dmu = maxDiff(m.dipole, r.dipole);
  check(`${name}: energia, cariche, dipolo e gradiente`, r.converged && dE < 1e-9 && dq < 1e-6 && dmu < 1e-6 && dg < 1e-8,
    `E = ${r.energy.toFixed(9)} Eh (tblite ${m.energy.toFixed(9)}), |ΔE| ${dE.toExponential(1)}, |Δq| ${dq.toExponential(1)}, |Δμ| ${dmu.toExponential(1)}, |Δg| ${dg.toExponential(1)} Eh/bohr, ${r.iterations} iterazioni, ${ms.toFixed(0)} ms`);
}

// gradiente alle differenze finite per un sistema con occupazioni frazionarie (Fe–CO, quintetto)
{
  const m = ref.molecules['Fe(CO)'];
  const calc = new GFN2xTB({ etol: 1e-12, ptol: 1e-10 });
  const pos = m.pos.flat();
  const g = calc.compute(m.Z, pos, { uhf: m.uhf }).gradient;
  let err = 0;
  const h = 1e-4;
  for (const k of [0, 4, 8]) {
    const p = pos.slice(); p[k] += h;
    const ep = calc.compute(m.Z, p, { uhf: m.uhf, gradient: false }).energy;
    p[k] -= 2 * h;
    const em = calc.compute(m.Z, p, { uhf: m.uhf, gradient: false }).energy;
    err = Math.max(err, Math.abs((ep - em) / (2 * h) - g[k]));
  }
  check('Fe–CO: gradiente analitico = differenze finite', err < 1e-7, `max scarto ${err.toExponential(1)} Eh/bohr`);
}

// invarianza per traslazione e rotazione
{
  const m = ref.molecules['HCl_H2O'];
  const calc = new GFN2xTB();
  const e0 = calc.compute(m.Z, m.pos.flat(), { gradient: false }).energy;
  const c = Math.cos(0.7), s = Math.sin(0.7);
  const rot = m.pos.flatMap(([x, y, z]) => [c * x - s * y + 3.1, s * x + c * y - 1.2, z + 0.4]);
  const e1 = calc.compute(m.Z, rot, { gradient: false }).energy;
  check('HCl·H₂O: energia invariante per rototraslazione', Math.abs(e1 - e0) < 1e-9, `scarto ${Math.abs(e1 - e0).toExponential(1)} Eh`);
}

if (failures) { console.log(`\n${failures} verifiche GFN2-xTB fallite`); process.exit(1); }
console.log('\nTutte le verifiche GFN2-xTB superate.');
