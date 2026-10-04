// GFN2-xTB in JavaScript contro tblite (tests/data/gfn2-reference.json, generato da tools/xtb/reference_gfn2.py):
// energia totale, cariche di Mulliken, dipolo e gradiente analitico per molecole neutre, ioni, un tripletto,
// metalli di transizione ed elementi pesanti fino al piombo; più un controllo del gradiente alle differenze finite.

import fs from 'fs';
import { GFN2xTB } from '../js/chem/xtb/gfn2.js';
import { makeGFN2Provider } from '../js/chem/xtb/provider.js';
import { computeGrid } from '../js/chem/densityWorker.js';

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

// campo elettrico uniforme: energia e dipolo contro tblite, gradiente contro le differenze finite
for (const [name, f] of Object.entries(ref.field)) {
  const m = ref.molecules[name];
  const calc = new GFN2xTB({ field: f.field, etol: 1e-12, ptol: 1e-10 });
  const pos = m.pos.flat(), opt = { charge: m.charge, uhf: m.uhf };
  const r = calc.compute(m.Z, pos, opt);
  let err = 0;
  for (const k of [1, 2, 4]) {
    const h = 1e-4, p = pos.slice();
    p[k] += h; const ep = calc.compute(m.Z, p, { ...opt, gradient: false }).energy;
    p[k] -= 2 * h; const em = calc.compute(m.Z, p, { ...opt, gradient: false }).energy;
    err = Math.max(err, Math.abs((ep - em) / (2 * h) - r.gradient[k]));
  }
  const dE = Math.abs(r.energy - f.energy), dmu = maxDiff(f.dipole, r.dipole);
  check(`${name} nel campo [${f.field}] au: energia, dipolo, gradiente`, dE < 1e-9 && dmu < 1e-6 && err < 1e-7,
    `|ΔE| ${dE.toExponential(1)} Eh, |Δμ| ${dmu.toExponential(1)} e·bohr, gradiente − differenze finite ${err.toExponential(1)} Eh/bohr`);
}

// solvente implicito ALPB (acqua): energia, cariche e gradiente contro tblite
for (const [name, a] of Object.entries(ref.alpb)) {
  const m = ref.molecules[name];
  const r = new GFN2xTB({ solvent: 'water' }).compute(m.Z, m.pos.flat(), { charge: m.charge, uhf: m.uhf });
  const dE = Math.abs(r.energy - a.energy), dq = maxDiff(a.charges, r.charges), dg = maxDiff(a.gradient.flat(), r.gradient);
  check(`${name} in acqua (ALPB): energia, cariche, gradiente`, r.converged && dE < 1e-9 && dq < 1e-6 && dg < 1e-8,
    `ΔG_solv = ${(r.parts.solvation * 627.509).toFixed(2)} kcal/mol, |ΔE| ${dE.toExponential(1)}, |Δq| ${dq.toExponential(1)}, |Δg| ${dg.toExponential(1)} Eh/bohr`);
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

// nuvola elettronica disegnata dalla sandbox: valenza GFN2 + core dell'atomo isolato = tutti gli elettroni
{
  const p = makeGFN2Provider();
  const pos = Float64Array.from([0, 0, 0, 0.757, 0.587, 0, -0.757, 0.587, 0].map((v, i) => v + [0.0371, 0.0613, 0.0517][i % 3]));
  p.compute([8, 1, 1], pos, new Float64Array(9), {});
  const w = p.wavefunction();
  const integ = (opts) => { const g = computeGrid({ Z: w.Z, pos: w.pos, box: 8, res: 121, P: w.P, ...opts }); const dv = (g.step / 0.52917721090) ** 3; return g.values.reduce((s, v) => s + v * v ** (opts.what === 'orbital' ? 1 : 0) * dv, 0); };
  const ne = integ({ mode: 'xtb', what: 'density' }), npro = integ({ mode: 'promolecular', what: 'density' });
  const norm = integ({ mode: 'xtb', what: 'orbital', orb: w.homo.c });
  check('H₂O: densità GFN2 + core e promolecolare integrano a 10 elettroni, HOMO normalizzato',
    Math.abs(ne - 10) < 0.2 && Math.abs(npro - 10) < 0.2 && Math.abs(norm - 1) < 1e-2, `∫ρ = ${ne.toFixed(3)} (GFN2), ${npro.toFixed(3)} (atomi isolati), ∫|ψ|² = ${norm.toFixed(4)}`);
}

if (failures) { console.log(`\n${failures} verifiche GFN2-xTB fallite`); process.exit(1); }
console.log('\nTutte le verifiche GFN2-xTB superate.');
