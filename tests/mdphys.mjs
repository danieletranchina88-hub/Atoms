// Verifiche degli strumenti di chimica fisica della dinamica molecolare: g(r), capacità termica, barostato.

import { Simulation } from '../js/chem/md.js';
import { cubePairDensity, RDF, MSD, HeatCapacity } from '../js/chem/mdAnalysis.js';
import { KB_EV } from '../js/chem/reactiveData.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};

// distanza fra punti uniformi in un cubo: normalizzata, ⟨r²⟩ = L²/2, ⟨r⟩ = 0,6617 L (costante di Robbins)
{
  const L = 20, n = 3000, rmax = Math.sqrt(3) * L;
  let s = 0, m1 = 0, m2 = 0;
  for (let i = 0; i < n; i++) { const r = (i + 0.5) / n * rmax; const p = cubePairDensity(r, L) * rmax / n; s += p; m1 += p * r; m2 += p * r * r; }
  check('densità delle distanze nel cubo', Math.abs(s - 1) < 1e-4 && Math.abs(m2 / (L * L) - 0.5) < 1e-4 && Math.abs(m1 / L - 0.66170718) < 1e-4, `∫p = ${s.toFixed(5)}, ⟨r⟩/L = ${(m1 / L).toFixed(5)}, ⟨r²⟩/L² = ${(m2 / L / L).toFixed(5)}`);
}
// punti casuali (gas ideale): g(r) = 1
{
  let seed = 7;
  const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const L = 20, N = 300;
  const rdf = new RDF(10, 20);
  rdf.reset([18, 18]);
  for (let k = 0; k < 40; k++) rdf.accumulate({ N, Z: new Array(N).fill(18), pos: Array.from({ length: 3 * N }, () => (rnd() - 0.5) * L), box: L, fragOf: null });
  const g = rdf.result().g;
  const dev = Math.max(...g.slice(4).map(v => Math.abs(v - 1)));
  check('g(r) del gas ideale = 1', dev < 0.06, `scarto massimo ${dev.toFixed(3)}`);
}
// argon gassoso a 300 K: C_V = 3/2 N k_B (gas monoatomico), g(r) → 1 a grande distanza
{
  const sim = new Simulation({ box: 42, T: 300, dt: 2, seed: 11 });
  sim.addMolecule({ Z: [18], pos: new Float64Array(3) }, 120, 300);
  const cv = new HeatCapacity(), msd = new MSD();
  for (let i = 0; i < 4000; i++) sim.step();
  msd.reset(sim);
  for (let i = 0; i < 40000; i++) { sim.step(); if (sim.stepCount % 50 === 0) { cv.sample(sim, 'k'); msd.sample(sim); } }
  const c = cv.result(300);
  const r = c.Cv / (sim.N * KB_EV);
  check('C_V dell\'argon gassoso dalle fluttuazioni', c.stationary && Math.abs(r - 1.5) < 3 * c.err / (sim.N * KB_EV) + 0.05, `${r.toFixed(3)} ± ${(c.err / (sim.N * KB_EV)).toFixed(3)} N k_B (atteso 1,5)`);
  const m = msd.result(sim.box);
  check('gas rarefatto: MSD saturato dalle pareti, nessun D', m.D === null && Math.abs(m.points.at(-1)[1] / m.saturation - 1) < 0.2, `MSD finale ${m.points.at(-1)[1].toFixed(0)} Å² (L²/2 = ${m.saturation.toFixed(0)})`);
}
// barostato: il volume si porta al valore del gas ideale V = N k_B T / P₀
{
  const sim = new Simulation({ box: 42, T: 300, dt: 2, seed: 5 });
  sim.addMolecule({ Z: [18], pos: new Float64Array(3) }, 120, 300);
  sim.barostat = { on: true, P0: 20, tau: 10000 };
  for (let i = 0; i < 100000; i++) sim.step();
  const Lid = Math.cbrt(120 * 1.380649e-23 * 300 / 20e5) * 1e10;
  check('barostato: volume del gas ideale a 20 bar', Math.abs(sim.box / Lid - 1) < 0.03, `L = ${sim.box.toFixed(2)} Å (gas ideale ${Lid.toFixed(2)} Å)`);
}

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nChimica fisica della dinamica molecolare: tutte le verifiche superate');
