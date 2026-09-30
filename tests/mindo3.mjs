// Verifiche del metodo quantistico MINDO/3 della sandbox contro i valori di riferimento di PySCF
// (pyscf-semiempirical 0.1.1: tests/test_mindo3.py e i controlli in rmindo3_grad.py / umindo3_grad.py).

import { Mindo3, mindo3Supports } from '../js/chem/mindo3.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};
const run = (Z, xyz, opts) => {
  const m = new Mindo3({ Tel: 0, conv: 1e-10, ...opts });
  const F = new Float64Array(3 * Z.length);
  const r = m.compute(Z, Float64Array.from(xyz.flat()), F);
  return { r, F, m };
};

// calori di formazione (kcal/mol) in geometrie fisse
{
  const { r } = run([8, 1, 1], [[0, 0, 0], [1, 0, 0], [0, 1, 0]], { nalpha: 4, nbeta: 4 });
  check('H₂O (RHF): ΔfH', Math.abs(r.Hf - -48.82621264564841) < 1e-4, `${r.Hf.toFixed(6)} kcal/mol (PySCF −48,826213)`);
  const c = run([6, 1, 1, 1, 1], [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]], { nalpha: 4, nbeta: 4 }).r;
  check('CH₄ deformato (RHF): ΔfH', Math.abs(c.Hf - 75.76019731515225) < 1e-4, `${c.Hf.toFixed(6)} kcal/mol (PySCF 75,760197)`);
  const oh = run([8, 1], [[0, 0, 0], [1, 0, 0]], { nalpha: 4, nbeta: 3 }).r;
  check('OH radicale (UHF): ΔfH', Math.abs(oh.Hf - 18.08247965492137) < 1e-4, `${oh.Hf.toFixed(6)} kcal/mol (PySCF 18,082480)`);
}
// energia totale dell'acqua all'equilibrio, singoletto e tripletto
const water = [[0, 0, 0], [0, -0.757, 0.587], [0, 0.757, 0.587]];
{
  const s = run([8, 1, 1], water, { nalpha: 4, nbeta: 4 }).r;
  check('H₂O: energia totale (RHF)', Math.abs(s.E - -341.50046431149383) < 1e-5, `${s.E.toFixed(6)} eV (PySCF −341,500464)`);
  const t = run([8, 1, 1], water, { nalpha: 5, nbeta: 3 }).r;
  check('H₂O tripletto: energia totale (UHF)', Math.abs(t.E - -336.25080977434175) < 1e-5, `${t.E.toFixed(6)} eV (PySCF −336,250810)`);
}
// gradiente analitico contro differenze finite
for (const [name, Z, xyz, opts] of [
  ['H₂O', [8, 1, 1], [[0.02, -0.01, 0.03], [0.1, -0.75, 0.6], [-0.05, 0.77, 0.55]], { nalpha: 4, nbeta: 4 }],
  ['HCN + H', [6, 7, 1, 1], [[0, 0, 0], [1.18, 0.05, 0], [-1.1, 0.1, 0.2], [0.3, 1.9, -0.4]], {}],
  ['O₂ + CH₄ (spin libero)', [8, 8, 6, 1, 1, 1, 1], [[0, 0, 0], [1.25, 0.1, 0], [2.9, 1.2, 0.3], [3.6, 0.5, 0.6], [3.3, 2.1, 0.8], [2.3, 1.5, -0.6], [2.2, 0.9, 1.1]], { Tel: 1000 }],
]) {
  const pos = xyz.flat();
  const { F, m } = run(Z, xyz, opts);
  let err = 0, fmax = 0;
  const h = 1e-4;
  for (let k = 0; k < pos.length; k++) {
    const p1 = Float64Array.from(pos), p2 = Float64Array.from(pos);
    p1[k] += h; p2[k] -= h;
    const e1 = m.compute(Z, p1, null).E, e2 = m.compute(Z, p2, null).E;
    const fd = -(e1 - e2) / (2 * h);
    err = Math.max(err, Math.abs(fd - F[k]));
    fmax = Math.max(fmax, Math.abs(F[k]));
  }
  check(`forze analitiche ${name}`, err < 2e-4 * Math.max(1, fmax), `scarto massimo ${err.toExponential(1)} eV/Å (|F| max ${fmax.toFixed(2)})`);
}
// stato di spin libero: O₂ tripletto, radicale metile doppietto
{
  const o2 = new Mindo3({ Tel: 300 }).compute([8, 8], Float64Array.from([0, 0, 0, 1.21, 0, 0]));
  check('O₂: lo stato fondamentale è un tripletto', Math.abs(Math.abs(o2.Sz) - 1) < 0.05, `S_z = ${o2.Sz.toFixed(3)}`);
  const ch3 = new Mindo3({ Tel: 300 }).compute([6, 1, 1, 1], Float64Array.from([0, 0, 0, 1.08, 0, 0, -0.54, 0.935, 0, -0.54, -0.935, 0]));
  check('CH₃: radicale (doppietto)', Math.abs(Math.abs(ch3.Sz) - 0.5) < 0.05, `S_z = ${ch3.Sz.toFixed(3)}`);
}
// orbitali e densità per il disegno: l'HOMO (base ortogonalizzata alla Löwdin) è normalizzato sulla griglia
{
  const { makeMindo3Provider } = await import('../js/chem/mindo3.js');
  const { computeGrid } = await import('../js/chem/densityWorker.js');
  const Z = [8, 1, 1], pos = Float64Array.from(water.flat());
  const p = makeMindo3Provider();
  p.compute(Z, pos, new Float64Array(9));
  const w = p.wavefunction();
  const g = computeGrid({ mode: 'sto', what: 'orbital', Z, pos, box: 6, res: 81, first: w.first, orb: w.homo.c });
  let s2 = 0;
  for (const v of g.values) s2 += v * v;
  s2 *= (g.step / 0.52917721090) ** 3;
  check('HOMO dell\'acqua normalizzato sulla griglia', Math.abs(s2 - 1) < 0.01, `∫|ψ|² = ${s2.toFixed(4)}, ε = ${w.homo.e.toFixed(2)} eV`);
}
// reazione vera nella dinamica quantistica: due radicali metile formano il legame C–C dell'etano
{
  const { Simulation } = await import('../js/chem/md.js');
  const { makeMindo3Provider } = await import('../js/chem/mindo3.js');
  const sim = new Simulation({ box: 12, T: 300, dt: 0.4, seed: 1 });
  const ch3 = (x, s) => [[6, x, 0, 0], [1, x + s * 0.35, 1.03, 0], [1, x + s * 0.35, -0.51, 0.89], [1, x + s * 0.35, -0.51, -0.89]];
  const at = [...ch3(-1.2, -1), ...ch3(1.2, 1)];
  sim.addAtoms(at.map(([Z, x, y, z]) => ({ Z, pos: [x, y, z] })));
  // spin opposti sui due carboni (coppia di singoletto): solo così si forma il legame
  sim.provider = makeMindo3Provider({ spinGuess: [0.5, 0, 0, 0, -0.5, 0, 0, 0] });
  sim.res = null;
  for (let i = 0; i < 600; i++) sim.step();
  const cc = sim.res.bonds.find(b => (b.i === 0 && b.j === 4) || (b.i === 4 && b.j === 0));
  const d = Math.hypot(sim.pos[0] - sim.pos[12], sim.pos[1] - sim.pos[13], sim.pos[2] - sim.pos[14]);
  // con spin paralleli (tripletto) i due radicali si respingono
  const sim2 = new Simulation({ box: 12, T: 300, dt: 0.4, seed: 1 });
  sim2.addAtoms(at.map(([Z, x, y, z]) => ({ Z, pos: [x, y, z] })));
  sim2.provider = makeMindo3Provider({ spinGuess: [0.5, 0, 0, 0, 0.5, 0, 0, 0] });
  sim2.res = null;
  for (let i = 0; i < 300; i++) sim2.step();
  const d2 = Math.hypot(sim2.pos[0] - sim2.pos[12], sim2.pos[1] - sim2.pos[13], sim2.pos[2] - sim2.pos[14]);
  check('radicali con spin paralleli non si legano (tripletto)', d2 > 3.2 && Math.abs(sim2.provider.info.Sz - 1) < 0.05, `C···C ${d2.toFixed(2)} Å, S_z = ${sim2.provider.info.Sz.toFixed(2)}`);
  check('2 CH₃ → C₂H₆ nella dinamica quantistica', cc && cc.n > 0.85 && d < 1.7, `ordine di legame C–C ${cc ? cc.n.toFixed(2) : 0}, distanza ${d.toFixed(2)} Å dopo ${(sim.time / 1000).toFixed(2)} ps, S_z = ${sim.provider.info.Sz.toFixed(2)}`);
}
check('elementi e coppie senza parametri rifiutati', !mindo3Supports([17, 8]).ok && !mindo3Supports([18]).ok && mindo3Supports([6, 1, 8, 7]).ok, 'Cl–O e Ar esclusi');

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nMINDO/3: tutte le verifiche superate');
