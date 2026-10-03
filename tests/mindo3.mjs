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
// ioni: calori di formazione contro PySCF (pyscf-semiempirical 0.1.1, MINDO/3 con mol.charge e mol.spin)
{
  const IONS = [
    ['NH₄⁺', [7, 1, 1, 1, 1], [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]], 1, 0, 162.4212838320509],
    ['OH⁻', [8, 1], [[0, 0, 0], [0, 0, 0.97]], -1, 0, -6.426754158389485],
    ['H₃O⁺', [8, 1, 1, 1], [[0, 0, 0.1], [0.95, 0, -0.2], [-0.47, 0.82, -0.2], [-0.47, -0.82, -0.2]], 1, 0, 146.8992323668399],
    ['NO₂⁺', [7, 8, 8], [[0, 0, 0], [0, 0, 1.15], [0, 0, -1.15]], 1, 0, 184.0780338150944],
    ['CH₃⁺', [6, 1, 1, 1], [[0, 0, 0], [1.09, 0, 0], [-0.545, 0.944, 0], [-0.545, -0.944, 0]], 1, 0, 260.3304317477141],
    ['H₂O⁺ (doppietto)', [8, 1, 1], [[0, 0, 0], [0, -0.757, 0.587], [0, 0.757, 0.587]], 1, 1, 231.2039852511698],
  ];
  let worst = 0, name = '';
  for (const [nm, Z, xyz, charge, spin, ref] of IONS) {
    const nel = Z.reduce((s, z) => s + ({ 1: 1, 6: 4, 7: 5, 8: 6 })[z], 0) - charge;
    const nb = (nel - spin) / 2;
    const r = new Mindo3({ Tel: 0, conv: 1e-10, charge, nalpha: nb + spin, nbeta: nb }).compute(Z, Float64Array.from(xyz.flat()), null);
    if (Math.abs(r.Hf - ref) > Math.abs(worst)) { worst = r.Hf - ref; name = nm; }
  }
  check('ioni: ΔfH come PySCF', Math.abs(worst) < 1e-3, `scarto massimo ${worst.toExponential(1)} kcal/mol (${name}), 6 ioni`);
}
// campo elettrico uniforme: −∂E/∂E = μ, forze esatte, forza netta qE su uno ione
{
  const Z = [8, 1, 1], pos = Float64Array.from(water.flat());
  const run = (field, p = pos, F = null, charge = 0, z = Z) => new Mindo3({ Tel: 0, conv: 1e-12, field, charge }).compute(z, p, F);
  const mu = run(null).dipole, h = 1e-4;
  let err = 0;
  for (let c = 0; c < 3; c++) { const a = [0, 0, 0], b = [0, 0, 0]; a[c] = h; b[c] = -h; err = Math.max(err, Math.abs(-(run(a).E - run(b).E) / (2 * h) - mu[c])); }
  check('campo: dipolo = −∂E/∂E', err < 1e-6, `μ(H₂O) = ${(Math.hypot(...mu) * 4.80320).toFixed(2)} D, scarto ${err.toExponential(1)} e·Å`);
  const fld = [0.3, -0.2, 0.5], F = new Float64Array(9);
  run(fld, pos, F);
  let ferr = 0;
  for (let k = 0; k < 9; k++) { const a = pos.slice(), b = pos.slice(); a[k] += h; b[k] -= h; ferr = Math.max(ferr, Math.abs(-(run(fld, a).E - run(fld, b).E) / (2 * h) - F[k])); }
  const F2 = new Float64Array(6); run(fld, Float64Array.from([0, 0, 0, 0, 0, 0.97]), F2, -1, [8, 1]);
  const net = [F2[0] + F2[3], F2[1] + F2[4], F2[2] + F2[5]];
  check('campo: forze e forza qE sullo ione', ferr < 1e-5 && net.every((v, c) => Math.abs(v + fld[c]) < 1e-9), `forze ${ferr.toExponential(1)} eV/Å; OH⁻: F = (${net.map(v => v.toFixed(3)).join('; ')}) eV/Å`);
}
// gruppi non accoppiati: nessun elettrone passa fra molecole lontane senza sovrapposizione degli orbitali
{
  const mk = (opts, formal) => { const m = new Mindo3({ Tel: 300, conv: 1e-9, ...opts }); m.formal = formal; return m; };
  const Zw = [8, 1, 1, 8, 1, 1], pw = Float64Array.from([0, 0, -7.5, 0, -0.757, -6.9, 0, 0.757, -6.9, 0, 0, 7.5, 0, -0.757, 8.1, 0, 0.757, 8.1]);
  const rw = mk({ field: [0, 0, 2] }, [0, 0, 0, 0, 0, 0]).compute(Zw, pw, null);
  const Zi = [6, 1, 1, 1, 7, 8, 8], pi = [-6, 0, 0, -6, 1.09, 0.3, -6, -0.6, 0.95, -6, -0.6, -0.95, 6, 0.1, 0, 6, 0, 1.15, 6, 0, -1.15];
  const formal = [-1, 0, 0, 0, 1, 0, 0];
  const ri = mk({}, formal).compute(Zi, Float64Array.from(pi), null);
  const F = new Float64Array(21); mk({ conv: 1e-11 }, formal).compute(Zi, Float64Array.from(pi), F);
  let ferr = 0;
  for (let k = 0; k < 21; k++) { const a = pi.slice(), b = pi.slice(); a[k] += 1e-4; b[k] -= 1e-4; ferr = Math.max(ferr, Math.abs(-(mk({ conv: 1e-11 }, formal).compute(Zi, Float64Array.from(a)).E - mk({ conv: 1e-11 }, formal).compute(Zi, Float64Array.from(b)).E) / 2e-4 - F[k])); }
  const qw = rw.q[0] + rw.q[1] + rw.q[2], qc = ri.q[0] + ri.q[1] + ri.q[2] + ri.q[3];
  check('molecole lontane: niente trasferimento di elettroni a distanza', Math.abs(qw) < 0.01 && Math.abs(qc + 1) < 0.01 && ferr < 1e-5,
    `H₂O a 15 Å in 2 V/Å: q = ${qw.toFixed(3)}; CH₃⁻ ··· NO₂⁺ a 12 Å: q(CH₃) = ${qc.toFixed(3)}; forze ${ferr.toExponential(1)}`);
}
check('elementi e coppie senza parametri rifiutati', !mindo3Supports([17, 8]).ok && !mindo3Supports([18]).ok && mindo3Supports([6, 1, 8, 7]).ok, 'Cl–O e Ar esclusi');

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nMINDO/3: tutte le verifiche superate');
