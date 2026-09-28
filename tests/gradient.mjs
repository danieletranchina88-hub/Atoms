// Gradiente analitico HF confrontato con PySCF 2.14 (stesse basi del Basis Set Exchange).
import { runHF, ANGSTROM_TO_BOHR as A } from '../js/chem/hf.js';
import { hfGradient } from '../js/chem/gradient.js';

const mol = (list) => list.map(([Z, x, y, z]) => ({ Z, xyz: [x * A, y * A, z * A] }));
const REF = {
  h2o: [[-0.0018849083102358555, 0.004328414438077903, 0.018757339682538454], [0.00013577412612964662, 0.007887948276082302, -0.007702352718102157], [0.0017491341841061603, -0.012216362714160045, -0.011054986964432745]],
  o2: [[-0.000959628875318197, 1.1809903337522812e-16, 0.023173118081173172], [0.0009596288753179749, -1.1814504605032085e-16, -0.023173118081173172]],
};
const cases = [
  ['h2o', mol([[8, 0, 0, 0.1173], [1, 0, 0.7572, -0.4692], [1, 0.1, -0.7572, -0.4692]]), { basis: '6-31G*' }],
  ['o2', mol([[8, 0, 0, 0.6037], [8, 0.05, 0, -0.6037]]), { basis: '6-31G', multiplicity: 3 }],
];
let ok = true;
for (const [name, atoms, opts] of cases) {
  const r = runHF(atoms, opts);
  const t0 = performance.now();
  const g = hfGradient(r);
  let err = 0;
  g.forEach((row, i) => row.forEach((v, k) => { err = Math.max(err, Math.abs(v - REF[name][i][k])); }));
  if (err > 1e-7) ok = false;
  console.log(name, 'errore massimo', err.toExponential(2), `${(performance.now() - t0).toFixed(0)} ms`);
}
console.log(ok ? 'OK' : 'FAIL');
process.exit(ok ? 0 : 1);
