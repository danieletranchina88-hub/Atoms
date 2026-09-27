// Proprietà dalla funzione d'onda, confrontate con PySCF 2.14 (stesse basi del Basis Set Exchange).
import { runHF, ANGSTROM_TO_BOHR as A } from '../js/chem/hf.js';
import { populationAnalysis, dipoleMoment, mp2Energy, electrostaticPotential } from '../js/chem/properties.js';

let ok = true;
const check = (label, got, want, tol) => {
  const pass = Math.abs(got - want) <= tol;
  if (!pass) ok = false;
  console.log(`${pass ? 'ok ' : 'NO '} ${label}: ${got.toFixed(9)} (PySCF ${want.toFixed(9)})`);
};
const mol = (list) => list.map(([Z, x, y, z]) => ({ Z, xyz: [x * A, y * A, z * A] }));
const water = runHF(mol([[8, 0, 0, 0.1173], [1, 0, 0.7572, -0.4692], [1, 0, -0.7572, -0.4692]]), { basis: '6-31G*' });
const pop = populationAnalysis(water);
check('carica di Mulliken O (H2O, 6-31G*)', pop.mulliken[0], -0.8662277513906318, 1e-6);
check('ordine di legame di Mayer O–H', pop.bondOrder[0][1], 0.7863824616802716, 1e-6);
check('momento di dipolo (D)', dipoleMoment(water).debye, 2.2260254922783984, 1e-6);
check('correlazione MP2 (Ha)', mp2Energy(water).energy, -0.1861195365804197, 1e-7);
const esp = electrostaticPotential(water, new Float64Array([0, 0, 3.0, 1.5, 1.0, -2.0, 0.2, 2.5, 0.4]));
[-0.0866131343, 0.0642049927, 0.0630497961].forEach((v, i) => check(`potenziale elettrostatico, punto ${i + 1}`, esp[i], v, 1e-8));
const o2 = runHF(mol([[8, 0, 0, 0.6037], [8, 0, 0, -0.6037]]), { basis: '6-31G', multiplicity: 3 });
check('correlazione UMP2 O2 tripletto (Ha)', mp2Energy(o2).energy, -0.2376320667350552, 1e-7);
console.log(ok ? 'OK' : 'FAIL');
process.exit(ok ? 0 : 1);
