// Confronto Hartree–Fock con valori di riferimento calcolati con PySCF 2.14
// (stessi dati di base del Basis Set Exchange, funzioni d cartesiane).
import { runHF, ANGSTROM_TO_BOHR } from '../js/chem/hf.js';

const A = ANGSTROM_TO_BOHR;
const mol = (list) => list.map(([Z, x, y, z]) => ({ Z, xyz: [x * A, y * A, z * A] }));
const water = mol([[8, 0, 0, 0.1173], [1, 0, 0.7572, -0.4692], [1, 0, -0.7572, -0.4692]]);
const o2 = mol([[8, 0, 0, 0.6037], [8, 0, 0, -0.6037]]);
const nh3 = mol([[7, 0, 0, 0.1127], [1, 0, 0.9377, -0.2630], [1, 0.8121, -0.4689, -0.2630], [1, -0.8121, -0.4689, -0.2630]]);
const REF = {
  'H2O STO-3G': -74.96302316286221,
  'H2O 6-31G*': -76.01050499531273,
  'O2 triplet UHF 6-31G': -149.54557863160665,
  'NH3 6-31G**': -56.19532257494224,
};
const cases = [
  ['H2O STO-3G', water, { basis: 'STO-3G' }],
  ['H2O 6-31G*', water, { basis: '6-31G*' }],
  ['O2 triplet UHF 6-31G', o2, { basis: '6-31G', multiplicity: 3 }],
  ['NH3 6-31G**', nh3, { basis: '6-31G**' }],
];
let ok = true;
for (const [name, atoms, opts] of cases) {
  const t0 = performance.now();
  const r = runHF(atoms, opts);
  const ref = REF[name];
  const diff = ref !== undefined ? r.energy - ref : NaN;
  if (!(Math.abs(diff) < 1e-9)) ok = false;
  console.log(`${name.padEnd(24)} nbf=${String(r.n).padStart(3)} it=${r.iterations} E=${r.energy.toFixed(10)} ref=${ref} Δ=${diff.toExponential(2)} ${r.S2 !== null ? `<S²>=${r.S2.toFixed(4)}` : ''} ${(performance.now() - t0).toFixed(0)} ms`);
}
console.log(ok ? 'OK' : 'FAIL');
process.exit(ok ? 0 : 1);
