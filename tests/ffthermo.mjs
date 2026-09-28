// Verifica termochimica del campo di forze reattivo della sandbox: entalpie di atomizzazione a 298 K
// calcolate (−E del minimo, gli atomi isolati hanno energia zero) contro la legge di Hess sui dati sperimentali
// ΔatH = Σ ΔfH(atomi) − ΔfH(molecola). Fonti: CODATA Key Values for Thermodynamics (Cox, Wagman, Medvedev 1989)
// per atomi e molecole inorganiche; tabelle NBS riportate da Atkins, Physical Chemistry, tab. 2.5, per le organiche.

import { ReactiveFF } from '../js/chem/reactive.js';
import { parseSmiles } from '../js/chem/smiles.js';
import { embedMolecule } from '../js/chem/embed.js';
import { analyzeStructure } from '../js/chem/structure.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};

const ATOM_DFH = { 1: 217.998, 6: 716.68, 7: 472.68, 8: 249.18, 9: 79.38, 16: 277.17, 17: 121.301 };
const MOLECULES = [
  ['H₂', '[H][H]', 0], ['O₂', 'O=O', 0], ['N₂', 'N#N', 0], ['Cl₂', 'ClCl', 0], ['F₂', 'FF', 0],
  ['H₂O', 'O', -241.826], ['NH₃', 'N', -45.94], ['HCl', 'Cl', -92.31], ['HF', 'F', -273.30], ['H₂S', 'S', -20.6],
  ['CO₂', 'O=C=O', -393.51], ['CO', '[C-]#[O+]', -110.53],
  ['CH₄', 'C', -74.81], ['C₂H₆', 'CC', -84.68], ['C₃H₈', 'CCC', -103.85], ['C₂H₄', 'C=C', 52.26], ['C₂H₂', 'C#C', 226.73],
  ['C₆H₆', 'c1ccccc1', 82.93], ['CH₃OH', 'CO', -200.66], ['C₂H₅OH', 'CCO', -235.10], ['HCHO', 'C=O', -108.57],
  ['HCN', 'C#N', 135.1], ['CH₃NH₂', 'CN', -22.97],
];

function geometry(smiles) {
  const g = parseSmiles(smiles);
  const an = analyzeStructure({ atoms: g.atoms, bonds: g.bonds, charge: 0 });
  const at = embedMolecule(g, an.atoms.map(a => a.lonePairs + a.radical), { seeds: 2, iterations: 1500 });
  return { Z: at.map(a => a.Z), pos: Float64Array.from(at.flatMap(a => a.xyz.map(v => v * 0.52917721090))) };
}
function minimize(Z, pos) {
  const ff = new ReactiveFF(); const F = new Float64Array(pos.length), v = new Float64Array(pos.length);
  let dt = 0.05, alpha = 0.1, npos = 0, res;
  for (let s = 0; s < 6000; s++) {
    res = ff.compute(Z, pos, F);
    let P = 0, vn = 0, fn = 0;
    for (let k = 0; k < pos.length; k++) { P += F[k] * v[k]; vn += v[k] * v[k]; fn += F[k] * F[k]; }
    if (Math.sqrt(fn) < 1e-5) break;
    vn = Math.sqrt(vn); fn = Math.sqrt(fn);
    for (let k = 0; k < pos.length; k++) v[k] = (1 - alpha) * v[k] + alpha * F[k] / fn * vn;
    if (P > 0) { if (++npos > 5) { dt = Math.min(dt * 1.1, 0.3); alpha *= 0.99; } } else { npos = 0; dt *= 0.5; v.fill(0); alpha = 0.1; }
    for (let k = 0; k < pos.length; k++) { v[k] += dt * F[k]; pos[k] += dt * v[k]; }
  }
  return res;
}

const rows = [];
for (const [name, smi, dfh] of MOLECULES) {
  const g = geometry(smi);
  const r = minimize(g.Z, g.pos);
  const exp = g.Z.reduce((a, z) => a + ATOM_DFH[z], 0) - dfh;
  const calc = -r.E * 96.48533212;
  rows.push({ name, exp, calc, dev: (calc - exp) / exp * 100 });
}
console.log('molecola    ΔatH sper.  campo reattivo  scarto');
for (const x of rows) console.log(`${x.name.padEnd(10)} ${x.exp.toFixed(1).padStart(9)} ${x.calc.toFixed(1).padStart(12)} ${x.dev.toFixed(1).padStart(8)} %`);
const regular = rows.filter(x => !['CO', 'CO₂', 'N₂'].includes(x.name));
const mad = regular.reduce((a, x) => a + Math.abs(x.dev), 0) / regular.length;
check('scarto medio (legami ordinari)', mad < 4, `${mad.toFixed(1)} % su ${regular.length} molecole`);
check('nessuno scarto oltre il 7 % sui legami ordinari', regular.every(x => Math.abs(x.dev) < 7), regular.filter(x => Math.abs(x.dev) >= 7).map(x => x.name).join(', ') || 'tutti entro il 7 %');
// limiti noti, controllati perché non peggiorino: legami multipli con cariche formali o molto corti
const lim = (n) => rows.find(x => x.name === n).dev;
check('limiti noti: N₂, CO₂, CO', lim('N₂') > -15 && lim('CO₂') > -22 && lim('CO') > -45, `N₂ ${lim('N₂').toFixed(1)} %, CO₂ ${lim('CO₂').toFixed(1)} %, CO ${lim('CO').toFixed(1)} %`);

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nTermochimica del campo reattivo: tutte le verifiche superate');
