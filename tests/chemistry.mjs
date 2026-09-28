// Verifiche della chimica "classica": SMILES, strutture di Lewis, numeri di ossidazione, VSEPR, simmetria.
import { parseSmiles, writeSmiles, hillFormula } from '../js/chem/smiles.js';
import { analyzeStructure } from '../js/chem/structure.js';
import { embedMolecule } from '../js/chem/embed.js';
import { analyzeSymmetry } from '../js/chem/symmetry.js';
import { ISOTOPE_MASS } from '../js/chem/vibrations.js';

let ok = true;
const check = (label, got, want) => {
  const pass = JSON.stringify(got) === JSON.stringify(want);
  if (!pass) ok = false;
  console.log(`${pass ? 'ok ' : 'NO '} ${label}: ${JSON.stringify(got)}${pass ? '' : ` (atteso ${JSON.stringify(want)})`}`);
};
const analyze = (s) => { const g = parseSmiles(s); return { g, a: analyzeStructure(g) }; };
const heavy = (r, key) => r.a.atoms.filter(x => x.Z > 1).map(x => x[key]);

check('formula acido acetico', hillFormula(parseSmiles('CC(=O)O').atoms), 'C₂H₄O₂');
check('formula ammoniaca', hillFormula(parseSmiles('N').atoms), 'NH₃');
check('SMILES scritto e riletto (benzene)', hillFormula(parseSmiles(writeSmiles(parseSmiles('c1ccccc1'))).atoms), 'C₆H₆');
check('risonanza benzene', analyze('c1ccccc1').a.resonanceCount, 2);
check('risonanza nitrato', analyze('[O-][N+](=O)[O-]').a.resonanceCount, 3);
check('risonanza carbonato', analyze('C(=O)([O-])[O-]').a.resonanceCount, 3);
check('ossidazione solfato', heavy(analyze('[O-]S(=O)(=O)[O-]'), 'oxidation'), [-2, 6, -2, -2, -2]);
check('ossidazione perossido', heavy(analyze('OO'), 'oxidation'), [-1, -1]);
check('ossidazione acido acetico', heavy(analyze('CC(=O)O'), 'oxidation'), [-3, 3, -2, -2]);
check('cariche formali CO', heavy(analyze('[C-]#[O+]'), 'formalCharge'), [-1, 1]);
check('VSEPR H2O', analyze('O').a.atoms[0].vsepr.molecularGeometry, 'angolare (piegata)');
check('VSEPR NH3', analyze('N').a.atoms[0].vsepr.molecularGeometry, 'piramidale trigonale');
check('VSEPR SF4', analyze('FS(F)(F)F').a.atoms[1].vsepr.molecularGeometry, 'a cavalletto (altalena)');
check('VSEPR ClF3', analyze('FCl(F)F').a.atoms[1].vsepr.molecularGeometry, 'a T');
check('VSEPR XeF4', analyze('F[Xe](F)(F)F').a.atoms[1].vsepr.molecularGeometry, 'quadrato planare');
check('ibridazione C etilene', analyze('C=C').a.atoms[0].hybridization, 'sp²');
check('ibridazione C acetilene', analyze('C#C').a.atoms[0].hybridization, 'sp');
check('legami σ/π acetonitrile', [analyze('CC#N').a.sigmaBonds, analyze('CC#N').a.piBonds], [5, 2]);

// geometrie ideali
const benzene = [];
for (let k = 0; k < 6; k++) {
  const t = k * Math.PI / 3;
  benzene.push({ Z: 6, xyz: [2.64 * Math.cos(t), 2.64 * Math.sin(t), 0] }, { Z: 1, xyz: [4.69 * Math.cos(t), 4.69 * Math.sin(t), 0] });
}
const sf6 = [{ Z: 16, xyz: [0, 0, 0] }, ...[[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]].map(v => ({ Z: 9, xyz: v.map(c => c * 2.95) }))];
for (const [name, geo, pg, sigma] of [['benzene', benzene, 'D6h', 12], ['SF6', sf6, 'Oh', 24]]) {
  const sym = analyzeSymmetry(geo, geo.map(x => ISOTOPE_MASS[x.Z]));
  check(`gruppo puntuale ${name} (geometria ideale)`, [sym.pointGroup, sym.sigma], [pg, sigma]);
}
for (const [s, pg, sigma] of [['O', 'C2v', 2], ['N', 'C3v', 3], ['C', 'Td', 12], ['C=C', 'D2h', 4], ['O=C=O', 'D∞h', 2]]) {
  const g = parseSmiles(s);
  const a = analyzeStructure(g);
  const geo = embedMolecule(g, a.atoms.map(x => x.lonePairs + x.radical), { seeds: 4, iterations: 4000 });
  const sym = analyzeSymmetry(geo, geo.map(x => ISOTOPE_MASS[x.Z]));
  check(`gruppo puntuale ${s}`, [sym.pointGroup, sym.sigma], [pg, sigma]);
}
console.log(ok ? 'OK' : 'FAIL');
process.exit(ok ? 0 : 1);
