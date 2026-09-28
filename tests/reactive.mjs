// Verifiche del campo di forze reattivo della sandbox:
//  1. forze = −gradiente dell'energia (differenze finite) su configurazioni casuali;
//  2. geometrie di equilibrio (lunghezze di legame, angoli, ordini di legame);
//  3. energie di reazione (legge di Hess con le energie medie di legame).

import { ReactiveFF } from '../js/chem/reactive.js';
import { parseSmiles } from '../js/chem/smiles.js';
import { embedMolecule } from '../js/chem/embed.js';
import { analyzeStructure } from '../js/chem/structure.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};

function random(seed) {
  let s = seed;
  return () => { s = (s * 1103515245 + 12345) % 2147483648; return s / 2147483648; };
}

// 1. gradiente
const R = random(11);
for (const Z of [[8, 1, 1, 1], [6, 1, 1, 1, 8, 8], [6, 6, 6, 1, 11, 17, 17], [7, 7, 18, 1, 1]]) {
  const pos = Float64Array.from({ length: 3 * Z.length }, () => (R() - 0.5) * 2.6);
  const F = new Float64Array(pos.length);
  new ReactiveFF().compute(Z, pos, F);
  let err = 0, fmax = 0;
  const h = 1e-5, tmp = new Float64Array(pos.length);
  for (let k = 0; k < pos.length; k++) {
    const p1 = pos.slice(), p2 = pos.slice();
    p1[k] += h; p2[k] -= h;
    const fd = -(new ReactiveFF().compute(Z, p1, tmp).E - new ReactiveFF().compute(Z, p2, tmp).E) / (2 * h);
    err = Math.max(err, Math.abs(fd - F[k]) / Math.max(1, Math.abs(F[k])));
    fmax = Math.max(fmax, Math.abs(F[k]));
  }
  check(`forze analitiche ${Z.join(',')}`, err < 1e-5, `errore relativo ${err.toExponential(1)} (|F| max ${fmax.toFixed(1)} eV/Å)`);
}

// 2. geometrie di equilibrio
function geometry(smiles) {
  const g = parseSmiles(smiles);
  const an = analyzeStructure({ atoms: g.atoms, bonds: g.bonds, charge: 0 });
  const at = embedMolecule(g, an.atoms.map(a => a.lonePairs + a.radical), { seeds: 2, iterations: 1500 });
  return { Z: at.map(a => a.Z), pos: Float64Array.from(at.flatMap(a => a.xyz.map(v => v * 0.52917721090))) };
}
function minimize(Z, pos) {
  const ff = new ReactiveFF();
  const F = new Float64Array(pos.length), v = new Float64Array(pos.length);
  let dt = 0.05, alpha = 0.1, npos = 0, res;
  for (let s = 0; s < 5000; s++) {
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
const dist = (p, i, j) => Math.hypot(p[3 * i] - p[3 * j], p[3 * i + 1] - p[3 * j + 1], p[3 * i + 2] - p[3 * j + 2]);
const angle = (p, i, j, k) => {
  const a = [0, 1, 2].map(c => p[3 * i + c] - p[3 * j + c]), b = [0, 1, 2].map(c => p[3 * k + c] - p[3 * j + c]);
  return Math.acos((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / Math.hypot(...a) / Math.hypot(...b)) * 180 / Math.PI;
};
const E = {};
const mols = {};
for (const smi of ['[H][H]', 'O=O', 'N#N', 'O', 'C', 'O=C=O', 'N', 'C=C', 'C#C', 'CC', 'c1ccccc1', 'Cl', 'ClCl']) {
  const g = geometry(smi);
  const r = minimize(g.Z, g.pos);
  E[smi] = r.E;
  mols[smi] = { ...g, res: r };
}
const bondOrder = (smi, i, j) => mols[smi].res.bonds.find(b => (b.i === i && b.j === j) || (b.i === j && b.j === i))?.n ?? 0;
const cases = [
  ['H–H', dist(mols['[H][H]'].pos, 0, 1), 0.741, 0.03],
  ['O=O', dist(mols['O=O'].pos, 0, 1), 1.21, 0.04],
  ['N≡N', dist(mols['N#N'].pos, 0, 1), 1.10, 0.04],
  ['O–H (acqua)', dist(mols.O.pos, 0, 1), 0.96, 0.03],
  ['H–O–H', angle(mols.O.pos, 1, 0, 2), 104.5, 5],
  ['H–C–H (metano)', angle(mols.C.pos, 1, 0, 2), 109.47, 1],
  ['O=C=O', angle(mols['O=C=O'].pos, 0, 1, 2), 180, 3],
  ['C–C (etano)', dist(mols.CC.pos, 0, 1), 1.54, 0.05],
  ['C=C (etilene)', dist(mols['C=C'].pos, 0, 1), 1.34, 0.04],
  ['C≡C (acetilene)', dist(mols['C#C'].pos, 0, 1), 1.20, 0.04],
];
for (const [name, v, ref, tol] of cases) check(name, Math.abs(v - ref) <= tol, `${v.toFixed(3)} (rif. ${ref})`);
check('ordini di legame 1, 2, 3', Math.abs(bondOrder('CC', 0, 1) - 1) < 0.15 && Math.abs(bondOrder('C=C', 0, 1) - 2) < 0.15 && Math.abs(bondOrder('C#C', 0, 1) - 3) < 0.15,
  `${bondOrder('CC', 0, 1).toFixed(2)}, ${bondOrder('C=C', 0, 1).toFixed(2)}, ${bondOrder('C#C', 0, 1).toFixed(2)}`);
check('benzene aromatico (n ≈ 1,5)', Math.abs(bondOrder('c1ccccc1', 0, 1) - 1.5) < 0.1, bondOrder('c1ccccc1', 0, 1).toFixed(2));
const qO = mols.O.res.q[0];
check('carica parziale di O nell\'acqua negativa', qO < -0.3 && qO > -1, qO.toFixed(2));

// legame a idrogeno: dimero d'acqua
{
  const w = mols.O;
  const E1 = E.O;
  const u = [0, 1, 2].map(c => w.pos[3 + c] - w.pos[c]);
  const L = Math.hypot(...u);
  const pos = new Float64Array(18);
  pos.set(w.pos);
  for (let k = 0; k < 3; k++) for (let c = 0; c < 3; c++) pos[9 + 3 * k + c] = w.pos[3 * k + c] + 3.0 * u[c] / L;
  const r = minimize([...w.Z, ...w.Z], pos);
  const bind = (r.E - 2 * E1) * 96.485;
  const oo = dist(pos, 0, 3);
  check('dimero d\'acqua (legame a idrogeno)', bind < -10 && bind > -30 && Math.abs(oo - 2.9) < 0.15, `${bind.toFixed(1)} kJ/mol, O···O ${oo.toFixed(2)} Å (sper. ≈ −21 kJ/mol, 2,91 Å)`);
}

// 3. energie di reazione (kJ/mol)
const KJ = 96.485;
const rx = [
  ['2 H₂ + O₂ → 2 H₂O', 2 * E.O - 2 * E['[H][H]'] - E['O=O'], -484, 60],
  ['H₂ + Cl₂ → 2 HCl', 2 * E.Cl - E['[H][H]'] - E.ClCl, -184, 30],
  ['C₂H₄ + H₂ → C₂H₆', E.CC - E['C=C'] - E['[H][H]'], -137, 40],
];
for (const [name, dE, ref, tol] of rx) check(name, Math.abs(dE * KJ - ref) <= tol, `${(dE * KJ).toFixed(0)} kJ/mol (sper. ΔH ${ref})`);

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nCampo di forze reattivo: tutte le verifiche superate');
