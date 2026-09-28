// Verifiche del solutore di equilibri in soluzione acquosa (becher) con casi da manuale.

import { Beaker, equilibrate } from '../js/chem/aqueous.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};
const run = (steps) => { const b = new Beaker(); let e = null; for (const [r, v] of steps) e = b.add(r, v); return { b, s: b.summary(), e }; };

// acido debole: pH = ½(pKa − log C) ≈ 2,88 per CH₃COOH 0,1 M
{
  const { s } = run([['HAc', 10], ['H2O', 90]]);
  check('pH di CH₃COOH 0,1 M', Math.abs(s.pH - 2.88) < 0.03, s.pH.toFixed(3));
}
// tampone equimolare: pH ≈ pKa (correzione di attività di Davies inclusa)
{
  const { s } = run([['HAc', 10], ['NaAc', 10], ['H2O', 80]]);
  check('tampone acetico equimolare', Math.abs(s.pH - 4.65) < 0.08, `pH ${s.pH.toFixed(3)} (pKa 4,76 − correzione di attività)`);
}
// base debole: NH₃ 0,1 M → pH ≈ 11,13
{
  const { s } = run([['NH3', 10], ['H2O', 90]]);
  check('pH di NH₃ 0,1 M', Math.abs(s.pH - 11.13) < 0.04, s.pH.toFixed(3));
}
// acido forte diluito e acqua pura
{
  const { s } = run([['HCl', 1], ['H2O', 99]]);
  check('pH di HCl 0,01 M', Math.abs(s.pH - 2.04) < 0.03, s.pH.toFixed(3));
  const w = equilibrate({}, 0.1);
  check('acqua pura: pH 7', Math.abs(w.pH - 7) < 1e-6, w.pH.toFixed(6));
}
// neutralizzazione: ΔH = −55,8 kJ/mol, ΔT = q/(m c_p)
{
  const { s } = run([['HCl', 25], ['NaOH', 25]]);
  check('pH al punto equivalente HCl/NaOH', Math.abs(s.pH - 7) < 0.02, s.pH.toFixed(3));
  const dT = s.T - 298.15;
  check('riscaldamento per neutralizzazione', Math.abs(s.heat - 0.025 * 55.8) < 0.02 && Math.abs(dT - 6.67) < 0.1, `q = ${s.heat.toFixed(3)} kJ, ΔT = ${dT.toFixed(2)} K (atteso 1,395 kJ, 6,67 K)`);
}
// prodotto di solubilità: AgCl in acqua pura, s = √Ksp (con γ ≈ 1)
{
  const eq = equilibrate({ Ag: 1e-3, Cl: 1e-3, Na: 0, NO3: 1e-3, K: 1e-3 }, 1);
  const ag = eq.species['Ag+'];
  const solid = eq.solids['AgCl(s)'] ?? 0;
  check('AgCl precipita da Ag⁺ 1 mM + Cl⁻ 1 mM', solid > 0.98e-3, `solido ${(solid * 1000).toFixed(4)} mmol, [Ag⁺] = ${ag.toExponential(2)} M (√Ksp = 1,35·10⁻⁵)`);
  check('solubilità di AgCl ≈ √Ksp', Math.abs(Math.log10(ag) - Math.log10(1.35e-5)) < 0.05, ag.toExponential(3));
}
// complessazione: AgCl si scioglie in ammoniaca
{
  const { s } = run([['NaCl', 5], ['H2O', 45], ['AgNO3', 10], ['NH3', 40]]);
  const agcl = s.solids.find(x => x.id === 'AgCl(s)');
  const complex = s.species.find(x => x.id === 'Ag(NH3)2+')?.c ?? 0;
  check('AgCl + NH₃ → [Ag(NH₃)₂]⁺', !agcl && complex > 0.009, `[Ag(NH₃)₂⁺] = ${complex.toExponential(2)} M`);
}
// anfoterismo: Al(OH)₃ precipita e si ridiscioglie in eccesso di base
{
  const a = run([['AlCl3', 20], ['H2O', 30], ['NaOH', 5]]).s;
  const b = run([['AlCl3', 20], ['H2O', 30], ['NaOH', 12]]).s;
  check('Al(OH)₃ anfotero', a.solids.some(x => x.id === 'Al(OH)3(s)') && !b.solids.some(x => x.id === 'Al(OH)3(s)'), `con 5 mL di NaOH: ${a.solids.map(x => x.label).join(', ')}; con 12 mL: ${b.solids.map(x => x.label).join(', ') || 'nessun precipitato'}`);
}
// redox: Zn + 2 H⁺ → Zn²⁺ + H₂; Cu + 2 Ag⁺ → Cu²⁺ + 2 Ag
{
  const { s } = run([['HCl', 40], ['Zn(s)', 0.5]]);
  const nZn = 0.5 / 65.38;
  check('Zn + HCl libera H₂', Math.abs(s.gasH2 - nZn) / nZn < 1e-3, `${(s.gasH2 * 1000).toFixed(3)} mmol di H₂ (Zn: ${(nZn * 1000).toFixed(3)} mmol)`);
  const c = run([['AgNO3', 50], ['Cu(s)', 0.5]]).s;
  const ag = c.metals.find(m => m.M === 'Ag');
  check('Cu + 2 Ag⁺ → Cu²⁺ + 2 Ag', ag && Math.abs(ag.n - 0.005) < 1e-5, `${ag ? (ag.n * 1000).toFixed(4) : 0} mmol di Ag depositato (atteso 5,000)`);
  const d = run([['HCl', 40], ['Cu(s)', 0.5]]).s;
  check('il rame non reagisce con HCl (E° > 0)', d.gasH2 === 0, `H₂ = ${d.gasH2}`);
}
// gas: CaCO₃ + 2 HCl → CaCl₂ + H₂O + CO₂
{
  const { s } = run([['HCl', 30], ['CaCO3(s)', 1]]);
  const n = 1 / 100.09;
  check('marmo + acido: CO₂', s.gasCO2 > 0.7 * n && s.gasCO2 < n, `${(s.gasCO2 * 1000).toFixed(2)} mmol di CO₂ liberati (resto disciolto fino alla solubilità di Henry)`);
}
// colore: rame tetraamminico blu, cromato giallo
{
  const cu = run([['CuSO4', 10], ['H2O', 40], ['NH3', 30]]).s.color;
  const cr = run([['K2CrO4', 10], ['H2O', 40]]).s.color;
  check('[Cu(NH₃)₄]²⁺ è blu', cu[2] > cu[0] + 80 && cu[2] > cu[1] + 80, `rgb(${cu})`);
  check('CrO₄²⁻ è giallo', cr[0] > 200 && cr[1] > 200 && cr[2] < 80, `rgb(${cr})`);
}

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nEquilibri in soluzione: tutte le verifiche superate');
