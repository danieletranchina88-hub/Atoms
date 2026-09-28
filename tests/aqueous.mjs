// Verifiche del solutore di equilibri in soluzione acquosa (becher, database MINTEQ v4) con casi da manuale.
// Valori di riferimento: Harris, Quantitative Chemical Analysis; CRC Handbook (Ksp, ΔsolH); Bandura e Lvov,
// J. Phys. Chem. Ref. Data 35, 15 (2006) per il prodotto ionico dell'acqua.

import { Beaker, equilibrate } from '../js/chem/aqueous.js';
import { CATIONS, ANIONS, saltRecipe, solutionOf, solidOf, SPECIAL, specialRecipe, METAL_SOLIDS, metalRecipe, MINERAL_SOLIDS, mineralRecipe, WATER, beakerPresets } from '../js/chem/beakerReagents.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};
const ci = (sp) => CATIONS.findIndex(x => x[0] === sp);
const ai = (sp) => ANIONS.findIndex(x => x[0] === sp);
const sol = (cat, an, c) => solutionOf(saltRecipe(ci(cat), ai(an)), c);
const NH3 = specialRecipe(SPECIAL.find(x => x.id === 'NH3'), 1);
const NH3c = specialRecipe(SPECIAL.find(x => x.id === 'NH3'), 6);
const metal = (id) => metalRecipe(METAL_SOLIDS.find(m => m.id === id));
const mineral = (ph) => mineralRecipe(MINERAL_SOLIDS.find(m => m.phase === ph));
const R = {
  HCl: sol('H+', 'Cl-', 1), NaOH: sol('Na+', 'OH-', 1), HAc: sol('H+', 'Acetate-', 1), NaAc: sol('Na+', 'Acetate-', 1),
  NaCl: sol('Na+', 'Cl-', 1), AgNO3: sol('Ag+', 'NO3-', 0.1), CuSO4: sol('Cu+2', 'SO4-2', 0.5), AlCl3: sol('Al+3', 'Cl-', 0.1),
  FeSO4: sol('Fe+2', 'SO4-2', 0.1), KMnO4: sol('K+', 'MnO4-', 0.02), H2SO4: sol('H+', 'SO4-2', 0.5), K2CrO4: sol('K+', 'CrO4-2', 0.1),
  FeNO33: sol('Fe+3', 'NO3-', 0.1), FeNO32: sol('Fe+2', 'NO3-', 0.1), HNO3: sol('H+', 'NO3-', 1), NH3, H2O: WATER,
};
const run = (steps) => { const b = new Beaker(); let e = null; for (const [r, v] of steps) e = b.add(r, v); return { b, s: b.summary(), e }; };
const solid = (s, n) => s.solids.find(x => x.n === n);
const conc = (s, n) => s.species.find(x => x.n === n)?.c ?? 0;

// acqua pura e dipendenza dalla temperatura: pKw(25 °C) = 13,995, pKw(60 °C) = 13,017 (Bandura e Lvov 2006)
{
  const w = equilibrate({}, 0.1);
  check('acqua pura a 25 °C: pH 7,00', Math.abs(w.pH - 7) < 0.01, w.pH.toFixed(4));
  const h = equilibrate({}, 0.1, 333.15);
  check('acqua pura a 60 °C: pH = pKw/2 ≈ 6,51', Math.abs(h.pH - 6.51) < 0.03, h.pH.toFixed(3));
}
// acido debole: pH ≈ 2,88 per CH₃COOH 0,1 M
{
  const { s } = run([[R.HAc, 10], [R.H2O, 90]]);
  check('pH di CH₃COOH 0,1 M', Math.abs(s.pH - 2.88) < 0.03, s.pH.toFixed(3));
}
// tampone equimolare: pH = pKa + log(γA/γHA) ≈ 4,63 a I = 0,1
{
  const { s } = run([[R.HAc, 10], [R.NaAc, 10], [R.H2O, 80]]);
  check('tampone acetico equimolare', Math.abs(s.pH - 4.65) < 0.08, `pH ${s.pH.toFixed(3)} (pKa 4,76 − correzione di attività)`);
}
// base debole: NH₃ 0,1 M → pH ≈ 11,13
{
  const { s } = run([[R.NH3, 10], [R.H2O, 90]]);
  check('pH di NH₃ 0,1 M', Math.abs(s.pH - 11.13) < 0.04, s.pH.toFixed(3));
}
// acido forte diluito
{
  const { s } = run([[R.HCl, 1], [R.H2O, 99]]);
  check('pH di HCl 0,01 M', Math.abs(s.pH - 2.04) < 0.03, s.pH.toFixed(3));
}
// neutralizzazione: ΔrH = −55,8 kJ/mol, ΔT = q/(m c_p); il pH neutro a 31,7 °C è pKw(T)/2 ≈ 6,86
{
  const { s } = run([[R.HCl, 25], [R.NaOH, 25]]);
  const dT = s.T - 298.15;
  check('riscaldamento per neutralizzazione', Math.abs(s.heat - 0.025 * 55.8) < 0.03 && Math.abs(dT - 6.67) < 0.15, `q = ${s.heat.toFixed(3)} kJ, ΔT = ${dT.toFixed(2)} K (atteso ≈ 1,40 kJ, 6,7 K)`);
  check('pH neutro alla temperatura raggiunta', Math.abs(s.pH - 6.86) < 0.03, `pH ${s.pH.toFixed(3)} a ${(s.T - 273.15).toFixed(1)} °C`);
}
// prodotto di solubilità: AgCl (log Ksp = −9,75)
{
  const { s } = run([[R.NaCl, 5], [R.H2O, 45], [R.AgNO3, 10]]);
  const agcl = solid(s, 'Chlorargyrite') ?? s.solids.find(x => /AgCl/.test(x.label));
  check('AgCl precipita', agcl && Math.abs(agcl.amt - 1e-3) < 2e-5, `${agcl ? (agcl.amt * 1000).toFixed(4) : 0} mmol (atteso 1,000)`);
  const { s: s2 } = run([[R.NaCl, 5], [R.H2O, 45], [R.AgNO3, 10], [R.NH3, 40]]);
  check('AgCl si scioglie in NH₃ come [Ag(NH₃)₂]⁺', !s2.solids.length && conc(s2, 'Ag(NH3)2+') > 0.009, `[Ag(NH₃)₂⁺] = ${conc(s2, 'Ag(NH3)2+').toExponential(2)} M`);
}
// anfoterismo di Al(OH)₃
{
  const a = run([[R.AlCl3, 20], [R.H2O, 30], [R.NaOH, 5]]).s;
  const b = run([[R.AlCl3, 20], [R.H2O, 30], [R.NaOH, 12]]).s;
  check('Al(OH)₃ anfotero', a.solids.length === 1 && !b.solids.length && conc(b, 'Al(OH)4-') > 0.03, `5 mL NaOH: ${a.solids.map(x => x.label).join(', ')}; 12 mL: ${b.solids.map(x => x.label).join(', ') || 'nessun precipitato'}`);
}
// redox: metalli e acidi (scala dei potenziali standard)
{
  const nZn = 0.5 / 65.38;
  const z = run([[R.HCl, 40], [metal('Zn'), 0.5]]).s;
  check('Zn + 2 H⁺ → Zn²⁺ + H₂', Math.abs(z.gasH2 - nZn) / nZn < 1e-3, `${(z.gasH2 * 1000).toFixed(3)} mmol di H₂ (Zn: ${(nZn * 1000).toFixed(3)} mmol)`);
  const d = run([[R.HCl, 40], [metal('Cu'), 0.5]]).s;
  check('il rame non reagisce con HCl (E° > 0)', d.gasH2 < 1e-9 && solid(d, 'Cumetal')?.amt > 0.0078, `H₂ = ${d.gasH2.toExponential(1)} mol`);
  const c = run([[R.AgNO3, 50], [metal('Cu'), 0.5]]).s;
  const ag = solid(c, 'Agmetal');
  check('Cu + 2 Ag⁺ → Cu²⁺ + 2 Ag', ag && Math.abs(ag.amt - 0.005) < 1e-5, `${ag ? (ag.amt * 1000).toFixed(4) : 0} mmol di Ag depositato (atteso 5,000)`);
}
// titolazione redox: 5 Fe²⁺ + MnO₄⁻ + 8 H⁺ → 5 Fe³⁺ + Mn²⁺ + 4 H₂O
{
  const { s } = run([[R.FeSO4, 25], [R.H2SO4, 10], [R.KMnO4, 5]]);
  const mn = s.species.filter(x => /^Mn/.test(x.n) && !/MnO4/.test(x.n)).reduce((a, x) => a + x.c, 0) * s.V / 1000;
  check('permanganato ridotto a Mn²⁺ dal ferro(II)', Math.abs(mn - 1e-4) < 1e-6 && conc(s, 'MnO4-') < 1e-12, `Mn(II) ${(mn * 1e3).toFixed(4)} mmol (atteso 0,1000)`);
  check('Eh definito dalla coppia Fe³⁺/Fe²⁺', s.Eh !== null && s.Eh > 0.6 && s.Eh < 0.8, `Eh = ${s.Eh?.toFixed(3)} V`);
}
// Nernst: Fe³⁺/Fe²⁺ equimolari in acido forte → Eh ≈ potenziale formale (E° = 0,771 V; E°' = 0,732 V in HClO₄ 1 M, Harris)
{
  const { s } = run([[R.FeNO33, 10], [R.FeNO32, 10], [R.HNO3, 10]]);
  check('Nernst per Fe³⁺/Fe²⁺', s.Eh !== null && Math.abs(s.Eh - 0.732) < 0.03, `Eh = ${s.Eh?.toFixed(3)} V (E°' ≈ 0,73 V)`);
  const { s: s2 } = run([[R.CuSO4, 10], [R.H2O, 40]]);
  check('senza coppia redox l\'Eh non è definito', s2.Eh === null, `Eh = ${s2.Eh}`);
}
// gas: CaCO₃ + 2 HCl → CaCl₂ + H₂O + CO₂
{
  const { s } = run([[R.HCl, 30], [mineral('Calcite'), 1]]);
  const n = 1 / 100.09;
  check('marmo + acido: CO₂', s.gasCO2 > 0.7 * n && s.gasCO2 < n, `${(s.gasCO2 * 1000).toFixed(2)} mmol di CO₂ liberati (il resto resta disciolto secondo Henry)`);
}
// colori: tetraamminorame blu, cromato giallo, dicromato arancio
{
  const cu = run([[R.CuSO4, 10], [R.H2O, 40], [NH3c, 20]]).s;
  check('[Cu(NH₃)₄]²⁺ è blu', cu.color[2] > cu.color[0] + 80 && cu.color[2] > cu.color[1] + 80 && !cu.solids.length, `rgb(${cu.color})`);
  const cr = run([[R.K2CrO4, 10], [R.H2O, 40]]).s.color;
  check('CrO₄²⁻ è giallo', cr[0] > 200 && cr[1] > 200 && cr[2] < 80, `rgb(${cr})`);
}
// solido sciolto in acqua: NH₄NO₃ raffredda (ΔsolH = +25,7 kJ/mol, CRC)
{
  const salt = solidOf(saltRecipe(ci('NH4+'), ai('NO3-')));
  const { s } = run([[R.H2O, 50], [salt, 4]]);
  check('NH₄NO₃ si scioglie assorbendo calore', s.T < 298.15 - 1.5, `ΔT = ${(s.T - 298.15).toFixed(2)} K`);
}
// libertà del banco: ogni combinazione catione × anione dà un equilibrio finito
{
  let bad = [], n = 0;
  for (let i = 0; i < CATIONS.length; i++) for (let j = 0; j < ANIONS.length; j++) {
    const b = new Beaker();
    b.add(R.H2O, 40);
    b.add(solutionOf(saltRecipe(i, j), 0.1), 10);
    const s = b.summary();
    n++;
    if (!Number.isFinite(s.pH) || !Number.isFinite(s.T) || !s.converged) bad.push(saltRecipe(i, j).formula);
  }
  check('tutte le combinazioni catione × anione convergono', bad.length === 0, `${n - bad.length}/${n}${bad.length ? `; non convergono: ${bad.slice(0, 12).join(', ')}` : ''}`);
}
// tutte le esperienze guidate girano fino in fondo
{
  const bad = [];
  for (const p of beakerPresets()) {
    const b = new Beaker();
    for (const [r, v] of p.steps) b.add(r, v);
    if (p.next) for (let k = 0; k < 5; k++) b.add(p.next, p.amount);
    const s = b.summary();
    if (!Number.isFinite(s.pH) || !s.converged) bad.push(p.id);
  }
  check('esperienze guidate', bad.length === 0, bad.join(', ') || 'tutte convergono');
}

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nEquilibri in soluzione: tutte le verifiche superate');
