import assert from 'node:assert/strict';
import {
  parseFormula,
  parseEquation,
  checkBalance,
  solveBalance,
  equationText,
  molarMass,
  stoichiometry,
  oxidation,
  redoxLedger,
  ionicFormula,
  dilution,
  isotopeMean,
  significantFigures,
  empiricalFormula,
  rutherford,
  chromatography,
  measurementStatistics,
  waterPhase,
  measuredCalculation,
  NA,
} from '../js/course/chemistry.js';
import {
  REACTIONS,
  COMPOUNDS,
  UNITS,
  CATIONS,
  ANIONS,
  QUESTIONS,
} from '../js/course/content.js';
const near = (a, b, t = 1e-10) =>
  assert.ok(Math.abs(a - b) <= t, `${a} differs from ${b}`);
const eq = (s) => parseEquation(s);
const balanced = (s) => {
  const d = eq(s);
  return solveBalance(d.left, d.right);
};

assert.deepEqual(parseFormula('Ca(OH)2').atoms, { Ca: 1, O: 2, H: 2 });
assert.deepEqual(parseFormula('Al2(SO4)3').atoms, { Al: 2, S: 3, O: 12 });
assert.deepEqual(parseFormula('CuSO4·5H2O').atoms, {
  Cu: 1,
  S: 1,
  O: 9,
  H: 10,
});
assert.deepEqual(parseFormula('[Fe(CN)6]^3-').atoms, { Fe: 1, C: 6, N: 6 });
assert.equal(parseFormula('SO₄²⁻').charge, -2);
assert.equal(parseFormula('Fe3+').charge, 3);
assert.equal(parseFormula('NH4+').atoms.H, 4);
assert.equal(parseFormula('NH4+').charge, 1);
assert.equal(parseFormula('H2O(l)').phase, 'l');
assert.deepEqual(parseFormula('e^-').atoms, {});
assert.equal(parseFormula('e^-').charge, -1);
for (const formula of [
  '',
  'NaClgarbage',
  'Ca(OH',
  'Ca()2',
  'H0',
  'Xx2',
  '(H2]',
  '0H2',
  'H2..O',
  'CuSO4·',
  'H2O^0-',
  '<img>',
])
  assert.throws(() => parseFormula(formula));
assert.throws(() => molarMass('Tc'));
near(molarMass('H2O'), 18.015, 1e-12);
near(molarMass('Ca(OH)2'), 74.092, 1e-10);
assert.throws(() => molarMass('e^-'));

let d = balanced('Fe + O2 -> Fe2O3');
assert.deepEqual(
  d.left.map((s) => s.coeff),
  [4, 3],
);
assert.deepEqual(
  d.right.map((s) => s.coeff),
  [2],
);
d = balanced('MnO4^- + Fe^2+ + H^+ -> Mn^2+ + Fe^3+ + H2O');
assert.deepEqual(
  d.left.map((s) => s.coeff),
  [1, 5, 8],
);
assert.deepEqual(
  d.right.map((s) => s.coeff),
  [1, 5, 4],
);
d = balanced('MnO4^- + H^+ + e^- -> Mn^2+ + H2O');
assert.deepEqual(
  d.left.map((s) => s.coeff),
  [1, 8, 5],
);
assert.deepEqual(
  d.right.map((s) => s.coeff),
  [1, 4],
);
assert.ok(checkBalance(d.left, d.right).balanced);
for (const formula of [
  'H2 -> H2O',
  'H2 + O2 -> H2O + H2O2',
  'H2 + O2 + N2 -> H2O',
]) {
  const q = eq(formula);
  assert.throws(() => solveBalance(q.left, q.right));
}
const bad = eq('H2 + O2 -> H2O');
assert.equal(checkBalance(bad.left, bad.right).balanced, false);
bad.left[0].coeff = 0;
assert.throws(() => checkBalance(bad.left, bad.right));
assert.throws(() => parseEquation('0 H2 + O2 -> H2O'));

for (const r of REACTIONS) {
  const initial = eq(r.equation),
    solved = solveBalance(initial.left, initial.right);
  assert.ok(checkBalance(solved.left, solved.right).balanced, r.id);
  assert.ok(
    checkBalance(...Object.values(eq(equationText(solved.left, solved.right))))
      .balanced,
    r.id,
  );
  if (r.id !== 'emati')
    assert.ok(checkBalance(initial.left, initial.right).balanced, r.id);
  if (r.net) {
    const net = eq(r.net);
    assert.ok(checkBalance(net.left, net.right).balanced, r.id + ' net');
  }
}
for (const c of COMPOUNDS) assert.ok(molarMass(c.formula) > 0);

const water = eq('2 H2 + O2 -> 2 H2O');
let st = stoichiometry(water.left, water.right, [
  { value: 3, unit: 'mol' },
  { value: 2, unit: 'mol' },
]);
near(st.maximum, 1.5);
assert.deepEqual(st.limiting, [0]);
near(st.remaining[1], 0.5);
near(st.products[0].mol, 3);
st = stoichiometry(
  water.left,
  water.right,
  [
    { value: 3, unit: 'mol' },
    { value: 2, unit: 'mol' },
  ],
  0.5,
  75,
);
near(st.remaining[0], 1.5);
near(st.products[0].mol, 1.5);
near(st.products[0].isolatedGrams, 3 * 18.015 * 0.75);
near(
  3 * molarMass('H2') + 2 * molarMass('O2'),
  st.remaining[0] * molarMass('H2') +
    st.remaining[1] * molarMass('O2') +
    st.products[0].grams,
  1e-10,
);
const grams = stoichiometry(water.left, water.right, [
  { value: 3 * molarMass('H2'), unit: 'g' },
  { value: 2 * molarMass('O2'), unit: 'g' },
]);
near(grams.maximum, 1.5);
for (let i = 0; i < 25; i++) {
  const a = i * 0.71,
    b = (25 - i) * 0.43,
    f = i / 24,
    result = stoichiometry(
      water.left,
      water.right,
      [
        { value: a, unit: 'mol' },
        { value: b, unit: 'mol' },
      ],
      f,
    );
  near(
    a * molarMass('H2') + b * molarMass('O2'),
    result.remaining[0] * molarMass('H2') +
      result.remaining[1] * molarMass('O2') +
      result.products[0].grams,
    1e-9,
  );
  assert.ok(result.remaining.every((x) => x >= 0));
}
assert.throws(() =>
  stoichiometry(water.left, water.right, [
    { value: -1, unit: 'mol' },
    { value: 1, unit: 'mol' },
  ]),
);
assert.throws(() =>
  stoichiometry(
    water.left,
    water.right,
    [
      { value: 1, unit: 'mol' },
      { value: 1, unit: 'mol' },
    ],
    1,
    101,
  ),
);

assert.equal(oxidation('H2O2').states.O, -1);
assert.equal(oxidation('NaH').states.H, -1);
assert.equal(oxidation('OF2').states.O, 2);
assert.equal(oxidation('KO2').states.O, -0.5);
assert.equal(oxidation('SO4^2-').states.S, 6);
near(oxidation('Fe3O4').states.Fe, 8 / 3);
assert.ok(oxidation('C2H5NO2').unresolved);
for (const [rxId, expected] of [
  ['termite', 6],
  ['metano', 8],
  ['solfito', 4],
  ['perossido', 2],
  ['cloro', 1],
  ['permanganato', 5],
  ['dicromato', 6],
  ['argento', 0],
]) {
  const r = REACTIONS.find((x) => x.id === rxId),
    reaction = eq(r.equation),
    ledger = redoxLedger(reaction.left, reaction.right);
  near(ledger.lost, expected);
  near(ledger.gained, expected);
  assert.ok(ledger.balanced, rxId);
}
const sulfur = eq('2 SO2 + O2 -> 2 SO3');
const oxygen = redoxLedger(sulfur.left, sulfur.right).changes.find(
  (x) => x.el === 'O',
);
assert.equal(oxygen.atoms, 2);
assert.equal(oxygen.from, 0);
assert.ok(redoxLedger(d.left, d.right).external);

for (const c of CATIONS)
  for (const a of ANIONS) {
    const out = ionicFormula(c, a),
      formula = parseFormula(out.formula);
    assert.equal(formula.charge, 0);
    near(
      out.cationCount * parseFormula(c).charge +
        out.anionCount * parseFormula(a).charge,
      0,
    );
  }
assert.equal(ionicFormula('Al^3+', 'SO4^2-').formula, 'Al2(SO4)3');
assert.equal(ionicFormula('NH4^+', 'SO4^2-').formula, '(NH4)2SO4');
assert.throws(() => ionicFormula('Cl^-', 'Na^+'));
near(dilution(2, 50, 200).finalConcentration, 0.5);
near(dilution(2, 50, 200).mol, 0.1);
assert.throws(() => dilution(2, 200, 50));
assert.throws(() => dilution(2, 0, 100));
near(isotopeMean(35, 37, 77.5), 35.45);
near(isotopeMean(35, 37, 100), 35);
near(isotopeMean(35, 37, 0), 37);
assert.throws(() => isotopeMean(35, 37, 101));
assert.equal(significantFigures('0.00450'), 3);
assert.equal(significantFigures('1.20e3'), 3);
assert.equal(significantFigures('1200.'), 4);
assert.equal(
  empiricalFormula([
    { el: 'C', mass: 40 },
    { el: 'H', mass: 6.71 },
    { el: 'O', mass: 53.29 },
  ]).formula,
  'CH2O',
);
assert.throws(() =>
  empiricalFormula([
    { el: 'C', mass: 0 },
    { el: 'H', mass: 2 },
  ]),
);
assert.equal(NA, 6.02214076e23);

near(chromatography(1, 1, 0.75).Rf, 0.5);
near(chromatography(4, 1, 0.75).fractionAlongPaper, 0.15);
near(chromatography(0, 1, 1).Rf, 1);
assert.throws(() => chromatography(-1));
assert.equal(waterPhase(-10), 'solido');
assert.equal(waterPhase(0), 'solido + liquido');
assert.equal(waterPhase(20), 'liquido');
assert.equal(waterPhase(100), 'liquido + gas');
assert.equal(waterPhase(110), 'gas');
assert.throws(() => waterPhase(-274));
const stats = measurementStatistics([98, 100, 102]);
near(stats.mean, 100);
near(stats.sd, 2);
assert.throws(() => measurementStatistics([1]));
assert.equal(measuredCalculation('1.20', '2.3', 'multiply').reported, '2.8');
assert.equal(measuredCalculation('12.11', '18.0', 'add').reported, '30.1');
assert.equal(measuredCalculation('1.20e3', '12', 'add').reported, '1210');
assert.throws(() => measuredCalculation('1', '0', 'divide'));
assert.throws(() =>
  stoichiometry(water.left, water.right, [
    { value: 1, unit: 'ml' },
    { value: 1, unit: 'mol' },
  ]),
);

// Independent invariants of the Coulomb path: energy and angular momentum.
for (const b of [0, 20, 200]) {
  const Z = 79,
    E = 5,
    model = rutherford({ Z, E, b, points: 1001 }),
    v0 = Math.sqrt(2 * E),
    u = b === 0 ? 2.8 : 3.2;
  const T = ((2 * model.a) / v0) * (model.eccentricity * Math.sinh(u) + u),
    k = 2 * Z * 1.439964548;
  near(model.theta, b === 0 ? Math.PI : 2 * Math.atan(k / (2 * E * b)), 1e-12);
  near(model.closest, model.a + Math.hypot(model.a, b), 1e-10);
  near(model.path[0].t, 0, 1e-12);
  near(model.path.at(-1).t, 1, 1e-12);
  for (let i = 10; i < 991; i += 11) {
    const p = model.path[i],
      a = model.path[i - 1],
      c = model.path[i + 1],
      dt = (c.t - a.t) * T,
      vx = (c.x - a.x) / dt,
      vy = (c.y - a.y) / dt;
    near((vx * vx + vy * vy) / 2 + k / Math.hypot(p.x, p.y), E, 2e-5);
    near(p.x * vy - p.y * vx, -b * v0, Math.max(1, b * v0) * 5e-6);
  }
}
assert.throws(() => rutherford({ E: 0 }));
assert.throws(() => rutherford({ b: -1 }));
for (const q of QUESTIONS) assert.ok(UNITS.some((u) => u.id === q.unit));
console.log(
  'Course chemistry: parsing, exact balancing, charge, mass conservation, limiting reagent, redox ledger, nomenclature ions, SI calculations and Coulomb trajectory invariants passed.',
);
