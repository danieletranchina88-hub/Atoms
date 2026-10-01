// Reagenti del becher: cationi e anioni combinabili liberamente, acidi, basi, metalli, indicatori, solidi.
// Ogni ricetta dice quali specie del database si mettono in soluzione (o quale fase solida si aggiunge).

import { molarMass, phaseByName } from './aqueous.js';
import { SOLUTION_ENTHALPY } from './aqueousExtra.js';

// [specie del database, formula del frammento, carica, nome italiano]
export const CATIONS = [
  ['H+', 'H', 1, 'idrogeno (acido)'],
  ['Li+', 'Li', 1, 'litio'], ['Na+', 'Na', 1, 'sodio'], ['K+', 'K', 1, 'potassio'], ['NH4+', 'NH4', 1, 'ammonio'],
  ['Mg+2', 'Mg', 2, 'magnesio'], ['Ca+2', 'Ca', 2, 'calcio'], ['Sr+2', 'Sr', 2, 'stronzio'], ['Ba+2', 'Ba', 2, 'bario'],
  ['Al+3', 'Al', 3, 'alluminio'], ['Cr+3', 'Cr', 3, 'cromo(III)'], ['Mn+2', 'Mn', 2, 'manganese(II)'],
  ['Fe+2', 'Fe', 2, 'ferro(II)'], ['Fe+3', 'Fe', 3, 'ferro(III)'], ['Co+2', 'Co', 2, 'cobalto(II)'], ['Ni+2', 'Ni', 2, 'nichel(II)'],
  ['Cu+2', 'Cu', 2, 'rame(II)'], ['Zn+2', 'Zn', 2, 'zinco'], ['Cd+2', 'Cd', 2, 'cadmio'], ['Ag+', 'Ag', 1, 'argento'], ['Pb+2', 'Pb', 2, 'piombo(II)'],
];
export const ANIONS = [
  ['OH-', 'OH', -1, 'idrossido'], ['F-', 'F', -1, 'fluoruro'], ['Cl-', 'Cl', -1, 'cloruro'], ['Br-', 'Br', -1, 'bromuro'], ['I-', 'I', -1, 'ioduro'],
  ['NO3-', 'NO3', -1, 'nitrato'], ['SO4-2', 'SO4', -2, 'solfato'], ['HSO4-', 'HSO4', -1, 'idrogenosolfato'],
  ['CO3-2', 'CO3', -2, 'carbonato'], ['HCO3-', 'HCO3', -1, 'idrogenocarbonato'],
  ['PO4-3', 'PO4', -3, 'fosfato'], ['HPO4-2', 'HPO4', -2, 'idrogenofosfato'], ['H2PO4-', 'H2PO4', -1, 'diidrogenofosfato'],
  ['S-2', 'S', -2, 'solfuro'], ['HS-', 'HS', -1, 'idrogenosolfuro'],
  ['CrO4-2', 'CrO4', -2, 'cromato'], ['Cr2O7-2', 'Cr2O7', -2, 'dicromato'], ['MnO4-', 'MnO4', -1, 'permanganato'],
  ['SCN-', 'SCN', -1, 'tiocianato'],
  ['Acetate-', 'CH3COO', -1, 'acetato'], ['Formate-', 'HCOO', -1, 'formiato'], ['Citrate-3', 'C6H5O7', -3, 'citrato'],
  ['Tartarate-2', 'C4H4O6', -2, 'tartrato'], ['Phthalate-2', 'C8H4O4', -2, 'ftalato'], ['Benzoate-', 'C6H5COO', -1, 'benzoato'],
  ['Salicylate-2', 'C7H4O3', -2, 'salicilato'], ['Glycine-', 'NH2CH2COO', -1, 'glicinato'], ['Edta-4', 'C10H12N2O8', -4, 'EDTA'],
];

const gcd = (a, b) => (b ? gcd(b, a % b) : a);
const SUB = '₀₁₂₃₄₅₆₇₈₉';
export const subscript = (f) => f.replace(/(\d+)/g, d => d.split('').map(x => SUB[+x]).join(''));
const needsParens = (frag) => /[A-Z].*[A-Z]|\d/.test(frag.replace(/^[A-Z][a-z]?$/, ''));

/** Sale neutro da catione e anione: formula, massa molare e specie in soluzione. */
export function saltRecipe(catIdx, anIdx) {
  const [cs, cf, cz, cname] = CATIONS[catIdx];
  const [as, af, az, aname] = ANIONS[anIdx];
  const g = gcd(cz, -az);
  const nc = -az / g, na = cz / g;
  let formula;
  if (cs === 'H+' && as === 'Acetate-') formula = 'CH3COOH';
  else if (cs === 'H+' && as === 'Formate-') formula = 'HCOOH';
  else if (cs === 'H+' && as === 'Benzoate-') formula = 'C6H5COOH';
  else {
    const part = (frag, n) => (n === 1 ? frag : needsParens(frag) ? `(${frag})${n}` : `${frag}${n}`);
    formula = part(cf, nc) + part(af, na);
  }
  const name = cs === 'H+' ? acidName(as, aname) : `${aname} di ${cname}`;
  const dsol = SOLUTION_ENTHALPY[formula] ?? null;
  return {
    label: subscript(formula), formula, name,
    species: [[cs, nc], [as, na]],
    M: molarMass(formula.replace('CH3COO', 'C2H3O2').replace('HCOO', 'CHO2').replace('C6H5COO', 'C7H5O2').replace('NH2CH2COO', 'C2H4NO2')),
    dsol,
  };
}

function acidName(as, aname) {
  const special = {
    'Cl-': 'acido cloridrico', 'Br-': 'acido bromidrico', 'I-': 'acido iodidrico', 'F-': 'acido fluoridrico', 'NO3-': 'acido nitrico',
    'SO4-2': 'acido solforico', 'PO4-3': 'acido fosforico', 'CO3-2': 'acido carbonico', 'Acetate-': 'acido acetico',
    'Formate-': 'acido formico', 'Citrate-3': 'acido citrico', 'S-2': 'acido solfidrico', 'OH-': 'acqua', 'CrO4-2': 'acido cromico',
    'Tartarate-2': 'acido tartarico', 'Phthalate-2': 'acido ftalico', 'Benzoate-': 'acido benzoico', 'Salicylate-2': 'acido salicilico',
    'Glycine-': 'glicina', 'Edta-4': 'acido etilendiamminotetraacetico', 'SCN-': 'acido tiocianico', 'MnO4-': 'acido permanganico',
  };
  return special[as] ?? `acido (${aname})`;
}

export const solutionOf = (salt, c) => ({ ...salt, kind: 'solution', c, label: `${salt.label} ${c.toLocaleString('it-IT')} M` });
export const solidOf = (salt) => ({ ...salt, kind: 'solid' });

// reagenti molecolari e speciali
export const SPECIAL = [
  { id: 'NH3', label: 'NH₃', name: 'ammoniaca', species: [['NH3', 1]], M: 17.031, c: 1 },
  { id: 'en', label: 'en', name: 'etilendiammina', species: [['Ethylenediamine', 1]], M: 60.10, c: 1 },
  { id: 'KHP', label: 'KHC₈H₄O₄', name: 'ftalato acido di potassio (standard primario)', species: [['K+', 1], ['H+', 1], ['Phthalate-2', 1]], M: 204.22, c: 0.1 },
  { id: 'Na2H2EDTA', label: 'Na₂H₂EDTA', name: 'EDTA disodico', species: [['Na+', 2], ['H+', 2], ['Edta-4', 1]], M: 336.21, c: 0.05 },
  { id: 'PP', label: 'fenolftaleina', name: 'fenolftaleina (indicatore)', species: [['HPP', 1]], M: 318.3, c: 0.003, indicator: true },
  { id: 'BTB', label: 'blu di bromotimolo', name: 'blu di bromotimolo (indicatore)', species: [['HBTB', 1]], M: 624.4, c: 0.0016, indicator: true },
  { id: 'MO', label: 'metilarancio', name: 'metilarancio (indicatore)', species: [['HMO', 1]], M: 327.3, c: 0.0015, indicator: true },
];

export const METAL_SOLIDS = [
  { id: 'Mg', phase: 'Mgmetal', name: 'magnesio (nastro)' },
  { id: 'Zn', phase: 'Znmetal', name: 'zinco (granuli)' },
  { id: 'Fe', phase: 'Femetal', name: 'ferro (limatura)' },
  { id: 'Pb', phase: 'Pbmetal', name: 'piombo' },
  { id: 'Cd', phase: 'Cdmetal(alpha)', name: 'cadmio' },
  { id: 'Cu', phase: 'Cumetal', name: 'rame (filo)' },
  { id: 'Ag', phase: 'Agmetal', name: 'argento' },
];
export const MINERAL_SOLIDS = [
  { phase: 'Malachite', name: 'malachite (Cu₂CO₃(OH)₂)' },
  { phase: 'Calcite', name: 'marmo, calcite (CaCO₃)' },
  { phase: 'Gypsum', name: 'gesso (CaSO₄·2H₂O)' },
  { phase: 'Pyrolusite', name: 'biossido di manganese (MnO₂)' },
  { phase: 'Fluorite', name: 'fluorite (CaF₂)' },
  { phase: 'Portlandite', name: 'calce spenta (Ca(OH)₂)' },
];

export function metalRecipe(m) {
  const p = phaseByName(m.phase);
  return { label: m.id, name: m.name, kind: 'metal', phase: m.phase, M: molarMass(p.f), formula: m.id };
}
export function mineralRecipe(m) {
  const p = phaseByName(m.phase);
  return { label: m.name, name: m.name, kind: 'phase', phase: m.phase, M: molarMass(p.f), formula: p.f };
}
export function specialRecipe(s, c = s.c) {
  return { label: `${s.label}${s.indicator ? '' : ` ${c.toLocaleString('it-IT')} M`}`, name: s.name, kind: 'solution', species: s.species, c, M: s.M, formula: s.label };
}
export const WATER = { label: 'H₂O', name: 'acqua distillata', kind: 'water', species: [], c: 0 };

// scaffale rapido: [catione, anione, concentrazione]
const idx = (list, sp) => list.findIndex(x => x[0] === sp);
export function quickShelf() {
  const s = (cat, an, c) => solutionOf(saltRecipe(idx(CATIONS, cat), idx(ANIONS, an)), c);
  return [
    ['Acidi', [s('H+', 'Cl-', 1), s('H+', 'SO4-2', 0.5), s('H+', 'NO3-', 1), s('H+', 'Acetate-', 1), s('H+', 'PO4-3', 0.5)]],
    ['Basi', [s('Na+', 'OH-', 1), specialRecipe(SPECIAL[0], 1), specialRecipe(SPECIAL[0], 6), s('Na+', 'CO3-2', 0.5), s('Na+', 'HCO3-', 0.5)]],
    ['Sali', [s('Na+', 'Cl-', 1), s('Na+', 'Acetate-', 1), s('NH4+', 'Cl-', 1), s('K+', 'I-', 0.1), s('K+', 'SCN-', 0.1), s('Na+', 'SO4-2', 0.1), s('K+', 'CrO4-2', 0.1)]],
    ['Ioni metallici', [s('Ag+', 'NO3-', 0.1), s('Pb+2', 'NO3-', 0.1), s('Ba+2', 'Cl-', 0.1), s('Ca+2', 'Cl-', 0.5), s('Cu+2', 'SO4-2', 0.5), s('Cu+2', 'Cl-', 0.1), s('Fe+3', 'Cl-', 0.1), s('Fe+2', 'SO4-2', 0.1), s('Al+3', 'Cl-', 0.1), s('Ni+2', 'SO4-2', 0.5), s('Zn+2', 'SO4-2', 0.5)]],
    ['Ossidanti', [s('K+', 'MnO4-', 0.02), s('K+', 'Cr2O7-2', 0.0167)]],
    ['Indicatori', SPECIAL.filter(x => x.indicator).map(x => specialRecipe(x))],
    ['Metalli e solidi', [...METAL_SOLIDS.map(metalRecipe), ...MINERAL_SOLIDS.map(mineralRecipe)]],
  ];
}

// esperienze: sequenze di [ricetta, quantità]
export function beakerPresets() {
  const q = (cat, an, c) => solutionOf(saltRecipe(idx(CATIONS, cat), idx(ANIONS, an)), c);
  const PP = specialRecipe(SPECIAL.find(x => x.id === 'PP'));
  const BTB = specialRecipe(SPECIAL.find(x => x.id === 'BTB'));
  const NH3 = specialRecipe(SPECIAL[0], 1);
  const NH3c = specialRecipe(SPECIAL[0], 6);
  const metal = (id) => metalRecipe(METAL_SOLIDS.find(m => m.id === id));
  return [
    { id: 'neutral', name: 'Neutralizzazione HCl + NaOH', steps: [[q('H+', 'Cl-', 1), 25], [PP, 0.3]], next: q('Na+', 'OH-', 1), amount: 1,
      text: 'H⁺ + OH⁻ → H₂O libera 55,8 kJ/mol (ΔrH di MINTEQ): la soluzione si scalda. Versa NaOH a piccole dosi: al punto equivalente la fenolftaleina diventa rosa.' },
    { id: 'permanganometry', name: 'Permanganometria: Fe²⁺ + MnO₄⁻', steps: [[q('Fe+2', 'SO4-2', 0.1), 25], [q('H+', 'SO4-2', 0.5), 10]], next: q('K+', 'MnO4-', 0.02), amount: 1,
      text: 'MnO₄⁻ + 5 Fe²⁺ + 8 H⁺ → Mn²⁺ + 5 Fe³⁺ + 4 H₂O. Il permanganato si decolora finché resta Fe²⁺; al punto equivalente (25 mL) la soluzione diventa viola. Il potenziale redox Eh salta come il pH in una titolazione acido-base.' },
    { id: 'silver', name: 'Alogenuri d\'argento e ammoniaca', steps: [[q('Na+', 'Cl-', 1), 5], [WATER, 45], [q('Ag+', 'NO3-', 0.1), 10]], next: NH3, amount: 5,
      text: 'Ag⁺ + Cl⁻ → AgCl(s), bianco. L\'ammoniaca scioglie il precipitato formando [Ag(NH₃)₂]⁺; lo ioduro invece riprecipita AgI, molto meno solubile.' },
    { id: 'copper', name: 'Rame: sale basico e complesso tetraamminico', steps: [[q('Cu+2', 'SO4-2', 0.5), 10], [WATER, 40]], next: NH3c, amount: 0.5,
      text: 'Poca ammoniaca fa precipitare un solfato basico di rame celeste (brochantite). In eccesso di NH₃ il precipitato si scioglie e si forma [Cu(NH₃)₄]²⁺, blu intenso. Serve ammoniaca 6 M: con NH₃ 1 M l\'ammoniaca libera resta troppo diluita e, secondo Ksp e β₄ del database, parte dell\'idrossido non si scioglie.' },
    { id: 'iron', name: 'Ferro(III) e tiocianato', steps: [[WATER, 50], [q('H+', 'NO3-', 1), 1], [q('Fe+3', 'Cl-', 0.1), 0.5], [q('K+', 'SCN-', 0.1), 0.5]], next: q('Na+', 'OH-', 1), amount: 1,
      text: 'Fe³⁺ + SCN⁻ ⇌ FeSCN²⁺, rosso (log K° = 2,85). Con NaOH precipita Fe(OH)₃ e il colore scompare: principio di Le Châtelier.' },
    { id: 'zinc', name: 'Zinco e acido: idrogeno', steps: [[q('H+', 'Cl-', 1), 40], [metal('Zn'), 0.5]], next: metal('Zn'), amount: 0.5,
      text: 'Zn + 2 H⁺ → Zn²⁺ + H₂↑. Lo zinco si scioglie finché il potenziale della soluzione resta sotto quello dell\'idrogeno (H₂ a fugacità unitaria; velocità non calcolata).' },
    { id: 'cu-ag', name: 'Rame in nitrato d\'argento', steps: [[q('Ag+', 'NO3-', 0.1), 50], [metal('Cu'), 0.5]], next: metal('Cu'), amount: 0.5,
      text: 'Cu + 2 Ag⁺ → Cu²⁺ + 2 Ag (E° = +0,46 V): si deposita argento e la soluzione diventa azzurra.' },
    { id: 'buffer', name: 'Tampone acetico', steps: [[q('H+', 'Acetate-', 1), 20], [q('Na+', 'Acetate-', 1), 20], [BTB, 0.5]], next: q('H+', 'Cl-', 1), amount: 1,
      text: 'pH ≈ pKa + log([A⁻]/[HA]). Aggiungi HCl o NaOH: il pH cambia poco finché resta una scorta di entrambe le forme.' },
    { id: 'marble', name: 'Marmo e acido: CO₂', steps: [[q('H+', 'Cl-', 1), 30], [mineralRecipe(MINERAL_SOLIDS.find(m=>m.phase==='Calcite')), 1]], next: q('H+', 'Cl-', 1), amount: 5,
      text: 'CaCO₃ + 2 H⁺ → Ca²⁺ + H₂O + CO₂↑: il gas si libera quando la sua pressione di equilibrio supera 1 atm (legge di Henry).' },
    { id: 'amphoteric', name: 'Idrossidi anfoteri', steps: [[q('Al+3', 'Cl-', 0.1), 20], [WATER, 30]], next: q('Na+', 'OH-', 1), amount: 1,
      text: 'NaOH fa precipitare Al(OH)₃ bianco, che in eccesso di base si ridiscioglie come Al(OH)₄⁻.' },
    { id: 'chromate', name: 'Cromato e dicromato', steps: [[q('K+', 'CrO4-2', 0.1), 10], [WATER, 40]], next: q('H+', 'Cl-', 1), amount: 2,
      text: '2 CrO₄²⁻ (giallo) + 2 H⁺ ⇌ Cr₂O₇²⁻ (arancione) + H₂O.' },
    { id: 'hardwater', name: 'Acqua dura e addolcimento', steps: [[q('Ca+2', 'Cl-', 0.5), 2], [q('Na+', 'HCO3-', 0.5), 4], [WATER, 94]], next: q('Na+', 'CO3-2', 0.5), amount: 1,
      text: 'Ca²⁺ e HCO₃⁻ sono gli ioni dell\'acqua dura. Scaldando (piastra a 80 °C) il carbonato di calcio precipita perché la sua solubilità diminuisce con la temperatura; il carbonato di sodio lo fa precipitare a freddo.' },
    { id: 'cu-hcl', name: 'Cu + HCl: controllo senza ossidante', steps: [[WATER, 90], [q('H+', 'Cl-', 1), 10], [metal('Cu'), 0.1]], next: q('H+', 'Cl-', 1), amount: 1, thermostat: true,
      text: 'Il rame non sposta quantità macroscopiche di H₂ da HCl diluito privo di ossidanti. Il modello non include corrosione da ossigeno atmosferico.' },
    { id: 'cu-carbonate', name: 'CuCl₂ + bicarbonato: carbonato basico', steps: [[WATER, 80], [q('Cu+2', 'Cl-', 0.1), 10]], next: q('Na+', 'HCO3-', 0.5), amount: 5, thermostat: true,
      text: 'Aggiungi NaHCO₃: confronta le fasi nel database, poi apri il reticolo della malachite Cu₂CO₃(OH)₂. Il prodotto basico non è un bicarbonato. Nucleazione e precursori amorfi non sono simulati.' },
    { id: 'malachite-acid', name: 'Dissoluzione della malachite', steps: [[WATER, 100], [mineralRecipe(MINERAL_SOLIDS.find(m=>m.phase==='Malachite')), 0.1]], next: q('H+', 'Cl-', 1), amount: 1, thermostat: true,
      text: 'Aggiungi acido: protonazione del carbonato e dissoluzione. Il rame resta Cu(II): questa trasformazione non è una redox.' },
    { id: 'barium', name: 'Ba²⁺ + SO₄²⁻: precipitazione selettiva', steps: [[WATER, 80], [q('Ba+2', 'Cl-', 0.1), 10]], next: q('Na+', 'SO4-2', 0.1), amount: 5,
      text: 'Aggiungi solfato: le quantità delle fasi solide seguono i prodotti delle attività e i bilanci di materia.' },
    { id: 'cementation', name: 'Zn + Cu²⁺: deposizione di rame', steps: [[WATER, 80], [q('Cu+2', 'Cl-', 0.1), 10]], next: metal('Zn'), amount: 0.1,
      text: 'Zn + Cu²⁺ → Zn²⁺ + Cu. Osserva il deposito metallico nel bilancio redox. Morfologia e velocità non sono previste.' },
    { id: 'empty', name: 'Becher vuoto', steps: [], text: 'Componi qualunque soluzione con il catione e l\'anione che vuoi, alla concentrazione che vuoi.' },
  ];
}

