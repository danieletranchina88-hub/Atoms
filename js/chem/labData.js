// Dati sperimentali per il laboratorio di chimica fisica (CRC Handbook of Chemistry and Physics).

// Costanti di van der Waals: a in L²·bar/mol², b in L/mol.
export const VDW_GASES = [
  { id: 'He', name: 'Elio', formula: 'He', a: 0.0346, b: 0.0238 },
  { id: 'H2', name: 'Idrogeno', formula: 'H₂', a: 0.2476, b: 0.02661 },
  { id: 'N2', name: 'Azoto', formula: 'N₂', a: 1.370, b: 0.0387 },
  { id: 'O2', name: 'Ossigeno', formula: 'O₂', a: 1.382, b: 0.03186 },
  { id: 'Ar', name: 'Argon', formula: 'Ar', a: 1.355, b: 0.03201 },
  { id: 'CH4', name: 'Metano', formula: 'CH₄', a: 2.283, b: 0.04278 },
  { id: 'CO2', name: 'Anidride carbonica', formula: 'CO₂', a: 3.640, b: 0.04267 },
  { id: 'NH3', name: 'Ammoniaca', formula: 'NH₃', a: 4.225, b: 0.0371 },
  { id: 'H2O', name: 'Acqua', formula: 'H₂O', a: 5.536, b: 0.03049 },
  { id: 'Cl2', name: 'Cloro', formula: 'Cl₂', a: 6.579, b: 0.05622 },
];

// Acidi (pKa successive a 25 °C). Le basi deboli sono descritte dal loro acido coniugato.
export const ACIDS = [
  { id: 'HCl', name: 'Acido cloridrico (forte)', formula: 'HCl', pKa: [-7] },
  { id: 'acetic', name: 'Acido acetico', formula: 'CH₃COOH', pKa: [4.76] },
  { id: 'formic', name: 'Acido formico', formula: 'HCOOH', pKa: [3.75] },
  { id: 'HF', name: 'Acido fluoridrico', formula: 'HF', pKa: [3.17] },
  { id: 'HCN', name: 'Acido cianidrico', formula: 'HCN', pKa: [9.21] },
  { id: 'HClO', name: 'Acido ipocloroso', formula: 'HClO', pKa: [7.53] },
  { id: 'benzoic', name: 'Acido benzoico', formula: 'C₆H₅COOH', pKa: [4.20] },
  { id: 'carbonic', name: 'Acido carbonico', formula: 'H₂CO₃', pKa: [6.35, 10.33] },
  { id: 'oxalic', name: 'Acido ossalico', formula: 'H₂C₂O₄', pKa: [1.25, 4.27] },
  { id: 'phosphoric', name: 'Acido fosforico', formula: 'H₃PO₄', pKa: [2.15, 7.20, 12.35] },
  { id: 'citric', name: 'Acido citrico', formula: 'C₆H₈O₇', pKa: [3.13, 4.76, 6.40] },
];
export const BASES = [
  { id: 'NaOH', name: 'Idrossido di sodio (forte)', formula: 'NaOH', pKaConj: 20 },
  { id: 'NH3', name: 'Ammoniaca', formula: 'NH₃', pKaConj: 9.25 },
  { id: 'methylamine', name: 'Metilammina', formula: 'CH₃NH₂', pKaConj: 10.64 },
  { id: 'pyridine', name: 'Piridina', formula: 'C₅H₅N', pKaConj: 5.23 },
];

// Indicatori acido-base: intervallo di viraggio.
export const INDICATORS = [
  { name: 'Metilarancio', from: 3.1, to: 4.4, colors: ['#e0402a', '#f2c22e'] },
  { name: 'Rosso metile', from: 4.4, to: 6.2, colors: ['#e0402a', '#f2d22e'] },
  { name: 'Blu di bromotimolo', from: 6.0, to: 7.6, colors: ['#e8d02a', '#2f6fd6'] },
  { name: 'Fenolftaleina', from: 8.2, to: 10.0, colors: ['#f5f5f5', '#d93a9c'] },
];

// Potenziali standard di riduzione a 25 °C (V rispetto all'elettrodo standard a idrogeno).
export const HALF_REACTIONS = [
  { id: 'Li', eq: 'Li⁺ + e⁻ → Li', ox: 'Li⁺', red: 'Li', n: 1, E: -3.04, metal: true },
  { id: 'K', eq: 'K⁺ + e⁻ → K', ox: 'K⁺', red: 'K', n: 1, E: -2.93, metal: true },
  { id: 'Ca', eq: 'Ca²⁺ + 2e⁻ → Ca', ox: 'Ca²⁺', red: 'Ca', n: 2, E: -2.87, metal: true },
  { id: 'Na', eq: 'Na⁺ + e⁻ → Na', ox: 'Na⁺', red: 'Na', n: 1, E: -2.71, metal: true },
  { id: 'Mg', eq: 'Mg²⁺ + 2e⁻ → Mg', ox: 'Mg²⁺', red: 'Mg', n: 2, E: -2.37, metal: true },
  { id: 'Al', eq: 'Al³⁺ + 3e⁻ → Al', ox: 'Al³⁺', red: 'Al', n: 3, E: -1.66, metal: true },
  { id: 'Zn', eq: 'Zn²⁺ + 2e⁻ → Zn', ox: 'Zn²⁺', red: 'Zn', n: 2, E: -0.76, metal: true },
  { id: 'Fe', eq: 'Fe²⁺ + 2e⁻ → Fe', ox: 'Fe²⁺', red: 'Fe', n: 2, E: -0.44, metal: true },
  { id: 'Ni', eq: 'Ni²⁺ + 2e⁻ → Ni', ox: 'Ni²⁺', red: 'Ni', n: 2, E: -0.25, metal: true },
  { id: 'Sn', eq: 'Sn²⁺ + 2e⁻ → Sn', ox: 'Sn²⁺', red: 'Sn', n: 2, E: -0.14, metal: true },
  { id: 'Pb', eq: 'Pb²⁺ + 2e⁻ → Pb', ox: 'Pb²⁺', red: 'Pb', n: 2, E: -0.13, metal: true },
  { id: 'H2', eq: '2H⁺ + 2e⁻ → H₂', ox: 'H⁺', red: 'H₂', n: 2, E: 0.00, metal: false },
  { id: 'Cu', eq: 'Cu²⁺ + 2e⁻ → Cu', ox: 'Cu²⁺', red: 'Cu', n: 2, E: 0.34, metal: true },
  { id: 'I2', eq: 'I₂ + 2e⁻ → 2I⁻', ox: 'I₂', red: 'I⁻', n: 2, E: 0.54, metal: false },
  { id: 'Fe3', eq: 'Fe³⁺ + e⁻ → Fe²⁺', ox: 'Fe³⁺', red: 'Fe²⁺', n: 1, E: 0.77, metal: false },
  { id: 'Ag', eq: 'Ag⁺ + e⁻ → Ag', ox: 'Ag⁺', red: 'Ag', n: 1, E: 0.80, metal: true },
  { id: 'Br2', eq: 'Br₂ + 2e⁻ → 2Br⁻', ox: 'Br₂', red: 'Br⁻', n: 2, E: 1.07, metal: false },
  { id: 'O2', eq: 'O₂ + 4H⁺ + 4e⁻ → 2H₂O', ox: 'O₂', red: 'H₂O', n: 4, E: 1.23, metal: false },
  { id: 'Cl2', eq: 'Cl₂ + 2e⁻ → 2Cl⁻', ox: 'Cl₂', red: 'Cl⁻', n: 2, E: 1.36, metal: false },
  { id: 'Au', eq: 'Au³⁺ + 3e⁻ → Au', ox: 'Au³⁺', red: 'Au', n: 3, E: 1.50, metal: true },
  { id: 'MnO4', eq: 'MnO₄⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O', ox: 'MnO₄⁻', red: 'Mn²⁺', n: 5, E: 1.51, metal: false },
  { id: 'F2', eq: 'F₂ + 2e⁻ → 2F⁻', ox: 'F₂', red: 'F⁻', n: 2, E: 2.87, metal: false },
];

export const R_GAS = 8.314462618;          // J/(mol·K)
export const R_LBAR = 0.08314462618;       // L·bar/(mol·K)
export const FARADAY = 96485.33212;        // C/mol
