// Dati che completano il database MINTEQ v4 (js/chem/aqueousDB.js), ciascuno con la sua fonte.

// --- specie aggiuntive (stessa convenzione di aqueousDB.js) -----------------------------------------------
export const EXTRA_COMPONENTS = ['SCN-', 'HPP', 'HBTB', 'HMO'];
export const EXTRA_SPECIES = [
  // tiocianato di ferro(III): log K a forza ionica zero, Inorg. Chim. Acta 469 (2018) 164, doi:10.1016/j.ica.2017.09.004
  { n: 'SCN-', z: -1, c: { 'SCN-': 1 }, w: 0, lk: 0, dh: 0, g: null, src: null },
  { n: 'FeSCN+2', z: 2, c: { 'Fe+3': 1, 'SCN-': 1 }, w: 0, lk: 2.85, dh: null, g: null, src: 'Inorg. Chim. Acta 2018 (log K₁° = 2,85 ± 0,08)' },
  { n: 'Fe(SCN)2+', z: 1, c: { 'Fe+3': 1, 'SCN-': 2 }, w: 0, lk: 4.36, dh: null, g: null, src: 'Inorg. Chim. Acta 2018 (log K₂° = 1,51 ± 0,13)' },
  // ammino-complessi mancanti in MINTEQ v4: database LLNL (thermo.com.V8.R6+, distribuito con PHREEQC; costanti da NIST SRD 46),
  // riscritti con NH₄⁺ come componente: log K = log β − n·9,244, ΔH = ΔH(LLNL) + n·52,0 kJ/mol (NH₄⁺ → NH₃ + H⁺, MINTEQ)
  { n: 'Cu(NH3)2+2', z: 2, c: { 'Cu+2': 1, 'NH4+': 2, 'H+': -2 }, w: 0, lk: 7.4512 - 2 * 9.244, dh: -45.1269 + 2 * 52, g: [4.5, 0], src: 'LLNL (log β₂ = 7,451)' },
  { n: 'Cu(NH3)3+2', z: 2, c: { 'Cu+2': 1, 'NH4+': 3, 'H+': -3 }, w: 0, lk: 10.2719 - 3 * 9.244, dh: -67.2779 + 3 * 52, g: [4.5, 0], src: 'LLNL (log β₃ = 10,272)' },
  // β₄: LibreTexts, Reference Tables E4 (K_f = 1,1·10¹³); ΔH non disponibile
  { n: 'Cu(NH3)4+2', z: 2, c: { 'Cu+2': 1, 'NH4+': 4, 'H+': -4 }, w: 0, lk: 13.04 - 4 * 9.244, dh: null, g: [4.5, 0], src: 'LibreTexts E4 (K_f = 1,1·10¹³)' },
  { n: 'Zn(NH3)4+2', z: 2, c: { 'Zn+2': 1, 'NH4+': 4, 'H+': -4 }, w: 0, lk: 8.3738 - 4 * 9.244, dh: -54.9027 + 4 * 52, g: [4.5, 0], src: 'LLNL (log β₄ = 8,374)' },
  { n: 'Ni(NH3)6+2', z: 2, c: { 'Ni+2': 1, 'NH4+': 6, 'H+': -6 }, w: 0, lk: 8.7344 - 6 * 9.244, dh: -88.0436 + 6 * 52, g: [4.5, 0], src: 'LLNL (log β₆ = 8,734)' },
  // indicatori acido-base: forma acida HIn (componente) e forma basica In⁻.
  // pKa: fenolftaleina 9,4; blu di bromotimolo 7,1; metilarancio 3,46 (Harris, Quantitative Chemical Analysis, tab. 11-4)
  { n: 'HPP', z: 0, c: { HPP: 1 }, w: 0, lk: 0, dh: null, g: null, src: null },
  { n: 'PP-', z: -1, c: { HPP: 1, 'H+': -1 }, w: 0, lk: -9.4, dh: null, g: null, src: 'Harris, QCA (pKa 9,4)' },
  { n: 'HBTB', z: 0, c: { HBTB: 1 }, w: 0, lk: 0, dh: null, g: null, src: null },
  { n: 'BTB-', z: -1, c: { HBTB: 1, 'H+': -1 }, w: 0, lk: -7.1, dh: null, g: null, src: 'Harris, QCA (pKa 7,1)' },
  { n: 'HMO', z: 0, c: { HMO: 1 }, w: 0, lk: 0, dh: null, g: null, src: null },
  { n: 'MO-', z: -1, c: { HMO: 1, 'H+': -1 }, w: 0, lk: -3.46, dh: null, g: null, src: 'Harris, QCA (pKa 3,46)' },
];

// metalli che MINTEQ non contiene: log K di M(s) = Mⁿ⁺ + n e⁻ dal potenziale standard, log K = −nE°/0,05916
// (E° da Bard, Parsons, Jordan, Standard Potentials in Aqueous Solution, IUPAC 1985)
export const EXTRA_PHASES = [
  { n: 'Mgmetal', f: 'Mg', c: { 'Mg+2': 1, 'e-': 2 }, w: 0, lk: +(2 * 2.372 / 0.05916).toFixed(3), dh: -466.85, src: 'Bard 1985: E°(Mg²⁺/Mg) = −2,372 V; ΔfH°(Mg²⁺) NBS' },
  { n: 'Femetal', f: 'Fe', c: { 'Fe+3': 1, 'e-': 3 }, w: 0, lk: +((2 * 0.447 - 0.771) / 0.05916).toFixed(3), dh: -48.5, src: 'Bard 1985: E°(Fe²⁺/Fe) = −0,447 V, E°(Fe³⁺/Fe²⁺) = +0,771 V' },
];

// --- nomi leggibili ---------------------------------------------------------------------------------------
export const ORGANIC_ABBR = {
  Acetate: 'Ac', Citrate: 'Cit', Edta: 'EDTA', Formate: 'Form', Phthalate: 'Ftal', Glycine: 'Gly',
  Tartarate: 'Tart', Benzoate: 'Benz', Ethylenediamine: 'en', Salicylate: 'Sal',
};
export const ORGANIC_LEGEND = {
  Ac: 'acetato CH₃COO⁻', Cit: 'citrato C₆H₅O₇³⁻', EDTA: 'etilendiamminotetraacetato', Form: 'formiato HCOO⁻',
  Ftal: 'ftalato C₆H₄(COO)₂²⁻', Gly: 'glicinato NH₂CH₂COO⁻', Tart: 'tartrato', Benz: 'benzoato C₆H₅COO⁻',
  en: 'etilendiammina NH₂CH₂CH₂NH₂', Sal: 'salicilato',
};

// fasi: nome italiano e colore osservato del solido (CRC Handbook, "Physical Constants of Inorganic Compounds")
export const PHASE_INFO = {
  Cerargyrite: ['AgCl', 'cloruro d\'argento', '#f2f2ee'], Bromyrite: ['AgBr', 'bromuro d\'argento', '#eee6bf'],
  Iodyrite: ['AgI', 'ioduro d\'argento', '#f0dc6c'], Ag2CrO4: ['Ag₂CrO₄', 'cromato d\'argento', '#9c2f25'],
  Ag2CO3: ['Ag₂CO₃', 'carbonato d\'argento', '#ece5b8'], Ag3PO4: ['Ag₃PO₄', 'fosfato d\'argento', '#e8d44d'],
  Ag2O: ['Ag₂O', 'ossido d\'argento', '#3b2e22'], Ag2SO4: ['Ag₂SO₄', 'solfato d\'argento', '#f5f5f5'],
  Acanthite: ['Ag₂S', 'solfuro d\'argento', '#1c1c1c'],
  Cotunnite: ['PbCl₂', 'cloruro di piombo', '#f5f5f5'], PbI2: ['PbI₂', 'ioduro di piombo', '#ffd21f'],
  PbBr2: ['PbBr₂', 'bromuro di piombo', '#f5f5f5'], Anglesite: ['PbSO₄', 'solfato di piombo', '#f5f5f5'],
  PbCrO4: ['PbCrO₄', 'cromato di piombo (giallo di cromo)', '#f7c21a'], Cerussite: ['PbCO₃', 'carbonato di piombo', '#f5f5f5'],
  'Pb(OH)2': ['Pb(OH)₂', 'idrossido di piombo', '#f0f0f0'], Galena: ['PbS', 'solfuro di piombo', '#222222'],
  PbF2: ['PbF₂', 'fluoruro di piombo', '#f5f5f5'], 'Pb3(PO4)2': ['Pb₃(PO₄)₂', 'fosfato di piombo', '#f5f5f5'],
  Barite: ['BaSO₄', 'solfato di bario', '#fbfbfb'], Witherite: ['BaCO₃', 'carbonato di bario', '#f8f8f8'],
  BaCrO4: ['BaCrO₄', 'cromato di bario', '#f3dd3c'], BaF2: ['BaF₂', 'fluoruro di bario', '#f8f8f8'], BaHPO4: ['BaHPO₄', 'idrogenofosfato di bario', '#f8f8f8'],
  Calcite: ['CaCO₃', 'carbonato di calcio', '#fafafa'], Gypsum: ['CaSO₄·2H₂O', 'gesso', '#fafafa'],
  Fluorite: ['CaF₂', 'fluoruro di calcio', '#f6f6f6'], Portlandite: ['Ca(OH)₂', 'idrossido di calcio', '#fafafa'],
  'Ca3(PO4)2(beta)': ['Ca₃(PO₄)₂', 'fosfato di calcio', '#fafafa'], CaCrO4: ['CaCrO₄', 'cromato di calcio', '#f3e04a'],
  Celestite: ['SrSO₄', 'solfato di stronzio', '#fafafa'], Strontianite: ['SrCO₃', 'carbonato di stronzio', '#fafafa'], SrF2: ['SrF₂', 'fluoruro di stronzio', '#fafafa'],
  Brucite: ['Mg(OH)₂', 'idrossido di magnesio', '#fafafa'], Nesquehonite: ['MgCO₃·3H₂O', 'carbonato di magnesio', '#fafafa'], MgF2: ['MgF₂', 'fluoruro di magnesio', '#fafafa'],
  'Al(OH)3(am)': ['Al(OH)₃', 'idrossido di alluminio (gelatinoso)', '#eeeeee'],
  Ferrihydrite: ['Fe(OH)₃', 'idrossido di ferro(III)', '#9b4a1c'], 'Fe(OH)2': ['Fe(OH)₂', 'idrossido di ferro(II)', '#8fa36a'],
  Siderite: ['FeCO₃', 'carbonato di ferro(II)', '#c9c2a4'], 'FeS(ppt)': ['FeS', 'solfuro di ferro(II)', '#1e1e1e'], Strengite: ['FePO₄·2H₂O', 'fosfato di ferro(III)', '#efe6d0'],
  'Cu(OH)2': ['Cu(OH)₂', 'idrossido di rame(II)', '#5e9fd8'], Malachite: ['Cu₂(OH)₂CO₃', 'carbonato basico di rame', '#3d9a73'],
  Brochantite: ['Cu₄(OH)₆SO₄', 'solfato basico di rame', '#6fbfa8'], Covellite: ['CuS', 'solfuro di rame(II)', '#1f2230'],
  Nantokite: ['CuCl', 'cloruro di rame(I)', '#f2f2f2'], CuI: ['CuI', 'ioduro di rame(I)', '#efe9e0'], 'Cu3(PO4)2': ['Cu₃(PO₄)₂', 'fosfato di rame', '#7fb6d6'],
  Cuprite: ['Cu₂O', 'ossido di rame(I)', '#b23a1e'],
  'Zn(OH)2(am)': ['Zn(OH)₂', 'idrossido di zinco', '#f4f4f4'], ZnCO3: ['ZnCO₃', 'carbonato di zinco', '#f4f4f4'],
  'ZnS(am)': ['ZnS', 'solfuro di zinco', '#f7f7f0'], 'Zn3(PO4)2': ['Zn₃(PO₄)₂', 'fosfato di zinco', '#f4f4f4'],
  'Ni(OH)2': ['Ni(OH)₂', 'idrossido di nichel', '#79b86d'], NiCO3: ['NiCO₃', 'carbonato di nichel', '#8fcf7f'],
  'NiS(alpha)': ['NiS', 'solfuro di nichel', '#1c1c1c'], 'Ni3(PO4)2': ['Ni₃(PO₄)₂', 'fosfato di nichel', '#a3d99a'],
  'Co(OH)2': ['Co(OH)₂', 'idrossido di cobalto', '#d98fa0'], CoCO3: ['CoCO₃', 'carbonato di cobalto', '#d98fb0'],
  'CoS(alpha)': ['CoS', 'solfuro di cobalto', '#1c1c1c'], 'Co3(PO4)2': ['Co₃(PO₄)₂', 'fosfato di cobalto', '#b76ea8'],
  'Cd(OH)2(am)': ['Cd(OH)₂', 'idrossido di cadmio', '#f5f5f5'], Otavite: ['CdCO₃', 'carbonato di cadmio', '#f5f5f5'], Greenockite: ['CdS', 'solfuro di cadmio', '#f2b705'],
  Pyrochroite: ['Mn(OH)₂', 'idrossido di manganese(II)', '#efe6e0'], Rhodochrosite: ['MnCO₃', 'carbonato di manganese', '#f0cfd0'],
  'MnS(pnk)': ['MnS', 'solfuro di manganese', '#e8a9a0'], Pyrolusite: ['MnO₂', 'diossido di manganese', '#2b2522'], Birnessite: ['MnO₂ (birnessite)', 'ossido di manganese', '#3a2d24'],
  'Cr(OH)3(am)': ['Cr(OH)₃', 'idrossido di cromo(III)', '#5f8f5a'],
  Agmetal: ['Ag', 'argento metallico', '#b9bec6'], Cumetal: ['Cu', 'rame metallico', '#c46f3c'], Pbmetal: ['Pb', 'piombo metallico', '#5d6168'],
  Znmetal: ['Zn', 'zinco metallico', '#9aa3ad'], 'Cdmetal(alpha)': ['Cd', 'cadmio metallico', '#8a929a'],
  Mgmetal: ['Mg', 'magnesio metallico', '#c9ccd1'], Femetal: ['Fe', 'ferro metallico', '#6d6f72'],
  'CO2(g)': ['CO₂(g)', 'anidride carbonica gassosa', '#ffffff'],
};

// --- spettri di assorbimento: [λ_max nm, ε_max L mol⁻¹ cm⁻¹, larghezza a metà altezza cm⁻¹] ---------------------
// Valori da letteratura: FeSCN²⁺ (Frank & Oswalt, JACS 1947); MnO₄⁻, CrO₄²⁻, Cr₂O₇²⁻ (Skoog, Principles of Instrumental
// Analysis); ioni dei metalli di transizione (Lever, Inorganic Electronic Spectroscopy). Larghezze tipiche delle bande d–d
// e di trasferimento di carica: servono solo a calcolare il colore, non entrano nell'equilibrio.
export const SPECTRA = {
  'Cu+2': [[810, 12, 6500]],
  'CuOH+': [[780, 15, 6500]],
  'CuNH3+2': [[745, 20, 6500]], 'Cu(NH3)2+2': [[680, 30, 6500]], 'Cu(NH3)3+2': [[630, 42, 6000]], 'Cu(NH3)4+2': [[600, 52, 5800]],
  'Cu(Edta)-2': [[735, 90, 6000]], 'CuEdta-2': [[735, 90, 6000]],
  'Ni+2': [[395, 5.2, 3500], [720, 2.1, 3500]],
  'Ni(NH3)2+2': [[380, 5, 3500], [640, 3, 3500]], 'Ni(NH3)6+2': [[355, 6, 3500], [570, 6, 3500]],
  'Co+2': [[510, 4.8, 3000]],
  'Fe+3': [[240, 4200, 9000]], 'FeOH+2': [[297, 2000, 9000]], 'Fe(OH)2+': [[310, 1500, 9000]],
  'FeCl+2': [[336, 2100, 8000]], 'FeCl2+': [[350, 2500, 8000]],
  'FeSCN+2': [[447, 4700, 6500]], 'Fe(SCN)2+': [[460, 5000, 6500]],
  'Fe+2': [[960, 1.1, 5000]],
  'MnO4-': [[525, 2400, 1700], [545, 2400, 1700], [311, 1800, 3000]],
  'MnO4-2': [[610, 1700, 3500], [440, 1000, 3500]],
  'Mn+2': [[400, 0.03, 2000]],
  'CrO4-2': [[372, 4800, 4200]], 'HCrO4-': [[350, 1600, 4500], [440, 150, 4000]],
  'Cr2O7-2': [[350, 3000, 4500], [445, 370, 5500]],
  'Cr+3': [[575, 13, 3500], [410, 15, 3500]], 'CrOH+2': [[580, 15, 3500], [415, 17, 3500]], 'Cr(OH)2+': [[585, 16, 3500], [420, 18, 3500]],
  'PP-': [[552, 30000, 1900], [374, 9000, 5000]],
  HBTB: [[432, 16000, 4200]], 'BTB-': [[617, 37000, 3500]],
  HMO: [[507, 40000, 3200]], 'MO-': [[464, 27000, 4200]],
};

// entalpie di soluzione a diluizione infinita (kJ/mol) per i sali solidi che MINTEQ non tratta come fasi
// (CRC Handbook, "Enthalpy of Solution of Electrolytes")
// (CRC Handbook, 90ª ed., "Enthalpy of Solution of Electrolytes"). Per gli altri sali solidi il calore di dissoluzione
// non è noto al programma e viene dichiarato come tale.
export const SOLUTION_ENTHALPY = { NaCl: 3.87, NH4NO3: 25.69, NaOH: -44.50, KOH: -57.61 };
