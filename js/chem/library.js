// Libreria di molecole con dati sperimentali di confronto
// (NIST CCCBDB, NIST Chemistry WebBook, CRC Handbook). Distanze in Å, angoli in gradi,
// momenti di dipolo in debye, frequenze fondamentali in cm⁻¹.

export const CATEGORIES = [
  { id: 'diatomic', label: 'Biatomiche' },
  { id: 'hydride', label: 'Idruri' },
  { id: 'multiple', label: 'Legami multipli' },
  { id: 'organic', label: 'Organiche' },
  { id: 'aromatic', label: 'Aromatiche' },
  { id: 'vsepr', label: 'VSEPR e ipervalenti' },
  { id: 'oxide', label: 'Ossidi e ossoacidi' },
  { id: 'ion', label: 'Ioni' },
  { id: 'radical', label: 'Radicali' },
];

export const MOLECULES = [
  // --- biatomiche ---
  { id: 'H2', name: 'Idrogeno', smiles: '[H][H]', cat: 'diatomic', exp: { r: { 'H-H': 0.741 }, freq: [4161] } },
  { id: 'N2', name: 'Azoto', smiles: 'N#N', cat: 'diatomic', exp: { r: { 'N-N': 1.098 }, freq: [2330] } },
  { id: 'O2', name: 'Ossigeno (tripletto)', smiles: 'O=O', multiplicity: 3, cat: 'diatomic', exp: { r: { 'O-O': 1.208 }, freq: [1556] },
    note: 'La struttura di Lewis prevede tutti gli elettroni appaiati, ma O₂ è paramagnetico: lo stato fondamentale è un tripletto (³Σg⁻) con due elettroni spaiati negli orbitali π*. Lo spiega solo la teoria degli orbitali molecolari.' },
  { id: 'F2', name: 'Fluoro', smiles: 'FF', cat: 'diatomic', exp: { r: { 'F-F': 1.412 }, freq: [892] } },
  { id: 'Cl2', name: 'Cloro', smiles: 'ClCl', cat: 'diatomic', exp: { r: { 'Cl-Cl': 1.988 }, freq: [560] } },
  { id: 'HF', name: 'Acido fluoridrico', smiles: 'F', cat: 'diatomic', exp: { r: { 'H-F': 0.917 }, dipole: 1.83, freq: [3962] } },
  { id: 'HCl', name: 'Acido cloridrico', smiles: 'Cl', cat: 'diatomic', exp: { r: { 'H-Cl': 1.275 }, dipole: 1.11, freq: [2886] } },
  { id: 'CO', name: 'Monossido di carbonio', smiles: '[C-]#[O+]', cat: 'diatomic', exp: { r: { 'C-O': 1.128 }, dipole: 0.11, freq: [2143] },
    note: 'Il dipolo sperimentale di CO è piccolissimo e con l\'estremo negativo sul carbonio, contro l\'intuizione basata sull\'elettronegatività. Hartree–Fock sbaglia il verso: serve la correlazione elettronica.' },
  { id: 'LiH', name: 'Idruro di litio', smiles: '[LiH]', cat: 'diatomic', exp: { r: { 'H-Li': 1.595 }, dipole: 5.88, freq: [1360] } },
  { id: 'LiF', name: 'Fluoruro di litio', smiles: '[Li]F', cat: 'diatomic', exp: { r: { 'F-Li': 1.564 }, dipole: 6.33, freq: [898] } },
  { id: 'NaCl', name: 'Cloruro di sodio (gas)', smiles: '[Na]Cl', cat: 'diatomic', exp: { r: { 'Cl-Na': 2.361 }, dipole: 9.00, freq: [366] },
    note: 'Legame ionico: Δχ = 2,23. Il trasferimento di carica dal sodio al cloro si vede nelle cariche parziali e nel grande momento di dipolo.' },

  // --- idruri ---
  { id: 'H2O', name: 'Acqua', smiles: 'O', cat: 'hydride', exp: { r: { 'H-O': 0.958 }, angle: { 'H-O-H': 104.5 }, dipole: 1.85, freq: [1595, 3657, 3756] } },
  { id: 'NH3', name: 'Ammoniaca', smiles: 'N', cat: 'hydride', exp: { r: { 'H-N': 1.012 }, angle: { 'H-N-H': 106.7 }, dipole: 1.47, freq: [950, 1627, 3337, 3444] } },
  { id: 'CH4', name: 'Metano', smiles: 'C', cat: 'hydride', exp: { r: { 'C-H': 1.087 }, angle: { 'H-C-H': 109.47 }, dipole: 0, freq: [1306, 1534, 2917, 3019] } },
  { id: 'H2S', name: 'Solfuro di idrogeno', smiles: 'S', cat: 'hydride', exp: { r: { 'H-S': 1.336 }, angle: { 'H-S-H': 92.1 }, dipole: 0.97, freq: [1183, 2615, 2626] } },
  { id: 'PH3', name: 'Fosfina', smiles: 'P', cat: 'hydride', exp: { r: { 'H-P': 1.420 }, angle: { 'H-P-H': 93.3 }, dipole: 0.57 } },
  { id: 'SiH4', name: 'Silano', smiles: '[SiH4]', cat: 'hydride', exp: { r: { 'H-Si': 1.480 }, dipole: 0 } },
  { id: 'BH3', name: 'Borano', smiles: '[BH3]', cat: 'hydride', exp: { r: { 'B-H': 1.190 }, angle: { 'H-B-H': 120 } },
    note: 'Il boro ha solo 6 elettroni di valenza attorno: ottetto incompleto. Per questo BH₃ è un acido di Lewis e in natura dimerizza a diborano (B₂H₆).' },
  { id: 'BeH2', name: 'Idruro di berillio', smiles: '[BeH2]', cat: 'hydride', exp: { r: { 'Be-H': 1.326 }, angle: { 'H-Be-H': 180 } } },

  // --- legami multipli ---
  { id: 'C2H6', formula: 'C₂H₆', name: 'Etano', smiles: 'CC', cat: 'multiple', exp: { r: { 'C-C': 1.535, 'C-H': 1.094 }, dipole: 0 } },
  { id: 'C2H4', formula: 'C₂H₄', name: 'Etilene (etene)', smiles: 'C=C', cat: 'multiple', exp: { r: { 'C-C': 1.339, 'C-H': 1.086 }, angle: { 'H-C-H': 117.4 }, dipole: 0 } },
  { id: 'C2H2', formula: 'C₂H₂', name: 'Acetilene (etino)', smiles: 'C#C', cat: 'multiple', exp: { r: { 'C-C': 1.203, 'C-H': 1.063 }, dipole: 0 } },
  { id: 'HCN', formula: 'HCN', name: 'Acido cianidrico', smiles: 'C#N', cat: 'multiple', exp: { r: { 'C-H': 1.066, 'C-N': 1.153 }, dipole: 2.98, freq: [712, 2097, 3311] } },
  { id: 'HNC', formula: 'HNC', name: 'Isocianuro di idrogeno', smiles: '[C-]#[NH+]', cat: 'multiple', exp: { r: { 'H-N': 0.995, 'C-N': 1.169 }, dipole: 3.05 } },
  { id: 'CO2', name: 'Anidride carbonica', smiles: 'O=C=O', cat: 'multiple', exp: { r: { 'C-O': 1.162 }, angle: { 'O-C-O': 180 }, dipole: 0, freq: [667, 1333, 2349] } },
  { id: 'H2CO', formula: 'H₂CO', name: 'Formaldeide', smiles: 'C=O', cat: 'multiple', exp: { r: { 'C-O': 1.208, 'C-H': 1.116 }, angle: { 'H-C-H': 116.5 }, dipole: 2.33 } },
  { id: 'N2H4', name: 'Idrazina', smiles: 'NN', cat: 'multiple', exp: { r: { 'N-N': 1.447 }, dipole: 1.75 } },
  { id: 'N2H2', name: 'Diazene (trans)', smiles: 'N=N', cat: 'multiple', exp: { r: { 'N-N': 1.252 } } },

  // --- organiche ---
  { id: 'CH3OH', formula: 'CH₃OH', name: 'Metanolo', smiles: 'CO', cat: 'organic', exp: { r: { 'C-O': 1.427 }, dipole: 1.70 } },
  { id: 'C2H5OH', formula: 'C₂H₅OH', name: 'Etanolo', smiles: 'CCO', cat: 'organic', exp: { r: { 'C-O': 1.431 }, dipole: 1.69 } },
  { id: 'CH3NH2', formula: 'CH₃NH₂', name: 'Metilammina', smiles: 'CN', cat: 'organic', exp: { r: { 'C-N': 1.471 }, dipole: 1.31 } },
  { id: 'HCOOH', formula: 'HCOOH', name: 'Acido formico', smiles: 'OC=O', cat: 'organic', exp: { r: { 'C-O': 1.202 }, dipole: 1.41 } },
  { id: 'CH3COOH', formula: 'CH₃COOH', name: 'Acido acetico', smiles: 'CC(=O)O', cat: 'organic', exp: { dipole: 1.74 } },
  { id: 'acetone', formula: '(CH₃)₂CO', name: 'Acetone', smiles: 'CC(=O)C', cat: 'organic', exp: { r: { 'C-O': 1.213 }, dipole: 2.88 } },
  { id: 'urea', formula: 'CO(NH₂)₂', name: 'Urea', smiles: 'NC(=O)N', cat: 'organic', exp: { dipole: 3.83 } },
  { id: 'CH3Cl', formula: 'CH₃Cl', name: 'Clorometano', smiles: 'CCl', cat: 'organic', exp: { r: { 'C-Cl': 1.776 }, dipole: 1.87 } },
  { id: 'CH3CN', formula: 'CH₃CN', name: 'Acetonitrile', smiles: 'CC#N', cat: 'organic', exp: { r: { 'C-N': 1.157 }, dipole: 3.92 } },
  { id: 'propene', formula: 'CH₃CH=CH₂', name: 'Propene', smiles: 'CC=C', cat: 'organic', exp: { dipole: 0.37 } },
  { id: 'glycine', formula: 'NH₂CH₂COOH', name: 'Glicina', smiles: 'NCC(=O)O', cat: 'organic', exp: {} },
  { id: 'cyclopropane', name: 'Ciclopropano', smiles: 'C1CC1', cat: 'organic', exp: { r: { 'C-C': 1.510 } },
    note: 'Angoli C–C–C di 60°, molto lontani dai 109,5° dell\'ibridazione sp³: tensione d\'anello (legami "a banana").' },

  // --- aromatiche ---
  { id: 'benzene', name: 'Benzene', smiles: 'c1ccccc1', cat: 'aromatic', exp: { r: { 'C-C': 1.397, 'C-H': 1.084 }, dipole: 0 },
    note: 'Due strutture di Kekulé equivalenti: tutti i legami C–C hanno ordine 1,5 e la stessa lunghezza. Gli elettroni π sono delocalizzati sull\'anello.' },
  { id: 'pyridine', name: 'Piridina', smiles: 'c1ccncc1', cat: 'aromatic', exp: { dipole: 2.22 } },
  { id: 'pyrrole', name: 'Pirrolo', smiles: 'c1cc[nH]c1', cat: 'aromatic', exp: { dipole: 1.74 } },
  { id: 'furan', name: 'Furano', smiles: 'c1ccoc1', cat: 'aromatic', exp: { dipole: 0.66 } },

  // --- VSEPR e ipervalenti ---
  { id: 'BF3', name: 'Trifluoruro di boro', smiles: 'FB(F)F', cat: 'vsepr', exp: { r: { 'B-F': 1.307 }, angle: { 'F-B-F': 120 }, dipole: 0 } },
  { id: 'CF4', name: 'Tetrafluorometano', smiles: 'FC(F)(F)F', cat: 'vsepr', exp: { r: { 'C-F': 1.315 }, dipole: 0 } },
  { id: 'PF5', name: 'Pentafluoruro di fosforo', smiles: 'FP(F)(F)(F)F', cat: 'vsepr', exp: { dipole: 0 } },
  { id: 'SF4', name: 'Tetrafluoruro di zolfo', smiles: 'FS(F)(F)F', cat: 'vsepr', exp: { dipole: 0.63 } },
  { id: 'ClF3', name: 'Trifluoruro di cloro', smiles: 'FCl(F)F', cat: 'vsepr', exp: { dipole: 0.6 } },
  { id: 'SF6', name: 'Esafluoruro di zolfo', smiles: 'FS(F)(F)(F)(F)F', cat: 'vsepr', exp: { r: { 'F-S': 1.561 }, dipole: 0 } },
  { id: 'BrF5', name: 'Pentafluoruro di bromo', smiles: 'FBr(F)(F)(F)F', cat: 'vsepr', exp: { dipole: 1.51 } },
  { id: 'KrF2', name: 'Difluoruro di kripton', smiles: 'F[Kr]F', cat: 'vsepr', exp: { r: { 'F-Kr': 1.875 }, dipole: 0 },
    note: 'Un composto di un gas nobile: AX₂E₃, lineare. Le tre coppie solitarie stanno sul piano equatoriale della bipiramide trigonale.' },

  // --- ossidi e ossoacidi ---
  { id: 'O3', name: 'Ozono', smiles: '[O-][O+]=O', cat: 'oxide', exp: { r: { 'O-O': 1.278 }, angle: { 'O-O-O': 116.8 }, dipole: 0.53 } },
  { id: 'SO2', name: 'Anidride solforosa', smiles: 'O=S=O', cat: 'oxide', exp: { r: { 'O-S': 1.431 }, angle: { 'O-S-O': 119.3 }, dipole: 1.63 } },
  { id: 'SO3', name: 'Anidride solforica', smiles: 'O=S(=O)=O', cat: 'oxide', exp: { r: { 'O-S': 1.420 }, dipole: 0 } },
  { id: 'H2O2', name: 'Perossido di idrogeno', smiles: 'OO', cat: 'oxide', exp: { r: { 'O-O': 1.475 }, dipole: 1.57 },
    note: 'Nel perossido ogni ossigeno ha numero di ossidazione −1: il legame O–O è tra atomi uguali e i suoi elettroni si dividono a metà.' },
  { id: 'N2O4', name: 'Tetrossido di diazoto', smiles: '[O-][N+](=O)[N+]([O-])=O', cat: 'oxide', exp: { r: { 'N-N': 1.782 } } },
  { id: 'HNO3', formula: 'HNO₃', name: 'Acido nitrico', smiles: 'O[N+](=O)[O-]', cat: 'oxide', exp: { dipole: 2.17 } },
  { id: 'H2SO4', name: 'Acido solforico', smiles: 'OS(=O)(=O)O', cat: 'oxide', exp: { dipole: 2.73 } },

  // --- ioni ---
  { id: 'Hplus', formula: 'H⁺', name: 'Protone', smiles: '[H+]', cat: 'ion', exp: {} },
  { id: 'NH4', formula: 'NH₄⁺', name: 'Ione ammonio', smiles: '[NH4+]', cat: 'ion', exp: {} },
  { id: 'H3O', formula: 'H₃O⁺', name: 'Ione ossonio', smiles: '[OH3+]', cat: 'ion', exp: {} },
  { id: 'OH-', formula: 'OH⁻', name: 'Ione idrossido', smiles: '[OH-]', cat: 'ion', exp: {} },
  { id: 'CN-', formula: 'CN⁻', name: 'Ione cianuro', smiles: '[C-]#N', cat: 'ion', exp: {} },
  { id: 'NO3-', formula: 'NO₃⁻', name: 'Ione nitrato', smiles: '[O-][N+](=O)[O-]', cat: 'ion', exp: {},
    note: 'Tre strutture di risonanza equivalenti: i tre legami N–O hanno ordine medio 4/3 e sono identici.' },
  { id: 'CO3', formula: 'CO₃²⁻', name: 'Ione carbonato', smiles: 'C(=O)([O-])[O-]', cat: 'ion', exp: {},
    note: 'Nel vuoto il dianione è instabile rispetto alla perdita di un elettrone: esiste stabilizzato da acqua o controioni. Il calcolo in fase gassosa mostra comunque la struttura con tre legami equivalenti.' },
  { id: 'SO4', formula: 'SO₄²⁻', name: 'Ione solfato', smiles: '[O-]S(=O)(=O)[O-]', cat: 'ion', exp: {} },

  // --- radicali ---
  { id: 'NO', name: 'Monossido di azoto', smiles: '[N]=O', cat: 'radical', exp: { r: { 'N-O': 1.151 }, dipole: 0.16, freq: [1876] } },
  { id: 'NO2', name: 'Diossido di azoto', smiles: '[O][N]=O', cat: 'radical', exp: { r: { 'N-O': 1.194 }, angle: { 'O-N-O': 134.1 }, dipole: 0.32 } },
  { id: 'CH3', formula: 'CH₃', name: 'Radicale metile', smiles: '[CH3]', cat: 'radical', exp: { r: { 'C-H': 1.079 }, angle: { 'H-C-H': 120 } } },
  { id: 'OH', formula: 'OH', name: 'Radicale ossidrile', smiles: '[OH]', cat: 'radical', exp: { r: { 'H-O': 0.970 }, dipole: 1.66 } },
];

/**
 * Reazioni in fase gassosa con entalpia standard sperimentale a 298,15 K
 * (da entalpie di formazione NIST/JANAF, kJ/mol).
 */
export const REACTIONS = [
  { id: 'water', name: 'Sintesi dell\'acqua', left: [[2, 'H2'], [1, 'O2']], right: [[2, 'H2O']], dH: -483.6 },
  { id: 'haber', name: 'Sintesi dell\'ammoniaca (Haber–Bosch)', left: [[1, 'N2'], [3, 'H2']], right: [[2, 'NH3']], dH: -91.8 },
  { id: 'methane', name: 'Combustione del metano', left: [[1, 'CH4'], [2, 'O2']], right: [[1, 'CO2'], [2, 'H2O']], dH: -802.3 },
  { id: 'hcl', name: 'Sintesi dell\'acido cloridrico', left: [[1, 'H2'], [1, 'Cl2']], right: [[2, 'HCl']], dH: -184.6 },
  { id: 'hf', name: 'Sintesi dell\'acido fluoridrico', left: [[1, 'H2'], [1, 'F2']], right: [[2, 'HF']], dH: -546.6 },
  { id: 'hydrog', name: 'Idrogenazione dell\'etilene', left: [[1, 'C2H4'], [1, 'H2']], right: [[1, 'C2H6']], dH: -136.3 },
  { id: 'hydrog2', name: 'Idrogenazione dell\'acetilene', left: [[1, 'C2H2'], [2, 'H2']], right: [[1, 'C2H6']], dH: -311.4 },
  { id: 'wgs', name: 'Reazione del gas d\'acqua (water-gas shift)', left: [[1, 'CO'], [1, 'H2O']], right: [[1, 'CO2'], [1, 'H2']], dH: -41.2 },
  { id: 'n2o4', name: 'Dimerizzazione di NO₂', left: [[2, 'NO2']], right: [[1, 'N2O4']], dH: -57.2 },
  { id: 'ozone', name: 'Decomposizione dell\'ozono', left: [[2, 'O3']], right: [[3, 'O2']], dH: -285.4 },
  { id: 'hcn', name: 'Isomerizzazione HCN → HNC', left: [[1, 'HCN']], right: [[1, 'HNC']], dH: 62 },
  { id: 'pa_nh3', name: 'Affinità protonica dell\'ammoniaca', left: [[1, 'NH3'], [1, 'Hplus']], right: [[1, 'NH4']], dH: -853.6 },
  { id: 'pa_h2o', name: 'Affinità protonica dell\'acqua', left: [[1, 'H2O'], [1, 'Hplus']], right: [[1, 'H3O']], dH: -691.0 },
  { id: 'chlorination', name: 'Clorurazione del metano', left: [[1, 'CH4'], [1, 'Cl2']], right: [[1, 'CH3Cl'], [1, 'HCl']], dH: -99.6 },
  { id: 'no_ox', name: 'Ossidazione di NO', left: [[2, 'NO'], [1, 'O2']], right: [[2, 'NO2']], dH: -114.1 },
  { id: 'so2_ox', name: 'Ossidazione di SO₂ (processo di contatto)', left: [[2, 'SO2'], [1, 'O2']], right: [[2, 'SO3']], dH: -197.8 },
];

export function moleculeById(id) {
  return MOLECULES.find(m => m.id === id);
}

/** Formula da mostrare: quella convenzionale se indicata, altrimenti null (si usa la formula bruta). */
export function displayFormula(id) {
  return moleculeById(id)?.formula ?? null;
}
