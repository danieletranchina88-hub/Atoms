// Tavola periodica: dati sperimentali di riferimento per i 118 elementi.
// Masse atomiche: IUPAC (valori abbreviati; [x] = numero di massa dell'isotopo più stabile).
// Energie di prima ionizzazione: NIST Atomic Spectra Database, in eV (null = non misurata).

const SYMBOLS = [
  'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
  'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
  'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
  'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
  'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
  'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
  'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
  'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
];

const NAMES = [
  'Idrogeno', 'Elio', 'Litio', 'Berillio', 'Boro', 'Carbonio', 'Azoto', 'Ossigeno', 'Fluoro', 'Neon',
  'Sodio', 'Magnesio', 'Alluminio', 'Silicio', 'Fosforo', 'Zolfo', 'Cloro', 'Argon', 'Potassio', 'Calcio',
  'Scandio', 'Titanio', 'Vanadio', 'Cromo', 'Manganese', 'Ferro', 'Cobalto', 'Nichel', 'Rame', 'Zinco',
  'Gallio', 'Germanio', 'Arsenico', 'Selenio', 'Bromo', 'Kripton', 'Rubidio', 'Stronzio', 'Ittrio', 'Zirconio',
  'Niobio', 'Molibdeno', 'Tecnezio', 'Rutenio', 'Rodio', 'Palladio', 'Argento', 'Cadmio', 'Indio', 'Stagno',
  'Antimonio', 'Tellurio', 'Iodio', 'Xeno', 'Cesio', 'Bario', 'Lantanio', 'Cerio', 'Praseodimio', 'Neodimio',
  'Promezio', 'Samario', 'Europio', 'Gadolinio', 'Terbio', 'Disprosio', 'Olmio', 'Erbio', 'Tulio', 'Itterbio',
  'Lutezio', 'Afnio', 'Tantalio', 'Tungsteno', 'Renio', 'Osmio', 'Iridio', 'Platino', 'Oro', 'Mercurio',
  'Tallio', 'Piombo', 'Bismuto', 'Polonio', 'Astato', 'Radon', 'Francio', 'Radio', 'Attinio', 'Torio',
  'Protoattinio', 'Uranio', 'Nettunio', 'Plutonio', 'Americio', 'Curio', 'Berkelio', 'Californio', 'Einsteinio', 'Fermio',
  'Mendelevio', 'Nobelio', 'Laurenzio', 'Rutherfordio', 'Dubnio', 'Seaborgio', 'Bohrio', 'Hassio', 'Meitnerio', 'Darmstadtio',
  'Roentgenio', 'Copernicio', 'Nihonio', 'Flerovio', 'Moscovio', 'Livermorio', 'Tennessinio', 'Oganesson',
];

// Numeri come stringhe per conservare le cifre significative; tra parentesi quadre gli elementi senza isotopi stabili.
const MASSES = [
  '1.008', '4.0026', '6.94', '9.0122', '10.81', '12.011', '14.007', '15.999', '18.998', '20.180',
  '22.990', '24.305', '26.982', '28.085', '30.974', '32.06', '35.45', '39.95', '39.098', '40.078',
  '44.956', '47.867', '50.942', '51.996', '54.938', '55.845', '58.933', '58.693', '63.546', '65.38',
  '69.723', '72.630', '74.922', '78.971', '79.904', '83.798', '85.468', '87.62', '88.906', '91.224',
  '92.906', '95.95', '[97]', '101.07', '102.91', '106.42', '107.87', '112.41', '114.82', '118.71',
  '121.76', '127.60', '126.90', '131.29', '132.91', '137.33', '138.91', '140.12', '140.91', '144.24',
  '[145]', '150.36', '151.96', '157.25', '158.93', '162.50', '164.93', '167.26', '168.93', '173.05',
  '174.97', '178.49', '180.95', '183.84', '186.21', '190.23', '192.22', '195.08', '196.97', '200.59',
  '204.38', '207.2', '208.98', '[209]', '[210]', '[222]', '[223]', '[226]', '[227]', '232.04',
  '231.04', '238.03', '[237]', '[244]', '[243]', '[247]', '[247]', '[251]', '[252]', '[257]',
  '[258]', '[259]', '[266]', '[267]', '[268]', '[269]', '[270]', '[269]', '[278]', '[281]',
  '[282]', '[285]', '[286]', '[289]', '[290]', '[293]', '[294]', '[294]',
];

const IONIZATION_EV = [
  13.598, 24.587, 5.392, 9.323, 8.298, 11.260, 14.534, 13.618, 17.423, 21.565,
  5.139, 7.646, 5.986, 8.152, 10.487, 10.360, 12.968, 15.760, 4.341, 6.113,
  6.561, 6.828, 6.746, 6.767, 7.434, 7.902, 7.881, 7.640, 7.726, 9.394,
  5.999, 7.899, 9.789, 9.752, 11.814, 14.000, 4.177, 5.695, 6.217, 6.634,
  6.759, 7.092, 7.119, 7.361, 7.459, 8.337, 7.576, 8.994, 5.786, 7.344,
  8.608, 9.010, 10.451, 12.130, 3.894, 5.212, 5.577, 5.539, 5.473, 5.525,
  5.582, 5.644, 5.670, 6.150, 5.864, 5.939, 6.022, 6.108, 6.184, 6.254,
  5.426, 6.825, 7.550, 7.864, 7.834, 8.438, 8.967, 8.959, 9.226, 10.438,
  6.108, 7.417, 7.286, 8.414, 9.318, 10.749, 4.073, 5.278, 5.380, 6.307,
  5.89, 6.194, 6.266, 6.026, 5.974, 5.991, 6.198, 6.282, 6.368, 6.50,
  6.58, 6.626, 4.96, null, null, null, null, null, null, null,
  null, null, null, null, null, null, null, null,
];

export const CATEGORY_LABELS = {
  'alkali': 'Metallo alcalino',
  'alkaline': 'Metallo alcalino-terroso',
  'transition': 'Metallo di transizione',
  'post-transition': 'Metallo del blocco p',
  'metalloid': 'Semimetallo',
  'nonmetal': 'Non metallo',
  'halogen': 'Alogeno',
  'noble': 'Gas nobile',
  'lanthanide': 'Lantanide',
  'actinide': 'Attinide',
  'unknown': 'Proprietà chimiche non note',
};

function categoryOf(Z) {
  if ([2, 10, 18, 36, 54, 86].includes(Z)) return 'noble';
  if ([9, 17, 35, 53, 85].includes(Z)) return 'halogen';
  if ([3, 11, 19, 37, 55, 87].includes(Z)) return 'alkali';
  if ([4, 12, 20, 38, 56, 88].includes(Z)) return 'alkaline';
  if ([1, 6, 7, 8, 15, 16, 34].includes(Z)) return 'nonmetal';
  if ([5, 14, 32, 33, 51, 52].includes(Z)) return 'metalloid';
  if (Z >= 57 && Z <= 71) return 'lanthanide';
  if (Z >= 89 && Z <= 103) return 'actinide';
  if (Z >= 109) return 'unknown';
  if ([13, 31, 49, 50, 81, 82, 83, 84].includes(Z)) return 'post-transition';
  return 'transition';
}

// Periodo e gruppo IUPAC (1–18). Lantanidi e attinidi: gruppo null, posti nelle righe del blocco f.
const PERIOD_STARTS = [1, 3, 11, 19, 37, 55, 87, 119];

function periodOf(Z) {
  let p = 0;
  while (Z >= PERIOD_STARTS[p + 1]) p++;
  return p + 1;
}

function groupOf(Z) {
  if (Z === 1) return 1;
  if (Z === 2) return 18;
  const period = periodOf(Z);
  const offset = Z - PERIOD_STARTS[period - 1]; // 0-based posizione nel periodo
  if (period <= 3) return offset < 2 ? offset + 1 : offset + 11;
  if (period <= 5) return offset + 1;
  // Periodi 6 e 7: 2 (s) + 15 (f, con La/Ac) + 9 (d rimanenti) + 6 (p)
  if (offset < 2) return offset + 1;
  if (offset < 17) return null; // blocco f (La–Lu, Ac–Lr)
  return offset - 17 + 4;
}

export const ELEMENTS = SYMBOLS.map((symbol, i) => {
  const Z = i + 1;
  return {
    Z,
    symbol,
    name: NAMES[i],
    mass: MASSES[i],
    ionizationEV: IONIZATION_EV[i],
    period: periodOf(Z),
    group: groupOf(Z),
    category: categoryOf(Z),
  };
});

export function element(Z) {
  return ELEMENTS[Z - 1];
}
