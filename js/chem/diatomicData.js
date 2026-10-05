// Costanti sperimentali dello stato fondamentale delle molecole biatomiche usate in "Come nasce un legame".
//
// ωe, ωexe, Be, re: K. P. Huber, G. Herzberg, "Molecular Spectra and Molecular Structure IV. Constants of Diatomic
// Molecules" (Van Nostrand Reinhold, 1979), letti dal NIST Chemistry WebBook (SRD 69, "Constants of diatomic
// molecules", https://webbook.nist.gov/cgi/cbook.cgi?ID=<CAS>&Mask=1000), riga dello stato X.
// D0 (energia per rompere la molecola dal livello vibrazionale v = 0) dalle fonti indicate per ciascuna molecola.
// De, la profondità della buca di potenziale (quella che calcola la chimica quantistica, senza energia di punto zero),
// si ottiene come De = D0 + G(0), con G(0) = ωe/2 − ωexe/4.
// Masse isotopiche (u): AME2020 / NIST "Atomic Weights and Isotopic Compositions"; le costanti si riferiscono
// all'isotopologo più abbondante (¹H, ⁷Li, ¹²C, ¹⁴N, ¹⁶O, ¹⁹F, ³⁵Cl).

export const CM1_EV = 1 / 8065.543937;           // 1 cm⁻¹ in eV (CODATA 2018)
export const EV_HARTREE = 1 / 27.211386245988;

export const ISOTOPE_MASS = { 1: 1.00782503223, 2: 4.00260325413, 3: 7.0160034366, 6: 12, 7: 14.00307400443, 8: 15.99491461957, 9: 18.99840316273, 17: 34.968852682 };

const HH = 'Huber e Herzberg (1979), via NIST WebBook';

/**
 * Z: numeri atomici [A, B]; charge, mult: carica e molteplicità dello stato fondamentale; atoms: molteplicità degli
 * atomi separati (stati fondamentali), atomCharges le loro cariche (H₂⁺ → H + H⁺); state: termine spettroscopico; exp: costanti sperimentali (cm⁻¹, Å, eV).
 */
export const DIATOMICS = [
  {
    id: 'H2+', name: 'H₂⁺', Z: [1, 1], charge: 1, mult: 2, atoms: [2, 1], atomCharges: [0, 1], state: 'X ²Σg⁺',
    blurb: 'Il legame più semplice che esista: un solo elettrone condiviso da due protoni. Hartree–Fock con un elettrone è esatto nella base.',
    exp: {
      we: 2321.7, wexe: 66.2, Be: 29.8, re: 1.052,
      D0: 2.6507, D0src: 'Huber e Herzberg (1979)',
      note: 'Per H₂⁺ l\'equazione di Schrödinger (Born–Oppenheimer) si risolve esattamente: Re = 1,997 bohr, De = 0,10263 Eh = 2,793 eV (Bates, Ledsham e Stewart, Phil. Trans. R. Soc. A 246, 215, 1953).',
    },
  },
  {
    id: 'H2', name: 'H₂', Z: [1, 1], charge: 0, mult: 1, atoms: [2, 2], state: 'X ¹Σg⁺',
    blurb: 'Due elettroni, due protoni: il legame covalente di Heitler e London (1927). Qui la soluzione esatta nella base (interazione di configurazioni completa) si confronta con Hartree–Fock.',
    exp: {
      we: 4401.21, wexe: 121.33, Be: 60.853, re: 0.74144,
      D0: 36118.0696 * CM1_EV, D0src: 'J. Liu et al., J. Chem. Phys. 130, 174306 (2009): D0 = 36 118,0696 cm⁻¹',
      note: 'Il calcolo teorico più accurato (Kołos e Wolniewicz, J. Chem. Phys. 49, 404, 1968; Pachucki 2010) dà De = 38 292,9 cm⁻¹ = 4,7477 eV.',
    },
  },
  {
    id: 'He2', name: 'He₂', Z: [2, 2], charge: 0, mult: 1, atoms: [1, 1], state: 'X ¹Σg⁺',
    scanMin: 1.1, blurb: 'Il legame che non si forma: σg² σu², ordine di legame zero. Resta solo un\'attrazione di van der Waals di circa 1 meV, che Hartree–Fock non può descrivere.',
    exp: {
      re: 2.9683, De: 10.956 * 8.617333262e-5, D0src: 'R. A. Aziz, A. R. Janzen, M. R. Moldover, Phys. Rev. Lett. 74, 1586 (1995): ε/k = 10,956 K, rm = 2,9683 Å',
      note: 'La buca è così poco profonda che (⁴He)₂ ha un solo livello legato, a circa 1 mK sotto la soglia (Grisenti et al., Phys. Rev. Lett. 85, 2284, 2000).',
      noMorse: true,
    },
  },
  {
    id: 'LiH', name: 'LiH', Z: [3, 1], charge: 0, mult: 1, atoms: [2, 2], state: 'X ¹Σ⁺',
    blurb: 'Legame polare: l\'elettrone 2s del litio si sposta in buona parte verso l\'idrogeno (Li^δ+ H^δ−).',
    exp: {
      we: 1405.65, wexe: 23.20, Be: 7.5131, re: 1.5957,
      D0: (34492.5 - 14903.9) * CM1_EV, D0src: 'NIST WebBook (nota allo stato A): limite di dissociazione 34 492,5 cm⁻¹ meno l\'eccitazione Li 2s→2p, 14 903,9 cm⁻¹ (NIST ASD)',
    },
  },
  {
    id: 'Li2', name: 'Li₂', Z: [3, 3], charge: 0, mult: 1, atoms: [2, 2], state: 'X ¹Σg⁺',
    blurb: 'Un legame σ fatto da due elettroni 2s molto diffusi: lungo (2,67 Å) e debole.',
    exp: {
      we: 351.43, wexe: 2.610, Be: 0.67264, re: 2.6729,
      De: 8516.78 * CM1_EV, D0src: 'R. J. Le Roy et al., J. Chem. Phys. 131, 204309 (2009): De = 8516,78 cm⁻¹',
    },
  },
  {
    id: 'N2', name: 'N₂', Z: [7, 7], charge: 0, mult: 1, atoms: [4, 4], state: 'X ¹Σg⁺',
    blurb: 'Il triplo legame (un σ e due π) fra due atomi con tre elettroni spaiati ciascuno: uno dei legami più forti della chimica.',
    exp: { we: 2358.57, wexe: 14.324, Be: 1.99824, re: 1.09768, D0: 9.759, D0src: HH },
  },
  {
    id: 'O2', name: 'O₂', Z: [8, 8], charge: 0, mult: 3, atoms: [3, 3], state: 'X ³Σg⁻',
    blurb: 'La molecola paramagnetica: i due elettroni negli orbitali π* degeneri restano spaiati (regola di Hund), come previsto da Mulliken.',
    exp: { we: 1580.19, wexe: 11.98, Be: 1.4377, re: 1.20752, D0: 5.116, D0src: HH },
  },
  {
    id: 'F2', name: 'F₂', Z: [9, 9], charge: 0, mult: 1, atoms: [2, 2], state: 'X ¹Σg⁺',
    blurb: 'Un legame σ singolo e debole: Hartree–Fock non lo trova nemmeno legato, è tutto merito della correlazione elettronica.',
    exp: { we: 916.64, wexe: 11.236, Be: 0.89019, re: 1.41193, D0: 1.602, D0src: HH },
  },
  {
    id: 'HF', name: 'HF', Z: [9, 1], charge: 0, mult: 1, atoms: [2, 2], state: 'X ¹Σ⁺',
    blurb: 'Legame covalente polare: il fluoro attira la coppia di legame (H^δ+ F^δ−).',
    exp: {
      we: 4138.32, wexe: 89.88, Be: 20.9557, re: 0.91680,
      D0: 47333 * CM1_EV, D0src: 'G. di Lonardo, A. E. Douglas, Can. J. Phys. 51, 434 (1973): 47 333 ± 60 cm⁻¹ (NIST WebBook)',
    },
  },
  {
    id: 'CO', name: 'CO', Z: [6, 8], charge: 0, mult: 1, atoms: [3, 3], state: 'X ¹Σ⁺',
    blurb: 'Isoelettronico con N₂, ma eteronucleare: il legame più forte che si conosca fra due atomi neutri.',
    exp: { we: 2169.81358, wexe: 13.28831, Be: 1.93128087, re: 1.128323, D0: 11.092, D0src: HH },
  },
  {
    id: 'HCl', name: 'HCl', Z: [17, 1], charge: 0, mult: 1, atoms: [2, 2], state: 'X ¹Σ⁺',
    blurb: 'Legame σ fra l\'1s dell\'idrogeno e un 3p del cloro.',
    exp: { we: 2990.946, wexe: 52.8186, Be: 10.59341, re: 1.27455, D0: 4.4336, D0src: HH },
  },
  {
    id: 'Cl2', name: 'Cl₂', Z: [17, 17], charge: 0, mult: 1, atoms: [2, 2], state: 'X ¹Σg⁺',
    blurb: 'Legame σ fra due orbitali 3p: più lungo e più debole di quelli della seconda riga.',
    exp: {
      we: 559.7, wexe: 2.67, Be: 0.2439, re: 1.987,
      D0: 19997.14 * CM1_EV, D0src: 'R. J. Le Roy, J. Chem. Phys. 59, 556 (1973), riportato dal NIST WebBook: D0 = 19 997,14 cm⁻¹',
    },
  },
];

/** Energia di punto zero G(0) = ωe/2 − ωexe/4 (cm⁻¹ → eV) e profondità della buca De (eV). */
export function experimentalWell(m) {
  const e = m.exp;
  const zpe = e.we ? (e.we / 2 - (e.wexe ?? 0) / 4) * CM1_EV : 0;
  const De = e.De ?? (e.D0 + zpe);
  return { De, D0: e.D0 ?? De - zpe, zpe, re: e.re, we: e.we };
}

/** Massa ridotta (u). */
export function reducedMass(Z) {
  const [a, b] = Z.map(z => ISOTOPE_MASS[z]);
  return a * b / (a + b);
}

/**
 * Potenziale di Morse costruito dalle costanti sperimentali: V(R) = De (1 − e^{−β(R−re)})² − De (zero agli atomi
 * separati), con β = ωe √(μ / 2De) (ωe come pulsazione). R in Å, V in eV.
 */
export function morseCurve(m) {
  const { De, re, we } = experimentalWell(m);
  if (m.exp.noMorse || !we) return null;
  const mu = reducedMass(m.Z) * 1.66053906660e-27;           // kg
  const omega = 2 * Math.PI * 2.99792458e10 * we;           // rad/s
  const DeJ = De * 1.602176634e-19;
  const beta = omega * Math.sqrt(mu / (2 * DeJ)) * 1e-10;   // Å⁻¹
  return { De, re, beta, V: (R) => De * (1 - Math.exp(-beta * (R - re))) ** 2 - De };
}
