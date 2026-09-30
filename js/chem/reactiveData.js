// Parametri del potenziale didattico del progetto. Combina dati medi, stime e costanti
// empiriche: non è un campo di forze pubblicato e validato per reazioni generali.
//
//  • valenza (numero di legami covalenti normali) ed elettroni di valenza: regola dell'ottetto;
//  • raggi covalenti di Pyykkö (Chem. Eur. J. 2009) quando manca la lunghezza di legame misurata;
//  • elettronegatività χ e durezza J del metodo QEq/UFF (Rappé e Goddard, J. Phys. Chem. 1991;
//    Rappé et al., JACS 1992): χ ≈ (I + A)/2, J ≈ I − A;
//  • parametri di van der Waals x (distanza di minimo) e D (profondità) del campo UFF;
//  • energie medie di legame (kJ/mol) e lunghezze (Å) delle tabelle di chimica generale (CRC Handbook).

export const RX_ELEMENTS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 35, 36, 53, 54];

// Z: [valenza, elettroni di valenza, raggio covalente Å, χ eV, J eV, x_vdW Å, D_vdW kcal/mol]
export const ATOM_PARAMS = {
  1: [1, 1, 0.32, 4.528, 13.890, 2.886, 0.044],
  2: [0, 2, 0.46, 9.660, 29.840, 2.362, 0.056],
  3: [1, 1, 1.33, 3.006, 4.772, 2.451, 0.025],
  4: [2, 2, 1.02, 4.877, 8.886, 2.745, 0.085],
  5: [3, 3, 0.85, 5.110, 9.500, 4.083, 0.180],
  6: [4, 4, 0.75, 5.343, 10.126, 3.851, 0.105],
  7: [3, 5, 0.71, 6.899, 11.760, 3.660, 0.069],
  8: [2, 6, 0.63, 8.741, 13.364, 3.500, 0.060],
  9: [1, 7, 0.64, 10.874, 14.948, 3.364, 0.050],
  10: [0, 8, 0.67, 11.040, 21.100, 3.243, 0.042],
  11: [1, 1, 1.55, 2.843, 4.592, 2.983, 0.030],
  12: [2, 2, 1.39, 3.951, 7.386, 3.021, 0.111],
  13: [3, 3, 1.26, 4.060, 7.180, 4.499, 0.505],
  14: [4, 4, 1.16, 4.168, 6.974, 4.295, 0.402],
  15: [3, 5, 1.11, 5.463, 8.000, 4.147, 0.305],
  16: [2, 6, 1.03, 6.928, 8.972, 4.035, 0.274],
  17: [1, 7, 0.99, 8.564, 9.892, 3.947, 0.227],
  18: [0, 8, 0.96, 9.465, 12.710, 3.868, 0.185],
  19: [1, 1, 1.96, 2.421, 3.840, 3.812, 0.035],
  20: [2, 2, 1.71, 3.231, 5.760, 3.399, 0.238],
  35: [1, 7, 1.14, 7.790, 8.850, 4.189, 0.251],
  36: [0, 8, 1.17, 8.505, 11.430, 4.141, 0.220],
  53: [1, 7, 1.33, 6.822, 7.524, 4.500, 0.339],
  54: [0, 8, 1.31, 7.595, 9.950, 4.404, 0.332],
};

// Legami: "Z1-Z2" (Z1 ≤ Z2) → [D singolo, r singolo, D doppio, r doppio, D triplo, r triplo]
// (kJ/mol, Å; null = dato non disponibile). Energie medie di legame e lunghezze tipiche.
export const BOND_TABLE = {
  '1-1': [436, 0.741],
  '3-3': [105, 2.67],
  '5-5': [293, 1.70],
  '6-6': [346, 1.54, 614, 1.34, 839, 1.20],
  '7-7': [167, 1.45, 418, 1.25, 945, 1.10],
  '8-8': [142, 1.48, 498, 1.21],
  '9-9': [159, 1.42],
  '11-11': [74, 3.08],
  '12-12': [8, 3.89],
  '13-13': [133, 2.70],
  '14-14': [226, 2.35],
  '15-15': [201, 2.21, null, null, 489, 1.89],
  '16-16': [266, 2.05, 425, 1.89],
  '17-17': [242, 1.99],
  '19-19': [55, 3.92],
  '20-20': [15, 4.28],
  '35-35': [193, 2.28],
  '53-53': [151, 2.67],
  // idrogeno
  '1-3': [243, 1.60], '1-4': [226, 1.34], '1-5': [389, 1.19], '1-6': [413, 1.09], '1-7': [391, 1.01],
  '1-8': [463, 0.96], '1-9': [567, 0.92], '1-11': [186, 1.89], '1-12': [127, 1.73], '1-13': [285, 1.65],
  '1-14': [318, 1.48], '1-15': [322, 1.42], '1-16': [363, 1.34], '1-17': [431, 1.27], '1-19': [182, 2.24],
  '1-20': [168, 2.00], '1-35': [366, 1.41], '1-53': [299, 1.61],
  // carbonio
  '5-6': [356, 1.56], '6-7': [305, 1.47, 615, 1.28, 891, 1.16], '6-8': [358, 1.43, 745, 1.21, 1072, 1.13],
  '6-9': [485, 1.35], '6-14': [318, 1.87], '6-15': [264, 1.84], '6-16': [272, 1.82, 573, 1.56],
  '6-17': [339, 1.77], '6-35': [285, 1.94], '6-53': [213, 2.14],
  // azoto
  '5-7': [389, 1.58], '7-8': [201, 1.40, 607, 1.21], '7-9': [283, 1.37], '7-14': [355, 1.74],
  '7-17': [200, 1.75],
  // ossigeno
  '3-8': [341, 1.60], '4-8': [437, 1.33], '5-8': [536, 1.37], '8-9': [190, 1.42], '8-11': [256, 1.92],
  '8-12': [363, 1.75], '8-13': [512, 1.62], '8-14': [452, 1.63], '8-15': [335, 1.63, 544, 1.48],
  '8-16': [265, 1.57, 522, 1.48], '8-17': [218, 1.70], '8-19': [278, 2.17], '8-20': [402, 1.82],
  '8-35': [201, 1.85],
  // fluoro
  '3-9': [577, 1.56], '4-9': [573, 1.40], '5-9': [613, 1.31], '9-11': [477, 1.93], '9-12': [462, 1.75],
  '9-13': [664, 1.65], '9-14': [565, 1.56], '9-15': [490, 1.54], '9-16': [327, 1.56], '9-17': [253, 1.63],
  '9-19': [489, 2.17], '9-20': [529, 1.97],
  // cloro
  '3-17': [469, 2.02], '4-17': [388, 1.80], '5-17': [456, 1.75], '11-17': [412, 2.36], '12-17': [312, 2.18],
  '13-17': [420, 2.06], '14-17': [381, 2.02], '15-17': [326, 2.04], '16-17': [255, 2.01],
  '17-19': [433, 2.67], '17-20': [409, 2.44], '17-35': [218, 2.14], '17-53': [211, 2.32],
  // bromo e iodio
  '3-35': [418, 2.17], '11-35': [363, 2.50], '19-35': [380, 2.82], '35-53': [179, 2.47],
  '3-53': [345, 2.39], '11-53': [304, 2.71], '19-53': [325, 3.05],
};

// Energie dei legami omonucleari usate nella formula di Pauling quando manca un dato eteronucleare.
const HOMONUCLEAR_FALLBACK = { 4: 59, 10: 0 };

export const KJ_PER_EV = 96.48533212;
export const KE = 14.39964548;          // e²/4πε₀ in eV·Å
export const KB_EV = 8.617333262e-5;    // eV/K
export const S_BRENNER = 1.3;           // rapporto tra le pendenze repulsiva e attrattiva (Brenner)
export const DEFAULT_SHIFT = 0.31;      // relazione di Pauling r(n) = r(1) − c ln n (Å)

export function pairKey(Za, Zb) {
  return Za <= Zb ? `${Za}-${Zb}` : `${Zb}-${Za}`;
}

/**
 * Dati sperimentali di una coppia: energia e lunghezza del legame singolo; esponente p
 * dell'energia E(n) = D₁ nᵖ e coefficiente c di Pauling r(n) = r₁ − c ln n, ricavati
 * con i minimi quadrati dai legami doppi e tripli quando sono tabulati.
 */
export function bondData(Za, Zb, paulingChi) {
  const row = BOND_TABLE[pairKey(Za, Zb)];
  const pa = ATOM_PARAMS[Za], pb = ATOM_PARAMS[Zb];
  let D1, r1, p = 1, c = DEFAULT_SHIFT, source = 'tabella';
  if (row) {
    [D1, r1] = row;
    let sp = 0, sc = 0, sw = 0;
    for (let k = 1; k <= 2; k++) {
      const Dn = row[2 * k], rn = row[2 * k + 1];
      if (Dn == null) continue;
      const L = Math.log(k + 1);
      sp += Math.log(Dn / D1) * L;
      sc += (r1 - rn) * L;
      sw += L * L;
    }
    if (sw > 0) { p = sp / sw; c = sc / sw; }
  } else {
    // formula di Pauling: D(A–B) = ½[D(A–A) + D(B–B)] + 96,5 (Δχ)² kJ/mol
    const homo = (Z) => { const value = BOND_TABLE[pairKey(Z, Z)]?.[0] ?? HOMONUCLEAR_FALLBACK[Z]; if (value === undefined) throw new Error(`Dati di legame mancanti per Z=${Z}.`); return value; };
    if (!Number.isFinite(paulingChi[Za]) || !Number.isFinite(paulingChi[Zb])) throw new Error('Elettronegatività non disponibile.');
    const dchi = paulingChi[Za] - paulingChi[Zb];
    D1 = 0.5 * (homo(Za) + homo(Zb)) + 96.5 * dchi * dchi;
    r1 = pa[2] + pb[2];
    source = 'stima empirica di Pauling (non misura)';
  }
  return { D1, r1, p, c, source };
}
