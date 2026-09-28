// Configurazioni elettroniche dello stato fondamentale, regole di Slater e regole di Hund.

export const L_LETTERS = ['s', 'p', 'd', 'f', 'g', 'h', 'i'];

// Ordine di riempimento di Madelung (regola n + l, a parità di n + l prima n minore).
export const MADELUNG_ORDER = (() => {
  const list = [];
  for (let n = 1; n <= 7; n++) {
    for (let l = 0; l < n && l <= 3; l++) list.push({ n, l });
  }
  list.sort((a, b) => (a.n + a.l) - (b.n + b.l) || a.n - b.n);
  // Oltre 7p (Z > 118) non serve.
  return list.filter(s => (s.n + s.l) <= 8);
})();

// Eccezioni sperimentali alla regola di Madelung (NIST ASD), come variazioni rispetto alla configurazione di Madelung.
// Chiave: Z. Valore: occupazioni che sostituiscono quelle di Madelung.
const EXCEPTIONS = {
  24: { '3d': 5, '4s': 1 },            // Cr
  29: { '3d': 10, '4s': 1 },           // Cu
  41: { '4d': 4, '5s': 1 },            // Nb
  42: { '4d': 5, '5s': 1 },            // Mo
  44: { '4d': 7, '5s': 1 },            // Ru
  45: { '4d': 8, '5s': 1 },            // Rh
  46: { '4d': 10, '5s': 0 },           // Pd
  47: { '4d': 10, '5s': 1 },           // Ag
  57: { '4f': 0, '5d': 1 },            // La
  58: { '4f': 1, '5d': 1 },            // Ce
  64: { '4f': 7, '5d': 1 },            // Gd
  78: { '5d': 9, '6s': 1 },            // Pt
  79: { '5d': 10, '6s': 1 },           // Au
  89: { '5f': 0, '6d': 1 },            // Ac
  90: { '5f': 0, '6d': 2 },            // Th
  91: { '5f': 2, '6d': 1 },            // Pa
  92: { '5f': 3, '6d': 1 },            // U
  93: { '5f': 4, '6d': 1 },            // Np
  96: { '5f': 7, '6d': 1 },            // Cm
  103: { '6d': 0, '7p': 1 },           // Lr
};

export function subshellLabel(n, l) {
  return `${n}${L_LETTERS[l]}`;
}

/**
 * Configurazione dello stato fondamentale.
 * Restituisce un array di sottolivelli {n, l, occ} ordinati per n e poi per l (notazione "per shell").
 */
export function groundStateConfiguration(Z) {
  const occ = new Map();
  let remaining = Z;
  for (const { n, l } of MADELUNG_ORDER) {
    if (remaining <= 0) break;
    const cap = 2 * (2 * l + 1);
    const q = Math.min(cap, remaining);
    occ.set(subshellLabel(n, l), q);
    remaining -= q;
  }
  const exc = EXCEPTIONS[Z];
  if (exc) {
    for (const [key, q] of Object.entries(exc)) occ.set(key, q);
  }
  const list = [];
  for (const [key, q] of occ) {
    if (q <= 0) continue;
    const n = parseInt(key, 10);
    const l = L_LETTERS.indexOf(key.slice(String(n).length));
    list.push({ n, l, occ: q });
  }
  list.sort((a, b) => a.n - b.n || a.l - b.l);
  const total = list.reduce((s, x) => s + x.occ, 0);
  if (total !== Z) throw new Error(`Configurazione incoerente per Z=${Z}: ${total} elettroni`);
  return list;
}

export function isMadelungException(Z) {
  return Z in EXCEPTIONS;
}

// Nuclei dei gas nobili per la notazione abbreviata.
const NOBLE_CORES = [
  { Z: 86, symbol: 'Rn' }, { Z: 54, symbol: 'Xe' }, { Z: 36, symbol: 'Kr' },
  { Z: 18, symbol: 'Ar' }, { Z: 10, symbol: 'Ne' }, { Z: 2, symbol: 'He' },
];

/** Notazione abbreviata con nucleo di gas nobile, es. [Ar] 3d⁵ 4s¹. */
export function configurationParts(Z, config) {
  for (const core of NOBLE_CORES) {
    if (core.Z >= Z) continue;
    const coreCfg = groundStateConfiguration(core.Z);
    const coreSet = new Map(coreCfg.map(s => [subshellLabel(s.n, s.l), s.occ]));
    const isSuperset = coreCfg.every(s => config.some(c => c.n === s.n && c.l === s.l && c.occ >= s.occ));
    if (!isSuperset) continue;
    const rest = config.filter(c => coreSet.get(subshellLabel(c.n, c.l)) !== c.occ);
    return { core: core.symbol, subshells: sortByFilling(rest) };
  }
  return { core: null, subshells: sortByFilling(config) };
}

// Dopo il nucleo di gas nobile si ordina per n e poi per l, come nelle tabelle NIST: [Xe] 4f¹⁴ 5d¹⁰ 6s¹.
function sortByFilling(list) {
  return [...list].sort((a, b) => a.n - b.n || a.l - b.l);
}

// ---------------------------------------------------------------------------
// Regole di Slater (1930): carica nucleare efficace Z_eff = Z − σ.
// Gruppi: [1s] [2s,2p] [3s,3p] [3d] [4s,4p] [4d] [4f] [5s,5p] [5d] [5f] [6s,6p] [6d] [7s,7p]
// ---------------------------------------------------------------------------

function slaterGroupKey(n, l) {
  // s e p condividono il gruppo; d e f formano gruppi separati, posti dopo ns,np dello stesso n.
  return n * 10 + (l <= 1 ? 0 : l - 1);
}

export function slaterZeff(Z, config, n, l) {
  const key = slaterGroupKey(n, l);
  let sigma = 0;
  for (const s of config) {
    const k = slaterGroupKey(s.n, s.l);
    let count = s.occ;
    if (k === key) {
      count -= 1; // escludi l'elettrone considerato
      if (count <= 0) continue;
      sigma += count * (n === 1 ? 0.30 : 0.35);
    } else if (k < key) {
      if (l <= 1) {
        // elettrone s o p: guscio n-1 → 0.85, gusci più interni → 1.00; nd, nf con stesso n stanno "a destra"
        if (s.n === n - 1) sigma += count * 0.85;
        else if (s.n < n - 1) sigma += count * 1.00;
      } else {
        // elettrone d o f: tutti i gruppi a sinistra schermano completamente
        sigma += count * 1.00;
      }
    }
  }
  return Z - sigma;
}

/** Numero quantico principale efficace di Slater n*. */
export function slaterNStar(n) {
  return [0, 1, 2, 3, 3.7, 4.0, 4.2, 4.4][n] ?? n;
}

// ---------------------------------------------------------------------------
// Regola di Hund: disposizione degli elettroni nelle caselle m_l.
// ---------------------------------------------------------------------------

/** Restituisce per ogni m_l (da +l a −l) il numero di elettroni con spin su/giù. */
export function hundBoxes(l, occ) {
  const m = 2 * l + 1;
  const boxes = [];
  for (let i = 0; i < m; i++) boxes.push({ ml: l - i, up: 0, down: 0 });
  for (let e = 0; e < occ; e++) {
    if (e < m) boxes[e].up = 1;
    else boxes[e - m].down = 1;
  }
  return boxes;
}

const L_TERM = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'Q'];

/**
 * Termine spettroscopico dello stato fondamentale previsto dalle regole di Hund
 * (accoppiamento Russell–Saunders, gusci aperti combinati con la stessa regola).
 */
export function hundTerm(config) {
  let S2 = 0; // 2S
  let L = 0;
  let parity = 0;
  let openShells = 0;
  let moreThanHalf = false;
  for (const s of config) {
    const cap = 2 * (2 * s.l + 1);
    parity += s.l * s.occ;
    if (s.occ === cap) continue;
    openShells++;
    const boxes = hundBoxes(s.l, s.occ);
    let ms2 = 0;
    let ml = 0;
    for (const b of boxes) {
      ms2 += b.up - b.down;
      ml += b.ml * (b.up + b.down);
    }
    S2 += ms2;
    L += ml;
    if (s.occ > cap / 2) moreThanHalf = true;
  }
  const S = S2 / 2;
  const J = moreThanHalf ? L + S : Math.abs(L - S);
  const odd = parity % 2 === 1;
  return {
    multiplicity: S2 + 1,
    L,
    letter: L_TERM[L] ?? `L=${L}`,
    J,
    Jlabel: Number.isInteger(J) ? String(J) : `${Math.round(2 * J)}/2`,
    odd,
    openShells,
    unpaired: S2,
  };
}
