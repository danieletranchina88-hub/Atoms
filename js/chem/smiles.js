// Lettore SMILES (Weininger, 1988): converte una stringa come "CC(=O)O" in un grafo molecolare
// con ordini di legame, cariche formali e idrogeni espliciti. Gli anelli aromatici vengono
// "kekulizzati" (assegnazione di legami singoli e doppi alternati).

import { ELEMENTS } from '../physics/elements.js';

const SYMBOL_TO_Z = new Map(ELEMENTS.map(e => [e.symbol, e.Z]));
const ORGANIC = ['Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I'];
const AROMATIC = ['b', 'c', 'n', 'o', 'p', 's'];
const NORMAL_VALENCE = { 5: [3], 6: [4], 7: [3, 5], 8: [2], 15: [3, 5], 16: [2, 4, 6], 9: [1], 17: [1], 35: [1], 53: [1] };

export class SmilesError extends Error {}

/**
 * @returns {{ atoms: [{Z, charge, aromatic, hExplicit}], bonds: [{a, b, order, aromatic}] }}
 * con gli idrogeni aggiunti come atomi espliciti.
 */
export function parseSmiles(smiles) {
  const s = smiles.trim();
  if (!s) throw new SmilesError('Stringa SMILES vuota');
  const atoms = [];
  const bonds = [];
  const stack = [];
  const rings = new Map();
  let prev = -1;
  let pendingBond = null;
  let i = 0;

  const addAtom = (atom) => {
    atoms.push(atom);
    const idx = atoms.length - 1;
    if (prev >= 0) addBond(prev, idx, pendingBond);
    pendingBond = null;
    prev = idx;
  };
  const addBond = (a, b, sym) => {
    let order = 1;
    let aromatic = false;
    if (sym === '=') order = 2;
    else if (sym === '#') order = 3;
    else if (sym === '$') order = 4;
    else if (sym === ':') aromatic = true;
    else if (sym === null && atoms[a].aromatic && atoms[b].aromatic) aromatic = true;
    if (bonds.some(x => (x.a === a && x.b === b) || (x.a === b && x.b === a))) throw new SmilesError('Legame duplicato');
    bonds.push({ a, b, order, aromatic });
  };

  while (i < s.length) {
    const ch = s[i];
    if (ch === '(') { stack.push(prev); i++; continue; }
    if (ch === ')') {
      if (!stack.length) throw new SmilesError('Parentesi chiusa senza apertura');
      prev = stack.pop(); i++; continue;
    }
    if ('-=#$:/\\'.includes(ch)) {
      pendingBond = (ch === '/' || ch === '\\' || ch === '-') ? '-' : ch;
      i++;
      continue;
    }
    if (ch === '.') { prev = -1; pendingBond = null; i++; continue; }
    if (/[0-9%]/.test(ch)) {
      let num;
      if (ch === '%') { num = parseInt(s.slice(i + 1, i + 3), 10); i += 3; } else { num = parseInt(ch, 10); i++; }
      if (prev < 0) throw new SmilesError('Chiusura di anello senza atomo');
      if (rings.has(num)) {
        const r = rings.get(num);
        rings.delete(num);
        addBond(r.atom, prev, pendingBond ?? r.bond);
      } else {
        rings.set(num, { atom: prev, bond: pendingBond });
      }
      pendingBond = null;
      continue;
    }
    if (ch === '[') {
      const end = s.indexOf(']', i);
      if (end < 0) throw new SmilesError('Parentesi quadra non chiusa');
      addAtom(parseBracket(s.slice(i + 1, end)));
      i = end + 1;
      continue;
    }
    // atomo del sottoinsieme organico
    const two = s.slice(i, i + 2);
    let sym = null;
    if (ORGANIC.includes(two)) sym = two;
    else if (ORGANIC.includes(ch)) sym = ch;
    else if (AROMATIC.includes(ch)) sym = ch;
    if (!sym) throw new SmilesError(`Simbolo non riconosciuto: "${ch}" in posizione ${i + 1}`);
    const aromatic = sym === sym.toLowerCase();
    const Z = SYMBOL_TO_Z.get(aromatic ? sym.toUpperCase() : sym);
    addAtom({ Z, charge: 0, aromatic, hExplicit: null, organic: true });
    i += sym.length;
  }
  if (stack.length) throw new SmilesError('Parentesi aperta non chiusa');
  if (rings.size) throw new SmilesError(`Anello ${[...rings.keys()][0]} non chiuso`);

  kekulize(atoms, bonds);
  return addHydrogens(atoms, bonds);
}

function parseBracket(txt) {
  const m = /^(\d+)?([A-Z][a-z]?|[bcnops]|se|as)(@{0,2})(H\d*)?([+-]+\d*|[+-]\d+)?(:\d+)?$/.exec(txt);
  if (!m) throw new SmilesError(`Atomo tra parentesi non valido: [${txt}]`);
  const aromatic = m[2] === m[2].toLowerCase();
  const sym = aromatic ? m[2][0].toUpperCase() + m[2].slice(1) : m[2];
  const Z = SYMBOL_TO_Z.get(sym);
  if (!Z) throw new SmilesError(`Elemento sconosciuto: ${m[2]}`);
  let h = 0;
  if (m[4]) h = m[4].length > 1 ? parseInt(m[4].slice(1), 10) : 1;
  let charge = 0;
  if (m[5]) {
    const sign = m[5][0] === '+' ? 1 : -1;
    const rest = m[5].slice(1);
    if (/^\d+$/.test(rest)) charge = sign * parseInt(rest, 10);
    else charge = sign * m[5].length;
  }
  return { Z, charge, aromatic, hExplicit: h, organic: false };
}

/** Assegna legami doppi agli anelli aromatici (accoppiamento perfetto per backtracking). */
function kekulize(atoms, bonds) {
  const arom = bonds.filter(b => b.aromatic);
  if (!arom.length) return;
  // Un atomo aromatico "richiede" un doppio legame se la sua valenza libera è 1:
  // valenza − (legami non aromatici + legami aromatici + idrogeni) = 1
  const needs = atoms.map((a, idx) => {
    if (!a.aromatic) return false;
    let valence;
    if (a.Z === 6) valence = 4 - Math.abs(a.charge);
    else if (a.Z === 5) valence = 3 + a.charge;
    else valence = (NORMAL_VALENCE[a.Z] ?? [3])[0] + a.charge;
    let used = 0;
    let aromaticBonds = 0;
    for (const b of bonds) {
      if (b.a !== idx && b.b !== idx) continue;
      if (b.aromatic) aromaticBonds++; else used += b.order;
    }
    let h = a.hExplicit ?? 0;
    if (a.organic && a.hExplicit === null && a.Z === 6 && used + aromaticBonds < 3) h = 1;
    return valence - used - aromaticBonds - h === 1;
  });
  const adj = atoms.map(() => []);
  arom.forEach((b, k) => { adj[b.a].push(k); adj[b.b].push(k); });
  const matched = new Array(atoms.length).fill(false);
  const order = atoms.map((_, i) => i).filter(i => needs[i]);
  const solve = (pos) => {
    while (pos < order.length && matched[order[pos]]) pos++;
    if (pos === order.length) return true;
    const a = order[pos];
    for (const k of adj[a]) {
      const b = arom[k];
      const other = b.a === a ? b.b : b.a;
      if (!needs[other] || matched[other]) continue;
      b.order = 2; matched[a] = matched[other] = true;
      if (solve(pos + 1)) return true;
      b.order = 1; matched[a] = matched[other] = false;
    }
    return false;
  };
  for (const b of arom) b.order = 1;
  if (!solve(0)) throw new SmilesError('Impossibile assegnare una struttura di Kekulé al sistema aromatico');
}

function addHydrogens(atoms, bonds) {
  const out = atoms.map(a => ({ Z: a.Z, charge: a.charge, aromatic: a.aromatic }));
  const outBonds = bonds.map(b => ({ a: b.a, b: b.b, order: b.order, aromatic: b.aromatic }));
  atoms.forEach((a, idx) => {
    let h = a.hExplicit ?? 0;
    if (a.organic && a.hExplicit === null) {
      const sum = bonds.reduce((s, b) => s + (b.a === idx || b.b === idx ? b.order : 0), 0);
      const vals = NORMAL_VALENCE[a.Z] ?? [0];
      const target = vals.find(v => v >= sum) ?? sum;
      h = Math.max(0, target - sum);
    }
    for (let k = 0; k < h; k++) {
      out.push({ Z: 1, charge: 0, aromatic: false });
      outBonds.push({ a: idx, b: out.length - 1, order: 1, aromatic: false });
    }
  });
  return { atoms: out, bonds: outBonds };
}

// Sequenza di elettronegatività della nomenclatura IUPAC (Red Book 2005, Tabella VI):
// nelle formule dei composti inorganici l'elemento che compare prima in questa lista si scrive per ultimo.
const IUPAC_SEQUENCE = ['F', 'Cl', 'Br', 'I', 'At', 'O', 'S', 'Se', 'Te', 'Po', 'H', 'N', 'P', 'As', 'Sb', 'Bi', 'C', 'Si', 'Ge', 'Sn', 'Pb',
  'B', 'Al', 'Ga', 'In', 'Tl', 'Zn', 'Cd', 'Hg', 'Cu', 'Ag', 'Au', 'Ni', 'Pd', 'Pt', 'Co', 'Rh', 'Ir', 'Fe', 'Ru', 'Os', 'Mn', 'Tc', 'Re',
  'Cr', 'Mo', 'W', 'V', 'Nb', 'Ta', 'Ti', 'Zr', 'Hf', 'Sc', 'Y', 'La', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr',
  'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'];

/**
 * Formula bruta: ordine di Hill (C, H, poi alfabetico) per i composti organici,
 * ordine di elettronegatività IUPAC per gli altri (NH₃, H₂O, SF₆, NaCl…).
 */
export function hillFormula(atoms, charge = 0) {
  const counts = new Map();
  for (const a of atoms) {
    const sym = ELEMENTS[a.Z - 1].symbol;
    counts.set(sym, (counts.get(sym) ?? 0) + 1);
  }
  const keys = [...counts.keys()];
  let order;
  if (counts.has('C') && counts.has('H')) order = ['C', 'H', ...keys.filter(k => k !== 'C' && k !== 'H').sort()];
  else {
    const rank = (k) => { const i = IUPAC_SEQUENCE.indexOf(k); return i < 0 ? 50 : i; };
    order = keys.sort((a, b) => rank(b) - rank(a));
  }
  const sub = (n) => String(n).split('').map(d => '₀₁₂₃₄₅₆₇₈₉'[+d]).join('');
  let f = order.map(k => k + (counts.get(k) > 1 ? sub(counts.get(k)) : '')).join('');
  if (charge) {
    const sup = (n) => String(n).split('').map(d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d]).join('');
    f += (Math.abs(charge) > 1 ? sup(Math.abs(charge)) : '') + (charge > 0 ? '⁺' : '⁻');
  }
  return f;
}
