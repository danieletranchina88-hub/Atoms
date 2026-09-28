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

// ---------------------------------------------------------------------------
// Scrittura SMILES e valenze standard (per l'editor della molecola)
// ---------------------------------------------------------------------------

const STANDARD_VALENCE = { 1: 1, 5: 3, 6: 4, 7: 3, 8: 2, 9: 1, 14: 4, 15: 3, 16: 2, 17: 1, 35: 1, 53: 1, 3: 1, 4: 2, 11: 1, 12: 2, 13: 3 };

/** Numero di legami "normale" per un atomo con la carica data (N⁺ ha 4 legami, O⁻ uno solo, …). */
export function standardValence(Z, charge = 0) {
  const v = STANDARD_VALENCE[Z];
  if (v === undefined) return null;
  if (Z === 7 || Z === 8 || Z === 15 || Z === 16) return v + charge;   // gruppi 15–16: la carica + aggiunge un legame
  if (Z === 5 || Z === 13) return v - charge;                          // boro: B⁻ ha 4 legami
  if (Z === 6) return v - Math.abs(charge);                            // carbocationi e carbanioni: 3 legami
  return v - Math.abs(charge);
}

/** Aggiunge o toglie idrogeni perché ogni atomo pesante abbia la valenza standard. */
export function fillHydrogens(graph) {
  const atoms = graph.atoms.map(a => ({ ...a }));
  let bonds = graph.bonds.map(b => ({ ...b }));
  const heavy = atoms.map((a, i) => i).filter(i => atoms[i].Z !== 1);
  // rimuove gli H esistenti legati ad atomi pesanti, poi li ricrea
  const removeH = new Set();
  bonds.forEach(b => {
    const [x, y] = [b.a, b.b];
    if (atoms[x].Z === 1 && atoms[y].Z !== 1) removeH.add(x);
    if (atoms[y].Z === 1 && atoms[x].Z !== 1) removeH.add(y);
  });
  const keep = atoms.map((a, i) => !removeH.has(i));
  const newIndex = [];
  const out = [];
  atoms.forEach((a, i) => { if (keep[i]) { newIndex[i] = out.length; out.push(a); } });
  bonds = bonds.filter(b => keep[b.a] && keep[b.b]).map(b => ({ ...b, a: newIndex[b.a], b: newIndex[b.b] }));
  const n = out.length;
  for (let i = 0; i < n; i++) {
    if (out[i].Z === 1) continue;
    const target = standardValence(out[i].Z, out[i].charge ?? 0);
    if (target === null) continue;
    const used = bonds.reduce((s, b) => s + (b.a === i || b.b === i ? b.order : 0), 0);
    for (let k = 0; k < target - used; k++) {
      out.push({ Z: 1, charge: 0, aromatic: false });
      bonds.push({ a: i, b: out.length - 1, order: 1, aromatic: false });
    }
  }
  void heavy;
  return { atoms: out, bonds };
}

/** Stringa SMILES (non canonica) da un grafo con idrogeni espliciti. */
export function writeSmiles(graph) {
  const { atoms, bonds } = graph;
  const n = atoms.length;
  const adj = Array.from({ length: n }, () => []);
  bonds.forEach(b => { adj[b.a].push({ j: b.b, order: b.order }); adj[b.b].push({ j: b.a, order: b.order }); });
  const isH = (i) => atoms[i].Z === 1 && adj[i].length === 1 && atoms[adj[i][0].j].Z !== 1 && !(atoms[i].charge);
  const heavy = [...Array(n).keys()].filter(i => !isH(i));
  if (!heavy.length) return '';
  const hCount = (i) => adj[i].filter(e => isH(e.j)).length;
  const bondSym = (o) => (o === 2 ? '=' : o === 3 ? '#' : '');
  const ORGANIC_SET = new Set([5, 6, 7, 8, 9, 15, 16, 17, 35, 53]);
  const atomText = (i) => {
    const a = atoms[i];
    const symb = ELEMENTS[a.Z - 1].symbol;
    const h = hCount(i);
    const charge = a.charge ?? 0;
    const heavyOrder = adj[i].filter(e => !isH(e.j)).reduce((s, e) => s + e.order, 0);
    const implicitOk = ORGANIC_SET.has(a.Z) && charge === 0 && (() => {
      const vals = NORMAL_VALENCE[a.Z] ?? [];
      const target = vals.find(v => v >= heavyOrder);
      return target !== undefined && target - heavyOrder === h;
    })();
    if (implicitOk) return symb;
    const hs = h ? (h > 1 ? `H${h}` : 'H') : '';
    const cs = charge ? (charge > 0 ? (charge > 1 ? `+${charge}` : '+') : (charge < -1 ? `-${-charge}` : '-')) : '';
    return `[${symb}${hs}${cs}]`;
  };
  const visited = new Array(n).fill(false);
  const ringLabels = new Map(); // "i-j" → numero
  let ringCounter = 1;
  // prima passata: trova i legami di chiusura d'anello con una DFS
  const closures = new Map(); // atomo → [{ partner, order, label }]
  const parent = new Array(n).fill(-1);
  const seen = new Array(n).fill(false);
  const dfs1 = (i) => {
    seen[i] = true;
    for (const e of adj[i]) {
      if (isH(e.j)) continue;
      if (e.j === parent[i]) continue;
      if (seen[e.j]) {
        const key = i < e.j ? `${i}-${e.j}` : `${e.j}-${i}`;
        if (!ringLabels.has(key)) {
          const label = ringCounter++;
          ringLabels.set(key, label);
          (closures.get(i) ?? closures.set(i, []).get(i)).push({ partner: e.j, order: e.order, label });
          (closures.get(e.j) ?? closures.set(e.j, []).get(e.j)).push({ partner: i, order: e.order, label });
        }
        continue;
      }
      parent[e.j] = i;
      dfs1(e.j);
    }
  };
  const parts = [];
  for (const start of heavy) {
    if (seen[start]) continue;
    dfs1(start);
    const write = (i) => {
      visited[i] = true;
      let s = atomText(i);
      for (const c of closures.get(i) ?? []) {
        const lab = c.label < 10 ? String(c.label) : `%${c.label}`;
        s += (visited[c.partner] ? bondSym(c.order) : '') + lab;
      }
      const children = adj[i].filter(e => !isH(e.j) && parent[e.j] === i && !visited[e.j]);
      children.forEach((e, k) => {
        const sub = bondSym(e.order) + write(e.j);
        s += k < children.length - 1 ? `(${sub})` : sub;
      });
      return s;
    };
    parts.push(write(start));
  }
  return parts.join('.');
}
