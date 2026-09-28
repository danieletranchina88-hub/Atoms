// Analisi chimica "classica" di una molecola: strutture di Lewis e risonanza, cariche formali,
// numeri di ossidazione, teoria VSEPR, ibridazione, polarità dei legami.

import { PAULING, valenceElectrons, covalentRadius } from './elementData.js';
import { ELEMENTS } from '../physics/elements.js';

const sym = (Z) => ELEMENTS[Z - 1].symbol;
const chi = (Z) => PAULING[Z] ?? 2.0;
const period = (Z) => ELEMENTS[Z - 1].period;

/** Legami dalla geometria: distanza < 1,25 × somma dei raggi covalenti singoli. */
export function connectivityFromGeometry(atoms) {
  const bonds = [];
  const BOHR_PM = 52.917721090;
  for (let i = 0; i < atoms.length; i++) {
    for (let j = 0; j < i; j++) {
      const a = atoms[i].xyz, b = atoms[j].xyz;
      const d = Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]) * BOHR_PM;
      if (d < 1.25 * (covalentRadius(atoms[i].Z) + covalentRadius(atoms[j].Z))) bonds.push({ a: j, b: i, order: 1 });
    }
  }
  return bonds;
}

function neighborsOf(n, bonds) {
  const nb = Array.from({ length: n }, () => []);
  bonds.forEach((b, k) => { nb[b.a].push({ atom: b.b, bond: k }); nb[b.b].push({ atom: b.a, bond: k }); });
  return nb;
}

/** Numero massimo di legami (somma degli ordini) per un atomo. */
function maxBondSum(Z) {
  if (Z === 1) return 1;
  if (Z === 2) return 0;
  if (period(Z) === 2) return Z === 5 ? 4 : Z === 9 ? 1 : Z === 8 ? 3 : 4;
  return 6;
}

/**
 * Distribuisce gli elettroni non di legame e valuta la struttura.
 * orders: ordini di legame; restituisce { lone (elettroni non condivisi per atomo), fc, score } o null.
 */
function evaluate(atoms, bonds, nb, orders, totalCharge) {
  const n = atoms.length;
  const V = atoms.map(a => valenceElectrons(a.Z));
  const bondSum = new Array(n).fill(0);
  bonds.forEach((b, k) => { bondSum[b.a] += orders[k]; bondSum[b.b] += orders[k]; });
  for (let i = 0; i < n; i++) if (bondSum[i] > maxBondSum(atoms[i].Z) && atoms[i].Z <= 10) return null;
  let remaining = V.reduce((s, v) => s + v, 0) - totalCharge - 2 * orders.reduce((s, o) => s + o, 0);
  if (remaining < 0) return null;
  const lone = new Array(n).fill(0);
  const target = atoms.map(a => (a.Z <= 2 ? 2 : 8));
  // completa gli ottetti, prima gli atomi più elettronegativi
  const byChi = [...Array(n).keys()].sort((a, b) => chi(atoms[b].Z) - chi(atoms[a].Z) || nb[a].length - nb[b].length);
  for (const i of byChi) {
    const need = Math.max(0, target[i] - 2 * bondSum[i]);
    const give = Math.min(need, remaining);
    lone[i] += give;
    remaining -= give;
  }
  // elettroni in eccesso: ottetto espanso sugli atomi del 3° periodo e oltre (prima quelli centrali)
  if (remaining > 0) {
    const expandable = [...Array(n).keys()].filter(i => period(atoms[i].Z) >= 3 && atoms[i].Z > 2)
      .sort((a, b) => nb[b].length - nb[a].length || chi(atoms[a].Z) - chi(atoms[b].Z));
    for (const i of expandable) {
      const room = 12 - (2 * bondSum[i] + lone[i]);
      const give = Math.min(room, remaining);
      if (give > 0) { lone[i] += give; remaining -= give; }
    }
  }
  if (remaining > 0) return null;
  const fc = atoms.map((a, i) => V[i] - lone[i] - bondSum[i]);
  // punteggio: ottetti incompleti, cariche formali, cariche "sbagliate", ottetti espansi
  let score = 0;
  for (let i = 0; i < n; i++) {
    const Z = atoms[i].Z;
    const e = 2 * bondSum[i] + lone[i];
    if (Z > 2) {
      const deficit = Math.max(0, 8 - e);
      const tolerant = Z === 5 || Z === 4 || Z === 13 || Z === 3 || Z === 11 || Z === 12;
      score += deficit * (tolerant ? 1.5 : 12);
      if (e > 8) score += (period(Z) >= 3 ? 1 : 1000) * (e - 8) / 2;
    } else if (Z === 1 && e !== 2) score += 50;
    score += Math.abs(fc[i]) * 8;
    if (fc[i] < 0) score += (4 - chi(Z)) * 1.5;
    if (fc[i] > 0) score += chi(Z) * 1.5;
  }
  return { lone, fc, score, orders: [...orders] };
}

/**
 * Cerca le strutture di Lewis migliori per una data connettività (risonanza inclusa).
 * Varia gli ordini dei legami tra atomi non terminali-idrogeno/fluoro.
 */
export function lewisStructures(atoms, bonds, totalCharge, { maxVariable = 16 } = {}) {
  const n = atoms.length;
  const nb = neighborsOf(n, bonds);
  const variable = bonds.map((b, k) => k).filter(k => {
    const za = atoms[bonds[k].a].Z, zb = atoms[bonds[k].b].Z;
    return za !== 1 && zb !== 1 && za !== 9 && zb !== 9 && za !== 2 && zb !== 2;
  });
  if (variable.length > maxVariable) return null;
  const orders = bonds.map(() => 1);
  const results = [];
  const bondSum = new Array(n).fill(0);
  bonds.forEach(b => { bondSum[b.a]++; bondSum[b.b]++; });
  const totalV = atoms.reduce((s, a) => s + valenceElectrons(a.Z), 0) - totalCharge;
  const rec = (p, usedPairs) => {
    if (2 * usedPairs > totalV) return;
    if (p === variable.length) {
      const r = evaluate(atoms, bonds, nb, orders, totalCharge);
      if (r) results.push(r);
      return;
    }
    const k = variable[p];
    const { a, b } = bonds[k];
    for (let o = 1; o <= 3; o++) {
      if (o > 1) {
        if (bondSum[a] + 1 > maxBondSum(atoms[a].Z) || bondSum[b] + 1 > maxBondSum(atoms[b].Z)) break;
        bondSum[a]++; bondSum[b]++;
      }
      orders[k] = o;
      rec(p + 1, usedPairs + o - 1);
    }
    bondSum[a] -= orders[k] - 1; bondSum[b] -= orders[k] - 1;
    orders[k] = 1;
  };
  rec(0, bonds.length);
  if (!results.length) return [];
  results.sort((x, y) => x.score - y.score);
  const best = results[0].score;
  return results.filter(r => r.score < best + 1e-6);
}

/** Numeri di ossidazione: gli elettroni di ogni legame vanno all'atomo più elettronegativo. */
export function oxidationNumbers(atoms, bonds, orders, lone) {
  const ox = atoms.map((a) => valenceElectrons(a.Z));
  atoms.forEach((a, i) => { ox[i] -= lone[i]; });
  bonds.forEach((b, k) => {
    const e = 2 * orders[k];
    const ca = chi(atoms[b.a].Z), cb = chi(atoms[b.b].Z);
    if (atoms[b.a].Z === atoms[b.b].Z || Math.abs(ca - cb) < 1e-9) { ox[b.a] -= e / 2; ox[b.b] -= e / 2; }
    else if (ca > cb) ox[b.a] -= e;
    else ox[b.b] -= e;
  });
  return ox;
}

const VSEPR = {
  1: ['—', '—', '—'],
  2: { 0: ['lineare', 'lineare', '180°'] },
  3: { 0: ['trigonale planare', 'trigonale planare', '120°'], 1: ['trigonale planare', 'angolare (piegata)', '< 120°'] },
  4: {
    0: ['tetraedrica', 'tetraedrica', '109,5°'],
    1: ['tetraedrica', 'piramidale trigonale', '< 109,5°'],
    2: ['tetraedrica', 'angolare (piegata)', '< 109,5°'],
    3: ['tetraedrica', 'lineare', '—'],
  },
  5: {
    0: ['bipiramidale trigonale', 'bipiramidale trigonale', '90°, 120°'],
    1: ['bipiramidale trigonale', 'a cavalletto (altalena)', '< 90°, < 120°'],
    2: ['bipiramidale trigonale', 'a T', '< 90°'],
    3: ['bipiramidale trigonale', 'lineare', '180°'],
  },
  6: {
    0: ['ottaedrica', 'ottaedrica', '90°'],
    1: ['ottaedrica', 'piramidale quadrata', '< 90°'],
    2: ['ottaedrica', 'quadrato planare', '90°'],
  },
  7: { 0: ['bipiramidale pentagonale', 'bipiramidale pentagonale', '72°, 90°'] },
};
const HYBRID = { 2: 'sp', 3: 'sp²', 4: 'sp³', 5: 'sp³d', 6: 'sp³d²', 7: 'sp³d³' };

/**
 * Analisi completa di una molecola.
 * mol: { atoms: [{Z, charge?}], bonds: [{a,b,order}], charge }
 */
export function analyzeStructure(mol) {
  const { atoms, bonds } = mol;
  const n = atoms.length;
  const totalCharge = mol.charge ?? atoms.reduce((s, a) => s + (a.charge ?? 0), 0);
  const nb = neighborsOf(n, bonds);
  // struttura di partenza (dallo SMILES) e ricerca delle strutture di risonanza
  let given = null;
  if (bonds.every(b => b.order)) {
    const V = atoms.map(a => valenceElectrons(a.Z));
    const bondSum = new Array(n).fill(0);
    bonds.forEach(b => { bondSum[b.a] += b.order; bondSum[b.b] += b.order; });
    const lone = atoms.map((a, i) => V[i] - (a.charge ?? 0) - bondSum[i]);
    if (lone.every(x => x >= 0)) given = { orders: bonds.map(b => b.order), lone, fc: atoms.map(a => a.charge ?? 0) };
  }
  const found = lewisStructures(atoms, bonds, totalCharge);
  let chosen = given;
  let resonance = found ?? [];
  if (!chosen && resonance.length) chosen = resonance[0];
  if (chosen && resonance.length) {
    // se la struttura data è tra le migliori usiamo quella, altrimenti la migliore trovata
    const same = resonance.find(r => r.orders.every((o, k) => o === chosen.orders[k]));
    if (!same && given) {
      const g = evaluate(atoms, bonds, nb, given.orders, totalCharge);
      if (g && g.score > resonance[0].score + 1e-6) chosen = resonance[0];
    }
  }
  if (!chosen) chosen = { orders: bonds.map(b => b.order ?? 1), lone: new Array(n).fill(0), fc: new Array(n).fill(0) };
  const orders = chosen.orders;
  const lone = chosen.lone;
  const fc = chosen.fc ?? atoms.map((a, i) => valenceElectrons(a.Z) - lone[i] - bonds.reduce((s, b, k) => s + (b.a === i || b.b === i ? orders[k] : 0), 0));
  const ox = oxidationNumbers(atoms, bonds, orders, lone);
  // ordine di legame medio sulle strutture di risonanza
  const avgOrder = bonds.map((b, k) => (resonance.length ? resonance.reduce((s, r) => s + r.orders[k], 0) / resonance.length : orders[k]));

  const atomInfo = atoms.map((a, i) => {
    const nNb = nb[i].length;
    const lp = Math.floor(lone[i] / 2);
    const radical = lone[i] % 2;
    const domains = nNb + lp + radical;
    let vsepr = null;
    if (nNb >= 2) {
      const e = lp + radical;
      const row = VSEPR[domains]?.[e];
      vsepr = {
        notation: `AX${nNb > 1 ? subscript(nNb) : ''}${e ? `E${e > 1 ? subscript(e) : ''}` : ''}`,
        stericNumber: domains,
        electronGeometry: row?.[0] ?? '—',
        molecularGeometry: row?.[1] ?? '—',
        idealAngle: row?.[2] ?? '—',
      };
    }
    const bondElectrons = 2 * bonds.reduce((s, b, k) => s + (b.a === i || b.b === i ? orders[k] : 0), 0);
    return {
      index: i, Z: a.Z, symbol: sym(a.Z),
      chi: PAULING[a.Z],
      valence: valenceElectrons(a.Z),
      lonePairs: lp, radical,
      formalCharge: fc[i],
      oxidation: ox[i],
      electronsAround: bondElectrons + lone[i],
      neighbors: nNb,
      hybridization: nNb >= 2 || (nNb === 1 && lp > 0 && a.Z > 2) ? (HYBRID[domains] ?? '—') : '—',
      vsepr,
    };
  });

  const bondInfo = bonds.map((b, k) => {
    const za = atoms[b.a].Z, zb = atoms[b.b].Z;
    const dchi = Math.abs(chi(za) - chi(zb));
    const ionic = 1 - Math.exp(-dchi * dchi / 4); // Pauling
    const order = orders[k];
    return {
      a: b.a, b: b.b, order, averageOrder: avgOrder[k],
      sigma: 1, pi: order - 1,
      deltaChi: dchi,
      ionicCharacter: ionic,
      type: dchi < 0.4 ? 'covalente apolare' : dchi < 1.9 ? 'covalente polare' : 'ionico',
      // lunghezza prevista dai raggi covalenti di Pyykkö (pm)
      predictedLength: covalentRadius(za, Math.round(avgOrder[k])) + covalentRadius(zb, Math.round(avgOrder[k])),
    };
  });
  const sigma = bondInfo.length;
  const pi = bondInfo.reduce((s, b) => s + b.pi, 0);
  const totalValence = atoms.reduce((s, a) => s + valenceElectrons(a.Z), 0) - totalCharge;
  return {
    atoms: atomInfo, bonds: bondInfo,
    charge: totalCharge,
    valenceElectrons: totalValence,
    unpaired: lone.reduce((s, x) => s + (x % 2), 0),
    sigmaBonds: sigma, piBonds: pi,
    resonanceCount: Math.max(1, resonance.length),
    resonance: resonance.slice(0, 6),
    orders, lone, formalCharges: fc,
    octetViolations: atomInfo.filter(a => a.Z > 2 && a.electronsAround !== 8).map(a => ({ index: a.index, electrons: a.electronsAround })),
  };
}

function subscript(n) {
  return String(n).split('').map(d => '₀₁₂₃₄₅₆₇₈₉'[+d]).join('');
}
