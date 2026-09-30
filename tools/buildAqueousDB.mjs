// Converte il database termodinamico MINTEQ v4 (U.S. EPA MINTEQA2 v4, distribuito con PHREEQC dell'USGS)
// nei dati usati dal becher: specie in soluzione e fasi solide, riscritte in funzione delle componenti
// (ioni principali, H⁺ ed elettroni), con log K a 25 °C, entalpia di reazione, parametri di attività e fonti.
//
// Uso: node tools/buildAqueousDB.mjs percorso/minteq.v4.dat > js/chem/aqueousDB.js
// Il file originale: https://github.com/usgs-coupled/phreeqc3/blob/master/database/minteq.v4.dat

import { readFileSync } from 'node:fs';

const file = process.argv[2];
const text = readFileSync(file, 'utf8').split(/\r?\n/);

// componenti (specie principali). Le coppie redox Fe, Mn, Cr, Cu, Co sono accoppiate dagli elettroni;
// N(5)/N(−3), S(6)/S(−2) restano separate (in laboratorio non si convertono l'una nell'altra).
const BASIS = [
  'H+', 'e-', 'Na+', 'K+', 'Li+', 'Ca+2', 'Mg+2', 'Ba+2', 'Sr+2', 'Al+3', 'Zn+2', 'Ni+2', 'Cd+2', 'Pb+2', 'Ag+',
  'Cl-', 'Br-', 'I-', 'F-', 'PO4-3', 'CO3-2', 'SO4-2', 'NO3-', 'NH4+', 'HS-',
  'Fe+3', 'Mn+3', 'CrO4-2', 'Cu+2', 'Co+3',
  'Acetate-', 'Citrate-3', 'Edta-4', 'Formate-', 'Phthalate-2', 'Glycine-', 'Tartarate-2', 'Benzoate-', 'Ethylenediamine', 'Salicylate-2',
];
// master secondari ammessi (stati di ossidazione collegati dagli elettroni)
const SECONDARY = ['Fe+2', 'Mn+2', 'MnO4-', 'MnO4-2', 'Cr(OH)2+', 'Cu+', 'Co+2'];
const ALLOWED = new Set([...BASIS, ...SECONDARY, 'H2O']);

// fasi solide che precipitano davvero in provetta a temperatura ambiente (niente minerali che si formano
// solo in tempi geologici, come ematite, dolomite o apatite ben cristallizzata)
const PHASES_KEEP = [
  'Cerargyrite', 'Bromyrite', 'Iodyrite', 'Ag2CrO4', 'Ag2CO3', 'Ag3PO4', 'Ag2O', 'Ag2SO4', 'Acanthite',
  'Cotunnite', 'PbI2', 'PbBr2', 'Anglesite', 'PbCrO4', 'Cerussite', 'Pb(OH)2', 'Galena', 'PbF2', 'Pb3(PO4)2',
  'Barite', 'Witherite', 'BaCrO4', 'BaF2', 'BaHPO4',
  'Calcite', 'Gypsum', 'Fluorite', 'Portlandite', 'Ca3(PO4)2(beta)', 'CaCrO4',
  'Celestite', 'Strontianite', 'SrF2',
  'Brucite', 'Nesquehonite', 'MgF2',
  'Al(OH)3(am)',
  'Ferrihydrite', 'Fe(OH)2', 'Siderite', 'FeS(ppt)', 'Strengite',
  'Cu(OH)2', 'Malachite', 'Brochantite', 'Covellite', 'Nantokite', 'CuI', 'Cu3(PO4)2',
  'Zn(OH)2(am)', 'ZnS(am)',
  'Ni(OH)2', 'NiCO3', 'NiS(alpha)', 'Ni3(PO4)2',
  'Co(OH)2', 'CoCO3', 'CoS(alpha)', 'Co3(PO4)2',
  'Cd(OH)2(am)', 'Otavite', 'Greenockite',
  'Pyrochroite', 'Rhodochrosite', 'MnS(pnk)', 'Pyrolusite', 'Birnessite',
  'Cr(OH)3(am)',
  'Agmetal', 'Cumetal', 'Pbmetal', 'Znmetal', 'Cdmetal(alpha)',
  'CO2(g)',
];

const parseTerms = (side) => side.trim().split(/\s+\+\s+/).map(t => {
  const m = t.trim().match(/^([0-9]*\.?[0-9]+)?\s*(\S+)$/);
  return [m[2], m[1] ? parseFloat(m[1]) : 1];
});

const KEYWORDS = new Set(['SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES', 'PHASES', 'SURFACE_MASTER_SPECIES', 'SURFACE_SPECIES',
  'EXCHANGE_MASTER_SPECIES', 'EXCHANGE_SPECIES', 'END', 'LLNL_AQUEOUS_MODEL_PARAMETERS', 'RATES', 'SOLUTION', 'NAMED_EXPRESSIONS']);

function readBlocks() {
  const species = new Map();
  const phases = new Map();
  let section = null;
  let cur = null;
  let phaseName = null;
  const finish = () => { cur = null; };
  for (const raw of text) {
    const line = raw.replace(/\t/g, '    ');
    const trimmed = line.trim();
    if (KEYWORDS.has(trimmed) && !line.startsWith(' ')) { section = trimmed; finish(); phaseName = null; continue; }
    if (section === 'SOLUTION_SPECIES') {
      if (trimmed.includes('=') && !trimmed.startsWith('#') && !trimmed.startsWith('-')) {
        const [lhs, rhs] = trimmed.split('#')[0].split('=');
        const R = parseTerms(rhs);
        cur = { name: R[0][0], lhs: parseTerms(lhs), rhs: R, logk: 0, dh: null, an: null, gamma: null, src: {} };
        species.set(cur.name, cur);
        continue;
      }
    } else if (section === 'PHASES') {
      if (!trimmed.includes('=') && /^[A-Za-z0-9(]/.test(trimmed) && !line.startsWith(' ')) { phaseName = trimmed.split(/\s+/)[0]; cur = null; continue; }
      if (trimmed.includes('=') && phaseName && !trimmed.startsWith('#') && !trimmed.startsWith('-')) {
        const [lhs, rhs] = trimmed.split('#')[0].split('=');
        cur = { name: phaseName, lhs: parseTerms(lhs), rhs: parseTerms(rhs), logk: 0, dh: null, an: null, src: {} };
        phases.set(phaseName, cur);
        phaseName = null;
        continue;
      }
    } else continue;
    if (!cur) continue;
    let m;
    if ((m = trimmed.match(/^-?log_k\s+(\S+)/))) cur.logk = parseFloat(m[1]);
    else if ((m = trimmed.match(/^-?delta_h\s+(\S+)\s*(\S*)/i))) {
      let v = parseFloat(m[1]);
      const unit = (m[2] || 'kcal').toLowerCase();
      if (unit.startsWith('kcal')) v *= 4.184;
      cur.dh = v;
    } else if ((m = trimmed.match(/^-analytic(?:al_expression)?\s+(.+)/))) cur.an = m[1].split(/\s+/).map(Number).filter(x => !Number.isNaN(x));
    else if ((m = trimmed.match(/^-gamma\s+(\S+)\s+(\S+)/))) cur.gamma = [parseFloat(m[1]), parseFloat(m[2])];
    else if ((m = trimmed.match(/#\s*log K source:\s*(.*)$/))) cur.src.k = m[1].trim();
    else if ((m = trimmed.match(/#\s*Delta H source:\s*(.*)$/))) cur.src.h = m[1].trim();
  }
  return { species, phases };
}

const { species, phases } = readBlocks();

// composizione in funzione delle componenti, con log K e ΔH accumulati lungo la catena di definizioni
const memo = new Map();
function compose(name, depth = 0) {
  if (memo.has(name)) return memo.get(name);
  if (depth > 20) return null;
  if (BASIS.includes(name)) {
    const r = { c: { [name]: 1 }, w: 0, logk: 0, dh: 0, dhKnown: true };
    memo.set(name, r); return r;
  }
  if (!ALLOWED.has(name) && !species.has(name)) return null;
  const s = species.get(name);
  if (!s) return null;
  // ammessi solo reagenti che siano componenti, master secondari consentiti o acqua
  const terms = [...s.lhs.map(([n, k]) => [n, k]), ...s.rhs.slice(1).map(([n, k]) => [n, -k])];
  const out = { c: {}, w: 0, logk: s.logk, dh: s.dh ?? 0, dhKnown: s.dh !== null };
  for (const [n, k] of terms) {
    if (n === 'H2O') { out.w += k; continue; }
    if (!ALLOWED.has(n)) { memo.set(name, null); return null; }
    if (n === name) { memo.set(name, null); return null; }
    const sub = compose(n, depth + 1);
    if (!sub) { memo.set(name, null); return null; }
    for (const [c, v] of Object.entries(sub.c)) out.c[c] = (out.c[c] ?? 0) + k * v;
    out.w += k * sub.w;
    out.logk += k * sub.logk;
    out.dh += k * sub.dh;
    out.dhKnown = out.dhKnown && sub.dhKnown;
  }
  for (const c of Object.keys(out.c)) if (Math.abs(out.c[c]) < 1e-12) delete out.c[c];
  memo.set(name, out);
  return out;
}

function chargeOf(name) {
  if (name === 'e-') return -1;
  const m = name.match(/([+-])(\d*)$/);
  if (!m) return 0;
  return (m[1] === '+' ? 1 : -1) * (m[2] ? parseInt(m[2], 10) : 1);
}

const outSpecies = [];
for (const name of species.keys()) {
  if (name === 'e-' || name === 'H2O') continue;
  const r = compose(name);
  if (!r) continue;
  if (Object.keys(r.c).includes('e-') && !Object.keys(r.c).some(c => ['Fe+3', 'Mn+3', 'CrO4-2', 'Cu+2', 'Co+3'].includes(c))) continue;
  const s = species.get(name);
  outSpecies.push({
    n: name, z: chargeOf(name), c: r.c, w: +r.w.toFixed(6), lk: +r.logk.toFixed(4),
    dh: r.dhKnown ? +r.dh.toFixed(3) : null,
    an: s?.an?.length && !(BASIS.includes(name)) ? s.an : null,
    g: s?.gamma ?? null,
    src: [s?.src.k, s?.src.h].filter(Boolean).join(' / ') || null,
  });
}
// i master secondari e le componenti devono esserci
for (const b of BASIS) if (b !== 'e-' && !outSpecies.find(s => s.n === b)) {
  const s = species.get(b);
  outSpecies.push({ n: b, z: chargeOf(b), c: { [b]: 1 }, w: 0, lk: 0, dh: 0, an: null, g: s?.gamma ?? null, src: null });
}

const outPhases = [];
for (const name of PHASES_KEEP) {
  const p = phases.get(name);
  if (!p) { console.error(`fase mancante: ${name}`); continue; }
  // reazione di dissoluzione: formula (+ altri reagenti) = prodotti; si esprime in componenti
  const terms = [...p.rhs.map(([n, k]) => [n, k]), ...p.lhs.slice(1).map(([n, k]) => [n, -k])];
  const out = { c: {}, w: 0, logk: p.logk, dh: p.dh ?? 0, dhKnown: p.dh !== null };
  let ok = true;
  for (const [n, k] of terms) {
    if (n === 'H2O') { out.w += k; continue; }
    const sub = compose(n);
    if (!sub) { ok = false; break; }
    for (const [c, v] of Object.entries(sub.c)) out.c[c] = (out.c[c] ?? 0) + k * v;
    out.w += k * sub.w;
    out.logk -= k * sub.logk;
    out.dh -= k * sub.dh;
    out.dhKnown = out.dhKnown && sub.dhKnown;
  }
  if (!ok) { console.error(`fase non componibile: ${name}`); continue; }
  for (const c of Object.keys(out.c)) if (Math.abs(out.c[c]) < 1e-12) delete out.c[c];
  outPhases.push({
    n: name, f: p.lhs[0][0], c: out.c, w: +out.w.toFixed(6),
    lk: +out.logk.toFixed(4), dh: out.dhKnown ? +out.dh.toFixed(3) : null,
    an: p.an?.length ? p.an : null,
    src: [p.src.k, p.src.h].filter(Boolean).join(' / ') || null,
  });
}

process.stdout.write(`// Generato da tools/buildAqueousDB.mjs a partire da minteq.v4.dat (U.S. EPA MINTEQA2 v4, distribuito con
// PHREEQC, USGS: https://github.com/usgs-coupled/phreeqc3/blob/master/database/minteq.v4.dat). Non modificare a mano.
//
// Specie: composizione nelle componenti (c), acqua consumata (w), log K a 25 °C della formazione dalle componenti (lk),
// ΔrH in kJ/mol (dh, null se non disponibile), espressione analitica di log K(T) se presente (an: A1..A5),
// parametri di attività di Debye–Hückel esteso (g: [a Å, b]), fonte bibliografica dei dati (src, codici MINTEQ).
// Fasi: formula (f); log K e ΔH della dissoluzione (solido → componenti), w = acqua consumata quando il solido si forma.

export const DB_COMPONENTS = ${JSON.stringify(BASIS)};
export const DB_SPECIES = [
${outSpecies.map(s => '  ' + JSON.stringify(s)).join(',\n')}
];
export const DB_PHASES = [
${outPhases.map(s => '  ' + JSON.stringify(s)).join(',\n')}
];
`);
console.error(`specie: ${outSpecies.length}, fasi: ${outPhases.length}`);
