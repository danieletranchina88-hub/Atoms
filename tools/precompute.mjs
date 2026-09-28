// Precalcolo della libreria di molecole con il motore quantistico del progetto:
// geometria (VSEPR + distanze) → ottimizzazione HF → frequenze armoniche → termochimica → MP2.
//
// Uso:  node tools/precompute.mjs [id ...]        (senza argomenti: tutte le molecole mancanti)
//       node tools/precompute.mjs --merge          (scrive js/chem/libraryData.js dalla cache)
//       node tools/precompute.mjs --parallel 4     (lancia 4 processi)

import fs from 'node:fs';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { MOLECULES, BARRIERS } from '../js/chem/library.js';
import { parseSmiles } from '../js/chem/smiles.js';
import { analyzeStructure } from '../js/chem/structure.js';
import { embedMolecule } from '../js/chem/embed.js';
import { buildBasis } from '../js/chem/integrals.js';
import { runHF } from '../js/chem/hf.js';
import { optimizeGeometry, optimizeTransitionState } from '../js/chem/optimize.js';
import { harmonicFrequencies, thermochemistry, ISOTOPE_MASS } from '../js/chem/vibrations.js';
import { mp2Energy, dipoleMoment } from '../js/chem/properties.js';
import { analyzeSymmetry } from '../js/chem/symmetry.js';

const ROOT = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const CACHE = path.join(ROOT, 'tools', 'cache');
fs.mkdirSync(CACHE, { recursive: true });
const BOHR_ANG = 0.52917721090;

/** Base usata: 6-31G* fino a 90 funzioni, altrimenti 3-21G. */
export function chooseBasis(atoms) {
  const n = buildBasis(atoms, '6-31G*').nbf;
  return n <= 90 ? '6-31G*' : '3-21G';
}

function compute(mol) {
  const t0 = Date.now();
  const graph = parseSmiles(mol.smiles);
  const charge = graph.atoms.reduce((s, a) => s + a.charge, 0);
  const nel = graph.atoms.reduce((s, a) => s + a.Z, 0) - charge;
  const multiplicity = mol.multiplicity ?? (nel % 2 === 0 ? 1 : 2);
  const st = analyzeStructure({ atoms: graph.atoms, bonds: graph.bonds, charge });
  const log = (...a) => console.log(`[${mol.id}]`, ...a);
  if (nel === 0) {
    // H⁺: nessun elettrone, solo moto traslazionale
    const atoms = [{ Z: 1, xyz: [0, 0, 0] }];
    const th = thermochemistry(atoms, [], 0, { multiplicity: 1 });
    return { id: mol.id, basis: '6-31G*', charge, multiplicity, geometry: [[1, 0, 0, 0]], energy: 0, mp2: 0, freqs: [], thermo: pickThermo(th), pointGroup: 'K_h', sigma: 1, time: 0 };
  }
  let atoms = embedMolecule(graph, st.atoms.map(a => a.lonePairs + a.radical));
  const basis = chooseBasis(atoms);
  const opts = { basis, charge, multiplicity };
  log('base', basis, 'atomi', atoms.length);
  let opt = { atoms, energy: null, converged: true };
  if (atoms.length > 1) {
    // pre-ottimizzazione veloce in STO-3G, poi nella base finale
    const pre = optimizeGeometry(atoms, { ...opts, basis: 'STO-3G', maxSteps: 40 });
    opt = optimizeGeometry(pre.atoms, { ...opts, maxSteps: 80, onStep: s => log('passo', s.step, s.energy.toFixed(8), s.gmax.toExponential(2)) });
    atoms = opt.atoms;
  }
  const res = runHF(atoms, opts);
  const mp2 = mp2Energy(res);
  const mu = dipoleMoment(res);
  let freqs = [];
  let th;
  if (atoms.length > 1) {
    const f = harmonicFrequencies(atoms, { ...opts, onProgress: p => { if (Math.round(p * 20) !== Math.round((p - 0.001) * 20)) log('hessiana', Math.round(p * 100) + '%'); } });
    freqs = f.modes.map(m => ({ freq: +m.freq.toFixed(2), ir: +(m.ir ?? 0).toFixed(3), mu: +m.reducedMass.toFixed(4), d: Array.from(m.displacement, v => +v.toFixed(4)) }));
    th = thermochemistry(atoms, f.modes, res.energy, { multiplicity });
  } else {
    th = thermochemistry(atoms, [], res.energy, { multiplicity });
  }
  const sym = analyzeSymmetry(atoms, atoms.map(a => ISOTOPE_MASS[a.Z]));
  return {
    id: mol.id, basis, charge, multiplicity,
    geometry: atoms.map(a => [a.Z, ...a.xyz.map(v => +(v * BOHR_ANG).toFixed(6))]),
    energy: res.energy, mp2: mp2.energy, converged: res.converged && opt.converged,
    dipole: mu.debye, S2: res.S2,
    freqs, thermo: pickThermo(th), pointGroup: sym.pointGroup, sigma: sym.sigma,
    time: Date.now() - t0,
  };
}

/** Stato di transizione: ricerca P-RFO dalla geometria di partenza, poi frequenze (una immaginaria). */
function computeTS(b) {
  const t0 = Date.now();
  const log = (...a) => console.log(`[ts:${b.id}]`, ...a);
  let atoms = b.tsGuess.map(([Z, x, y, z]) => ({ Z, xyz: [x / BOHR_ANG, y / BOHR_ANG, z / BOHR_ANG] }));
  const basis = '6-31G*';
  const opts = { basis, charge: 0, multiplicity: b.multiplicity };
  const ts = optimizeTransitionState(atoms, { ...opts, onStep: s => log('passo', s.step, s.energy.toFixed(8), s.gmax.toExponential(2)) });
  atoms = ts.atoms;
  const res = runHF(atoms, opts);
  const mp2 = mp2Energy(res);
  const f = harmonicFrequencies(atoms, opts);
  const freqs = f.modes.map(m => ({ freq: +m.freq.toFixed(2), ir: +(m.ir ?? 0).toFixed(3), mu: +m.reducedMass.toFixed(4), d: Array.from(m.displacement, v => +v.toFixed(4)) }));
  const th = thermochemistry(atoms, f.modes, res.energy, { multiplicity: b.multiplicity });
  const sym = analyzeSymmetry(atoms, atoms.map(a => ISOTOPE_MASS[a.Z]));
  return {
    id: `ts-${b.id}`, basis, charge: 0, multiplicity: b.multiplicity, ts: true,
    geometry: atoms.map(a => [a.Z, ...a.xyz.map(v => +(v * BOHR_ANG).toFixed(6))]),
    energy: res.energy, mp2: mp2.energy, converged: ts.converged, S2: res.S2,
    imaginary: f.modes.filter(m => m.freq < 0).length,
    freqs, thermo: pickThermo(th), pointGroup: sym.pointGroup, sigma: sym.sigma,
    time: Date.now() - t0,
  };
}

/**
 * Se una geometria ha frequenze immaginarie non è un minimo ma un punto di sella:
 * si sposta la struttura lungo il modo immaginario e si riottimizza, finché tutte le frequenze sono reali.
 */
function refineMinimum(id) {
  const mol = MOLECULES.find(m => m.id === id);
  const d = JSON.parse(fs.readFileSync(path.join(CACHE, `${id}.json`), 'utf8'));
  const log = (...a) => console.log(`[${id}]`, ...a);
  let atoms = d.geometry.map(([Z, x, y, z]) => ({ Z, xyz: [x / BOHR_ANG, y / BOHR_ANG, z / BOHR_ANG] }));
  const opts = { basis: d.basis, charge: d.charge, multiplicity: d.multiplicity };
  let freqs = d.freqs;
  for (let round = 0; round < 4; round++) {
    const imag = freqs.filter(f => f.freq < 0);
    if (!imag.length) break;
    log('frequenze immaginarie', imag.map(f => f.freq.toFixed(0)).join(', '), '→ spostamento lungo il modo');
    const mode = imag[0].d;
    atoms = atoms.map((a, i) => ({ ...a, xyz: a.xyz.map((v, k) => v + 0.3 * mode[3 * i + k]) }));
    const opt = optimizeGeometry(atoms, { ...opts, maxSteps: 100, gmax: 1.5e-4, grms: 1e-4, onStep: s => log('passo', s.step, s.energy.toFixed(8), s.gmax.toExponential(2)) });
    atoms = opt.atoms;
    const f = harmonicFrequencies(atoms, opts);
    freqs = f.modes.map(m => ({ freq: +m.freq.toFixed(2), ir: +(m.ir ?? 0).toFixed(3), mu: +m.reducedMass.toFixed(4), d: Array.from(m.displacement, v => +v.toFixed(4)) }));
  }
  const res = runHF(atoms, opts);
  const th = thermochemistry(atoms, freqs.map(f => ({ freq: f.freq })), res.energy, { multiplicity: d.multiplicity });
  const sym = analyzeSymmetry(atoms, atoms.map(a => ISOTOPE_MASS[a.Z]));
  const out = {
    ...d, geometry: atoms.map(a => [a.Z, ...a.xyz.map(v => +(v * BOHR_ANG).toFixed(6))]),
    energy: res.energy, mp2: mp2Energy(res).energy, dipole: dipoleMoment(res).debye,
    freqs, thermo: pickThermo(th), pointGroup: sym.pointGroup, sigma: sym.sigma,
  };
  fs.writeFileSync(path.join(CACHE, `${id}.json`), JSON.stringify(out));
  log('fatto:', sym.pointGroup, 'E', res.energy.toFixed(8), 'immaginarie rimaste', freqs.filter(f => f.freq < 0).length);
  void mol;
}

function pickThermo(th) {
  return { zpe: th.zpe, Hcorr: th.Hcorr, Gcorr: th.Gcorr, S: th.S, Cv: th.Cv, T: th.T };
}

function merge() {
  const data = {};
  for (const m of [...MOLECULES, ...BARRIERS.map(b => ({ id: `ts-${b.id}` }))]) {
    const f = path.join(CACHE, `${m.id}.json`);
    if (!fs.existsSync(f)) continue;
    const d = JSON.parse(fs.readFileSync(f, 'utf8'));
    // gruppo puntuale e σ ricalcolati con il codice attuale
    if (d.geometry.length > 1) {
      const atoms = d.geometry.map(([Z, x, y, z]) => ({ Z, xyz: [x / BOHR_ANG, y / BOHR_ANG, z / BOHR_ANG] }));
      const sym = analyzeSymmetry(atoms, atoms.map(a => ISOTOPE_MASS[a.Z]));
      d.pointGroup = sym.pointGroup;
      d.sigma = sym.sigma;
    }
    data[m.id] = d;
  }
  const out = `// Generato da tools/precompute.mjs con il motore quantistico del progetto (HF + MP2).\n// Geometrie ottimizzate in Å, energie in hartree, frequenze armoniche in cm⁻¹ (non scalate).\nexport const LIBRARY_DATA = ${JSON.stringify(data)};\n`;
  fs.writeFileSync(path.join(ROOT, 'js', 'chem', 'libraryData.js'), out);
  console.log(`libraryData.js: ${Object.keys(data).length} molecole`);
}

const args = process.argv.slice(2);
if (args[0] === '--refine') {
  for (const id of args.slice(1)) refineMinimum(id);
  merge();
} else if (args[0] === '--merge') {
  merge();
} else if (args[0] === '--parallel') {
  const nproc = +args[1] || 4;
  const todo = MOLECULES.filter(m => !fs.existsSync(path.join(CACHE, `${m.id}.json`))).map(m => m.id);
  // ordina per dimensione (le più grandi per prime) e distribuisci
  const size = (id) => parseSmiles(MOLECULES.find(m => m.id === id).smiles).atoms.length;
  todo.sort((a, b) => size(b) - size(a));
  const groups = Array.from({ length: nproc }, () => []);
  todo.forEach((id, k) => groups[k % nproc].push(id));
  let running = 0;
  for (const g of groups) {
    if (!g.length) continue;
    running++;
    const p = spawn(process.execPath, [fileURLToPath(import.meta.url), ...g], { stdio: 'inherit' });
    p.on('exit', () => { if (--running === 0) merge(); });
  }
} else {
  const ids = args.length ? args : MOLECULES.map(m => m.id);
  for (const id of ids) {
    if (id.startsWith('ts-')) {
      const b = BARRIERS.find(x => `ts-${x.id}` === id);
      try {
        const r = computeTS(b);
        fs.writeFileSync(path.join(CACHE, `${id}.json`), JSON.stringify(r));
        console.log(`[${id}] fatto in ${(r.time / 1000).toFixed(1)} s, E = ${r.energy.toFixed(8)}, ${r.imaginary} frequenze immaginarie`);
      } catch (e) {
        console.error(`[${id}] ERRORE`, e.stack);
      }
      continue;
    }
    const mol = MOLECULES.find(m => m.id === id);
    if (!mol) { console.error('molecola sconosciuta', id); continue; }
    try {
      const r = compute(mol);
      fs.writeFileSync(path.join(CACHE, `${id}.json`), JSON.stringify(r));
      console.log(`[${id}] fatto in ${(r.time / 1000).toFixed(1)} s, E = ${r.energy.toFixed(8)}, ${r.pointGroup}`);
    } catch (e) {
      console.error(`[${id}] ERRORE`, e.stack);
    }
  }
}
