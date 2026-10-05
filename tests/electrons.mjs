// Viste elettroniche della sandbox: orbitali localizzati (Pipek–Mezey) e ibridazione, elettroni campionati dalla
// densità, densità di legame, griglie limitate alle molecole; più il diagonalizzatore per righe usato da GFN2-xTB.

import { GFN2xTB } from '../js/chem/xtb/gfn2.js';
import { makeGFN2Provider } from '../js/chem/xtb/provider.js';
import { localizedOrbitals, computeGrid, samplePoints, fitGrid } from '../js/chem/densityWorker.js';
import { eigh, eighRows } from '../js/chem/linalg.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};
const BOHR = 0.52917721090;

// geometrie sperimentali (Å)
const tet = (r) => [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(v => v.map(x => x * r / Math.sqrt(3)));
const ring = (r, z = 0) => Array.from({ length: 6 }, (_, k) => [r * Math.cos(k * Math.PI / 3), r * Math.sin(k * Math.PI / 3), z]);
const MOLS = {
  CH4: { Z: [6, 1, 1, 1, 1], xyz: [[0, 0, 0], ...tet(1.09)] },
  NH3: { Z: [7, 1, 1, 1], xyz: [[0, 0, 0.12], [0.94, 0, -0.27], [-0.47, 0.814, -0.27], [-0.47, -0.814, -0.27]] },
  H2O: { Z: [8, 1, 1], xyz: [[0, 0, 0], [0.757, 0.586, 0], [-0.757, 0.586, 0]] },
  C2H4: { Z: [6, 6, 1, 1, 1, 1], xyz: [[0, 0, 0.665], [0, 0, -0.665], [0, 0.92, 1.23], [0, -0.92, 1.23], [0, 0.92, -1.23], [0, -0.92, -1.23]] },
  C2H2: { Z: [6, 6, 1, 1], xyz: [[0, 0, 0.6], [0, 0, -0.6], [0, 0, 1.66], [0, 0, -1.66]] },
  CO2: { Z: [8, 6, 8], xyz: [[0, 0, -1.16], [0, 0, 0], [0, 0, 1.16]] },
  C6H6: { Z: [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1], xyz: [...ring(1.39), ...ring(2.48)] },
};
// atteso: [σ, π, doppietti], ibrido dell'atomo 0 (o 1 per CO2) e intervallo di n = p/s
const EXPECT = {
  CH4: { count: [4, 0, 0], atom: 0, label: 'sp³', n: [2.6, 3.5] },
  NH3: { count: [3, 0, 1], atom: 0, label: 'sp³', n: [2.6, 3.6] },
  H2O: { count: [2, 0, 2], atom: 0, label: 'sp³', n: [2.5, 3.7] },
  C2H4: { count: [5, 1, 0], atom: 0, label: 'sp²', n: [1.5, 2.5] },
  C2H2: { count: [3, 2, 0], atom: 0, label: 'sp', n: [0.65, 1.35] },
  CO2: { count: [2, 4, 2], atom: 1, label: 'sp', n: [0.65, 1.35] },
  C6H6: { count: [12, 3, 0], atom: 0, label: 'sp²', n: [1.5, 2.5] },
};

for (const [name, m] of Object.entries(MOLS)) {
  const p = makeGFN2Provider();
  const pos = Float64Array.from(m.xyz.flat());
  p.compute(m.Z, pos, new Float64Array(pos.length), {});
  const w = p.wavefunction({ occupied: true });
  const r = localizedOrbitals({ Z: w.Z, pos: w.pos, occ: w.occ });
  const cnt = ['sigma', 'pi', 'lone'].map(k => r.orbitals.filter(o => o.kind === k).length);
  const e = EXPECT[name], h = r.atoms[e.atom];
  const ok = cnt.every((c, k) => c === e.count[k]) && h?.label === e.label && h.n >= e.n[0] && h.n <= e.n[1] && r.converged;
  check(`${name}: ${e.count[0]} σ, ${e.count[1]} π, ${e.count[2]} doppietti; atomo ${e.atom + 1} ${e.label}`, ok,
    `trovati ${cnt[0]} σ, ${cnt[1]} π, ${cnt[2]} doppietti; ibrido ${h?.label} (n = ${h?.n.toFixed(2)}), ${r.sweeps} sweep, ${r.meshes.length} superfici`);
}

// la localizzazione è una rotazione unitaria: la densità 2 Σ c cᵀ non cambia
{
  const m = MOLS.C6H6, calc = new GFN2xTB();
  calc.compute(m.Z, Float64Array.from(m.xyz.flat(), v => v / BOHR), { gradient: false });
  const w = calc.last, n = w.n;
  let mo = 0; while (w.f[mo] > 1.5) mo++;
  const { pipekMezey } = await import('../js/chem/xtb/localize.js');
  const aoL = Int32Array.from({ length: n }, (_, mu) => w.basis.shells[w.basis.aoShell[mu]].l);
  const L = pipekMezey(w.Ct, mo, n, w.S, w.basis.aoAtom, aoL, m.Z.length, Float64Array.from(m.xyz.flat()));
  let err = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    let a = 0, b = 0;
    for (let k = 0; k < mo; k++) { a += 2 * w.Ct[k * n + i] * w.Ct[k * n + j]; b += 2 * L.C[k * n + i] * L.C[k * n + j]; }
    err = Math.max(err, Math.abs(a - b));
  }
  // ortonormalità degli orbitali localizzati nella metrica S
  let orth = 0;
  for (let a = 0; a < mo; a++) for (let b = 0; b <= a; b++) {
    let s = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) s += L.C[a * n + i] * w.S[i * n + j] * L.C[b * n + j];
    orth = Math.max(orth, Math.abs(s - (a === b ? 1 : 0)));
  }
  check('benzene: gli orbitali localizzati danno la stessa densità e restano ortonormali', err < 1e-10 && orth < 1e-10, `max |ΔP| ${err.toExponential(1)}, max |CᵀSC − 1| ${orth.toExponential(1)}`);
}

// elettroni campionati, densità di legame e griglia limitata alle molecole
{
  const m = MOLS.H2O, p = makeGFN2Provider();
  const pos = Float64Array.from(m.xyz.flat().map((v, i) => v + [0.031, 0.047, 0.023][i % 3]));
  p.compute(m.Z, pos, new Float64Array(9), {});
  const w = p.wavefunction();
  const fg = fitGrid(w.pos, { maxHalf: 99, step: 0.12, maxRes: 120 });
  const val = computeGrid({ mode: 'xtb', what: 'density', noCore: true, Z: w.Z, pos: w.pos, P: w.P, box: 0, ...fg });
  const s = samplePoints(val, w.pos, 20000, 11);
  let far = 0;
  for (let k = 0; k < s.atom.length; k++) if (Math.hypot(s.pts[3 * k], s.pts[3 * k + 1], s.pts[3 * k + 2]) > 4) far++;
  check('H₂O: 8 elettroni di valenza campionati, punti vicini ai nuclei', Math.abs(s.electrons - 8) < 0.05 && far < 50, `∫ρ_valenza = ${s.electrons.toFixed(3)}, ${far} punti su 20 000 oltre 4 Å dal loro atomo`);
  const dv = (fg.half * 2 / (fg.res - 1) / BOHR) ** 3;
  // solo valenza: il core 1s è così piccato che il suo integrale dipende dal passo della griglia
  const full = computeGrid({ mode: 'xtb', what: 'density', noCore: true, Z: w.Z, pos: w.pos, P: w.P, box: 10, res: 167 });
  const fit = computeGrid({ mode: 'xtb', what: 'density', noCore: true, Z: w.Z, pos: w.pos, P: w.P, box: 0, ...fg });
  const intg = (g) => g.values.reduce((a, v) => a + v, 0) * (g.step / BOHR) ** 3;
  check('griglia limitata alle molecole: stessa carica di valenza della griglia su tutta la scatola', Math.abs(intg(full) - intg(fit)) < 0.01,
    `${intg(fit).toFixed(3)} e (cubo di ${(2 * fg.half).toFixed(1)} Å, ${fg.res}³ punti) contro ${intg(full).toFixed(3)} e (scatola intera, 167³ punti)`);
  const d = computeGrid({ mode: 'xtb', what: 'deformation', Z: w.Z, pos: w.pos, P: w.P, box: 0, ...fg });
  const net = d.values.reduce((a, v) => a + v, 0) * dv;
  // nel legame covalente apolare di H₂ gli elettroni si accumulano a metà fra i nuclei
  const h2 = makeGFN2Provider();
  h2.compute([1, 1], Float64Array.from([0, 0, -0.37, 0, 0, 0.37]), new Float64Array(6), {});
  const w2 = h2.wavefunction();
  const vmid = computeGrid({ mode: 'xtb', what: 'deformation', Z: w2.Z, pos: w2.pos, P: w2.P, box: 0, center: [0, 0, 0], half: 0.01, res: 2 }).values[0];
  check('densità di legame: integrale nullo in H₂O, accumulo a metà del legame di H₂', Math.abs(net) < 0.02 && vmid > 0.03, `∫Δρ = ${net.toFixed(4)} e, Δρ(metà H–H) = ${vmid.toFixed(3)} e/bohr³`);
}

// diagonalizzatore per righe (usato in GFN2-xTB) contro quello per colonne
{
  const n = 90, A = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) { const v = Math.sin(i * 7.3 + j * 1.1) + (i === j ? i * 0.01 : 0); A[i * n + j] = A[j * n + i] = v; }
  const r1 = eigh(A, n), r2 = eighRows(A, n);
  let de = 0, res = 0;
  for (let k = 0; k < n; k++) {
    de = Math.max(de, Math.abs(r1.values[k] - r2.values[k]));
    for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += A[i * n + j] * r2.rows[k * n + j]; res = Math.max(res, Math.abs(s - r2.values[k] * r2.rows[k * n + i])); }
  }
  check('eighRows: stessi autovalori di eigh e residui A·v − λ·v nulli', de < 1e-12 && res < 1e-12, `max |Δλ| ${de.toExponential(1)}, residuo ${res.toExponential(1)}`);
}

if (failures) { console.log(`\n${failures} verifiche delle viste elettroniche fallite`); process.exit(1); }
console.log('\nViste elettroniche: tutte le verifiche superate.');
