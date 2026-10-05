// Worker per il disegno della materia "come è": densità elettronica e orbitali su una griglia 3D.
//
//  • funzione d'onda MINDO/3 (orbitali di Slater di valenza): ρ(r) = Σ P_μν φ_μ(r) φ_ν(r) per la valenza,
//    più la densità degli elettroni di core presa dall'atomo isolato calcolato (DFT-LDA del sito);
//  • funzione d'onda GFN2-xTB (base STO-nG di valenza): come MINDO/3, valenza più core dell'atomo isolato;
//  • funzione d'onda Hartree–Fock (base gaussiana, tutti gli elettroni): ρ(r) = Σ P_μν φ_μ φ_ν;
//  • campo classico: densità promolecolare, somma delle densità degli atomi isolati (DFT-LDA).
// Le griglie sono in Å, centrate nell'origine della scatola; la densità è in e/bohr³, gli orbitali in bohr^(−3/2).

import { computeAtom } from '../physics/atom.js';
import { makeGrid } from '../physics/scf.js';
import { basisOnAtom, nBasis } from './mindo3.js';
import { buildBasis, basisValues } from './integrals.js';
import { buildBasis as xtbBasis, shellValues, nSph, moleculeIntegrals } from './xtb/gto.js';
import { pipekMezey, hybridLabel, atomHybridization } from './xtb/localize.js';
import { GFN2 } from './xtb/gfn2Data.js';
import { marchingCubes } from '../render/marching.js';

const BOHR = 0.52917721090;
const RMAX = 6;        // Å: raggio delle tabelle radiali
const DR = 0.005;      // Å
const CUT = 4.5;       // Å: oltre, le funzioni di valenza sono trascurabili ai valori di isosuperficie usati

// densità radiali degli atomi isolati (totale e di core) su una griglia uniforme in Å
const radial = new Map();
// core = orbitali dell'atomo isolato più interni del guscio di valenza (MINDO/3), oppure tutti quelli che
// non sono shell della base GFN2 (per i metalli di transizione il 3d è valenza anche se n è minore)
function atomTables(Z, xtb = false) {
  const key = xtb ? -Z : Z;
  if (radial.has(key)) return radial.get(key);
  const xtbShells = xtb ? new Set(GFN2.elements[Z - 1].shells.map(sh => `${sh.n},${sh.l}`)) : null;
  const a = computeAtom(Z);
  const g = makeGrid(Z, { xmin: a.grid.xmin, h: a.grid.h });
  const nVal = Math.max(...a.orbitals.map(o => o.n));
  const nPts = Math.ceil(RMAX / DR) + 1;
  const total = new Float32Array(nPts), core = new Float32Array(nPts);
  // ρ(r) = Σ occ u²/(4π r²) sulla griglia logaritmica, poi interpolazione lineare
  const rhoT = new Float64Array(g.N), rhoC = new Float64Array(g.N);
  for (const o of a.orbitals) {
    for (let i = 0; i < g.N; i++) {
      const r = g.r[i];
      const v = o.occ * o.u[i] * o.u[i] / (4 * Math.PI * r * r);
      rhoT[i] += v;
      if (xtb ? !xtbShells.has(`${o.n},${o.l}`) : o.n < nVal) rhoC[i] += v;
    }
  }
  let j = 0;
  for (let k = 0; k < nPts; k++) {
    const r = Math.max(k * DR / BOHR, g.r[0]);
    while (j < g.N - 2 && g.r[j + 1] < r) j++;
    const t = Math.min(1, Math.max(0, (r - g.r[j]) / (g.r[j + 1] - g.r[j])));
    total[k] = rhoT[j] + t * (rhoT[j + 1] - rhoT[j]);
    core[k] = rhoC[j] + t * (rhoC[j + 1] - rhoC[j]);
  }
  const tab = { total, core };
  radial.set(key, tab);
  return tab;
}
const lookup = (arr, r) => { const x = r / DR; const k = x | 0; if (k >= arr.length - 1) return 0; const t = x - k; return arr[k] + t * (arr[k + 1] - arr[k]); };

/** Elenco di celle per trovare in fretta gli atomi vicini a un punto (cubo di lato 2·half centrato in center). */
function cellList(pos, N, half, center = [0, 0, 0]) {
  const nc = Math.max(1, Math.ceil(2 * half / CUT));
  const cells = new Map();
  const cellOf = (x, ax) => Math.min(nc - 1, Math.max(0, Math.floor((x - center[ax] + half) / CUT)));
  for (let A = 0; A < N; A++) {
    const key = cellOf(pos[3 * A], 0) + nc * (cellOf(pos[3 * A + 1], 1) + nc * cellOf(pos[3 * A + 2], 2));
    if (!cells.has(key)) cells.set(key, []);
    cells.get(key).push(A);
  }
  return { nc, cells, cellOf };
}

function singleGrid({ mode, what, Z, pos, box, res, first, P, orb, basisName, atoms, Ps, center = [0, 0, 0], half = box / 2 + 1.2, noCore = false }) {
  if (what === 'spin') P = Ps;
  if (what === 'deformation') {
    if (mode !== 'xtb') throw new Error('La densità di legame è calcolata nella base di GFN2-xTB.');
    P = deformationMatrix(P, Z);
  }
  const isDensity = what === 'density' || what === 'spin' || what === 'deformation';
  const withCore = what === 'density' && !noCore;
  if (isDensity && mode !== 'promolecular' && !P) throw new Error('Matrice densità non disponibile.');
  const N = Z.length;
  const step = 2 * half / (res - 1);
  const values = new Float32Array(res * res * res);
  const cl = cellList(pos, N, half, center);
  const tabs = new Map();
  if (withCore && mode !== 'gauss' || mode === 'promolecular')
    for (const z of new Set(Z)) tabs.set(z, atomTables(z, mode === 'xtb'));
  // base GFN2: shell per atomo, indici delle funzioni
  let xb = null;
  if (mode === 'xtb') xb = xtbBasis(Z, GFN2.elements);
  // base gaussiana (HF): funzioni raggruppate per atomo
  let gb = null, phiG = null;
  if (mode === 'gauss') {
    gb = buildBasis(atoms, basisName);
    phiG = new Float64Array(gb.nbf);
  }
  const near = new Int32Array(N);
  const xphi = new Float64Array(9 * N), xidx = new Int32Array(9 * N);
  const phi = new Float64Array(4 * N);
  const idxs = new Int32Array(4 * N);
  for (let kz = 0; kz < res; kz++) {
    const z = center[2] - half + kz * step;
    const cz = cl.cellOf(z, 2);
    for (let ky = 0; ky < res; ky++) {
      const y = center[1] - half + ky * step;
      const cy = cl.cellOf(y, 1);
      for (let kx = 0; kx < res; kx++) {
        const x = center[0] - half + kx * step;
        const cx = cl.cellOf(x, 0);
        // atomi entro il raggio di taglio
        let nn = 0;
        for (let dz = -1; dz <= 1; dz++) for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
          const a = cx + dx, b = cy + dy, c = cz + dz;
          if (a < 0 || b < 0 || c < 0 || a >= cl.nc || b >= cl.nc || c >= cl.nc) continue;
          const list = cl.cells.get(a + cl.nc * (b + cl.nc * c));
          if (!list) continue;
          for (const A of list) {
            const ddx = x - pos[3 * A], ddy = y - pos[3 * A + 1], ddz = z - pos[3 * A + 2];
            if (ddx * ddx + ddy * ddy + ddz * ddz < CUT * CUT) near[nn++] = A;
          }
        }
        const gi = kx + res * (ky + res * kz);
        if (!nn) continue;
        if (mode === 'promolecular') {
          let rho = 0;
          for (let q = 0; q < nn; q++) {
            const A = near[q];
            rho += lookup(tabs.get(Z[A]).total, Math.hypot(x - pos[3 * A], y - pos[3 * A + 1], z - pos[3 * A + 2]));
          }
          values[gi] = rho;
          continue;
        }
        if (mode === 'xtb') {
          let m = 0;
          for (let q = 0; q < nn; q++) {
            const A = near[q];
            const dx = (x - pos[3 * A]) / BOHR, dy = (y - pos[3 * A + 1]) / BOHR, dz = (z - pos[3 * A + 2]) / BOHR;
            for (const k of xb.atomShells[A]) {
              const sh = xb.shells[k];
              shellValues(sh, dx, dy, dz, xphi, m);
              for (let c = 0; c < nSph(sh.l); c++) xidx[m + c] = sh.ao + c;
              m += nSph(sh.l);
            }
          }
          if (isDensity) {
            const nb = xb.nao;
            let rho = 0;
            for (let a = 0; a < m; a++) {
              const pa = xphi[a];
              if (pa === 0) continue;
              const row = xidx[a] * nb;
              let s = 0;
              for (let b = 0; b < m; b++) s += P[row + xidx[b]] * xphi[b];
              rho += pa * s;
            }
            if (withCore) for (let q = 0; q < nn; q++) {
              const A = near[q];
              rho += lookup(tabs.get(Z[A]).core, Math.hypot(x - pos[3 * A], y - pos[3 * A + 1], z - pos[3 * A + 2]));
            }
            values[gi] = rho;
          } else {
            let s = 0;
            for (let a = 0; a < m; a++) s += orb[xidx[a]] * xphi[a];
            values[gi] = s;
          }
          continue;
        }
        if (mode === 'sto') {
          // funzioni di valenza degli atomi vicini
          let m = 0;
          for (let q = 0; q < nn; q++) {
            const A = near[q];
            basisOnAtom(Z[A], x - pos[3 * A], y - pos[3 * A + 1], z - pos[3 * A + 2], phi, m);
            for (let k = 0; k < nBasis(Z[A]); k++) idxs[m + k] = first[A] + k;
            m += nBasis(Z[A]);
          }
          if (isDensity) {
            const nb = first[N];
            let rho = 0;
            for (let a = 0; a < m; a++) {
              const pa = phi[a];
              if (pa === 0) continue;
              const row = idxs[a] * nb;
              let s = 0;
              for (let b = 0; b < m; b++) s += P[row + idxs[b]] * phi[b];
              rho += pa * s;
            }
            if (withCore) for (let q = 0; q < nn; q++) {
              const A = near[q];
              rho += lookup(tabs.get(Z[A]).core, Math.hypot(x - pos[3 * A], y - pos[3 * A + 1], z - pos[3 * A + 2]));
            }
            values[gi] = rho;
          } else {
            let s = 0;
            for (let a = 0; a < m; a++) s += orb[idxs[a]] * phi[a];
            values[gi] = s;
          }
          continue;
        }
        // base gaussiana: tutte le funzioni (sistemi piccoli)
        basisValues(gb, x / BOHR, y / BOHR, z / BOHR, phiG);
        const nb = gb.nbf;
        if (isDensity) {
          let rho = 0;
          for (let a = 0; a < nb; a++) {
            const pa = phiG[a];
            if (Math.abs(pa) < 1e-12) continue;
            let s = 0;
            for (let b = 0; b < nb; b++) s += P[a * nb + b] * phiG[b];
            rho += pa * s;
          }
          values[gi] = rho;
        } else {
          let s = 0;
          for (let a = 0; a < nb; a++) s += orb[a] * phiG[a];
          values[gi] = s;
        }
      }
    }
  }
  return { values, res, half, step, center };
}

/**
 * Matrice P − P⁰, con P⁰ la densità degli atomi neutri sferici nella stessa base (occupazioni di riferimento della
 * shell distribuite in modo uguale sulle sue 2l+1 funzioni): la differenza è la "densità di deformazione".
 */
function deformationMatrix(P, Z) {
  const b = xtbBasis(Z, GFN2.elements), n = b.nao, D = Float64Array.from(P);
  for (const sh of b.shells) {
    const k = nSph(sh.l), occ = sh.ref.refocc / k;
    for (let m = 0; m < k; m++) D[(sh.ao + m) * n + sh.ao + m] -= occ;
  }
  return D;
}

/** Difference of two stored SCF densities on exactly the same laboratory grid. */
export function computeGrid(m) {
  if (m.what !== 'difference') return singleGrid(m);
  const r = m.reference;
  if (!r || r.mode !== m.mode || r.Z.join(',') !== m.Z.join(',')) throw new Error('Fissa un riferimento con gli stessi atomi e lo stesso modello.');
  // deform: differenza delle densità di deformazione (ρ − ρ_atomi): toglie la parte degli atomi che si sposta
  // rigidamente con i nuclei e lascia la ridistribuzione degli elettroni fra legami e doppietti
  const what = m.deform && m.mode === 'xtb' ? 'deformation' : 'density';
  const current = singleGrid({ ...m, what });
  const baseline = singleGrid({ ...r, box: m.box, res: m.res, center: m.center, half: m.half, noCore: m.noCore, what });
  for (let k = 0; k < current.values.length; k++) current.values[k] -= baseline.values[k];
  return current;
}

export function gridSlice(g, { axis = 'z', offset = 0 } = {}) {
  if (!['x', 'y', 'z'].includes(axis) || !Number.isFinite(offset)) throw new Error('Sezione non valida.');
  const ax = { x: 0, y: 1, z: 2 }[axis], c = g.center ?? [0, 0, 0];
  const k = Math.max(0, Math.min(g.res - 1, Math.round((offset - c[ax] + g.half) / g.step)));
  const values = new Float32Array(g.res * g.res);
  for (let v = 0; v < g.res; v++) for (let u = 0; u < g.res; u++) {
    const [x, y, z] = axis === 'x' ? [k, u, v] : axis === 'y' ? [u, k, v] : [u, v, k];
    values[u + g.res * v] = g.values[x + g.res * (y + g.res * z)];
  }
  return { values, res: g.res, half: g.half, axis, offset: c[ax] - g.half + k * g.step, center: c };
}

/**
 * Isosuperfici (marching cubes) della griglia, calcolate qui per non rallentare il disegno.
 * surfaces: [{ iso, sign, colorByAtom }]; colors: rgb per atomo (0–1) per colorare ogni vertice come l'atomo più vicino.
 */
export function computeSurfaces(g, surfaces, pos, colors) {
  const out = [];
  const N = pos.length / 3;
  const c = g.center ?? [0, 0, 0];
  const cl = colors ? cellList(pos, N, g.half, c) : null;
  for (const sf of surfaces) {
    const { positions, normals } = marchingCubes(g, sf.iso, sf.sign);
    for (let v = 0; v < positions.length; v += 3) { positions[v] += c[0]; positions[v + 1] += c[1]; positions[v + 2] += c[2]; }
    let col = null;
    if (sf.colorByAtom && colors && positions.length) {
      col = new Float32Array(positions.length);
      for (let v = 0; v < positions.length; v += 3) {
        const x = positions[v], y = positions[v + 1], z = positions[v + 2];
        const cx = cl.cellOf(x, 0), cy = cl.cellOf(y, 1), cz = cl.cellOf(z, 2);
        let best = -1, bd = Infinity;
        for (let dz = -1; dz <= 1; dz++) for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
          const list = cl.cells.get((cx + dx) + cl.nc * ((cy + dy) + cl.nc * (cz + dz)));
          if (!list || cx + dx < 0 || cy + dy < 0 || cz + dz < 0 || cx + dx >= cl.nc || cy + dy >= cl.nc || cz + dz >= cl.nc) continue;
          for (const A of list) {
            const d = (x - pos[3 * A]) ** 2 + (y - pos[3 * A + 1]) ** 2 + (z - pos[3 * A + 2]) ** 2;
            if (d < bd) { bd = d; best = A; }
          }
        }
        if (best < 0) { col[v] = col[v + 1] = col[v + 2] = 0.8; continue; }
        col[v] = colors[3 * best]; col[v + 1] = colors[3 * best + 1]; col[v + 2] = colors[3 * best + 2];
      }
    }
    out.push({ positions, normals, colors: col, iso: sf.iso, sign: sf.sign });
  }
  return out;
}

/**
 * Cubo che contiene gli atomi con un margine (Å), al passo voluto: evita di calcolare la densità nel vuoto della
 * scatola. Restituisce { center, half, res }.
 */
export function fitGrid(pos, { margin = 3.2, step = 0.24, maxHalf = Infinity, minRes = 24, maxRes = 72 } = {}) {
  const N = pos.length / 3, lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let A = 0; A < N; A++) for (let c = 0; c < 3; c++) { lo[c] = Math.min(lo[c], pos[3 * A + c]); hi[c] = Math.max(hi[c], pos[3 * A + c]); }
  const center = [0, 1, 2].map(c => (lo[c] + hi[c]) / 2);
  const half = Math.min(maxHalf, Math.max(...[0, 1, 2].map(c => (hi[c] - lo[c]) / 2)) + margin);
  const res = Math.max(minRes, Math.min(maxRes, Math.round(2 * half / step) + 1));
  return { center, half, res };
}

/**
 * Campionamento della densità (regola di Born): K punti estratti con probabilità ρ(r)·dV, ciascuno attaccato
 * all'atomo più vicino perché nel disegno segua il nucleo fra un calcolo e il successivo.
 */
export function samplePoints(g, pos, K, seed = 1) {
  const n = g.values.length, cum = new Float64Array(n);
  let tot = 0;
  for (let k = 0; k < n; k++) { tot += Math.max(0, g.values[k]); cum[k] = tot; }
  const dv = (g.step / BOHR) ** 3, electrons = tot * dv;
  const pts = new Float32Array(3 * K), atom = new Int32Array(K);
  if (!(tot > 0)) return { pts: new Float32Array(0), atom: new Int32Array(0), electrons: 0 };
  let s = seed >>> 0 || 1;
  const rnd = () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
  const N = pos.length / 3, res = g.res, c = g.center ?? [0, 0, 0];
  for (let p = 0; p < K; p++) {
    const u = rnd() * tot;
    let lo = 0, hi = n - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (cum[mid] < u) lo = mid + 1; else hi = mid; }
    const kx = lo % res, ky = ((lo / res) | 0) % res, kz = (lo / (res * res)) | 0;
    const x = c[0] - g.half + (kx + rnd() - 0.5) * g.step, y = c[1] - g.half + (ky + rnd() - 0.5) * g.step, z = c[2] - g.half + (kz + rnd() - 0.5) * g.step;
    let best = 0, bd = Infinity;
    for (let A = 0; A < N; A++) { const d = (x - pos[3 * A]) ** 2 + (y - pos[3 * A + 1]) ** 2 + (z - pos[3 * A + 2]) ** 2; if (d < bd) { bd = d; best = A; } }
    pts[3 * p] = x - pos[3 * best]; pts[3 * p + 1] = y - pos[3 * best + 1]; pts[3 * p + 2] = z - pos[3 * best + 2];
    atom[p] = best;
  }
  return { pts, atom, electrons };
}

/**
 * Orbitali localizzati di Pipek–Mezey dalla funzione d'onda GFN2 (orbitali doppiamente occupati) e loro
 * isosuperfici, ciascuna su una piccola griglia attorno ai suoi centri.
 */
export function localizedOrbitals({ Z, pos, occ, iso = 0.1, maxOrbitals = 80 }) {
  const b = xtbBasis(Z, GFN2.elements), n = b.nao;
  if (!occ || occ.n !== n) throw new Error('Orbitali occupati non disponibili per questo fotogramma.');
  const { S } = moleculeIntegrals(b, Float64Array.from(pos, v => v / BOHR));
  const aoL = Int32Array.from({ length: n }, (_, mu) => b.shells[b.aoShell[mu]].l);
  const L = pipekMezey(occ.Ct, occ.m, n, S, b.aoAtom, aoL, Z.length, pos);
  const ORDER = { sigma: 0, pi: 1, multi: 2, lone: 3 };
  const list = L.orbitals.map(o => ({
    kind: o.kind, label: o.label, index: o.index,
    centers: o.centers.map(c => ({ atom: c.atom, pop: c.pop })),
    hybrids: o.hybrids.map(h => ({ atom: h.atom, label: hybridLabel(h), s: h.s, p: h.p, d: h.d, pop: h.pop })),
    key: `${o.kind}:${o.centers.slice(0, 3).map(c => c.atom).sort((x, y) => x - y).join('-')}`,
  }));
  // chiavi uniche (due π fra gli stessi atomi): numero progressivo
  const seen = new Map();
  for (const o of list) { const k = seen.get(o.key) ?? 0; seen.set(o.key, k + 1); o.key += `#${k}`; }
  list.sort((x, y) => ORDER[x.kind] - ORDER[y.kind] || x.key.localeCompare(y.key));
  const meshes = [];
  for (const [rank, o] of list.slice(0, maxOrbitals).entries()) {
    const at = o.centers.filter(c => c.pop > 0.05).map(c => c.atom);
    const sub = Float64Array.from(at.flatMap(A => [pos[3 * A], pos[3 * A + 1], pos[3 * A + 2]]));
    const fg = fitGrid(sub, { margin: 2.6, step: 0.2, minRes: 20, maxRes: 44 });
    const orb = L.C.subarray(o.index * n, (o.index + 1) * n);
    const g = singleGrid({ mode: 'xtb', what: 'orbital', Z, pos, box: 0, res: fg.res, center: fg.center, half: fg.half, orb });
    for (const sf of computeSurfaces(g, [{ iso, sign: 1 }, { iso, sign: -1 }], pos, null)) {
      if (sf.positions.length) meshes.push({ rank, key: o.key, sign: sf.sign, positions: sf.positions, normals: sf.normals });
    }
  }
  const atoms = Z.map((_, A) => atomHybridization(L.orbitals, A));
  return { orbitals: list, meshes, atoms, sweeps: L.sweeps, converged: L.converged };
}

// nel browser è un Web Worker; in Node (test) si usano direttamente computeGrid e computeSurfaces
if (typeof self !== 'undefined') self.onmessage = (ev) => {
  const m = ev.data;
  try {
    if (m.what === 'lmo') {
      const r = localizedOrbitals(m);
      const transfer = r.meshes.flatMap(x => [x.positions.buffer, x.normals.buffer]);
      postMessage({ id: m.id, what: 'lmo', ...r }, transfer);
      return;
    }
    // griglia limitata alla zona delle molecole (con il margine della superficie a 0,002 e/bohr³)
    if (m.fit) Object.assign(m, fitGrid(m.pos, { maxHalf: m.box / 2 + 1.2 }));
    if (m.what === 'points') {
      const g = computeGrid({ ...m, what: 'density' });
      const r = samplePoints(g, m.pos, m.count ?? 20000, m.seed);
      postMessage({ id: m.id, what: 'points', ...r }, [r.pts.buffer, r.atom.buffer]);
      return;
    }
    const g = computeGrid(m);
    if (m.surfaces) {
      const meshes = computeSurfaces(g, m.surfaces, m.pos, m.colors);
      const transfer = meshes.flatMap(x => [x.positions.buffer, x.normals.buffer, ...(x.colors ? [x.colors.buffer] : [])]);
      const slice = m.slice ? gridSlice(g, m.slice) : null;
      if (slice) transfer.push(slice.values.buffer);
      postMessage({ id: m.id, what: m.what, meshes, slice }, transfer);
    } else postMessage({ id: m.id, what: m.what, ...g }, [g.values.buffer]);
  } catch (e) {
    postMessage({ id: m.id, error: e.message });
  }
};
