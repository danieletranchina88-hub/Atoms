// Worker per il disegno della materia "come è": densità elettronica e orbitali su una griglia 3D.
//
//  • funzione d'onda MINDO/3 (orbitali di Slater di valenza): ρ(r) = Σ P_μν φ_μ(r) φ_ν(r) per la valenza,
//    più la densità degli elettroni di core presa dall'atomo isolato calcolato (DFT-LDA del sito);
//  • funzione d'onda Hartree–Fock (base gaussiana, tutti gli elettroni): ρ(r) = Σ P_μν φ_μ φ_ν;
//  • campo classico: densità promolecolare, somma delle densità degli atomi isolati (DFT-LDA).
// Le griglie sono in Å, centrate nell'origine della scatola; la densità è in e/bohr³, gli orbitali in bohr^(−3/2).

import { computeAtom } from '../physics/atom.js';
import { makeGrid } from '../physics/scf.js';
import { basisOnAtom, nBasis } from './mindo3.js';
import { buildBasis, basisValues } from './integrals.js';
import { marchingCubes } from '../render/marching.js';

const BOHR = 0.52917721090;
const RMAX = 6;        // Å: raggio delle tabelle radiali
const DR = 0.005;      // Å
const CUT = 4.5;       // Å: oltre, le funzioni di valenza sono trascurabili ai valori di isosuperficie usati

// densità radiali degli atomi isolati (totale e di core) su una griglia uniforme in Å
const radial = new Map();
function atomTables(Z) {
  if (radial.has(Z)) return radial.get(Z);
  const a = computeAtom(Z);
  const g = makeGrid(Z);
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
      if (o.n < nVal) rhoC[i] += v;
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
  radial.set(Z, tab);
  return tab;
}
const lookup = (arr, r) => { const x = r / DR; const k = x | 0; if (k >= arr.length - 1) return 0; const t = x - k; return arr[k] + t * (arr[k + 1] - arr[k]); };

/** Elenco di celle per trovare in fretta gli atomi vicini a un punto. */
function cellList(pos, N, half) {
  const nc = Math.max(1, Math.ceil(2 * half / CUT));
  const cells = new Map();
  const cellOf = (x) => Math.min(nc - 1, Math.max(0, Math.floor((x + half) / CUT)));
  for (let A = 0; A < N; A++) {
    const key = cellOf(pos[3 * A]) + nc * (cellOf(pos[3 * A + 1]) + nc * cellOf(pos[3 * A + 2]));
    if (!cells.has(key)) cells.set(key, []);
    cells.get(key).push(A);
  }
  return { nc, cells, cellOf };
}

export function computeGrid({ mode, what, Z, pos, box, res, first, P, orb, basisName, atoms }) {
  const N = Z.length;
  const half = box / 2 + 1.2;
  const step = 2 * half / (res - 1);
  const values = new Float32Array(res * res * res);
  const cl = cellList(pos, N, half);
  const tabs = new Map();
  for (const z of new Set(Z)) tabs.set(z, atomTables(z));
  // base gaussiana (HF): funzioni raggruppate per atomo
  let gb = null, phiG = null;
  if (mode === 'gauss') {
    gb = buildBasis(atoms, basisName);
    phiG = new Float64Array(gb.nbf);
  }
  const near = new Int32Array(N);
  const phi = new Float64Array(4 * N);
  const idxs = new Int32Array(4 * N);
  for (let kz = 0; kz < res; kz++) {
    const z = -half + kz * step;
    const cz = cl.cellOf(z);
    for (let ky = 0; ky < res; ky++) {
      const y = -half + ky * step;
      const cy = cl.cellOf(y);
      for (let kx = 0; kx < res; kx++) {
        const x = -half + kx * step;
        const cx = cl.cellOf(x);
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
        if (mode === 'sto') {
          // funzioni di valenza degli atomi vicini
          let m = 0;
          for (let q = 0; q < nn; q++) {
            const A = near[q];
            basisOnAtom(Z[A], x - pos[3 * A], y - pos[3 * A + 1], z - pos[3 * A + 2], phi, m);
            for (let k = 0; k < nBasis(Z[A]); k++) idxs[m + k] = first[A] + k;
            m += nBasis(Z[A]);
          }
          if (what === 'density') {
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
            for (let q = 0; q < nn; q++) {
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
        if (what === 'density') {
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
  return { values, res, half, step };
}

/**
 * Isosuperfici (marching cubes) della griglia, calcolate qui per non rallentare il disegno.
 * surfaces: [{ iso, sign, colorByAtom }]; colors: rgb per atomo (0–1) per colorare ogni vertice come l'atomo più vicino.
 */
export function computeSurfaces(g, surfaces, pos, colors) {
  const out = [];
  const N = pos.length / 3;
  const cl = colors ? cellList(pos, N, g.half) : null;
  for (const sf of surfaces) {
    const { positions, normals } = marchingCubes(g, sf.iso, sf.sign);
    let col = null;
    if (sf.colorByAtom && colors && positions.length) {
      col = new Float32Array(positions.length);
      for (let v = 0; v < positions.length; v += 3) {
        const x = positions[v], y = positions[v + 1], z = positions[v + 2];
        const cx = cl.cellOf(x), cy = cl.cellOf(y), cz = cl.cellOf(z);
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

// nel browser è un Web Worker; in Node (test) si usano direttamente computeGrid e computeSurfaces
if (typeof self !== 'undefined') self.onmessage = (ev) => {
  const m = ev.data;
  try {
    const g = computeGrid(m);
    if (m.surfaces) {
      const meshes = computeSurfaces(g, m.surfaces, m.pos, m.colors);
      const transfer = meshes.flatMap(x => [x.positions.buffer, x.normals.buffer, ...(x.colors ? [x.colors.buffer] : [])]);
      postMessage({ id: m.id, what: m.what, meshes }, transfer);
    } else postMessage({ id: m.id, what: m.what, ...g }, [g.values.buffer]);
  } catch (e) {
    postMessage({ id: m.id, error: e.message });
  }
};
