// Marching cubes (Lorensen & Cline, 1987) per estrarre l'isosuperficie ψ = ±soglia.
// Le tabelle dei casi sono quelle di Paul Bourke, esportate da three.js.

import { edgeTable, triTable } from 'three/addons/objects/MarchingCubes.js';

const CORNERS = [
  [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
];
const EDGES = [
  [0, 1], [1, 2], [2, 3], [3, 0],
  [4, 5], [5, 6], [6, 7], [7, 4],
  [0, 4], [1, 5], [2, 6], [3, 7],
];

/**
 * @param {Float32Array} values campo scalare res³ (indice i + res·(j + res·k))
 * @param {number} sign +1 per il lobo positivo (ψ = +iso), −1 per il negativo (ψ = −iso)
 * @returns {{positions: Float32Array, normals: Float32Array}} in coordinate reali
 */
export function marchingCubes({ values, res, half, step }, iso, sign = 1) {
  const f = (i, j, k) => sign * values[i + res * (j + res * k)];
  const grad = (i, j, k, out) => {
    const i0 = Math.max(i - 1, 0), i1 = Math.min(i + 1, res - 1);
    const j0 = Math.max(j - 1, 0), j1 = Math.min(j + 1, res - 1);
    const k0 = Math.max(k - 1, 0), k1 = Math.min(k + 1, res - 1);
    out[0] = (f(i1, j, k) - f(i0, j, k)) / (i1 - i0);
    out[1] = (f(i, j1, k) - f(i, j0, k)) / (j1 - j0);
    out[2] = (f(i, j, k1) - f(i, j, k0)) / (k1 - k0);
  };

  const pos = [];
  const nor = [];
  const cv = new Float64Array(8);
  const cg = Array.from({ length: 8 }, () => [0, 0, 0]);
  const vertList = Array.from({ length: 12 }, () => [0, 0, 0, 0, 0, 0]);

  for (let k = 0; k < res - 1; k++) {
    for (let j = 0; j < res - 1; j++) {
      for (let i = 0; i < res - 1; i++) {
        let cubeIndex = 0;
        for (let c = 0; c < 8; c++) {
          const [di, dj, dk] = CORNERS[c];
          cv[c] = f(i + di, j + dj, k + dk);
          if (cv[c] < iso) cubeIndex |= 1 << c;
        }
        const edges = edgeTable[cubeIndex];
        if (edges === 0) continue;
        for (let c = 0; c < 8; c++) {
          const [di, dj, dk] = CORNERS[c];
          grad(i + di, j + dj, k + dk, cg[c]);
        }
        for (let e = 0; e < 12; e++) {
          if (!(edges & (1 << e))) continue;
          const [a, b] = EDGES[e];
          const va = cv[a];
          const vb = cv[b];
          const t = Math.abs(vb - va) > 1e-30 ? (iso - va) / (vb - va) : 0.5;
          const pa = CORNERS[a];
          const pb = CORNERS[b];
          const v = vertList[e];
          v[0] = -half + (i + pa[0] + t * (pb[0] - pa[0])) * step;
          v[1] = -half + (j + pa[1] + t * (pb[1] - pa[1])) * step;
          v[2] = -half + (k + pa[2] + t * (pb[2] - pa[2])) * step;
          // normale uscente: opposta al gradiente del campo (il campo cresce verso l'interno del lobo)
          const ga = cg[a];
          const gb = cg[b];
          let nx = -(ga[0] + t * (gb[0] - ga[0]));
          let ny = -(ga[1] + t * (gb[1] - ga[1]));
          let nz = -(ga[2] + t * (gb[2] - ga[2]));
          const len = Math.hypot(nx, ny, nz) || 1;
          v[3] = nx / len;
          v[4] = ny / len;
          v[5] = nz / len;
        }
        const base = cubeIndex * 16;
        for (let t = 0; triTable[base + t] !== -1; t += 3) {
          // ordine invertito così che le facce esterne siano orientate in senso antiorario
          for (const e of [triTable[base + t], triTable[base + t + 2], triTable[base + t + 1]]) {
            const v = vertList[e];
            pos.push(v[0], v[1], v[2]);
            nor.push(v[3], v[4], v[5]);
          }
        }
      }
    }
  }
  return { positions: new Float32Array(pos), normals: new Float32Array(nor) };
}
