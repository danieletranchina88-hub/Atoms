// Controlli analitici: orbitali dell'idrogeno e normalizzazione delle armoniche sferiche.
import { runSCF, makeGrid, solveRadial } from '../js/physics/scf.js';
import { RadialFunction, samplePoints, term, sampleGrid } from '../js/physics/wavefunction.js';
import { realYlm } from '../js/physics/harmonics.js';

let ok = true;
const check = (label, got, want, tol) => {
  const pass = Math.abs(got - want) <= tol;
  if (!pass) ok = false;
  console.log(`${pass ? 'ok ' : 'NO '} ${label}: ${got.toFixed(6)} (atteso ${want})`);
};

// 1) Atomo di idrogeno esatto: E_n = −1/(2n²), ⟨r⟩ = (3n² − l(l+1))/2
const grid = makeGrid(1);
const V = new Float64Array(grid.N).map((_, i) => -1 / grid.r[i]);
for (const [n, l] of [[1, 0], [2, 0], [2, 1], [3, 2], [4, 3], [6, 2], [7, 0]]) {
  const res = solveRadial(grid, V, n, l, undefined, 1);
  check(`H ${n},${l} energia`, res.e, -1 / (2 * n * n), 1e-7);
  let rAvg = 0;
  for (let i = 0; i < grid.N; i++) rAvg += res.u[i] ** 2 * grid.r[i] ** 2 * grid.h;
  check(`H ${n},${l} <r>`, rAvg, (3 * n * n - l * (l + 1)) / 2, 1e-3 * n * n);
  // R_10 analitico = 2 e^{-r}
  if (n === 1) {
    const R = new RadialFunction(grid, res.u, 0);
    check('R_10(1)', R.at(1), 2 * Math.exp(-1), 1e-5);
  }
}

// 2) Normalizzazione di Y_lm: ∫|Y|² dΩ = 1, ortogonalità
function integrateSphere(fn) {
  const nt = 200, np = 400;
  let s = 0;
  for (let i = 0; i < nt; i++) {
    const th = Math.PI * (i + 0.5) / nt;
    for (let j = 0; j < np; j++) {
      const ph = 2 * Math.PI * (j + 0.5) / np;
      s += fn(Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)) * Math.sin(th);
    }
  }
  return s * (Math.PI / nt) * (2 * Math.PI / np);
}
for (let l = 0; l <= 4; l++) {
  for (let m = -l; m <= l; m++) {
    const v = integrateSphere((x, y, z) => realYlm(l, m, x, y, z) ** 2);
    if (Math.abs(v - 1) > 1e-3) { ok = false; console.log(`NO norm Y${l}${m} = ${v}`); }
  }
}
check('<Y21|Y2-1>', integrateSphere((x, y, z) => realYlm(2, 1, x, y, z) * realYlm(2, -1, x, y, z)), 0, 1e-6);
// p_x ha il lobo positivo lungo +x
check('segno p_x in +x', Math.sign(realYlm(1, 1, 1, 0, 0)), 1, 0);
check('d_xy ∝ xy', realYlm(2, -2, 1, 1, 0) / realYlm(2, -2, 1, -1, 0), -1, 1e-12);
check('d_x2-y2 lungo x', Math.sign(realYlm(2, 2, 1, 0, 0)), 1, 0);

// 3) Campionamento: ⟨r⟩ dei punti di un 2p dell'idrogeno ≈ 5 a₀
const r21 = solveRadial(grid, V, 2, 1, undefined, 1);
const R21 = new RadialFunction(grid, r21.u, 1);
const pts = samplePoints([term(R21, 1, 0)], 40000, 3);
let mean = 0;
let zpos = 0;
for (let i = 0; i < pts.count; i++) {
  const [x, y, z] = [pts.positions[3 * i], pts.positions[3 * i + 1], pts.positions[3 * i + 2]];
  mean += Math.hypot(x, y, z);
  if ((z > 0) === (pts.signs[i] > 0)) zpos++;
}
check('Monte Carlo <r> 2p', mean / pts.count, 5, 0.08);
check('fase coerente con z per 2p_z', zpos / pts.count, 1, 1e-9);
// LCAO: la combinazione campiona senza errori
const s1 = new RadialFunction(grid, solveRadial(grid, V, 1, 0, undefined, 1).u, 0);
const mo = samplePoints([term(s1, 0, 0, 1, [0, 0, -0.7]), term(s1, 0, 0, -1, [0, 0, 0.7])], 5000, 5);
let wrong = 0;
for (let i = 0; i < mo.count; i++) if ((mo.positions[3 * i + 2] < 0) !== (mo.signs[i] > 0)) wrong++;
check('σ* 1s: segno opposto ai due lati del piano nodale', wrong, 0, 0);
const g = sampleGrid([term(R21, 1, 0)], 40, 20, 0.9);
console.log(`isovalore 90% per 2p_z: ${g.iso.toExponential(3)}`);

console.log(ok ? 'OK' : 'FAIL');
process.exit(ok ? 0 : 1);
