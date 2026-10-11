// Test della metadinamica well-tempered:
// 1. Verifica gradienti analitici contro differenze finite (tutti i tipi di CV)
// 2. Esplorazione riproducibile del doppio pozzo con una CV firmata
// 3. Ricostruzione FES contro una somma gaussiana analitica
// 4. Diagnostiche di convergenza
// 5. Test FES 2D
// 6. Test reset e riusabilita

import { CollectiveVariable, WellTemperedMetaDynamics, estimateSigma } from '../js/chem/metadynamics.js';
import { makeRandom } from '../js/physics/wavefunction.js';

let failures = 0;
let passes = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (ok) passes++; else failures++;
};

// ─── 1. Gradiente distanza: verifica contro differenze finite ───
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.15});
  const Z = [1, 1];
  const pos = new Float64Array([0, 0, 0, 1.5, 0.3, -0.2]);
  const grad = cv.gradient(Z, pos);

  const h = 1e-6;
  let maxErr = 0;
  for (let i = 0; i < 6; i++) {
    const p1 = pos.slice(); p1[i] += h;
    const p2 = pos.slice(); p2[i] -= h;
    const fd = (cv.evaluate(Z, p1) - cv.evaluate(Z, p2)) / (2 * h);
    maxErr = Math.max(maxErr, Math.abs(fd - grad[i]));
  }

  check('Gradiente distanza (diff. finite)', maxErr < 1e-8,
    `err max = ${maxErr.toExponential(2)}`);
}

// ─── 2. Gradiente angolo: verifica contro differenze finite ───
{
  const cv = new CollectiveVariable('angle', [0, 1, 2], {sigma: 0.2});
  const Z = [1, 1, 1];
  // Configurazione non degenere (angolo ~ 104.5 gradi come H2O)
  const pos = new Float64Array([0.96, 0, 0, 0, 0, 0, -0.24, 0.93, 0]);
  const grad = cv.gradient(Z, pos);

  const h = 1e-6;
  let maxErr = 0;
  for (let i = 0; i < 9; i++) {
    const p1 = pos.slice(); p1[i] += h;
    const p2 = pos.slice(); p2[i] -= h;
    const fd = (cv.evaluate(Z, p1) - cv.evaluate(Z, p2)) / (2 * h);
    maxErr = Math.max(maxErr, Math.abs(fd - grad[i]));
  }

  check('Gradiente angolo (diff. finite)', maxErr < 1e-7,
    `err max = ${maxErr.toExponential(2)}, angolo = ${(cv.evaluate(Z, pos) * 180 / Math.PI).toFixed(1)}°`);
}

// ─── 3. Gradiente coordinazione: verifica contro differenze finite ───
{
  const cv = new CollectiveVariable('coordination', [0], {r0: 2.0, n: 6, m: 12, sigma: 0.1});
  const Z = [8, 1, 1, 1];
  const pos = new Float64Array([0,0,0, 1.8,0,0, 0,1.8,0, 0,0,1.8]);
  const grad = cv.gradient(Z, pos);

  const h = 1e-5;
  let maxRelErr = 0;
  for (let i = 0; i < 12; i++) {
    const p1 = pos.slice(); p1[i] += h;
    const p2 = pos.slice(); p2[i] -= h;
    const fd = (cv.evaluate(Z, p1) - cv.evaluate(Z, p2)) / (2 * h);
    const relErr = Math.abs(fd - grad[i]) / Math.max(1e-6, Math.abs(fd));
    maxRelErr = Math.max(maxRelErr, relErr);
  }

  check('Gradiente coordinazione (diff. finite)', maxRelErr < 1e-4,
    `err rel max = ${maxRelErr.toExponential(2)}, CN = ${cv.evaluate(Z, pos).toFixed(3)}`);
}

// ─── 4. Conservazione del momento: somma gradienti = 0 ───
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.1});
  const Z = [6, 8];
  const pos = new Float64Array([0.5, -0.3, 0.1, -0.7, 0.4, 0.9]);
  const grad = cv.gradient(Z, pos);

  let sumX = 0, sumY = 0, sumZ = 0;
  for (let i = 0; i < Z.length; i++) {
    sumX += grad[3*i]; sumY += grad[3*i+1]; sumZ += grad[3*i+2];
  }

  check('Conservazione momento (distanza)', Math.abs(sumX) + Math.abs(sumY) + Math.abs(sumZ) < 1e-14,
    `|sum F| = ${(Math.abs(sumX) + Math.abs(sumY) + Math.abs(sumZ)).toExponential(1)}`);
}

{
  const cv = new CollectiveVariable('angle', [0, 1, 2], {sigma: 0.2});
  const Z = [1, 8, 1];
  const pos = new Float64Array([0.96, 0, 0, 0, 0, 0, -0.24, 0.93, 0]);
  const grad = cv.gradient(Z, pos);

  let sumX = 0, sumY = 0, sumZ = 0;
  for (let i = 0; i < Z.length; i++) {
    sumX += grad[3*i]; sumY += grad[3*i+1]; sumZ += grad[3*i+2];
  }

  check('Conservazione momento (angolo)', Math.abs(sumX) + Math.abs(sumY) + Math.abs(sumZ) < 1e-12,
    `|sum F| = ${(Math.abs(sumX) + Math.abs(sumY) + Math.abs(sumZ)).toExponential(1)}`);
}

// ─── 5. Esplorazione riproducibile del doppio pozzo ───
// A distance is |x| and can NEVER visit negative values. This toy system needs
// a signed Cartesian CV, defined only in the test (no sandbox model changes).
// Langevin BAOAB uses Gaussian noise, a fixed seed and kBT=0.1 Ha, consistent
// with a 1 Ha toy barrier. A 300 K bath cannot equilibrate that barrier here.
class SignedToyCoordinate {
  constructor(sigma) { this.sigma = sigma; }
  evaluate(Z, pos) { return pos[3]; }
  gradient(Z, pos) { const grad = new Float64Array(pos.length); grad[3] = 1; return grad; }
}
{
  const cv = new SignedToyCoordinate(0.2);
  const kB = 3.166808578e-6, kT = 0.1;
  const metad = new WellTemperedMetaDynamics([cv], {
    W0: 0.02, biasFactor: 15, T: kT / kB, depositionStride: 20,
  });
  const rand = makeRandom(7319);
  const gaussian = () => Math.sqrt(-2 * Math.log(Math.max(rand(), 1e-15))) * Math.cos(2 * Math.PI * rand());
  const dt = 0.025, gamma = 0.5;
  const decay = Math.exp(-gamma * dt);
  const noiseScale = Math.sqrt(kT * (1 - decay * decay));
  let x = 0.8, v = 0;
  const Z = [1, 1];
  for (let step = 0; step < 12000; step++) {
    const pos = new Float64Array([0, 0, 0, x, 0, 0]);
    const forces = new Float64Array(6);
    forces[3] = -4 * x * (x * x - 1);
    metad.computeBiasForces(Z, pos, forces);
    // B A O A B: recompute the force at the new coordinate without another hill.
    v += 0.5 * dt * forces[3];
    x += 0.5 * dt * v;
    v = decay * v + noiseScale * gaussian();
    x += 0.5 * dt * v;
    v += 0.5 * dt * (-4 * x * (x * x - 1) - metad.biasGradientCV([x])[0]);
    if (!Number.isFinite(x) || !Number.isFinite(v)) throw new Error('Toy integrator diverged');
  }
  const visited = metad.history.map(h => h.s[0]);
  const left = visited.filter(x => x < -0.5).length;
  const right = visited.filter(x => x > 0.5).length;
  check('Esplorazione doppio pozzo (CV firmata, seed fisso)', left > 50 && right > 50,
    `sinistro: ${left}, destro: ${right}, gaussiane: ${metad.history.length}`);

  // Independent finite difference oracle for the production bias force.
  const pos = new Float64Array([0, 0, 0, 0.37, 0, 0]);
  const h = 1e-5;
  metad.depositionStride = 1000000;
  const forces = new Float64Array(6);
  metad.computeBiasForces(Z, pos, forces);
  const fdForce = -(metad.evaluateBias([pos[3] + h]) - metad.evaluateBias([pos[3] - h])) / (2 * h);
  check('Forza di bias contro differenze finite', Math.abs(forces[3] - fdForce) < 1e-7,
    `errore = ${Math.abs(forces[3] - fdForce).toExponential(2)}`);
}

// ─── 6. FES 1D contro un bias analitico noto ───
// Finite-time exploration is not a proof of FES convergence. Verify the actual
// reconstruction formula and normalization with a symmetric two-hill oracle.
{
  const sigma = 0.2, weight = 0.08, gamma = 15;
  const metad = new WellTemperedMetaDynamics([new SignedToyCoordinate(sigma)], {
    W0: weight, biasFactor: gamma, T: 300,
  });
  metad.history = [-1, 1].map(x => ({ s: new Float64Array([x]), weight }));
  const fes = metad.reconstructFES(201, [-2, 2]);
  const raw = fes.s.map(x => -gamma / (gamma - 1) * weight * (
    Math.exp(-((x + 1) ** 2) / (2 * sigma * sigma)) +
    Math.exp(-((x - 1) ** 2) / (2 * sigma * sigma))
  ));
  const min = Math.min(...raw);
  const err = Math.max(...fes.F.map((F, i) => Math.abs(F - (raw[i] - min))));
  check('Ricostruzione FES analitica e normalizzazione', err < 1e-12 && Math.min(...fes.F) === 0,
    `errore massimo = ${err.toExponential(2)} Ha`);
  const asymmetry = Math.max(...fes.F.map((F, i) => Math.abs(F - fes.F[fes.F.length - 1 - i])));
  check('FES simmetrica con due pozzi', asymmetry < 1e-12 && fes.F[100] > fes.F[50] && fes.F[100] > fes.F[150],
    `asimmetria = ${asymmetry.toExponential(2)} Ha`);
}

// ─── 7. Diagnostiche di convergenza ───
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.1});
  const metad = new WellTemperedMetaDynamics([cv], {
    W0: 0.01, biasFactor: 10, T: 300, depositionStride: 1
  });

  // Riempi con gaussiane
  for (let i = 0; i < 200; i++) {
    const x = 1.0 + Math.sin(i * 0.1) * 0.5;
    const Z = [1, 1];
    const pos = new Float64Array([0, 0, 0, x, 0, 0]);
    const forces = new Float64Array(6);
    metad.computeBiasForces(Z, pos, forces);
  }

  const diag = metad.diagnostics();

  check('Diagnostiche convergenza',
    diag.nGaussians > 0 && diag.convergenceRatio < 1.0,
    `gauss=${diag.nGaussians}, ratio=${diag.convergenceRatio.toFixed(3)} (< 1 = convergente)`);
}

// ─── 8. FES 2D ───
{
  const cv1 = new CollectiveVariable('distance', [0, 1], {sigma: 0.15});
  const cv2 = new CollectiveVariable('distance', [0, 2], {sigma: 0.15});
  const metad = new WellTemperedMetaDynamics([cv1, cv2], {
    W0: 0.01, biasFactor: 10, T: 300, depositionStride: 1
  });

  // Accumula storia in 2D
  for (let i = 0; i < 300; i++) {
    const x1 = 1.0 + 0.5 * Math.sin(i * 0.05);
    const x2 = 1.2 + 0.3 * Math.cos(i * 0.07);
    const Z = [1, 1, 1];
    const pos = new Float64Array([0,0,0, x1,0,0, 0,x2,0]);
    const forces = new Float64Array(9);
    metad.computeBiasForces(Z, pos, forces);
  }

  const fes2d = metad.reconstructFES2D(20);

  const values = fes2d.F.flatMap(row => Array.from(row));
  check('FES 2D (struttura e valori finiti)', fes2d.s1.length === 20 && fes2d.s2.length === 20 && fes2d.F.length === 20 && fes2d.F.every(row => row.length === 20) && values.every(Number.isFinite) && Math.min(...values) === 0,
    `griglia ${fes2d.s1.length}x${fes2d.s2.length}, F range [0, ${Math.max(...values).toFixed(2)}]`);
}

// ─── 9. Reset ───
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.1});
  const metad = new WellTemperedMetaDynamics([cv], {W0: 0.01, biasFactor: 10, T: 300});

  const Z = [1, 1];
  const pos = new Float64Array([0, 0, 0, 1.5, 0, 0]);
  const forces = new Float64Array(6);

  for (let i = 0; i < 200; i++) metad.computeBiasForces(Z, pos, forces);
  const beforeReset = metad.history.length;

  metad.reset();

  check('Reset completo', beforeReset > 0 && metad.history.length === 0 && metad.timestep === 0,
    `prima: ${beforeReset} gaussiane, dopo: ${metad.history.length}`);
}

// ─── 10. estimateSigma ───
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.1});
  const Z = [1, 1];
  const pos = new Float64Array([0, 0, 0, 1.4, 0, 0]);

  const sigmas = estimateSigma([cv], Z, pos, 500);

  check('estimateSigma (stima ragionevole)', sigmas[0] > 0.001 && sigmas[0] < 0.1,
    `sigma stimato = ${sigmas[0].toFixed(4)} (atteso ~0.005-0.01)`);
}

// ─── Riepilogo ───
console.log(`\n${'═'.repeat(50)}`);
console.log(`Risultati: ${passes} superati, ${failures} falliti`);
console.log(`${'═'.repeat(50)}`);
process.exit(failures === 0 ? 0 : 1);
