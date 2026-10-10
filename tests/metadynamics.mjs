// Test della metadinamica well-tempered:
// 1. Verifica gradienti analitici contro differenze finite (tutti i tipi di CV)
// 2. Convergenza del bias su sistema modello (doppio pozzo 1D)
// 3. Ricostruzione della superficie di energia libera
// 4. Diagnostiche di convergenza
// 5. Test FES 2D
// 6. Test reset e riusabilita

import { CollectiveVariable, WellTemperedMetaDynamics, estimateSigma } from '../js/chem/metadynamics.js';

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

// ─── 5. Convergenza bias su doppio pozzo 1D ───
// V(x) = (x^2 - 1)^2: minimi a x = ±1, barriera a x = 0 (altezza = 1 Hartree)
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.15});
  const metad = new WellTemperedMetaDynamics([cv], {
    W0: 0.01,
    biasFactor: 15,
    T: 300,
    depositionStride: 1
  });

  let x = 0.8;
  let v = 0;
  const dt = 0.005;
  const mass = 1.0;
  const nsteps = 3000;

  for (let step = 0; step < nsteps; step++) {
    // Forza dal potenziale V(x) = (x^2 - 1)^2, F = -dV/dx = -4x(x^2-1)
    const force = -4 * x * (x * x - 1);

    const Z = [1, 1];
    const pos = new Float64Array([0, 0, 0, x, 0, 0]);
    const forces = new Float64Array(6);
    forces[3] = force;
    metad.computeBiasForces(Z, pos, forces);

    // Velocity Verlet
    v += 0.5 * dt * forces[3] / mass;
    x += v * dt;
    v += 0.5 * dt * forces[3] / mass;
    v *= 0.999; // debole dissipazione

    // Confina
    if (x < -2.5) { x = -2.5; v = Math.abs(v); }
    if (x > 2.5) { x = 2.5; v = -Math.abs(v); }
  }

  const visited = metad.history.map(h => h.s[0]);
  const leftWell = visited.filter(s => s < 0).length;
  const rightWell = visited.filter(s => s > 0).length;

  check('Esplorazione doppio pozzo (metaD)', leftWell > 50 && rightWell > 50,
    `sinistro: ${leftWell}, destro: ${rightWell}, gaussiane: ${metad.history.length}`);
}

// ─── 6. Ricostruzione FES 1D ───
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.12});
  const metad = new WellTemperedMetaDynamics([cv], {
    W0: 0.02,
    biasFactor: 20,
    T: 300,
    depositionStride: 1
  });

  let x = 1.0;
  let v = 0;
  const dt = 0.004;
  const nsteps = 6000;

  for (let step = 0; step < nsteps; step++) {
    const force = -4 * x * (x * x - 1);
    const Z = [1, 1];
    const pos = new Float64Array([0, 0, 0, x, 0, 0]);
    const forces = new Float64Array(6);
    forces[3] = force;
    metad.computeBiasForces(Z, pos, forces);

    v += 0.5 * dt * forces[3];
    x += v * dt;
    v += 0.5 * dt * forces[3];
    v *= 0.998;
    if (x < -2.5) { x = -2.5; v = Math.abs(v); }
    if (x > 2.5) { x = 2.5; v = -Math.abs(v); }
  }

  const fes = metad.reconstructFES(80);

  // Il FES deve avere un massimo (barriera) vicino a x = 0
  // e minimi vicino a x = ±1
  let barrier_val = -Infinity, barrier_pos = 0;
  for (let i = 0; i < fes.s.length; i++) {
    if (Math.abs(fes.s[i]) < 0.3 && fes.F[i] > barrier_val) {
      barrier_val = fes.F[i];
      barrier_pos = fes.s[i];
    }
  }

  // Energia della barriera del doppio pozzo = 1 Hartree
  // Con bias factor 20, ci aspettiamo una ricostruzione ragionevole
  check('Ricostruzione FES (barriera a x=0)', barrier_val > 0.3 && Math.abs(barrier_pos) < 0.3,
    `F(barriera) = ${barrier_val.toFixed(3)} Hartree a x = ${barrier_pos.toFixed(2)} (atteso ~1.0 a x=0)`);
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

  check('FES 2D (struttura)', fes2d.s1.length === 20 && fes2d.s2.length === 20 && fes2d.F.length === 20,
    `griglia ${fes2d.s1.length}x${fes2d.s2.length}, F range [0, ${Math.max(...fes2d.F.flat()).toFixed(2)}]`);
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
