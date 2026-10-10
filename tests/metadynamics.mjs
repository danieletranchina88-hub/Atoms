// Test della metadinamica well-tempered:
// 1. Verifica gradienti analitici contro differenze finite
// 2. Convergenza del bias su sistema modello (doppio pozzo)
// 3. Ricostruzione della superficie di energia libera

import { CollectiveVariable, WellTemperedMetaDynamics } from '../js/chem/metadynamics.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};

// 1. Gradiente distanza: verifica contro differenze finite
{
  const cv = new CollectiveVariable('distance', [0, 1]);
  const Z = [1, 1];
  const pos = new Float64Array([0, 0, 0, 1.5, 0.3, -0.2]);
  const grad = cv.gradient(Z, pos);
  
  const h = 1e-6;
  let maxErr = 0;
  for (let i = 0; i < 6; i++) {
    const p1 = pos.slice();
    const p2 = pos.slice();
    p1[i] += h;
    p2[i] -= h;
    const fd = (cv.evaluate(Z, p1) - cv.evaluate(Z, p2)) / (2 * h);
    const err = Math.abs(fd - grad[i]);
    maxErr = Math.max(maxErr, err);
  }
  
  check('Gradiente distanza (differenze finite)', maxErr < 1e-8, 
    `Errore massimo: ${maxErr.toExponential(2)}`);
}

// 2. Gradiente angolo: verifica contro differenze finite
{
  const cv = new CollectiveVariable('angle', [0, 1, 2]);
  const Z = [1, 1, 1];
  const pos = new Float64Array([1, 0, 0, 0, 0, 0, 0, 1, 0]); // angolo 90°
  const grad = cv.gradient(Z, pos);
  
  const h = 1e-6;
  let maxErr = 0;
  for (let i = 0; i < 9; i++) {
    const p1 = pos.slice();
    const p2 = pos.slice();
    p1[i] += h;
    p2[i] -= h;
    const fd = (cv.evaluate(Z, p1) - cv.evaluate(Z, p2)) / (2 * h);
    const err = Math.abs(fd - grad[i]);
    maxErr = Math.max(maxErr, err);
  }
  
  check('Gradiente angolo (differenze finite)', maxErr < 1e-7, 
    `Errore massimo: ${maxErr.toExponential(2)}`);
}

// 3. Gradiente coordinazione: verifica contro differenze finite
{
  const cv = new CollectiveVariable('coordination', [0], {r0: 2.0, n: 6, m: 12, sigma: 0.1});
  const Z = [1, 1, 1, 1];
  const pos = new Float64Array([0,0,0, 1.8,0,0, 0,1.8,0, 0,0,1.8]);
  const grad = cv.gradient(Z, pos);
  
  const h = 1e-5;
  let maxErr = 0;
  for (let i = 0; i < 12; i++) {
    const p1 = pos.slice();
    const p2 = pos.slice();
    p1[i] += h;
    p2[i] -= h;
    const fd = (cv.evaluate(Z, p1) - cv.evaluate(Z, p2)) / (2 * h);
    const err = Math.abs(fd - grad[i]);
    maxErr = Math.max(maxErr, err / Math.max(1, Math.abs(fd)));
  }
  
  check('Gradiente coordinazione (differenze finite)', maxErr < 1e-4, 
    `Errore relativo massimo: ${maxErr.toExponential(2)}`);
}

// 4. Convergenza bias su sistema modello
// Sistema 1D: V(x) = (x^2 - 1)^2 (doppio pozzo a x = ±1)
// La metadinamica deve riempire entrambi i minimi
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.15});
  const metad = new WellTemperedMetaDynamics([cv], {
    W0: 0.01,
    biasFactor: 15,
    T: 300,
    depositionStride: 1
  });
  
  // Simulazione semplice: particella in potenziale doppio pozzo
  let x = 0.5; // parte da un minimo
  const dt = 0.01;
  const nsteps = 2000;
  
  for (let step = 0; step < nsteps; step++) {
    // Forza dal potenziale V(x) = (x^2 - 1)^2
    const force = -4 * x * (x * x - 1);
    
    // Bias force
    const Z = [1, 1];
    const pos = new Float64Array([0, 0, 0, x, 0, 0]);
    const forces = new Float64Array(6);
    forces[3] = force;
    metad.computeBiasForces(Z, pos, forces);
    
    // Velocity Verlet semplice
    const vel = forces[3] * dt;
    x += vel * dt;
  }
  
  // Verifica che il bias abbia esplorato entrambi i minimi
  const visited = metad.history.map(h => h.s[0]);
  const leftWell = visited.filter(v => v < 0).length;
  const rightWell = visited.filter(v => v > 0).length;
  
  check('Esplorazione doppio pozzo (metaD)', 
    leftWell > 100 && rightWell > 100,
    `Visite: pozzo sinistro ${leftWell}, pozzo destro ${rightWell}`);
}

// 5. Ricostruzione FES
{
  const cv = new CollectiveVariable('distance', [0, 1], {sigma: 0.1});
  const metad = new WellTemperedMetaDynamics([cv], {
    W0: 0.02,
    biasFactor: 20,
    T: 300,
    depositionStride: 1
  });
  
  // Simulazione più lunga per accumulo bias
  let x = 1.0;
  const dt = 0.01;
  const nsteps = 5000;
  
  for (let step = 0; step < nsteps; step++) {
    const force = -4 * x * (x * x - 1);
    const Z = [1, 1];
    const pos = new Float64Array([0, 0, 0, x, 0, 0]);
    const forces = new Float64Array(6);
    forces[3] = force;
    metad.computeBiasForces(Z, pos, forces);
    
    const vel = forces[3] * dt;
    x += vel * dt;
    x = Math.max(-2, Math.min(2, x)); // confina
  }
  
  const fes = metad.reconstructFES(50);
  
  // Verifica che FES abbia due minimi vicino a x = ±1
  let min1 = -1, min2 = 1;
  let idx1 = fes.s.findIndex(s => s > -1.2 && s < -0.8);
  let idx2 = fes.s.findIndex(s => s > 0.8 && s < 1.2);
  
  let F_min1 = Infinity, F_min2 = Infinity;
  for (let i = Math.max(0, idx1-5); i < Math.min(fes.F.length, idx1+5); i++) {
    F_min1 = Math.min(F_min1, fes.F[i]);
  }
  for (let i = Math.max(0, idx2-5); i < Math.min(fes.F.length, idx2+5); i++) {
    F_min2 = Math.min(F_min2, fes.F[i]);
  }
  
  // F(barriera) > F(minimo) di almeno 0.5 kcal/mol
  const barrier_idx = fes.s.findIndex(s => s > -0.2 && s < 0.2);
  const F_barrier = barrier_idx >= 0 ? fes.F[barrier_idx] : Infinity;
  
  check('Ricostruzione FES doppio pozzo', 
    F_barrier > Math.max(F_min1, F_min2) + 0.5,
    `F(min1) = ${F_min1.toFixed(3)}, F(min2) = ${F_min2.toFixed(3)}, F(barriera) = ${F_barrier.toFixed(3)}`);
}

console.log(`\n${failures === 0 ? 'Tutti i test superati!' : `${failures} test falliti.`}`);
process.exit(failures === 0 ? 0 : 1);
