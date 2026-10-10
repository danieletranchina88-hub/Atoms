// Esempio di utilizzo della metadinamica well-tempered nella sandbox.
// Caso studio: dissociazione di H2O -> OH + H
//
// Questo esempio mostra come:
// 1. Definire variabili collettive appropriate per la reazione
// 2. Inizializzare il bias metadinamico
// 3. Integrare con il campo di forze reattivo
// 4. Ricostruire la superficie di energia libera
// 5. Identificare lo stato di transizione e la barriera di energia libera
//
// Riferimento teorico:
//   A. Barducci, G. Bussi, M. Parrinello, Phys. Rev. Lett. 100, 020603 (2008)
//
// Esecuzione:
//   node examples/metadynamics-water-dissociation.js

import { CollectiveVariable, WellTemperedMetaDynamics, estimateSigma } from '../js/chem/metadynamics.js';
import { ReactiveFF } from '../js/chem/reactive.js';

// ─── Configurazione iniziale: molecola di H2O ───
// Geometria di equilibrio sperimentale (r_OH = 0.957 A, angolo HOH = 104.5 gradi)
const atoms = [
  { Z: 8, mass: 15.999, pos: [0.000,  0.000,  0.000] },   // O
  { Z: 1, mass: 1.008,  pos: [0.957,  0.000,  0.000] },   // H1 (legame che si rompe)
  { Z: 1, mass: 1.008,  pos: [-0.240,  0.927,  0.000] }   // H2 (legame spettatore)
];

const N = atoms.length;
const Z = atoms.map(a => a.Z);
const masses = new Float64Array(N);
atoms.forEach((a, i) => { masses[i] = a.mass * 1822.888; }); // amu -> massa elettronica
const pos = new Float64Array(atoms.flatMap(a => a.pos));
const vel = new Float64Array(pos.length);

// Inizializza velocita con distribuzione di Maxwell-Boltzmann a T = 300 K
const kB = 3.166808578e-6; // Hartree/K
const T_init = 300;
let seed = 42;
const rng = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
const gauss = () => { const u1 = rng(), u2 = rng(); return Math.sqrt(-2*Math.log(u1+1e-10))*Math.cos(2*Math.PI*u2); };
for (let i = 0; i < N; i++) {
  const sigma_v = Math.sqrt(kB * T_init / masses[i]);
  vel[3*i] = gauss() * sigma_v;
  vel[3*i+1] = gauss() * sigma_v;
  vel[3*i+2] = gauss() * sigma_v;
}

// ─── Definizione delle variabili collettive ───
// CV1: distanza O-H1 (il legame che si rompe)
const cv_OH1 = new CollectiveVariable('distance', [0, 1], { sigma: 0.12 });

// Stima automatica di sigma (verifica che sia ragionevole)
const sigmas = estimateSigma([cv_OH1], Z, pos, 300);
console.log(`Sigma stimato per CV(O-H): ${sigmas[0].toFixed(4)} A (usato: ${cv_OH1.sigma})`);

// ─── Inizializzazione metadinamica ───
const metad = new WellTemperedMetaDynamics([cv_OH1], {
  W0: 0.003,              // Altezza iniziale gaussiana (Hartree) ~ 1.9 kcal/mol
  biasFactor: 12,         // gamma = T_bias/T
  T: T_init,              // Temperatura (K)
  depositionStride: 50    // Deposita gaussiana ogni 50 step MD
});

// ─── Parametri di simulazione ───
const dt = 20.0;          // timestep in unita atomiche (~0.48 fs)
const nsteps = 20000;     // passi totali
const ff = new ReactiveFF();
const forces = new Float64Array(pos.length);

// Array per tracciamento
const trajectory = [];
const cvHistory = [];
const energyHistory = [];

console.log('=== Dissociazione H2O -> OH + H (Metadinamica Well-Tempered) ===');
console.log(`Parametri: W0=${metad.W0} Ha, gamma=${metad.biasFactor}, T=${metad.T} K, stride=${metad.depositionStride}`);
console.log(`dt = ${dt} u.a. (~${(dt * 0.02419).toFixed(3)} fs), passi = ${nsteps}`);
console.log('');

// ─── Loop di dinamica molecolare con bias ───
for (let step = 0; step < nsteps; step++) {
  // Calcola forze dal campo reattivo
  const energy = ff.compute(Z, pos, forces);

  // Aggiungi forze di bias dalla metadinamica
  const biasEnergy = metad.computeBiasForces(Z, pos, forces);

  // Integrazione Velocity Verlet
  for (let i = 0; i < N; i++) {
    const inv_m = 1.0 / masses[i];
    for (let a = 0; a < 3; a++) {
      const idx = 3*i + a;
      vel[idx] += 0.5 * dt * forces[idx] * inv_m;
      pos[idx] += dt * vel[idx];
    }
  }

  // Ricalcola forze per seconda meta del Verlet
  ff.compute(Z, pos, forces);
  metad.computeBiasForces(Z, pos, forces);

  for (let i = 0; i < N; i++) {
    const inv_m = 1.0 / masses[i];
    for (let a = 0; a < 3; a++) {
      const idx = 3*i + a;
      vel[idx] += 0.5 * dt * forces[idx] * inv_m;
    }
  }

  // Termostato debole (riscalamento verso T_target)
  if (step % 100 === 0) {
    let KE = 0;
    for (let i = 0; i < N; i++) {
      KE += 0.5 * masses[i] * (vel[3*i]**2 + vel[3*i+1]**2 + vel[3*i+2]**2);
    }
    const T_curr = 2 * KE / (3 * N * kB);
    const scale = Math.sqrt(T_init / Math.max(1, T_curr));
    const lambda = 0.01; // accoppiamento debole
    const s = 1 + lambda * (scale - 1);
    for (let k = 0; k < vel.length; k++) vel[k] *= s;
  }

  // Tracciamento
  const cv_val = cv_OH1.evaluate(Z, pos);
  if (step % 100 === 0) {
    cvHistory.push({ step, s: cv_val, bias: biasEnergy });
    energyHistory.push({ step, E: energy, Ebias: biasEnergy });
  }

  // Salva traiettoria ogni 500 step
  if (step % 500 === 0) {
    trajectory.push({ step, pos: Float64Array.from(pos), cv: cv_val });
  }
}

// ─── Ricostruzione della FES ───
const fes = metad.reconstructFES(100, [0.5, 4.0]);

// ─── Analisi risultati ───
console.log(`Passi simulati: ${nsteps}`);
console.log(`Gaussiane depositate: ${metad.history.length}`);
console.log(`Tempo simulato: ${(nsteps * dt * 0.02419 / 1000).toFixed(2)} ps`);
console.log('');

const r_final = cv_OH1.evaluate(Z, pos);
const r_OH2 = new CollectiveVariable('distance', [0, 2]).evaluate(Z, pos);
console.log(`Distanze finali:`);
console.log(`  O-H1: ${r_final.toFixed(3)} A ${r_final > 2.0 ? '(DISSOCIATO)' : '(legato)'}`);
console.log(`  O-H2: ${r_OH2.toFixed(3)} A ${r_OH2 > 2.0 ? '(DISSOCIATO)' : '(legato)'}`);
console.log('');

// ─── FES: stato di transizione e barriera ───
let ts_idx = 0, F_max = -Infinity;
let reactant_idx = 0, F_reactant = Infinity;

for (let i = 0; i < fes.s.length; i++) {
  if (fes.F[i] > F_max) { F_max = fes.F[i]; ts_idx = i; }
  if (fes.s[i] < 1.2 && fes.F[i] < F_reactant) { F_reactant = fes.F[i]; reactant_idx = i; }
}

const barrier_Ha = F_max - F_reactant;
const barrier_kcal = barrier_Ha * 627.509;

console.log(`Superficie di Energia Libera:`);
console.log(`  Range CV: ${fes.s[0].toFixed(2)} - ${fes.s[fes.s.length-1].toFixed(2)} A`);
console.log(`  Reagente: r = ${fes.s[reactant_idx].toFixed(2)} A, F = ${F_reactant.toFixed(4)} Ha`);
console.log(`  Stato di transizione: r = ${fes.s[ts_idx].toFixed(2)} A, F = ${F_max.toFixed(4)} Ha`);
console.log(`  Barriera di energia libera: dG‡ = ${barrier_Ha.toFixed(4)} Ha = ${barrier_kcal.toFixed(1)} kcal/mol`);
console.log(`  (Sperimentale per H2O -> OH + H: ~118 kcal/mol; il campo reattivo da un valore approssimato)`);
console.log('');

// ─── Convergenza ───
const diag = metad.diagnostics();
console.log(`Convergenza bias:`);
console.log(`  Gaussiane totali: ${diag.nGaussians}`);
console.log(`  Altezza media iniziale: ${diag.avgWeightEarly.toFixed(5)} Ha`);
console.log(`  Altezza media recente: ${diag.avgWeightRecent.toFixed(5)} Ha`);
console.log(`  Rapporto convergenza: ${diag.convergenceRatio.toFixed(3)} (${diag.convergenceRatio < 0.5 ? 'BUONA' : 'in corso'})`);

export { trajectory, cvHistory, fes, metad };
