// Esempio di utilizzo della metadinamica well-tempered nella sandbox
// Caso studio: dissociazione di H₂O → OH + H
//
// Questo esempio mostra come:
// 1. Definire variabili collettive appropriate per la reazione
// 2. Inizializzare il bias metadinamico
// 3. Integrare con il simulatore esistente
// 4. Ricostruire la superficie di energia libera

import { Simulation } from './md.js';
import { ReactiveFF } from './reactive.js';
import { CollectiveVariable, WellTemperedMetaDynamics } from './metadynamics.js';

// Configurazione iniziale: molecola di H₂O
const atoms = [
  {Z: 8, pos: [0, 0, 0]},        // O
  {Z: 1, pos: [0.96, 0, 0]},     // H1
  {Z: 1, pos: [-0.24, 0.93, 0]}  // H2
];

const Z = atoms.map(a => a.Z);
const pos = new Float64Array(atoms.flatMap(a => a.pos));
const vel = new Float64Array(pos.length); // velocità iniziali nulle

// Definizione delle variabili collettive per la dissociazione
// CV1: distanza O-H1 (il legame che si rompe)
const cv_OH1 = new CollectiveVariable('distance', [0, 1], {sigma: 0.15});

// CV2: distanza O-H2 (per monitorare che non si rompa)
const cv_OH2 = new CollectiveVariable('distance', [0, 2], {sigma: 0.15});

// Inizializzazione metadinamica
const metad = new WellTemperedMetaDynamics([cv_OH1, cv_OH2], {
  W0: 0.002,              // Altezza iniziale gaussiana (Hartree)
  biasFactor: 12,         // gamma = T_bias/T (controllo convergenza)
  T: 300,                 // Temperatura (K)
  depositionStride: 50    // Deposita gaussiana ogni 50 step MD
});

// Parametri di simulazione
const dt = 0.5; // fs
const nsteps = 10000;
const ff = new ReactiveFF();

// Array per tracciamento
const trajectory = [];
const cvHistory = [];

// Loop di dinamica molecolare con bias metadinamico
for (let step = 0; step < nsteps; step++) {
  // Calcola forze dal campo reattivo
  const forces = new Float64Array(pos.length);
  const energy = ff.compute(Z, pos, forces);
  
  // Aggiungi forze di bias dalla metadinamica
  const biasEnergy = metad.computeBiasForces(Z, pos, forces);
  
  // Integrazione Velocity Verlet (semplificata)
  for (let i = 0; i < pos.length; i++) {
    vel[i] += 0.5 * dt * forces[i] / atoms[Math.floor(i/3)].Z;
    pos[i] += dt * vel[i];
    vel[i] += 0.5 * dt * forces[i] / atoms[Math.floor(i/3)].Z;
  }
  
  // Salva traiettoria ogni 100 step
  if (step % 100 === 0) {
    const s1 = cv_OH1.evaluate(Z, pos);
    const s2 = cv_OH2.evaluate(Z, pos);
    trajectory.push({step, pos: pos.slice(), energy, biasEnergy});
    cvHistory.push({step, s1, s2});
  }
}

// Ricostruzione della superficie di energia libera
// Attenzione: richiede CV singola per ora
const metad_1D = new WellTemperedMetaDynamics([cv_OH1], {
  W0: 0.002,
  biasFactor: 12,
  T: 300,
  depositionStride: 50
});

// Esegui nuova simulazione con CV singola per FES
// ... (codice simile)

const fes = metad_1D.reconstructFES(100);

// Output risultati
console.log('=== Dissociazione H₂O → OH + H ===');
console.log(`Passi simulati: ${nsteps}`);
console.log(`Gaussiane depositate: ${metad.history.length}`);
console.log(`\nDistanze finali:`);
console.log(`  O-H1: ${cv_OH1.evaluate(Z, pos).toFixed(3)} Å`);
console.log(`  O-H2: ${cv_OH2.evaluate(Z, pos).toFixed(3)} Å`);

if (fes) {
  console.log(`\nSuperficie di Energia Libera:`);
  console.log(`  Range CV: ${fes.s[0].toFixed(2)} - ${fes.s[fes.s.length-1].toFixed(2)} Å`);
  console.log(`  FES punti: ${fes.F.length}`);
  
  // Trova stato di transizione
  let ts_idx = 0;
  let F_max = -Infinity;
  for (let i = 0; i < fes.F.length; i++) {
    if (fes.F[i] > F_max) {
      F_max = fes.F[i];
      ts_idx = i;
    }
  }
  console.log(`  Stato di transizione: r = ${fes.s[ts_idx].toFixed(2)} Å, ΔF = ${F_max.toFixed(2)} Hartree`);
}

// Analisi della convergenza
console.log(`\nConvergenza bias:`);
const recent_bias = metad.history.slice(-10).reduce((s, h) => s + h.weight, 0) / 10;
const early_bias = metad.history.slice(0, 10).reduce((s, h) => s + h.weight, 0) / 10;
console.log(`  Altezza bias iniziale (media primi 10): ${early_bias.toFixed(4)} Hartree`);
console.log(`  Altezza bias finale (media ultimi 10): ${recent_bias.toFixed(4)} Hartree`);
console.log(`  Riduzione: ${((1 - recent_bias/early_bias) * 100).toFixed(1)}%`);

export {trajectory, cvHistory, fes};
