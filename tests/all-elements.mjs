// Verifica che il calcolo autoconsistente converga per tutti i 118 elementi e confronta
// l'energia di ionizzazione ΔSCF con quella sperimentale.
import { runSCF } from '../js/physics/scf.js';
import { groundStateConfiguration } from '../js/physics/configuration.js';
import { ionizationEnergy } from '../js/physics/atom.js';
import { ELEMENTS } from '../js/physics/elements.js';

let fails = 0;
let tmax = 0;
for (const el of ELEMENTS) {
  const t0 = performance.now();
  const res = runSCF(el.Z);
  const ion = ionizationEnergy(res);
  const dt = performance.now() - t0;
  tmax = Math.max(tmax, dt);
  const homo = res.orbitals.reduce((a, b) => (b.e > a.e ? b : a));
  const bad = !res.converged || !ion.converged || res.orbitals.some(o => !Number.isFinite(o.e));
  if (bad) fails++;
  console.log(`${String(el.Z).padStart(3)} ${el.symbol.padEnd(2)} it=${String(res.iterations).padStart(3)} ${res.converged ? 'ok ' : 'NO '} ${dt.toFixed(0).padStart(5)}ms  E=${res.energy.total.toFixed(4).padStart(13)}  IE=${ion.eV.toFixed(2)} eV (exp ${el.ionizationEV ?? '—'}) ${bad ? '<<<<' : ''}`);
}
console.log(`fallimenti: ${fails}, tempo massimo ${tmax.toFixed(0)} ms`);
process.exit(fails ? 1 : 0);
