// Atomi relativistici: confronto con i dati di riferimento NIST ScRLDA (Kotochigova, Levine, Shirley, Stiles,
// Clark, "Atomic reference data for electronic structure calculations", https://math.nist.gov/DFTdata/atomdata/):
// equazione di Koelling–Harmon, scambio-correlazione VWN con la correzione relativistica di MacDonald–Vosko,
// densità dalla grande componente. Energie in Hartree.
import { runSCF } from '../js/physics/scf.js';

const REF = {'1': {'total': -0.445668, 'e': {'1s': -0.233463}}, '3': {'total': -7.335231, 'e': {'1s': -1.878227, '2s': -0.105541}}, '8': {'total': -74.508648, 'e': {'1s': -18.758595, '2p': -0.338161}}, '18': {'total': -527.516861, 'e': {'1s': -114.090739, '3p': -0.381755}}, '29': {'total': -1650.858424, 'e': {'1s': -323.674291, '4s': -0.17804, '3d': -0.196202}}, '47': {'total': -5304.298884, 'e': {'1s': -925.674173, '5s': -0.172492, '4d': -0.282334, '3d': -13.184815, '4s': -3.482162}}, '55': {'total': -7769.373337, 'e': {'1s': -1307.794166, '6s': -0.08162, '4d': -2.731851, '5s': -0.986167, '3d': -26.023693, '4s': -8.051051}}, '79': {'total': -18963.123612, 'e': {'1s': -2948.072638, '6s': -0.221999, '5d': -0.262411, '4d': -11.962303, '5s': -3.976365, '3d': -80.917806, '4s': -26.660632}}, '80': {'total': -19572.039277, 'e': {'1s': -3035.505153, '6s': -0.260552, '5d': -0.371655, '4d': -12.800533, '5s': -4.33481, '3d': -84.195815, '4s': -28.086349}}, '82': {'total': -20827.188652, 'e': {'1s': -3215.862876, '6p': -0.136111, '5d': -0.784774, '6s': -0.450774, '4d': -14.75623, '5s': -5.282191, '3d': -91.172486, '4s': -31.282373}}, '86': {'total': -23492.950469, 'e': {'1s': -3598.855827, '6p': -0.289891, '5d': -1.708347, '6s': -0.811281, '4d': -19.059234, '5s': -7.383977, '3d': -106.090128, '4s': -38.353106}}, '92': {'total': -27899.641656, 'e': {'1s': -4234.557545, '7s': -0.160828, '5d': -3.638317, '6s': -1.757867, '6p': -0.856604, '4d': -26.70564, '5s': -11.350493, '3d': -131.090673, '4s': -51.001789, '5f': -0.13339}}};
const L = 'spdf';
let failures = 0;
for (const [Z, ref] of Object.entries(REF)) {
  const r = runSCF(+Z, undefined, { relativistic: true });
  const dE = r.energy.total - ref.total;
  let worst = 0, which = '';
  for (const o of r.orbitals) {
    const k = `${o.n}${L[o.l]}`;
    if (ref.e[k] === undefined) continue;
    const d = o.e - ref.e[k];
    if (Math.abs(d) > Math.abs(worst)) { worst = d; which = k; }
  }
  const ok = r.converged && Math.abs(dE) < 1e-3 && Math.abs(worst) < 1e-3;
  if (!ok) failures++;
  console.log(`${ok ? '✓' : '✗'} Z = ${Z}: E = ${r.energy.total.toFixed(6)} (NIST ${ref.total}), ΔE = ${dE.toExponential(1)}, Δε max = ${worst.toExponential(1)} (${which})`);
}
// limite non relativistico: c → ∞ riproduce l'LDA (stesso confronto di scf-reference)
const nr = runSCF(18);
if (Math.abs(nr.energy.total - -525.946195) > 2e-5) { failures++; console.log('✗ LDA non relativistica cambiata'); }
if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nAtomi relativistici (ScRLDA NIST): tutte le verifiche superate');
