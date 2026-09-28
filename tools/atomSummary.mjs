// Calcola per tutti i 118 elementi le grandezze usate per colorare la tavola periodica
// (energia di ionizzazione ΔSCF, raggio al 90% della carica, energia dell'orbitale più esterno).
import fs from 'node:fs';
import { computeAtom } from '../js/physics/atom.js';

const HARTREE_EV = 27.211386245988;
const BOHR_PM = 52.917721090;
const out = {};
for (let Z = 1; Z <= 118; Z++) {
  const a = computeAtom(Z);
  const homo = a.orbitals.reduce((x, y) => (y.e > x.e ? y : x));
  out[Z] = { ie: +a.ionization.eV.toFixed(3), r90: +(a.r90 * BOHR_PM).toFixed(1), homo: +(homo.e * HARTREE_EV).toFixed(3), rmax: +(homo.rMaxProb * BOHR_PM).toFixed(1) };
  if (Z % 20 === 0) console.log(Z);
}
fs.writeFileSync(new URL('../js/physics/atomSummary.js', import.meta.url),
  `// Generato da tools/atomSummary.mjs: grandezze calcolate (DFT-LDA) per i 118 elementi.\n// ie: energia di ionizzazione ΔSCF (eV); r90: raggio che contiene il 90% degli elettroni (pm);\n// homo: energia dell'orbitale occupato più alto (eV); rmax: suo raggio di massima probabilità (pm).\nexport const ATOM_SUMMARY = ${JSON.stringify(out)};\n`);
console.log('fatto');
