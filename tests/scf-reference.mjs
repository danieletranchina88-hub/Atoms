// Confronto con i dati di riferimento NIST (Kotochigova et al., LDA non relativistica).
import { runSCF } from '../js/physics/scf.js';

const REF = {
  1: { total: -0.445671, e: { '1s': -0.233471 } },
  2: { total: -2.834836, e: { '1s': -0.570425 } },
  10: { total: -128.233481, e: { '1s': -30.305855, '2s': -1.322809, '2p': -0.498034 } },
  18: { total: -525.946195, e: { '1s': -113.800134, '2p': -8.443439, '3p': -0.382330 } },
};
const L = 'spdf';
let ok = true;
for (const [Z, ref] of Object.entries(REF)) {
  const t0 = performance.now();
  const res = runSCF(+Z);
  const dt = performance.now() - t0;
  console.log(`Z=${Z} iter=${res.iterations} conv=${res.converged} t=${dt.toFixed(0)}ms Etot=${res.energy.total.toFixed(6)} (ref ${ref.total})`);
  if (Math.abs(res.energy.total - ref.total) > 2e-5) ok = false;
  for (const o of res.orbitals) {
    const key = `${o.n}${L[o.l]}`;
    console.log(`   ${key} e=${o.e.toFixed(6)} ref=${ref.e?.[key] ?? '—'}`);
    if (ref.e?.[key] !== undefined && Math.abs(o.e - ref.e[key]) > 2e-5) ok = false;
  }
}
console.log(ok ? 'OK' : 'FAIL');
process.exit(ok ? 0 : 1);
