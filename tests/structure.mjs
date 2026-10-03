// Verifiche dell'ordine strutturale (q₆ di Steinhardt), delle equazioni nette del becher e del passo della sandbox.

import { localOrder, phaseVerdict, Mobility } from '../js/chem/structureOrder.js';
import { PhaseMD } from '../js/chem/phaseMD.js';
import { Beaker } from '../js/chem/aqueous.js';
import { beakerPresets } from '../js/chem/beakerReagents.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};

// FCC perfetto: Q₆ = 0,57452 (Steinhardt et al. 1983), tutti gli atomi cristallini
{
  const md = new PhaseMD({ density: 0.95, temperature: 0.35 });
  const o = localOrder(md.x, { box: md.L, cutoff: 1.5 });
  check('FCC perfetto: Q₆ e atomi cristallini', Math.abs(o.Q6 - 0.57452) < 2e-4 && o.fraction.crystal === 1, `Q₆ = ${o.Q6.toFixed(5)}`);
  // invarianza per rotazione: un blocco ruotato senza condizioni periodiche resta cristallino all'interno
  const c = Math.cos(0.7), s = Math.sin(0.7), x = md.x.slice();
  for (let i = 0; i < md.n; i++) { const a = x[3 * i], b = x[3 * i + 1]; x[3 * i] = c * a - s * b; x[3 * i + 1] = s * a + c * b; }
  const r = localOrder(x, { cutoff: 1.5 });
  check('cristallo ruotato riconosciuto', r.fraction.crystal > 0.8, `${(100 * r.fraction.crystal).toFixed(0)}% cristallino (superficie esclusa)`);
  md.step(800);
  const t = localOrder(md.x, { box: md.L, cutoff: 1.5 });
  check('cristallo a T* = 0,35 resta cristallo', phaseVerdict(t.fraction).id === 'solido', `${(100 * t.fraction.crystal).toFixed(0)}% cristallino`);
}
{
  const liq = new PhaseMD({ density: 0.85, temperature: 1.5 });
  liq.step(1500);
  const o = localOrder(liq.x, { box: liq.L, cutoff: 1.5 });
  check('liquido a T* = 1,5: nessun ordine', o.fraction.crystal < 0.05 && phaseVerdict(o.fraction).id === 'liquido', `${(100 * o.fraction.crystal).toFixed(0)}% cristallino, Q₆ = ${o.Q6.toFixed(3)}`);
  const gas = new PhaseMD({ density: 0.04, temperature: 2 });
  gas.step(500);
  const g = localOrder(gas.x, { box: gas.L, cutoff: 1.5 });
  check('gas diluito', phaseVerdict(g.fraction).id === 'gas', `${(100 * g.fraction.vapor).toFixed(0)}% atomi isolati`);
}
// mobilità: le gabbie di vicini si rinnovano nel liquido, non nel cristallo
{
  const mob = (md) => { const m = new Mobility(); for (let k = 0; k < 4; k++) { m.sample(md.x, md.time, 1, 1, md.L); md.step(400); } return [m.sample(md.x, md.time, 1, 1, md.L), m.lostFraction]; };
  const [c, lc] = mob(new PhaseMD({ density: 0.95, temperature: 0.35 }));
  const [l, ll] = mob(new PhaseMD({ density: 0.85, temperature: 1.5 }));
  check('mobilità: cristallo fermo, liquido che scorre', c === false && l === true, `vicini persi in 3τ: ${(100 * lc).toFixed(0)}% e ${(100 * ll).toFixed(0)}%`);
}
// equazioni nette del becher: mostrate solo se bilanciate
{
  const want = { 'cu-carbonate': '2 Cu²⁺ + 4 HCO₃⁻ → 3 CO₂(aq) + Cu₂(OH)₂CO₃(s) + H₂O', 'malachite-acid': '4 H⁺ + Cu₂(OH)₂CO₃(s) → 2 Cu²⁺ + CO₂(aq) + 3 H₂O', barium: 'Ba²⁺ + SO₄²⁻ → BaSO₄(s)' };
  const got = {};
  for (const p of beakerPresets().filter(p => want[p.id])) {
    const b = new Beaker();
    if (p.thermostat) b.setThermostat(true, 298.15);
    for (const [r, n] of p.steps) b.add(r, n);
    got[p.id] = b.add(p.next, p.amount).equation?.text;
  }
  const bad = Object.keys(want).filter(k => got[k] !== want[k]);
  check('equazioni nette bilanciate', !bad.length, bad.map(k => `${k}: ${got[k]}`).join(' | ') || Object.values(got).join(' · '));
}

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nOrdine strutturale ed equazioni: tutte le verifiche superate');
