// Verifiche del motore cinetico: soluzioni esatte, un problema rigido di riferimento, dati di letteratura.

import { parseMechanism, simulate, rateConstant, orderAnalysis, eyring, balance, KINETICS_PRESETS } from '../js/chem/kinetics.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};
const preset = (id) => KINETICS_PRESETS.find(p => p.id === id);
const runPreset = (p) => { const m = parseMechanism(p.mech); return { m, r: simulate(m, p.c0, p.T, p.tEnd) }; };
const col = (r, s) => r.c.map(c => c[r.vars.indexOf(s)]);

// sintassi
{
  const m = parseMechanism('2 A -> B ; k=1 ord.A=1\nX -> Y ; t12=2 h\nA + C <=> D ; A=1e10 Ea=50 kr=3\ncostanti: C\nN2O5 -> 2 NO2 + 1/2 O2 ; k=1');
  check('lettura del meccanismo', !m.errors.length && m.reactions.length === 4 && m.fixed.has('C')
    && Math.abs(m.reactions[1].fwd.k - Math.LN2 / 7200) < 1e-15 && m.reactions[3].products[1][1] === 0.5 && m.reactions[0].orders.A === 1,
  m.errors.join('; ') || `${m.species.join(', ')}`);
  const bad = parseMechanism('A B C\nA -> B ; x=3\nA <=> B ; k=1');
  check('errori segnalati', bad.errors.length === 3, bad.errors.join(' | '));
  check('bilancio di massa e carica', balance(parseMechanism('C2H5Br + OH- -> C2H5OH + Br- ; k=1').reactions[0]).ok
    && !balance(parseMechanism('H2 + O2 -> H2O ; k=1').reactions[0]).ok, 'C₂H₅Br + OH⁻ bilanciata, H₂ + O₂ → H₂O no');
}
// primo ordine: Arrhenius di Atkins contro la costante misurata (ciclopropano a 500 °C: 6,71·10⁻⁴ s⁻¹)
{
  const p = preset('cyclopropane');
  const { m, r } = runPreset(p);
  const k = rateConstant(m.reactions[0].fwd, p.T);
  check('ciclopropano: Arrhenius riproduce k misurata', Math.abs(k / 6.71e-4 - 1) < 0.03, `k = ${k.toExponential(3)} s⁻¹ (Atkins: 6,71·10⁻⁴)`);
  const oa = orderAnalysis(r.t, col(r, 'C3H6cyclo'));
  check('primo ordine: t½ = ln 2 / k', oa.best.order === 1 && Math.abs(oa.t12 * k / Math.LN2 - 1) < 1e-4, `t½ = ${(oa.t12 / 60).toFixed(3)} min`);
  const e = eyring(m.reactions[0].fwd, p.T);
  check('parametri di Eyring', Math.abs(e.dH / 1000 - (272 - 8.314462618e-3 * p.T)) < 1e-6 && e.dS > 0, `ΔH‡ = ${(e.dH / 1000).toFixed(1)} kJ/mol, ΔS‡ = ${e.dS.toFixed(1)} J/(K mol)`);
}
// secondo ordine con concentrazioni uguali: 1/[A] = 1/[A]₀ + k t
{
  const p = preset('sn2');
  const { m, r } = runPreset(p);
  const k = rateConstant(m.reactions[0].fwd, p.T);
  let err = 0;
  r.t.forEach((t, i) => { const a = 1 / (1 / 0.1 + k * t); err = Math.max(err, Math.abs(r.c[i][0] - a) / a); });
  check('secondo ordine: soluzione esatta', err < 1e-4, `errore relativo massimo ${err.toExponential(1)}`);
  const oa = orderAnalysis(r.t, col(r, 'C2H5Br'));
  check('ordine apparente 2', oa.best.order === 2 && Math.abs(oa.best.kobs / k - 1) < 1e-3, `k = ${oa.best.kobs.toExponential(3)} M⁻¹ s⁻¹`);
}
// serie radioattiva: soluzione di Bateman
{
  const { r } = runPreset(preset('decay'));
  const k1 = Math.LN2 / (5.012 * 86400), k2 = Math.LN2 / (138.376 * 86400);
  let err = 0;
  r.t.forEach((t, i) => { err = Math.max(err, Math.abs(r.c[i][1] - k1 / (k2 - k1) * (Math.exp(-k1 * t) - Math.exp(-k2 * t)))); });
  check('²¹⁰Bi → ²¹⁰Po → ²⁰⁶Pb (Bateman)', err < 1e-5, `errore massimo ${err.toExponential(1)} (frazione di ²¹⁰Bi iniziale)`);
}
// problema rigido di Robertson (1966); riferimento a t = 40 da Hairer e Wanner, Solving ODEs II
{
  const m = parseMechanism('A -> B ; k=0.04\n2 B -> B + C ; k=3e7\nB + C -> A + C ; k=1e4');
  const r = simulate(m, { A: 1, B: 0, C: 0 }, 298.15, 40);
  const y = r.c[r.c.length - 1];
  const ref = [0.7158270687, 9.185534764e-6, 0.2841637457];
  const e = Math.max(...ref.map((v, i) => Math.abs(y[i] - v) / v));
  check('Robertson (rigido): valori di riferimento a t = 40', e < 1e-4 && r.stats.accepted < 10000, `errore relativo ${e.toExponential(1)}, ${r.stats.accepted} passi`);
}
// equilibrio: Q → K = k/k_r
{
  const m = parseMechanism('A + B <=> C ; k=2 kr=0.5');
  const r = simulate(m, { A: 1, B: 0.5, C: 0 }, 298.15, 50);
  const [a, b, c] = r.c[r.c.length - 1];
  check('equilibrio: Q = k/k_r', Math.abs(c / (a * b) / 4 - 1) < 1e-4 && Math.abs(a + c - 1) < 1e-9, `Q = ${(c / (a * b)).toFixed(5)}, K = 4`);
}
// ciclo di Chapman: stato stazionario [O₃] = [O₂] √(J₁ k₂[M] / (J₃ k₄))
{
  const p = preset('chapman');
  const { m, r } = runPreset(p);
  const k = m.reactions.map(x => rateConstant(x.fwd, p.T));
  const k2M = k[1] * p.c0.N2 + k[2] * p.c0.O2;
  const o3 = p.c0.O2 * Math.sqrt(k[0] * k2M / (k[3] * k[4]));
  const sim = col(r, 'O3').at(-1);
  check('Chapman: stato stazionario analitico', Math.abs(sim / o3 - 1) < 1e-3, `[O₃] = ${sim.toExponential(3)} (formula: ${o3.toExponential(3)}) molecole/cm³`);
}
// Oregonator: ciclo limite
{
  const { r } = runPreset(preset('oregonator'));
  const ce = col(r, 'Ce4+');
  let peaks = 0;
  for (let i = 1; i < ce.length - 1; i++) if (ce[i] > ce[i - 1] && ce[i] >= ce[i + 1] && ce[i] > 1e-4) peaks++;
  check('Oregonator oscilla', peaks >= 4, `${peaks} massimi di Ce(IV) in ${(preset('oregonator').tEnd / 60).toFixed(0)} min`);
}
// tutte le esperienze arrivano in fondo
{
  const bad = KINETICS_PRESETS.filter(p => { const { m, r } = runPreset(p); return m.errors.length || !r.finished || r.c.at(-1).some(v => !Number.isFinite(v)); });
  check('esperienze di cinetica', bad.length === 0, bad.map(p => p.id).join(', ') || `${KINETICS_PRESETS.length} completate`);
}

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nCinetica: tutte le verifiche superate');
