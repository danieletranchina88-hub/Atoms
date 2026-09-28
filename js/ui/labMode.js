// Modalità "Laboratorio": leggi della chimica fisica con grafici calcolati in tempo reale.

import { VDW_GASES, ACIDS, BASES, INDICATORS, HALF_REACTIONS, R_LBAR, R_GAS, FARADAY } from '../chem/labData.js';
import {
  vdwPressure, vdwCritical, vdwVolumes, maxwellConstruction, fractions, titrationPH,
  galvanicCell, integratedRate, halfLife, consecutive,
} from '../chem/labPhysics.js';
import { drawXY } from './chemCharts.js';

const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const nf = (v, d = 2) => (v === null || v === undefined || !Number.isFinite(v) ? '—' : Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d, useGrouping: Math.abs(v) >= 1e4 }));
const sci = (v, d = 2) => {
  if (!Number.isFinite(v)) return '—';
  if (v === 0) return '0';
  const e = Math.floor(Math.log10(Math.abs(v)));
  if (e >= -2 && e <= 4) return nf(v, d);
  const s = String(e).split('').map(c => ({ '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' }[c])).join('');
  return `${nf(v / Math.pow(10, e), d)}·10${s}`;
};

let active = false;
const S = {
  tool: 'gas',
  gas: { id: 'CO2', Tr: 0.92, V: null },
  acid: { analyte: 'acetic', Ca: 0.1, Va: 25, Ct: 0.1, indicator: 'Fenolftaleina' },
  cell: { a: 'Zn', b: 'Cu', cCat: 1, cAn: 1, T: 298.15 },
  kin: { order: 1, k: 0.1, A0: 1, k2: 0.05, Ea: 50, Afac: 1e10 },
};

const TOOLS = [
  { id: 'gas', label: 'Gas reali' },
  { id: 'acid', label: 'Acidi e basi' },
  { id: 'cell', label: 'Elettrochimica' },
  { id: 'kin', label: 'Cinetica' },
];

export function activateLab() {
  active = true;
  document.body.dataset.mode = 'lab';
  renderAll();
}
export function deactivateLab() { active = false; }
export function labRedraw() { if (active) drawAll(); }

function renderAll() {
  renderSide();
  drawAll();
}

function drawAll() {
  if (S.tool === 'gas') drawGas();
  else if (S.tool === 'acid') drawAcid();
  else if (S.tool === 'cell') drawCell();
  else drawKinetics();
}

function slider(id, label, min, max, step, value, out) {
  return `<div class="ctl"><label class="lbl" for="${id}">${label}</label>
    <div class="range-row"><input type="range" id="${id}" min="${min}" max="${max}" step="${step}" value="${value}"><output id="${id}-out">${out}</output></div></div>`;
}

function renderSide() {
  const tabs = `<div class="seg lab-tabs">${TOOLS.map(t => `<button type="button" data-tool="${t.id}" aria-pressed="${S.tool === t.id}">${t.label}</button>`).join('')}</div>`;
  let body = '';
  if (S.tool === 'gas') {
    const g = VDW_GASES.find(x => x.id === S.gas.id);
    const { Tc } = vdwCritical(g.a, g.b);
    body = `
      <label class="lbl" for="gas-sel">Gas</label>
      <select id="gas-sel">${VDW_GASES.map(x => `<option value="${x.id}" ${x.id === S.gas.id ? 'selected' : ''}>${x.formula} · ${x.name}</option>`).join('')}</select>
      ${slider('gas-T', 'Temperatura', 0.6, 2.0, 0.01, S.gas.Tr, `${nf(S.gas.Tr * Tc, 1)} K`)}
      <p class="hint">a = ${nf(g.a, 4)} L²·bar/mol², b = ${nf(g.b, 5)} L/mol. La temperatura è mostrata anche come frazione di T<sub>c</sub>: ${nf(S.gas.Tr, 2)}.</p>`;
  } else if (S.tool === 'acid') {
    const opts = [...ACIDS.map(a => `<option value="${a.id}" ${a.id === S.acid.analyte ? 'selected' : ''}>${a.formula} · ${a.name}</option>`),
      ...BASES.filter(b => b.id !== 'NaOH').map(b => `<option value="${b.id}" ${b.id === S.acid.analyte ? 'selected' : ''}>${b.formula} · ${b.name} (base)</option>`)].join('');
    body = `
      <label class="lbl" for="acid-sel">Sostanza da titolare</label>
      <select id="acid-sel">${opts}</select>
      ${slider('acid-C', 'Concentrazione', -3, 0, 0.05, Math.log10(S.acid.Ca), `${nf(S.acid.Ca, 3)} M`)}
      ${slider('acid-V', 'Volume iniziale', 5, 50, 1, S.acid.Va, `${S.acid.Va} mL`)}
      ${slider('acid-Ct', 'Concentrazione del titolante', -3, 0, 0.05, Math.log10(S.acid.Ct), `${nf(S.acid.Ct, 3)} M`)}
      <label class="lbl" for="acid-ind">Indicatore</label>
      <select id="acid-ind">${INDICATORS.map(i => `<option ${i.name === S.acid.indicator ? 'selected' : ''}>${i.name}</option>`).join('')}</select>
      <p class="hint">Titolante: ${isBase(S.acid.analyte) ? 'HCl (acido forte)' : 'NaOH (base forte)'}.</p>`;
  } else if (S.tool === 'cell') {
    const sel = (id, v) => `<select id="${id}">${HALF_REACTIONS.map(h => `<option value="${h.id}" ${h.id === v ? 'selected' : ''}>${h.eq}  (${h.E > 0 ? '+' : ''}${nf(h.E, 2)} V)</option>`).join('')}</select>`;
    body = `
      <label class="lbl" for="cell-a">Semireazione 1</label>${sel('cell-a', S.cell.a)}
      <label class="lbl" for="cell-b">Semireazione 2</label>${sel('cell-b', S.cell.b)}
      ${slider('cell-ccat', 'Concentrazione dello ione al catodo', -6, 0, 0.1, Math.log10(S.cell.cCat), `${sci(S.cell.cCat)} M`)}
      ${slider('cell-can', 'Concentrazione dello ione all\'anodo', -6, 0, 0.1, Math.log10(S.cell.cAn), `${sci(S.cell.cAn)} M`)}
      ${slider('cell-T', 'Temperatura', 273, 373, 1, S.cell.T, `${nf(S.cell.T, 0)} K`)}
      <p class="hint">Il catodo è la coppia con il potenziale di riduzione più alto. Le concentrazioni entrano nell'equazione di Nernst (attività ≈ concentrazioni).</p>`;
  } else {
    body = `
      <label class="lbl">Ordine di reazione</label>
      <div class="seg" id="kin-order">${[0, 1, 2].map(o => `<button type="button" data-v="${o}" aria-pressed="${S.kin.order === o}">ordine ${o}</button>`).join('')}</div>
      ${slider('kin-k', 'Costante di velocità k', -3, 1, 0.05, Math.log10(S.kin.k), sci(S.kin.k))}
      ${slider('kin-A0', 'Concentrazione iniziale [A]₀', -2, 1, 0.05, Math.log10(S.kin.A0), `${sci(S.kin.A0)} M`)}
      ${slider('kin-k2', 'k₂ per A → B → C', -3, 1, 0.05, Math.log10(S.kin.k2), sci(S.kin.k2))}
      ${slider('kin-Ea', 'Energia di attivazione Eₐ', 10, 200, 1, S.kin.Ea, `${S.kin.Ea} kJ/mol`)}
      <p class="hint">Tempo in secondi. Unità di k: M·s⁻¹ (ordine 0), s⁻¹ (ordine 1), M⁻¹·s⁻¹ (ordine 2).</p>`;
  }
  $('element-card').innerHTML = `<h3 class="side-h">Laboratorio di chimica fisica</h3>${tabs}<div class="lab-inputs">${body}</div>`;
  $('element-card').querySelectorAll('.lab-tabs button').forEach(b => b.addEventListener('click', () => { S.tool = b.dataset.tool; renderAll(); }));
  const on = (id, ev, fn) => $(id)?.addEventListener(ev, fn);
  const live = (id, fn, fmt) => on(id, 'input', (e) => { fn(+e.target.value); $(`${id}-out`).textContent = fmt(); drawAll(); });
  if (S.tool === 'gas') {
    on('gas-sel', 'change', (e) => { S.gas.id = e.target.value; renderAll(); });
    live('gas-T', (v) => { S.gas.Tr = v; }, () => { const g = VDW_GASES.find(x => x.id === S.gas.id); return `${nf(S.gas.Tr * vdwCritical(g.a, g.b).Tc, 1)} K`; });
  } else if (S.tool === 'acid') {
    on('acid-sel', 'change', (e) => { S.acid.analyte = e.target.value; renderAll(); });
    on('acid-ind', 'change', (e) => { S.acid.indicator = e.target.value; drawAll(); });
    live('acid-C', (v) => { S.acid.Ca = Math.pow(10, v); }, () => `${nf(S.acid.Ca, 3)} M`);
    live('acid-V', (v) => { S.acid.Va = v; }, () => `${S.acid.Va} mL`);
    live('acid-Ct', (v) => { S.acid.Ct = Math.pow(10, v); }, () => `${nf(S.acid.Ct, 3)} M`);
  } else if (S.tool === 'cell') {
    on('cell-a', 'change', (e) => { S.cell.a = e.target.value; drawAll(); });
    on('cell-b', 'change', (e) => { S.cell.b = e.target.value; drawAll(); });
    live('cell-ccat', (v) => { S.cell.cCat = Math.pow(10, v); }, () => `${sci(S.cell.cCat)} M`);
    live('cell-can', (v) => { S.cell.cAn = Math.pow(10, v); }, () => `${sci(S.cell.cAn)} M`);
    live('cell-T', (v) => { S.cell.T = v; }, () => `${nf(S.cell.T, 0)} K`);
  } else {
    $('kin-order').querySelectorAll('button').forEach(b => b.addEventListener('click', () => { S.kin.order = +b.dataset.v; renderAll(); }));
    live('kin-k', (v) => { S.kin.k = Math.pow(10, v); }, () => sci(S.kin.k));
    live('kin-A0', (v) => { S.kin.A0 = Math.pow(10, v); }, () => `${sci(S.kin.A0)} M`);
    live('kin-k2', (v) => { S.kin.k2 = Math.pow(10, v); }, () => sci(S.kin.k2));
    live('kin-Ea', (v) => { S.kin.Ea = v; }, () => `${S.kin.Ea} kJ/mol`);
  }
}

function isBase(id) { return BASES.some(b => b.id === id); }

function setTexts({ title, sub, note, c1, n1, c2, n2, c3, n3, results, analysis }) {
  $('viewport-title').innerHTML = `${title}<small>${sub ?? ''}</small>`;
  $('viewport-legend').innerHTML = '';
  $('viewport-note').innerHTML = note ?? '';
  $('radial-title').textContent = c1; $('radial-note').innerHTML = n1 ?? '';
  $('levels-title').textContent = c2; $('levels-note').innerHTML = n2 ?? '';
  $('slice-title').textContent = c3; $('slice-note').innerHTML = n3 ?? '';
  $('controls').innerHTML = '';
  $('panel-body').innerHTML = results;
  $('analysis').innerHTML = analysis;
}

// ---------------------------------------------------------------------------
// Gas reali
// ---------------------------------------------------------------------------

function drawGas() {
  const g = VDW_GASES.find(x => x.id === S.gas.id);
  const { Tc, Pc, Vc } = vdwCritical(g.a, g.b);
  const T = S.gas.Tr * Tc;
  const TB = g.a / (R_LBAR * g.b); // temperatura di Boyle
  const Vmin = g.b * 1.12, Vmax = 10 * Vc;
  const Vs = [];
  for (let i = 0; i <= 600; i++) Vs.push(Vmin * Math.pow(Vmax / Vmin, i / 600));
  const P = Vs.map(V => vdwPressure(V, T, g.a, g.b));
  const Pid = Vs.map(V => R_LBAR * T / V);
  const mx = maxwellConstruction(T, g.a, g.b);
  // cupola di coesistenza liquido-vapore (binodale)
  const domeL = [], domeG = [];
  for (let tr = 0.55; tr < 0.9999; tr += 0.01) {
    const m = maxwellConstruction(tr * Tc, g.a, g.b);
    if (m) { domeL.push([m.Vl, m.P]); domeG.push([m.Vg, m.P]); }
  }
  const dome = [...domeL, [Vc, Pc], ...domeG.reverse()];
  const Pplot = P.map((p, i) => (mx && Vs[i] > mx.Vl && Vs[i] < mx.Vg ? p : p));
  const series = [
    { xs: dome.map(d => d[0]), ys: dome.map(d => d[1]), color: cssVar('--muted'), dash: [3, 3], label: 'coesistenza liquido–vapore', fill: `color-mix(in srgb, ${cssVar('--phase-neg')} 10%, transparent)` },
    { xs: Vs, ys: Pid, color: cssVar('--phase-neg'), dash: [6, 4], label: 'gas ideale PV = RT' },
    { xs: Vs, ys: Pplot, color: cssVar('--accent'), width: 2.4, label: 'van der Waals' },
  ];
  if (mx) series.push({ xs: [mx.Vl, mx.Vg], ys: [mx.P, mx.P], color: cssVar('--text'), width: 2, label: 'Maxwell' });
  const main = $('lab-plot');
  drawXY(main, {
    series, xmin: 0, xmax: Vmax, ymin: 0, ymax: Math.max(2.2 * Pc, mx ? mx.P * 1.5 : 0),
    xlabel: 'volume molare V (L/mol)', ylabel: 'pressione P (bar)',
    points: [{ x: Vc, y: Pc, color: cssVar('--phase-pos'), label: 'punto critico' }],
  });
  // Z in funzione di P, a varie temperature
  const zSeries = [S.gas.Tr, 1.0, 1.5, TB / Tc].map((tr, k) => {
    const TT = tr * Tc;
    const xs = [], ys = [];
    for (let i = 600; i >= 0; i--) {
      const V = Vmin * Math.pow(200 * Vc / Vmin, i / 600);
      const p = vdwPressure(V, TT, g.a, g.b);
      if (p < 0 || p > 6 * Pc) continue;
      // lungo le isoterme sotto Tc si segue il ramo stabile (gas, poi liquido)
      xs.push(p); ys.push(p * V / (R_LBAR * TT));
    }
    return { xs, ys, color: [cssVar('--accent'), cssVar('--phase-pos'), cssVar('--block-d'), cssVar('--block-f')][k], label: k === 3 ? `T Boyle ${nf(TB, 0)} K` : `${nf(TT, 0)} K`, width: k === 0 ? 2.2 : 1.4 };
  });
  drawXY($('chart-radial'), { series: zSeries, xmin: 0, xmax: 5 * Pc, ymin: 0, ymax: 1.6, xlabel: 'P (bar)', ylabel: 'Z = PV/RT', hlines: [{ y: 1, label: 'gas ideale', color: cssVar('--muted') }] });
  // tensione di vapore
  const pv = [], lnp = [], invT = [];
  for (let tr = 0.55; tr < 0.999; tr += 0.005) {
    const m = maxwellConstruction(tr * Tc, g.a, g.b);
    if (m) { pv.push([tr * Tc, m.P]); lnp.push(Math.log(m.P)); invT.push(1000 / (tr * Tc)); }
  }
  drawXY($('chart-levels'), {
    series: [{ xs: pv.map(p => p[0]), ys: pv.map(p => p[1]), color: cssVar('--accent'), width: 2.2, label: 'P di vapore (van der Waals)' }],
    xmin: pv[0][0] * 0.95, xmax: Tc * 1.05, ymin: 0, ymax: Pc * 1.1, xlabel: 'T (K)', ylabel: 'P (bar)',
    points: [{ x: Tc, y: Pc, color: cssVar('--phase-pos'), label: 'C' }, ...(mx ? [{ x: T, y: mx.P, color: cssVar('--text'), label: `${nf(mx.P, 1)} bar` }] : [])],
  });
  // Clausius–Clapeyron: ln P contro 1/T, pendenza −ΔH_vap/R
  const k1 = Math.floor(lnp.length * 0.2), k2 = Math.floor(lnp.length * 0.5);
  const slope = (lnp[k2] - lnp[k1]) / (invT[k2] - invT[k1]);
  const dHvap = -slope * R_GAS; // J/mol·(1000) → kJ/mol perché 1/T in 1000/K
  drawXY($('chart-slice'), {
    series: [{ xs: invT, ys: lnp, color: cssVar('--accent'), width: 2.2, label: 'ln P_vap' }],
    xmin: Math.min(...invT) * 0.97, xmax: Math.max(...invT) * 1.02, ymin: Math.min(...lnp) - 0.3, ymax: Math.max(...lnp) + 0.3,
    xlabel: '1000/T (K⁻¹)', ylabel: 'ln (P/bar)',
  });
  const Vgas = vdwVolumes(1.01325, T, g.a, g.b);
  setTexts({
    title: `Isoterma di van der Waals · ${g.formula}`, sub: `T = ${nf(T, 1)} K (${nf(S.gas.Tr, 2)} T_c)`,
    note: `(P + a/V²)(V − b) = RT. Il termine a/V² tiene conto delle attrazioni tra molecole, b del loro volume proprio. ${mx ? 'Sotto la temperatura critica l\'isoterma ha un\'oscillazione non fisica: la costruzione di Maxwell (aree uguali) la sostituisce con il segmento orizzontale della transizione liquido–vapore.' : 'Sopra la temperatura critica non esiste distinzione tra liquido e gas (fluido supercritico).'}`,
    c1: 'Fattore di compressibilità', n1: `Z < 1: prevalgono le attrazioni; Z > 1: prevale il volume proprio. Alla temperatura di Boyle (${nf(TB, 0)} K) Z ≈ 1 fino a pressioni elevate.`,
    c2: 'Curva della tensione di vapore', n2: 'Pressione di coesistenza liquido–vapore dalla costruzione di Maxwell; termina nel punto critico C.',
    c3: 'Clausius–Clapeyron', n3: `ln P = −ΔH<sub>vap</sub>/RT + cost. Dalla pendenza: ΔH<sub>vap</sub> ≈ ${nf(dHvap, 1)} kJ/mol (nel modello di van der Waals).`,
    results: `
      <div><h3>${g.name}</h3>
      <dl class="info-list">
        <dt>Temperatura critica T<sub>c</sub> = 8a/27Rb</dt><dd>${nf(Tc, 1)} K</dd>
        <dt>Pressione critica P<sub>c</sub> = a/27b²</dt><dd>${nf(Pc, 1)} bar</dd>
        <dt>Volume critico V<sub>c</sub> = 3b</dt><dd>${nf(Vc * 1000, 1)} mL/mol</dd>
        <dt>Z critico</dt><dd>3/8 = 0,375</dd>
        <dt>Temperatura di Boyle a/Rb</dt><dd>${nf(TB, 0)} K</dd>
        ${mx ? `<dt>Tensione di vapore a ${nf(T, 0)} K</dt><dd>${nf(mx.P, 2)} bar</dd>
        <dt>V liquido / V vapore</dt><dd>${nf(mx.Vl * 1000, 1)} / ${nf(mx.Vg * 1000, 0)} mL/mol</dd>` : ''}
        <dt>V a 1 atm (van der Waals / ideale)</dt><dd>${Vgas.length ? nf(Vgas[Vgas.length - 1], 3) : '—'} / ${nf(R_LBAR * T / 1.01325, 3)} L/mol</dd>
      </dl></div>`,
    analysis: `<h2>Le leggi</h2><div class="analysis-grid">
      <article><h3>Equazione di van der Waals (1873)</h3><p class="eq">P = RT/(V − b) − a/V²</p><p class="desc">Al punto critico l'isoterma ha un flesso orizzontale: ∂P/∂V = ∂²P/∂V² = 0. Da qui T<sub>c</sub>, P<sub>c</sub>, V<sub>c</sub> in funzione di a e b. Per CO₂ il modello dà T<sub>c</sub> = ${nf(vdwCritical(3.64, 0.04267).Tc, 1)} K (sperimentale 304,1 K).</p></article>
      <article><h3>Costruzione di Maxwell</h3><p class="eq">∫<sub>V<sub>l</sub></sub><sup>V<sub>g</sub></sup> P dV = P<sub>sat</sub>(V<sub>g</sub> − V<sub>l</sub>)</p><p class="desc">Le due fasi coesistono quando hanno lo stesso potenziale chimico: le due aree tra l'isoterma e il segmento orizzontale sono uguali. Il modello è qualitativo: per CO₂ a 280 K dà ${nf(maxwellConstruction(280, 3.64, 0.04267)?.P ?? NaN, 1)} bar contro 41,6 bar misurati.</p></article>
    </div>`,
  });
}

// ---------------------------------------------------------------------------
// Acidi e basi
// ---------------------------------------------------------------------------

function drawAcid() {
  const acid = ACIDS.find(a => a.id === S.acid.analyte);
  const base = BASES.find(b => b.id === S.acid.analyte);
  const analyte = acid ? { kind: 'acid', pKa: acid.pKa } : { kind: 'base', pKaConj: base.pKaConj };
  const nProt = acid ? acid.pKa.length : 1;
  const Veq = S.acid.Ca * S.acid.Va / S.acid.Ct; // mL per protone
  const Vmax = Veq * (nProt + 0.6);
  const Vs = [], pHs = [];
  for (let i = 0; i <= 500; i++) {
    const V = Vmax * i / 500;
    Vs.push(V);
    pHs.push(titrationPH({ analyte, Canalyte: S.acid.Ca, Vanalyte: S.acid.Va, Ctitrant: S.acid.Ct, Vadded: V }));
  }
  const ind = INDICATORS.find(i => i.name === S.acid.indicator);
  const eqs = [];
  for (let k = 1; k <= nProt; k++) {
    const V = k * Veq;
    const pH = titrationPH({ analyte, Canalyte: S.acid.Ca, Vanalyte: S.acid.Va, Ctitrant: S.acid.Ct, Vadded: V });
    eqs.push({ V, pH });
  }
  const pH0 = pHs[0];
  const strong = acid?.id === 'HCl';
  const halfPoints = acid && !strong ? acid.pKa.map((pk, k) => ({ V: (k + 0.5) * Veq, pH: titrationPH({ analyte, Canalyte: S.acid.Ca, Vanalyte: S.acid.Va, Ctitrant: S.acid.Ct, Vadded: (k + 0.5) * Veq }) })) : [];
  drawXY($('lab-plot'), {
    series: [{ xs: Vs, ys: pHs, color: cssVar('--accent'), width: 2.6, label: 'pH calcolato' }],
    xmin: 0, xmax: Vmax, ymin: 0, ymax: 14,
    xlabel: `volume di ${acid ? 'NaOH' : 'HCl'} aggiunto (mL)`, ylabel: 'pH',
    hbands: [{ y0: ind.from, y1: ind.to, color: `color-mix(in srgb, ${ind.colors[1]} 22%, transparent)`, label: `${ind.name}: ${nf(ind.from, 1)}–${nf(ind.to, 1)}` }],
    vlines: eqs.map((e, k) => ({ x: e.V, label: `equivalenza ${nProt > 1 ? k + 1 : ''}`, color: cssVar('--phase-pos'), row: k % 2 })),
    points: [...eqs.map(e => ({ x: e.V, y: e.pH, color: cssVar('--phase-pos'), label: `pH ${nf(e.pH, 2)}` })),
      ...halfPoints.map((p, k) => ({ x: p.V, y: p.pH, color: cssVar('--phase-neg'), label: `pH = pKa${nProt > 1 ? k + 1 : ''}` }))],
  });
  // diagramma di distribuzione
  const pHaxis = [];
  for (let p = 0; p <= 14.001; p += 0.05) pHaxis.push(p);
  let distSeries;
  if (acid) {
    const names = [];
    const n = acid.pKa.length;
    for (let i = 0; i <= n; i++) names.push(i === 0 ? `H${n > 1 ? subN(n) : ''}A` : `${n - i > 0 ? `H${n - i > 1 ? subN(n - i) : ''}` : ''}A${supCharge(i)}`);
    const colors = [cssVar('--phase-neg'), cssVar('--accent'), cssVar('--block-f'), cssVar('--phase-pos')];
    distSeries = names.map((name, i) => ({ xs: pHaxis, ys: pHaxis.map(p => fractions(Math.pow(10, -p), acid.pKa)[i]), color: colors[i], label: name, width: 2 }));
  } else {
    const Ka = Math.pow(10, -base.pKaConj);
    distSeries = [
      { xs: pHaxis, ys: pHaxis.map(p => { const H = Math.pow(10, -p); return H / (H + Ka); }), color: cssVar('--phase-neg'), label: 'BH⁺', width: 2 },
      { xs: pHaxis, ys: pHaxis.map(p => { const H = Math.pow(10, -p); return Ka / (H + Ka); }), color: cssVar('--accent'), label: 'B', width: 2 },
    ];
  }
  drawXY($('chart-radial'), { series: distSeries, xmin: 0, xmax: 14, ymin: 0, ymax: 1.05, xlabel: 'pH', ylabel: 'frazione α' });
  // derivata dpH/dV
  const dxs = [], dys = [];
  for (let i = 1; i < Vs.length - 1; i++) { dxs.push(Vs[i]); dys.push((pHs[i + 1] - pHs[i - 1]) / (Vs[i + 1] - Vs[i - 1])); }
  drawXY($('chart-levels'), { series: [{ xs: dxs, ys: dys.map(Math.abs), color: cssVar('--accent'), width: 2, label: '|dpH/dV|' }], xmin: 0, xmax: Vmax, ymin: 0, ymax: Math.max(...dys.map(Math.abs)) * 1.1, xlabel: 'V (mL)', ylabel: 'pH/mL' });
  // capacità tampone β = dC_titolante/dpH
  const bxs = [], bys = [];
  for (let i = 1; i < Vs.length - 1; i++) {
    const dn = S.acid.Ct * (Vs[i + 1] - Vs[i - 1]) / (S.acid.Va + Vs[i]);
    const dp = Math.abs(pHs[i + 1] - pHs[i - 1]);
    if (dp > 1e-9) { bxs.push(Vs[i]); bys.push(dn / dp); }
  }
  drawXY($('chart-slice'), { series: [{ xs: bxs, ys: bys, color: cssVar('--block-f'), width: 2, label: 'β = dC/dpH' }], xmin: 0, xmax: Vmax, ymin: 0, ymax: Math.max(...bys) * 1.1, xlabel: 'V (mL)', ylabel: 'mol L⁻¹ pH⁻¹' });
  const goodInd = eqs.some(e => e.pH >= ind.from - 0.3 && e.pH <= ind.to + 0.3);
  const name = acid ? acid.name : base.name;
  setTexts({
    title: `Titolazione: ${acid ? acid.formula : base.formula}`, sub: `${nf(S.acid.Ca, 3)} M, ${S.acid.Va} mL con ${acid ? 'NaOH' : 'HCl'} ${nf(S.acid.Ct, 3)} M`,
    note: 'pH calcolato esattamente dal bilancio di carica e dalle costanti di equilibrio (anche l\'autoprotolisi dell\'acqua), senza le approssimazioni dei casi limite.',
    c1: 'Diagramma di distribuzione', n1: 'Frazione di ciascuna specie in funzione del pH. Le curve si incrociano a pH = pKa.',
    c2: 'Derivata della curva di titolazione', n2: 'Il massimo di dpH/dV individua il punto di equivalenza.',
    c3: 'Capacità tampone', n3: 'β è massima dove pH ≈ pKa: la soluzione resiste meglio all\'aggiunta di acido o base (zona tampone).',
    results: `
      <div><h3>${name}</h3>
      <dl class="info-list">
        ${acid ? acid.pKa.map((p, k) => `<dt>pKa${acid.pKa.length > 1 ? k + 1 : ''} (Ka)</dt><dd>${strong ? 'acido forte' : `${nf(p, 2)} (${sci(Math.pow(10, -p))})`}</dd>`).join('') : `<dt>pKb (pKa di BH⁺)</dt><dd>${nf(14 - base.pKaConj, 2)} (${nf(base.pKaConj, 2)})</dd>`}
        <dt>pH iniziale</dt><dd>${nf(pH0, 2)}</dd>
        ${eqs.map((e, k) => `<dt>Equivalenza ${nProt > 1 ? k + 1 : ''}</dt><dd>${nf(e.V, 2)} mL, pH ${nf(e.pH, 2)}</dd>`).join('')}
        ${halfPoints.map((p, k) => `<dt>Semi-equivalenza ${nProt > 1 ? k + 1 : ''}</dt><dd>pH ${nf(p.pH, 2)}</dd>`).join('')}
        <dt>${ind.name}</dt><dd>${goodInd ? 'adatto' : 'non adatto'}</dd>
      </dl>
      <p class="desc-muted" style="margin-top:6px">Un indicatore è adatto se il suo intervallo di viraggio cade nel salto di pH attorno all'equivalenza.</p></div>`,
    analysis: `<h2>Le leggi</h2><div class="analysis-grid">
      <article><h3>Equilibri e bilancio di carica</h3><p class="eq">K<sub>a</sub> = [H⁺][A⁻]/[HA],&nbsp;&nbsp; K<sub>w</sub> = [H⁺][OH⁻] = 10⁻¹⁴</p><p class="eq">[H⁺] + [Na⁺] = [OH⁻] + Σ i·α<sub>i</sub>·C<sub>A</sub></p><p class="desc">Frazioni delle specie: α<sub>i</sub> = [H⁺]<sup>n−i</sup> K<sub>1</sub>…K<sub>i</sub> / Σ<sub>j</sub>[H⁺]<sup>n−j</sup> K<sub>1</sub>…K<sub>j</sub>. L'equazione si risolve numericamente per [H⁺] a ogni volume aggiunto.</p></article>
      <article><h3>Henderson–Hasselbalch</h3><p class="eq">pH = pKa + log([A⁻]/[HA])</p><p class="desc">È un caso particolare del bilancio di carica, valido nella zona tampone. A metà titolazione [A⁻] = [HA] e quindi pH = pKa: il grafico lo mostra con i punti blu.</p></article>
    </div>`,
  });
}

function subN(n) { return String(n).split('').map(d => '₀₁₂₃₄₅₆₇₈₉'[+d]).join(''); }
function supCharge(i) { return (i > 1 ? '⁰¹²³⁴⁵⁶⁷⁸⁹'[i] : '') + '⁻'; }

// ---------------------------------------------------------------------------
// Elettrochimica
// ---------------------------------------------------------------------------

function drawCell() {
  const h1 = HALF_REACTIONS.find(h => h.id === S.cell.a);
  const h2 = HALF_REACTIONS.find(h => h.id === S.cell.b);
  const c = galvanicCell(h1, h2, { cOxCathode: S.cell.cCat, cOxAnode: S.cell.cAn, T: S.cell.T });
  drawCellDiagram($('lab-plot'), c);
  // Nernst: E contro log Q
  const lq = [], E = [];
  for (let x = -8; x <= 8; x += 0.1) { lq.push(x); E.push(c.E0 - R_GAS * S.cell.T / (c.n * FARADAY) * Math.log(Math.pow(10, x))); }
  drawXY($('chart-radial'), {
    series: [{ xs: lq, ys: E, color: cssVar('--accent'), width: 2.2, label: 'E = E° − (RT/nF) ln Q' }],
    xmin: -8, xmax: 8, ymin: Math.min(...E, 0) - 0.05, ymax: Math.max(...E, 0) + 0.05, xlabel: 'log₁₀ Q', ylabel: 'E (V)',
    points: [{ x: Math.log10(c.Q), y: c.E, color: cssVar('--phase-pos'), label: `${nf(c.E, 3)} V` }],
    hlines: [{ y: 0, label: 'equilibrio: E = 0, Q = K', color: cssVar('--muted') }],
  });
  // serie elettrochimica
  drawSeries($('chart-levels'), c);
  // scarica della pila: E in funzione della frazione di reazione
  const fx = [], fy = [];
  const ca0 = S.cell.cAn, cc0 = S.cell.cCat;
  for (let f = 0; f <= 0.999; f += 0.005) {
    const cc = cc0 * (1 - f);
    const ca = ca0 + cc0 * f * (c.cathode.n / c.anode.n);
    const Q = Math.pow(ca, c.n / c.anode.n) / Math.pow(cc, c.n / c.cathode.n);
    fx.push(f * 100); fy.push(c.E0 - R_GAS * S.cell.T / (c.n * FARADAY) * Math.log(Q));
  }
  drawXY($('chart-slice'), { series: [{ xs: fx, ys: fy, color: cssVar('--block-f'), width: 2.2, label: 'E durante la scarica' }], xmin: 0, xmax: 100, ymin: 0, ymax: Math.max(...fy) * 1.1, xlabel: 'ione del catodo consumato (%)', ylabel: 'E (V)' });
  const notation = `${c.anode.red} | ${c.anode.ox} || ${c.cathode.ox} | ${c.cathode.red}`;
  const log10K = c.lnK / Math.LN10;
  setTexts({
    title: `Pila galvanica`, sub: notation,
    note: 'All\'anodo (−) avviene l\'ossidazione, al catodo (+) la riduzione. Gli elettroni scorrono nel circuito esterno dall\'anodo al catodo; il ponte salino chiude il circuito con il moto degli ioni.',
    c1: 'Equazione di Nernst', n1: `A 25 °C: E = E° − (0,0592/n) log Q. La pila si scarica finché E = 0, cioè Q = K.`,
    c2: 'Serie elettrochimica', n2: 'Potenziali standard di riduzione. Più in alto: ossidanti più forti. Più in basso: riducenti più forti.',
    c3: 'Scarica della pila', n3: 'Mentre la pila lavora gli ioni del catodo si consumano e quelli dell\'anodo si accumulano: Q cresce ed E cala, lentamente fino alla fine.',
    results: `
      <div><h3>${notation}</h3>
      <dl class="info-list">
        <dt>Catodo (riduzione, +)</dt><dd>${c.cathode.eq}</dd>
        <dt>Anodo (ossidazione, −)</dt><dd>${c.anode.eq.replace('→', '←')}</dd>
        <dt>Elettroni scambiati n</dt><dd>${c.n}</dd>
        <dt>E° = E°<sub>cat</sub> − E°<sub>an</sub></dt><dd>${nf(c.E0, 3)} V</dd>
        <dt>E (Nernst)</dt><dd><b>${nf(c.E, 3)} V</b></dd>
        <dt>ΔG° = −nFE°</dt><dd>${nf(c.dG0, 1)} kJ/mol</dd>
        <dt>ΔG = −nFE</dt><dd>${nf(c.dG, 1)} kJ/mol</dd>
        <dt>K = e<sup>nFE°/RT</sup></dt><dd>10<sup>${nf(log10K, 1)}</sup></dd>
      </dl></div>`,
    analysis: `<h2>Le leggi</h2><div class="analysis-grid">
      <article><h3>Energia libera e lavoro elettrico</h3><p class="eq">ΔG = −nFE,&nbsp;&nbsp; F = 96 485 C/mol</p><p class="desc">Il lavoro elettrico massimo che una pila può compiere è uguale alla diminuzione di energia libera della reazione. E > 0 ⇔ ΔG < 0: reazione spontanea.</p></article>
      <article><h3>Equazione di Nernst (1889)</h3><p class="eq">E = E° − (RT/nF) ln Q,&nbsp;&nbsp; E° = (RT/nF) ln K</p><p class="desc">Viene da ΔG = ΔG° + RT ln Q. Nella pila di concentrazione (stessa coppia ai due elettrodi) E° = 0 e la tensione nasce solo dalla differenza di concentrazione.</p></article>
    </div>`,
  });
}

function drawCellDiagram(canvas, c) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  canvas.width = w * dpr; canvas.height = h * dpr;
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  const ink = cssVar('--text'), muted = cssVar('--muted'), accent = cssVar('--accent');
  const bw = Math.min(200, w * 0.28), bh = Math.min(200, h * 0.42);
  const y0 = h * 0.42;
  const xa = w * 0.25 - bw / 2, xc = w * 0.75 - bw / 2;
  const beaker = (x, label, ion, color) => {
    ctx.fillStyle = `color-mix(in srgb, ${color} 18%, transparent)`;
    ctx.fillRect(x, y0 + bh * 0.25, bw, bh * 0.75);
    ctx.strokeStyle = muted; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y0 + bh); ctx.lineTo(x + bw, y0 + bh); ctx.lineTo(x + bw, y0); ctx.stroke();
    ctx.fillStyle = muted; ctx.font = '12px "IBM Plex Sans", sans-serif'; ctx.textAlign = 'center';
    ctx.fillText(ion, x + bw / 2, y0 + bh * 0.85);
    ctx.fillStyle = ink; ctx.font = '600 13px "IBM Plex Sans", sans-serif';
    ctx.fillText(label, x + bw / 2, y0 + bh + 20);
  };
  beaker(xa, `Anodo (−): ossidazione`, `${c.anode.ox} ${sci(S.cell.cAn)} M`, cssVar('--phase-neg'));
  beaker(xc, `Catodo (+): riduzione`, `${c.cathode.ox} ${sci(S.cell.cCat)} M`, cssVar('--phase-pos'));
  // elettrodi
  const elec = (x, name) => {
    ctx.fillStyle = cssVar('--scene-axis');
    ctx.fillRect(x + bw / 2 - 9, y0 - 40, 18, bh * 0.7 + 40);
    ctx.fillStyle = ink; ctx.font = '600 12px "IBM Plex Mono", monospace'; ctx.textAlign = 'center';
    ctx.fillText(name, x + bw / 2, y0 - 46);
  };
  elec(xa, c.anode.metal ? c.anode.red : `Pt (${c.anode.red})`);
  elec(xc, c.cathode.metal ? c.cathode.red : `Pt (${c.cathode.red})`);
  // ponte salino
  ctx.strokeStyle = muted; ctx.lineWidth = 10; ctx.lineCap = 'round';
  ctx.beginPath(); ctx.moveTo(xa + bw * 0.8, y0 + bh * 0.45); ctx.lineTo(xa + bw * 0.8, y0 - 10); ctx.lineTo(xc + bw * 0.2, y0 - 10); ctx.lineTo(xc + bw * 0.2, y0 + bh * 0.45); ctx.stroke();
  ctx.lineCap = 'butt';
  ctx.fillStyle = muted; ctx.font = '11px "IBM Plex Sans", sans-serif'; ctx.textAlign = 'center';
  ctx.fillText('ponte salino (KNO₃)', w / 2, y0 + 8);
  // circuito esterno con voltmetro
  const top = y0 - 95;
  ctx.strokeStyle = ink; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(xa + bw / 2, y0 - 40); ctx.lineTo(xa + bw / 2, top); ctx.lineTo(w / 2 - 36, top);
  ctx.moveTo(w / 2 + 36, top); ctx.lineTo(xc + bw / 2, top); ctx.lineTo(xc + bw / 2, y0 - 40); ctx.stroke();
  ctx.beginPath(); ctx.arc(w / 2, top, 36, 0, 2 * Math.PI); ctx.stroke();
  ctx.fillStyle = accent; ctx.font = '700 18px "IBM Plex Mono", monospace'; ctx.textBaseline = 'middle';
  ctx.fillText(`${nf(c.E, 3)} V`, w / 2, top);
  ctx.textBaseline = 'alphabetic';
  // freccia degli elettroni
  ctx.fillStyle = accent; ctx.font = '12px "IBM Plex Sans", sans-serif';
  ctx.fillText('e⁻ →', w * 0.34, top - 8);
  ctx.fillText('e⁻ →', w * 0.66, top - 8);
}

function drawSeries(canvas, c) {
  const list = [...HALF_REACTIONS].sort((a, b) => b.E - a.E);
  const xs = list.map((_, i) => i);
  drawXY(canvas, {
    series: [{ xs, ys: list.map(h => h.E), color: cssVar('--muted'), width: 1 }],
    xmin: -0.5, xmax: list.length - 0.5, ymin: -3.3, ymax: 3.1, xlabel: 'semireazioni in ordine di E°', ylabel: 'E° (V)', legend: false, noXTicks: true,
    points: list.map((h, i) => ({ x: i, y: h.E, r: h === c.cathode || h === c.anode ? 5 : 3, color: h === c.cathode ? cssVar('--phase-pos') : h === c.anode ? cssVar('--phase-neg') : cssVar('--muted'), label: h === c.cathode || h === c.anode ? `${h.ox}/${h.red}` : '' })),
    hlines: [{ y: 0, label: 'SHE', color: cssVar('--muted') }],
  });
}

// ---------------------------------------------------------------------------
// Cinetica
// ---------------------------------------------------------------------------

function drawKinetics() {
  const { order, k, A0, k2, Ea } = S.kin;
  const t12 = halfLife(order, A0, k);
  const tmax = Math.min(8 * t12, order === 0 ? 1.3 * A0 / k : 8 * t12);
  const ts = [], As = [];
  for (let i = 0; i <= 400; i++) { const t = tmax * i / 400; ts.push(t); As.push(integratedRate(order, A0, k, t)); }
  drawXY($('lab-plot'), {
    series: [{ xs: ts, ys: As, color: cssVar('--accent'), width: 2.6, label: `[A](t), ordine ${order}` }],
    xmin: 0, xmax: tmax, ymin: 0, ymax: A0 * 1.05, xlabel: 't (s)', ylabel: '[A] (M)',
    vlines: [1, 2, 3].filter(n => n * t12 < tmax && (order === 1 || n === 1)).map(n => ({ x: n * t12, label: `${n > 1 ? n : ''}t½`, color: cssVar('--phase-neg') })),
  });
  // grafico linearizzato
  const lin = order === 0 ? As : order === 1 ? As.map(a => Math.log(a)) : As.map(a => 1 / a);
  const finite = lin.filter(Number.isFinite);
  drawXY($('chart-radial'), {
    series: [{ xs: ts, ys: lin, color: cssVar('--accent'), width: 2, label: order === 0 ? '[A]' : order === 1 ? 'ln [A]' : '1/[A]' }],
    xmin: 0, xmax: tmax, ymin: Math.min(...finite), ymax: Math.max(...finite), xlabel: 't (s)', ylabel: order === 0 ? '[A]' : order === 1 ? 'ln [A]' : '1/[A] (M⁻¹)',
  });
  // reazioni consecutive
  const tt = [], a = [], b = [], cc = [];
  const T2 = 6 / Math.min(k, k2);
  for (let i = 0; i <= 400; i++) { const t = T2 * i / 400; const r = consecutive(A0, k, k2, t); tt.push(t); a.push(r.A); b.push(r.B); cc.push(r.C); }
  const tBmax = Math.abs(k - k2) > 1e-12 ? Math.log(k / k2) / (k - k2) : 1 / k;
  drawXY($('chart-levels'), {
    series: [
      { xs: tt, ys: a, color: cssVar('--phase-neg'), label: 'A', width: 2 },
      { xs: tt, ys: b, color: cssVar('--accent'), label: 'B (intermedio)', width: 2 },
      { xs: tt, ys: cc, color: cssVar('--block-f'), label: 'C', width: 2 },
    ],
    xmin: 0, xmax: T2, ymin: 0, ymax: A0 * 1.05, xlabel: 't (s)', ylabel: 'M',
    vlines: [{ x: tBmax, label: 'max di B', color: cssVar('--accent') }],
  });
  // Arrhenius: k(T) con fattore pre-esponenziale tale che k(298 K) = k
  const Afac = k / Math.exp(-Ea * 1000 / (R_GAS * 298.15));
  const Ts = [], ks = [];
  for (let T = 250; T <= 500; T += 2) { Ts.push(T); ks.push(Afac * Math.exp(-Ea * 1000 / (R_GAS * T))); }
  drawXY($('chart-slice'), {
    series: [{ xs: Ts, ys: ks, color: cssVar('--accent'), width: 2, label: `k(T), Eₐ = ${Ea} kJ/mol` }],
    xmin: 250, xmax: 500, ymin: Math.max(Math.min(...ks), 1e-30), ymax: Math.max(...ks) * 2, xlabel: 'T (K)', ylabel: 'k', logY: true,
    vlines: [{ x: 298.15, label: '25 °C', color: cssVar('--muted') }, { x: 308.15, label: '+10 °C', color: cssVar('--muted'), row: 1 }],
  });
  const ratio10 = Math.exp(-Ea * 1000 / R_GAS * (1 / 308.15 - 1 / 298.15));
  const law = order === 0 ? '[A] = [A]₀ − kt' : order === 1 ? '[A] = [A]₀ e<sup>−kt</sup>' : '1/[A] = 1/[A]₀ + kt';
  const t90 = order === 0 ? 0.9 * A0 / k : order === 1 ? Math.log(10) / k : 9 / (k * A0);
  setTexts({
    title: `Legge cinetica di ordine ${order}`, sub: `v = k[A]${order === 1 ? '' : order === 0 ? '⁰' : '²'}`,
    note: `Soluzione esatta dell'equazione differenziale d[A]/dt = −k[A]${order === 1 ? '' : order === 0 ? '⁰' : '²'}: ${law}.`,
    c1: 'Grafico linearizzato', n1: `Una retta conferma l'ordine di reazione: ${order === 0 ? '[A]' : order === 1 ? 'ln [A]' : '1/[A]'} contro t ha pendenza ${order === 2 ? '+' : '−'}k.`,
    c2: 'Reazioni consecutive A → B → C', n2: `Soluzione di Bateman. L'intermedio B raggiunge il massimo a t = ln(k₁/k₂)/(k₁ − k₂) = ${sci(tBmax)} s.`,
    c3: 'Equazione di Arrhenius', n3: `k = A e<sup>−Eₐ/RT</sup>. Con Eₐ = ${Ea} kJ/mol, salendo da 25 °C a 35 °C la velocità si moltiplica per ${nf(ratio10, 2)}.`,
    results: `
      <div><h3>Risultati</h3>
      <dl class="info-list">
        <dt>Tempo di dimezzamento t½</dt><dd>${sci(t12)} s</dd>
        <dt>Formula di t½</dt><dd>${order === 0 ? '[A]₀/2k' : order === 1 ? 'ln 2/k' : '1/(k[A]₀)'}</dd>
        <dt>Tempo per consumare il 90%</dt><dd>${sci(t90)} s</dd>
        <dt>Velocità iniziale</dt><dd>${sci(k * Math.pow(A0, order))} M/s</dd>
        <dt>Fattore di Arrhenius A</dt><dd>${sci(Afac)}</dd>
      </dl>
      <p class="desc-muted" style="margin-top:6px">Solo per l'ordine 1 il tempo di dimezzamento non dipende dalla concentrazione iniziale (decadimento radioattivo).</p></div>`,
    analysis: `<h2>Le leggi</h2><div class="analysis-grid">
      <article><h3>Leggi cinetiche integrate</h3><p class="eq">ordine 0: [A] = [A]₀ − kt&nbsp;&nbsp; ordine 1: ln[A] = ln[A]₀ − kt&nbsp;&nbsp; ordine 2: 1/[A] = 1/[A]₀ + kt</p><p class="desc">Si ottengono integrando v = −d[A]/dt = k[A]ⁿ. L'ordine si determina sperimentalmente cercando quale grafico dà una retta.</p></article>
      <article><h3>Dalla teoria dello stato di transizione</h3><p class="eq">k = (k<sub>B</sub>T/h) e<sup>−ΔG‡/RT</sup></p><p class="desc">L'equazione di Eyring collega k alla barriera di energia libera. Nella modalità Reazioni la barriera viene calcolata dalla meccanica quantistica per alcune reazioni reali.</p></article>
    </div>`,
  });
}
