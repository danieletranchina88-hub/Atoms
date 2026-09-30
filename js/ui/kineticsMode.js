// Modalità "Cinetica": meccanismi scritti liberamente, integrati con un metodo implicito per sistemi rigidi.
// Si scelgono concentrazioni iniziali, temperatura e durata; si ottengono le curve c(t), l'ordine apparente
// dalle leggi integrate, il tempo di dimezzamento, il grafico di Arrhenius e i parametri di Eyring.

import {
  parseMechanism, simulate, rateConstant, orderAnalysis, eyring, balance, buildModel, KINETICS_PRESETS, R_GAS,
} from '../chem/kinetics.js';
import { drawXY } from './chemCharts.js';

const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const nf = (v, d = 2) => (v === null || v === undefined || !Number.isFinite(v) ? '—' : Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d, useGrouping: Math.abs(v) >= 1e4 }));
const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const sup = (n) => String(n).split('').map(c => SUP[c] ?? c).join('');
const sci = (v, d = 3) => {
  if (!Number.isFinite(v)) return '—';
  if (v === 0) return '0';
  const e = Math.floor(Math.log10(Math.abs(v)));
  if (e >= -2 && e <= 3) return nf(v, Math.max(0, d - 1 - Math.max(e, 0)) + Math.max(0, -e));
  return `${nf(v / Math.pow(10, e), d - 1)}·10${sup(e)}`;
};
const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const SUBS = '₀₁₂₃₄₅₆₇₈₉';
const SUPS = '⁰¹²³⁴⁵⁶⁷⁸⁹';
/** Formula leggibile: pedici, carica in apice; gli isotopi (Bi-210) restano come sono. */
function pretty(sp) {
  if (/-\d+$/.test(sp)) return esc(sp);
  let body = sp, charge = '';
  const m = sp.match(/^(.*?)(\d*)([+-])$/);
  if (m) {
    const mono = /^[A-Z][a-z]?$/.test(m[1]);
    body = mono ? m[1] : m[1] + m[2];
    const n = mono ? m[2] : '';
    charge = n.split('').map(d => SUPS[+d]).join('') + (m[3] === '+' ? '⁺' : '⁻');
  }
  return esc(body).replace(/([A-Za-z)\]])(\d+)/g, (x, a, d) => a + d.split('').map(c => SUBS[+c]).join('')) + charge;
}

const TIME_UNITS = [['s', 1], ['min', 60], ['h', 3600], ['d', 86400], ['anni', 365.25 * 86400]];
// tavolozza categoriale leggibile sia sul tema scuro sia su quello chiaro
const COLORS = ['#e8a33d', '#4aa3df', '#e0605e', '#5cb87a', '#a88ee0', '#d980b8', '#45b5b5', '#b8a040', '#8c9aa8', '#f08a4b'];

let active = false;
const K = {
  preset: KINETICS_PRESETS[0],
  text: KINETICS_PRESETS[0].mech,
  mech: null,
  c0: { ...KINETICS_PRESETS[0].c0 },
  T: KINETICS_PRESETS[0].T,
  tEnd: KINETICS_PRESETS[0].tEnd,
  logT: false, logY: false,
  focus: KINETICS_PRESETS[0].focus,
  result: null,
  hidden: new Set(),
  started: false,
};

export function activateKinetics() {
  active = true;
  document.body.dataset.mode = 'kinetics';
  $('busy').hidden = true;
  if (!K.started) { K.started = true; loadPreset(K.preset); }
  renderAll();
}
export function deactivateKinetics() { active = false; }
export function kineticsRedraw() { if (active) { drawMain(); drawCharts(); } }

function loadPreset(p) {
  K.preset = p;
  K.text = p.mech;
  K.c0 = { ...p.c0 };
  K.T = p.T;
  K.tEnd = p.tEnd;
  K.logT = !!p.logT; K.logY = !!p.logY;
  K.focus = p.focus;
  K.hidden = new Set();
  compile();
  run();
}

function compile() {
  K.mech = parseMechanism(K.text);
  for (const s of K.mech.species) if (K.c0[s] === undefined) K.c0[s] = 0;
  if (!K.mech.species.includes(K.focus)) K.focus = K.mech.species.find(s => !K.mech.fixed.has(s)) ?? null;
}

function run() {
  K.result = null;
  if (!K.mech || K.mech.errors.length || !K.mech.reactions.length) return;
  const t0 = performance.now();
  try {
    K.result = simulate(K.mech, K.c0, K.T, K.tEnd);
    K.result.ms = performance.now() - t0;
  } catch (err) {
    K.mech.errors.push(`errore nell'integrazione: ${err.message}`);
  }
}

// ---------------------------------------------------------------------------
// Pannelli
// ---------------------------------------------------------------------------

function renderAll() {
  if (!active) return;
  renderSide();
  renderControls();
  renderPanel();
  drawMain();
  drawCharts();
  renderAnalysis();
}

function renderSide() {
  const m = K.mech;
  $('element-card').innerHTML = `
    <h3 class="side-h">Meccanismo di reazione</h3>
    <label class="lbl" for="kn-preset">Esperienza</label>
    <select id="kn-preset">${KINETICS_PRESETS.map(p => `<option value="${p.id}" ${p.id === K.preset.id ? 'selected' : ''}>${esc(p.name)}</option>`).join('')}</select>
    <p class="mol-note">${K.preset.text}</p>
    <label class="lbl" for="kn-text">Reazioni (una per riga)</label>
    <textarea id="kn-text" class="mech-input" spellcheck="false" rows="12">${esc(K.text)}</textarea>
    <div class="btn-row"><button type="button" class="btn" id="kn-apply">Applica e simula</button><button type="button" class="btn" id="kn-reset">Ripristina</button></div>
    ${m?.errors.length ? `<p class="hint warn">${m.errors.map(esc).join('<br>')}</p>` : ''}
    <details class="hint"><summary>Sintassi</summary>
      <p><code>A + B -> C ; k=0.5</code> costante fissa<br>
      <code>A -> B ; A=1e13 Ea=100</code> Arrhenius, E<sub>a</sub> in kJ/mol (o <code>EaR=</code> in K)<br>
      <code>A + M -> B + M ; k0=6e-34 n=-2.4</code> k = k₀(T/300)ⁿ<br>
      <code>A <=> B ; k=2 kr=1</code> reversibile (anche <code>Ar=</code>, <code>Ear=</code>)<br>
      <code>X -> Y ; t12=3.5 h</code> tempo di dimezzamento (s, min, h, d, a)<br>
      <code>2 A -> B ; k=1 ord.A=1</code> ordine diverso dalla stechiometria<br>
      <code>costanti: M, H2O</code> specie a concentrazione fissa<br>
      <code>unità: molecole/cm3</code> per la fase gas (altrimenti mol/L)<br>
      Coefficienti frazionari ammessi (<code>1/2 O2</code>); il tempo è in secondi.</p>
    </details>`;
  $('kn-preset').addEventListener('change', (e) => { loadPreset(KINETICS_PRESETS.find(p => p.id === e.target.value)); renderAll(); });
  $('kn-apply').addEventListener('click', () => { K.text = $('kn-text').value; compile(); run(); renderAll(); });
  $('kn-reset').addEventListener('click', () => { loadPreset(K.preset); renderAll(); });
}

const fmtIn = (v) => (v === 0 ? '0' : Math.abs(v) >= 1e5 || Math.abs(v) < 1e-3 ? v.toExponential(4).replace(/\.?0+e/, 'e') : String(+v.toPrecision(6)));

function bestUnit(t) {
  let u = TIME_UNITS[0];
  for (const x of TIME_UNITS) if (t / x[1] >= 1) u = x;
  return u;
}

function renderControls() {
  const m = K.mech;
  const cu = m?.units === 'molec' ? 'molecole/cm³' : 'mol/L';
  const [un, uf] = bestUnit(K.tEnd);
  $('controls').innerHTML = `
    <div class="ctl"><span class="lbl">Concentrazioni iniziali (${cu})</span>
      <div class="kn-grid">${(m?.species ?? []).map(s => `<label title="${m.fixed.has(s) ? 'concentrazione costante' : ''}">${pretty(s)}${m.fixed.has(s) ? ' <span class="desc-muted">fissa</span>' : ''}</label><input type="text" inputmode="decimal" data-sp="${esc(s)}" value="${fmtIn(K.c0[s] ?? 0)}">`).join('')}</div>
    </div>
    <div class="ctl"><span class="lbl">Condizioni</span>
      <div class="kn-grid">
        <label for="kn-T">T (K)</label><input id="kn-T" type="text" inputmode="decimal" value="${+K.T.toFixed(2)}">
        <label for="kn-tend">durata</label><div class="range-row"><input id="kn-tend" type="text" inputmode="decimal" value="${+(K.tEnd / uf).toPrecision(6)}" style="width:90px"><select id="kn-tunit">${TIME_UNITS.map(([n, f]) => `<option value="${f}" ${n === un ? 'selected' : ''}>${n}</option>`).join('')}</select></div>
      </div>
      <p class="hint">${nf(K.T - 273.15, 2)} °C. Le costanti con A ed E<sub>a</sub> seguono Arrhenius; quelle date come k o t½ restano fisse.</p>
      <div class="btn-row"><button type="button" class="btn" id="kn-run">Simula</button></div>
    </div>
    <div class="ctl"><span class="lbl">Grafico</span>
      <div class="seg"><button type="button" id="kn-logt" aria-pressed="${K.logT}">tempo logaritmico</button><button type="button" id="kn-logy" aria-pressed="${K.logY}">concentrazioni logaritmiche</button></div>
    </div>`;
  const readAll = () => {
    $('controls').querySelectorAll('[data-sp]').forEach(inp => { const v = parseFloat(String(inp.value).replace(',', '.')); if (Number.isFinite(v) && v >= 0) K.c0[inp.dataset.sp] = v; });
    const T = parseFloat(String($('kn-T').value).replace(',', '.'));
    if (T > 0) K.T = T;
    const d = parseFloat(String($('kn-tend').value).replace(',', '.'));
    if (d > 0) K.tEnd = d * +$('kn-tunit').value;
  };
  $('kn-run').addEventListener('click', () => { readAll(); run(); renderAll(); });
  $('controls').querySelectorAll('input').forEach(i => i.addEventListener('keydown', (e) => { if (e.key === 'Enter') { readAll(); run(); renderAll(); } }));
  $('kn-logt').addEventListener('click', () => { K.logT = !K.logT; renderControls(); drawMain(); });
  $('kn-logy').addEventListener('click', () => { K.logY = !K.logY; renderControls(); drawMain(); });
}

function kUnits(order, units) {
  const n = order;
  if (Math.abs(n - 1) < 1e-9) return 's⁻¹';
  const p = 1 - n;
  const fmt = (x) => (Number.isInteger(x) ? String(x) : x.toFixed(2).replace('.', ','));
  if (units === 'molec') return `cm${sup(fmt(3 * (n - 1)))} molecole${sup(fmt(p))} s⁻¹`;
  return `M${sup(fmt(p))} s⁻¹`;
}

function stepRows() {
  const m = K.mech;
  const rows = [];
  for (const r of m.reactions) {
    const ord = r.reactants.map(([s, v]) => r.orders[s] ?? v).reduce((a, b) => a + b, 0) + Object.entries(r.orders).filter(([s]) => !r.reactants.some(x => x[0] === s)).reduce((a, [, o]) => a + o, 0);
    rows.push({ r, dir: '→', k: rateConstant(r.fwd, K.T), order: ord, spec: r.fwd });
    if (r.rev) rows.push({ r, dir: '←', k: rateConstant(r.rev, K.T), order: r.products.reduce((a, [, v]) => a + v, 0), spec: r.rev });
  }
  return rows;
}

function renderPanel() {
  const m = K.mech, res = K.result;
  $('viewport-title').innerHTML = `${esc(K.preset.name)}<small>${nf(K.T, 2)} K · ${fmtTime(K.tEnd)}</small>`;
  $('viewport-legend').innerHTML = '';
  $('viewport-note').innerHTML = res ? `Integrazione implicita (Rosenbrock ROS2, L-stabile): ${res.stats.accepted} passi accettati, ${res.stats.rejected} rifiutati, ${nf(res.ms, 0)} ms. Tolleranza relativa 10⁻⁶.${res.finished ? '' : ' <b>Integrazione interrotta prima della fine: troppi passi.</b>'}` : '';
  if (!m || m.errors.length || !res) { $('panel-body').innerHTML = '<p class="hint">Correggi il meccanismo per avviare la simulazione.</p>'; return; }
  const rows = stepRows();
  const focus = K.focus;
  const fi = res.vars.indexOf(focus);
  const oa = fi >= 0 ? orderAnalysis(res.t, res.c.map(c => c[fi])) : null;
  const last = res.c[res.c.length - 1];
  const cu = m.units === 'molec' ? 'molecole/cm³' : 'mol/L';
  $('panel-body').innerHTML = `
    <div><h3>Costanti cinetiche a ${nf(K.T, 1)} K</h3>
    <table class="data-table"><thead><tr><th>Passo</th><th class="num">k</th><th>unità</th></tr></thead><tbody>
    ${rows.map(x => `<tr title="${esc(specText(x.spec))}"><td class="eqn">${x.dir === '←' ? '← ' : ''}${esc(x.r.text)}</td><td class="num">${sci(x.k)}</td><td>${kUnits(x.order, m.units)}</td></tr>`).join('')}
    </tbody></table></div>
    <div><h3>Specie seguita</h3>
      <select id="kn-focus">${res.vars.map(s => `<option value="${esc(s)}" ${s === focus ? 'selected' : ''}>${esc(s)}</option>`).join('')}</select>
      ${oa ? `<dl class="info-list">
        <dt>Ordine apparente</dt><dd>${oa.best.order}</dd>
        <dt>k osservata (${oa.best.label})</dt><dd>${sci(oa.best.kobs)}</dd>
        <dt>Tempo di dimezzamento</dt><dd>${oa.t12 ? fmtTime(oa.t12) : '—'}</dd>
        ${oa.fits.map(f => `<dt class="desc-muted">R² ordine ${f.order}</dt><dd>${nf(f.fit.r2, 6)}</dd>`).join('')}
      </dl>` : '<p class="hint">La specie non si consuma abbastanza (servono almeno metà del reagente e alcuni punti) per stimare l\'ordine.</p>'}
    </div>
    <div><h3>Alla fine (${fmtTime(res.t[res.t.length - 1])})</h3>
    <table class="data-table"><thead><tr><th>Specie</th><th class="num">iniziale</th><th class="num">finale (${cu})</th></tr></thead><tbody>
    ${m.species.map(s => { const i = res.vars.indexOf(s); return `<tr><td>${pretty(s)}${m.fixed.has(s) ? ' <span class="desc-muted">fissa</span>' : ''}</td><td class="num">${sci(K.c0[s] ?? 0)}</td><td class="num">${sci(i >= 0 ? last[i] : K.c0[s])}</td></tr>`; }).join('')}
    </tbody></table></div>
    ${equilibriumHtml(last)}`;
  $('kn-focus').addEventListener('change', (e) => { K.focus = e.target.value; renderPanel(); drawCharts(); });
}

function specText(s) {
  if (!s) return '';
  if (s.kind === 'k') return 'costante data (indipendente da T)';
  if (s.kind === 't12') return `t½ = ${fmtTime(s.t12)} → k = ln 2 / t½`;
  if (s.kind === 'Tn') return `k = ${s.A} (T/300)^${s.n}`;
  return `k = A${s.n ? ` (T/300)^${s.n}` : ''} exp(−Ea/RT), A = ${s.A}, Ea = ${(s.EaR * R_GAS / 1000).toFixed(2)} kJ/mol`;
}

function equilibriumHtml(last) {
  const m = K.mech, res = K.result;
  const rev = m.reactions.filter(r => r.rev);
  if (!rev.length) return '';
  const val = (s) => { const i = res.vars.indexOf(s); return i >= 0 ? last[i] : K.c0[s] ?? 0; };
  return `<div><h3>Equilibri</h3><dl class="info-list">${rev.map(r => {
    const Kc = rateConstant(r.fwd, K.T) / rateConstant(r.rev, K.T);
    const Q = r.products.reduce((a, [s, v]) => a * Math.pow(val(s), v), 1) / r.reactants.reduce((a, [s, v]) => a * Math.pow(val(s), v), 1);
    return `<dt>${esc(r.text)}: K = k/k<sub>r</sub></dt><dd>${sci(Kc)}</dd><dt class="desc-muted">quoziente Q alla fine</dt><dd>${sci(Q)}</dd>`;
  }).join('')}</dl><p class="hint">All'equilibrio Q = K: il principio del bilancio dettagliato lega le due costanti cinetiche alla costante termodinamica.</p></div>`;
}

function fmtTime(t) {
  const [u, f] = bestUnit(t);
  return `${nf(t / f, t / f < 10 ? 2 : 1)} ${u}`;
}

// ---------------------------------------------------------------------------
// Grafici
// ---------------------------------------------------------------------------

function colorOf(i) { return COLORS[i % COLORS.length]; }

function drawMain() {
  const canvas = $('lab-plot');
  const res = K.result;
  const g = canvas.getContext('2d');
  if (!res) { g.clearRect(0, 0, canvas.width, canvas.height); return; }
  const [un, uf] = bestUnit(K.tEnd);
  const vars = res.vars;
  const series = [];
  let ymax = 0, ymin = Infinity;
  vars.forEach((s, i) => {
    if (K.hidden.has(s)) return;
    const ys = res.c.map(c => c[i]);
    for (const v of ys) { if (v > ymax) ymax = v; if (v > 0 && v < ymin) ymin = v; }
    series.push({ xs: res.t.map(t => t / uf), ys, color: colorOf(i), width: 2, label: s });
  });
  if (!series.length) { g.clearRect(0, 0, canvas.width, canvas.height); return; }
  const tmin = res.t.find(t => t > 0) / uf;
  const lo = K.logY ? Math.max(ymin, ymax * 1e-12) : 0;
  drawXY(canvas, {
    series, xmin: K.logT ? tmin : 0, xmax: K.tEnd / uf, ymin: lo, ymax: (K.logY ? ymax * 3 : ymax * 1.08) || 1,
    xlabel: `t (${un})`, ylabel: K.mech.units === 'molec' ? 'molecole/cm³' : 'c (mol/L)', logX: K.logT, logY: K.logY, legend: false, padLeft: 68,
  });
  $('viewport-legend').innerHTML = vars.map((s, i) => `<button type="button" class="legend-chip" data-sp="${esc(s)}" aria-pressed="${!K.hidden.has(s)}"><i class="swatch" style="background:${colorOf(i)}"></i>${pretty(s)}</button>`).join('');
  $('viewport-legend').querySelectorAll('[data-sp]').forEach(b => b.addEventListener('click', () => {
    const s = b.dataset.sp;
    if (K.hidden.has(s)) K.hidden.delete(s); else K.hidden.add(s);
    drawMain();
  }));
}

function drawCharts() {
  if (!active) return;
  const res = K.result;
  $('radial-title').textContent = 'Leggi cinetiche integrate';
  $('levels-title').textContent = 'Grafico di Arrhenius';
  $('slice-title').textContent = 'Velocità dei singoli passi';
  if (!res) return;
  const fi = res.vars.indexOf(K.focus);
  const oa = fi >= 0 ? orderAnalysis(res.t, res.c.map(c => c[fi])) : null;
  const [un, uf] = bestUnit(K.tEnd);
  if (oa) {
    const b = oa.best;
    const tr = [(c) => c, Math.log, (c) => 1 / c][b.order];
    const ylab = ['[A]', 'ln [A]', '1/[A]'][b.order];
    const xs = [], ys = [];
    for (let i = 0; i < res.t.length; i++) { const c = res.c[i][fi]; if (c < 0.05 * res.c[0][fi]) break; if (c > 0) { xs.push(res.t[i] / uf); ys.push(tr(c)); } }
    const fitY = xs.map(x => b.fit.a + b.fit.b * x * uf);
    const all = [...ys, ...fitY];
    const ymn = Math.min(...all), ymx = Math.max(...all);
    const pad = (ymx - ymn) * 0.08 || 1;
    drawXY($('chart-radial'), {
      series: [{ xs, ys, color: cssVar('--accent'), width: 2.4, label: `${ylab} simulato` }, { xs, ys: fitY, color: cssVar('--muted'), width: 1.4, dash: [5, 4], label: `retta, R² = ${nf(b.fit.r2, 5)}` }],
      xmin: 0, xmax: Math.max(...xs), ymin: ymn - pad, ymax: ymx + pad, xlabel: `t (${un})`, ylabel: ylab,
    });
    $('radial-note').innerHTML = `Per ${pretty(K.focus)}: la legge ${b.label} dà la retta migliore (ordine ${b.order}), fino al 90 % di conversione. R² degli altri ordini: ${oa.fits.filter(f => f !== b).map(f => `${f.order}: ${nf(f.fit.r2, 4)}`).join(', ')}.`;
  } else {
    drawXY($('chart-radial'), { series: [], xmin: 0, xmax: 1, ymin: 0, ymax: 1, xlabel: 't', ylabel: '' });
    $('radial-note').innerHTML = 'Scegli nel pannello una specie che si consuma per analizzare l\'ordine di reazione.';
  }
  // Arrhenius
  const rows = stepRows().filter(x => x.spec.kind === 'arr' && x.spec.EaR > 0);
  if (rows.length) {
    const Ts = [];
    const T0 = Math.max(150, K.T - 120), T1 = K.T + 150;
    for (let T = T0; T <= T1; T += (T1 - T0) / 60) Ts.push(T);
    const series = rows.slice(0, 6).map((x, i) => ({ xs: Ts.map(T => 1000 / T), ys: Ts.map(T => Math.log(rateConstant(x.spec, T))), color: colorOf(i), width: 2, label: `${x.dir === '←' ? '← ' : ''}${x.r.text.slice(0, 28)}` }));
    const ally = series.flatMap(s => s.ys).filter(Number.isFinite);
    drawXY($('chart-levels'), {
      series, xmin: 1000 / T1, xmax: 1000 / T0, ymin: Math.min(...ally) - 1, ymax: Math.max(...ally) + 1,
      xlabel: '1000 / T (K⁻¹)', ylabel: 'ln k', vlines: [{ x: 1000 / K.T, label: `${nf(K.T, 0)} K` }],
    });
    $('levels-note').innerHTML = `ln k = ln A − E<sub>a</sub>/(RT): una retta di pendenza −E<sub>a</sub>/R. ${rows.slice(0, 3).map(x => `${esc(x.r.text)}: E<sub>a</sub> = ${nf(x.spec.EaR * R_GAS / 1000, 1)} kJ/mol, k raddoppia ogni ${nf(doublingDT(x.spec, K.T), 1)} K`).join('; ')}.`;
  } else {
    drawXY($('chart-levels'), { series: [], xmin: 0, xmax: 1, ymin: 0, ymax: 1, xlabel: '1000 / T', ylabel: 'ln k' });
    $('levels-note').innerHTML = 'Nessun passo con parametri di Arrhenius (A ed E<sub>a</sub>): le costanti date come k o t½ non dipendono dalla temperatura.';
  }
  // velocità dei passi
  const { steps } = buildModel(K.mech, K.T);
  const fixedVal = new Map([...K.mech.fixed].map(s => [s, K.c0[s] ?? 0]));
  const vi = new Map(res.vars.map((s, i) => [s, i]));
  const rSeries = steps.slice(0, 8).map((st, j) => ({
    xs: res.t.map(t => t / uf),
    ys: res.c.map(c => st.ord.reduce((a, [s, o]) => a * Math.pow(Math.max(vi.has(s) ? c[vi.get(s)] : fixedVal.get(s) ?? 0, 0), o), st.k)),
    color: colorOf(j), width: 1.8, label: `passo ${j + 1}`,
  }));
  const rv = rSeries.flatMap(s => s.ys).filter(v => v > 0);
  if (rv.length) {
    const rmax = Math.max(...rv), rmin = Math.max(Math.min(...rv), rmax * 1e-14);
    drawXY($('chart-slice'), {
      series: rSeries, xmin: K.logT ? res.t.find(t => t > 0) / uf : 0, xmax: K.tEnd / uf, ymin: rmin, ymax: rmax * 3, logY: true, logX: K.logT,
      xlabel: `t (${un})`, ylabel: K.mech.units === 'molec' ? 'r (molecole cm⁻³ s⁻¹)' : 'r (M s⁻¹)',
    });
    $('slice-note').innerHTML = 'r = k Π cᵢ^oᵢ per ogni passo (le reazioni reversibili contano due passi). Il passo più lento di una sequenza determina la velocità complessiva; in uno stato stazionario le velocità di produzione e consumo di un intermedio si uguagliano.';
  }
}

function doublingDT(spec, T) {
  // T' tale che k(T') = 2 k(T): 1/T' = 1/T − ln 2 / (Ea/R)
  const Tp = 1 / (1 / T - Math.LN2 / spec.EaR);
  return Tp - T;
}

// ---------------------------------------------------------------------------
// Analisi
// ---------------------------------------------------------------------------

function renderAnalysis() {
  const m = K.mech;
  const ey = m ? m.reactions.map(r => {
    const uni = r.reactants.filter(([s]) => !m.fixed.has(s) && (r.orders[s] ?? 1) !== 0).length === 1 && r.reactants.reduce((a, [s, v]) => a + (r.orders[s] ?? v), 0) - r.reactants.filter(([s]) => m.fixed.has(s)).reduce((a, [s, v]) => a + (r.orders[s] ?? v), 0) === 1;
    return uni ? { r, e: eyring(r.fwd, K.T) } : null;
  }).filter(x => x?.e) : [];
  const bal = m ? m.reactions.map(r => ({ r, b: balance(r) })) : [];
  $('analysis').innerHTML = `
    <h2>Che cosa calcola la cinetica</h2>
    <div class="analysis-grid">
      <div><h3>Bilancio dei passi</h3>
      <table class="data-table"><thead><tr><th>Passo</th><th>Elementi e carica</th></tr></thead><tbody>
      ${bal.map(({ r, b }) => `<tr><td class="eqn">${esc(r.text)}</td><td>${b === null ? '<span class="desc-muted">non verificabile (nomi non chimici o isotopi)</span>' : b.ok ? 'bilanciato' : `<span class="warn">non bilanciato: ${b.bad.map(x => (x === 'q' ? 'carica' : x)).join(', ')}</span>`}</td></tr>`).join('')}
      </tbody></table></div>
      ${ey.length ? `<div><h3>Stato di transizione (Eyring)</h3>
      <p class="desc">Per i passi unimolecolari: k = (k<sub>B</sub>T/h) e<sup>ΔS‡/R</sup> e<sup>−ΔH‡/RT</sup>, con ΔH‡ = E<sub>a</sub> − RT.</p>
      <dl class="info-list">${ey.map(({ r, e }) => `<dt>${esc(r.text)}</dt><dd>ΔH‡ = ${nf(e.dH / 1000, 1)} kJ/mol · ΔS‡ = ${nf(e.dS, 1)} J/(K mol) · ΔG‡ = ${nf(e.dG / 1000, 1)} kJ/mol</dd>`).join('')}</dl>
      <p class="hint">ΔS‡ > 0: stato di transizione più "sciolto" dei reagenti (rottura di un legame); ΔS‡ < 0: più ordinato.</p></div>` : ''}
      <article><h3>Come viene calcolato</h3>
        <p>Ogni riga è un passo con velocità data dalla legge di azione di massa; le concentrazioni evolvono secondo il sistema di equazioni differenziali</p>
        <p class="eq">r<sub>j</sub> = k<sub>j</sub>(T) Π<sub>i</sub> c<sub>i</sub><sup>o<sub>ij</sub></sup>,&nbsp;&nbsp; dc<sub>i</sub>/dt = Σ<sub>j</sub> ν<sub>ij</sub> r<sub>j</sub>,&nbsp;&nbsp; k(T) = A (T/300)<sup>n</sup> e<sup>−E<sub>a</sub>/RT</sup></p>
        <p>I sistemi reali sono rigidi: specie che vivono millisecondi (O atomico, HBrO₂) accanto ad altre che cambiano in mesi. Un metodo esplicito richiederebbe passi minuscoli; qui si usa il metodo di Rosenbrock ROS2 (Verwer et al., SIAM J. Sci. Comput. 1999), implicito lineare e L-stabile, con jacobiana analitica e passo adattivo guidato dalla differenza con la soluzione del primo ordine incorporata.</p>
      </article>
      <article><h3>Fonti e limiti</h3>
        <p class="desc-muted">Parametri di Arrhenius: Atkins, <i>Physical Chemistry</i>, tabelle 22.1 e 22.4 (dati di Laidler, Pilling e Seakins, Nicholas, Frost e Pearson). Chimica dell'atmosfera: IUPAC, Atkinson et al., Atmos. Chem. Phys. 4, 1461 (2004); atmosfera standard USA 1976. Tempi di dimezzamento: NNDC/ENSDF. Oregonator: Field e Noyes (1974), costanti di Field e Försterling (1986). Costanti fisiche: CODATA 2018. Ogni simulazione è isoterma (il calore di reazione non cambia T) e a volume costante; le costanti di fotolisi (J) sono ordini di grandezza.</p>
      </article>
    </div>`;
}
