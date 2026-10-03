// Modalità "Fasi · MD": 256 atomi Lennard–Jones in unità ridotte, con lettura della fase dalla struttura locale.

import { PhaseMD, ARGON } from '../chem/phaseMD.js';
import { localOrder, phaseVerdict, Mobility } from '../chem/structureOrder.js';
import { drawAtoms, bindRotation } from './atomicCanvas.js';
import { drawXY } from './chemCharts.js';

const $ = (id) => document.getElementById(id);
const fmt = (v, n = 3) => Number.isFinite(v) ? Number(v).toLocaleString('it-IT', { maximumFractionDigits: n }) : '—';

// primo minimo di g(r) del fluido LJ vicino al punto triplo: delimita la prima sfera di coordinazione
const CUTOFF = 1.5;
const COLORS = { crystal: '#8b93ff', liquid: '#f0b429', vapor: '#9fd7ff' };
const PHASE_COLORS = { gas: '#7ec8ff', liquido: '#f0b429', solido: '#8b93ff', misto: '#e07a9a' };

const PRESETS = {
  solid: { name: 'Cristallo FCC', density: 0.95, temperature: 0.35 },
  melt: { name: 'Fusione: FCC riscaldato', density: 0.85, temperature: 1.5 },
  gas: { name: 'Gas diluito', density: 0.04, temperature: 2 },
  cool: { name: 'Raffreddamento del liquido', density: 0.95, temperature: 1.6 },
  coexist: { name: 'Liquido e vapore', density: 0.3, temperature: 0.75 },
};

let active = false, raf = 0, model = null, running = false, unbind = null, frames = 0, history = [], error = '';
let preset = 'solid', order = null, mobile = null, mobility = new Mobility();
const view = { yaw: 0.55, pitch: 0.3 };

export function activatePhase() {
  active = true;
  document.body.dataset.mode = 'phase';
  $('busy').hidden = true;
  if (!model) reset(); else render();
  unbind = bindRotation($('lab-plot'), view, draw);
  loop();
}
export function deactivatePhase() { active = false; running = false; cancelAnimationFrame(raf); unbind?.(); unbind = null; }
export function phaseRedraw() { if (active) { draw(); charts(); } }

function analyse() {
  order = localOrder(model.x, { box: model.L, cutoff: CUTOFF });
  mobile = mobility.sample(model.x, model.time, 1, 1, model.L);   // unità ridotte: τ = 1, σ = 1
  return order;
}

function reset() {
  model = new PhaseMD(PRESETS[preset]);
  running = false; error = ''; history = []; mobility = new Mobility(); mobile = null;
  analyse(); record(); render(); draw(); charts();
}

function record() {
  history.push({ ...model.stats(), crystal: order.fraction.crystal, Q6: order.Q6 });
  if (history.length > 800) history.shift();
}

function loop() {
  if (!active) return;
  if (running) {
    try {
      model.step(8);
      frames++;
      if (frames % 3 === 0) analyse();
      if (frames % 6 === 0) { record(); metrics(); }
      if (frames % 30 === 0) charts();
    } catch (e) { running = false; error = e.message; render(); }
  }
  draw();
  raf = requestAnimationFrame(loop);
}

function draw() {
  if (!model) return;
  if (!order) analyse();
  const atoms = Array.from({ length: model.n }, (_, i) => {
    const k = order.kind[i];
    return {
      element: 'Ar', radius: k === 'crystal' ? 0.3 : 0.26, color: COLORS[k],
      position: Array.from(model.x.slice(3 * i, 3 * i + 3), v => ((v % model.L) + model.L) % model.L),
    };
  });
  drawAtoms($('lab-plot'), atoms, { ...view, span: model.L * 1.15, unit: 'σ' });
  paintBanner(phaseVerdict(order.fraction, mobile));
}

function paintBanner(phase) {
  const canvas = $('lab-plot'), g = canvas.getContext('2d');
  const dpr = Math.min(devicePixelRatio || 1, 2);
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.font = '600 22px sans-serif';
  const w = g.measureText(phase.name).width + 28;
  g.fillStyle = 'rgba(10,16,24,.72)';
  g.beginPath(); g.roundRect(14, 52, w, 36, 10); g.fill();
  g.fillStyle = PHASE_COLORS[phase.id] ?? '#fff';
  g.fillText(phase.name, 28, 76);
  g.font = '12px sans-serif';
  const why = phase.why, ww = g.measureText(why).width + 12;
  g.fillStyle = 'rgba(10,16,24,.6)';
  g.fillRect(14, 92, ww, 18);
  g.fillStyle = '#e3e8ef';
  g.fillText(why, 20, 105);
}

function render() {
  $('element-card').innerHTML = `
    <h3 class="side-h">Sandbox · stati della materia</h3>
    <p class="hint">256 atomi · forze Lennard–Jones · condizioni periodiche</p>
    <label class="lbl" for="ph-preset">Preparazione iniziale</label>
    <select id="ph-preset">${Object.entries(PRESETS).map(([id, p]) => `<option value="${id}" ${preset === id ? 'selected' : ''}>${p.name}</option>`).join('')}</select>
    <p class="hint">${{
      cool: 'Avvia a T* = 1,6 finché il reticolo si scioglie, poi raffredda a 0,35. La nucleazione può richiedere tempo o non avvenire in questa piccola cella: il riconoscimento q₆ la segnala con qualunque orientazione del nuovo cristallo.',
      coexist: 'A densità intermedia e sotto il punto critico il fluido si separa in una goccia (o una lastra) di liquido e nel suo vapore.',
    }[preset] ?? 'La fase emerge dalle forze. I nomi indicano la preparazione, non impongono lo stato finale.'}</p>
    <div class="btn-row"><a class="btn" href="#sandbox">Reazioni molecolari</a><a class="btn" href="#becher">Soluzioni e precipitati</a></div>
    <h3 class="side-h">Esperimenti</h3>
    <p class="hint">Riscalda il cristallo e osserva sparire i legami solidi e i picchi di g(r). Raffredda il liquido e cerca la ricristallizzazione. Riduci la densità per esplorare gas e coesistenza. Spegni il termostato per verificare la conservazione dell'energia.</p>`;
  $('ph-preset').onchange = (e) => { preset = e.target.value; reset(); };

  $('controls').innerHTML = `
    <div class="btn-row"><button class="btn" id="ph-play">${running ? 'Pausa' : 'Avvia'}</button><button class="btn" id="ph-step">10 passi</button><button class="btn" id="ph-reset">Ripristina</button></div>
    <label class="lbl" for="ph-T">Temperatura del bagno T* = k<sub>B</sub>T/ε</label>
    <input id="ph-T" type="number" min="0.05" max="3" step="0.05" value="${model.target}">
    <div class="btn-row"><button class="btn" id="ph-cool">Raffredda · 0,35</button><button class="btn" id="ph-heat">Riscalda · 1,6</button></div>
    <label class="lbl"><input type="checkbox" id="ph-bath" ${model.thermostat ? 'checked' : ''}> Termostato Langevin (NVT)</label>
    <label class="lbl" for="ph-rho">Densità ρ* = Nσ³/V</label>
    <input type="number" id="ph-rho" min="0.02" max="1.1" step="0.01" value="${model.stats().density.toFixed(3)}">
    <button class="btn" id="ph-volume">Applica volume</button>
    <p class="hint">Δt* = 0,002 · seed 2026. Il lavoro di compressione è contabilizzato; il riferimento MSD riparte al cambio di volume.</p>
    <button class="btn" id="ph-csv">Esporta misure CSV</button>
    <p class="hint warn" id="ph-error">${error}</p>`;
  $('ph-play').onclick = () => { running = !running; $('ph-play').textContent = running ? 'Pausa' : 'Avvia'; };
  $('ph-step').onclick = () => {
    running = false;
    try { model.step(10); analyse(); record(); error = ''; } catch (e) { error = e.message; }
    render(); charts();
  };
  $('ph-reset').onclick = reset;
  $('ph-T').onchange = (e) => { const t = +e.target.value; if (Number.isFinite(t) && t >= 0.05 && t <= 3) model.target = t; else e.target.value = model.target; };
  $('ph-bath').onchange = (e) => { model.thermostat = e.target.checked; };
  $('ph-cool').onclick = () => { model.target = 0.35; model.thermostat = true; render(); };
  $('ph-heat').onclick = () => { model.target = 1.6; model.thermostat = true; render(); };
  $('ph-volume').onclick = () => {
    try { model.setDensity(+$('ph-rho').value); mobility = new Mobility(); analyse(); history = []; record(); error = ''; } catch (e) { error = e.message; }
    render(); charts();
  };
  $('ph-csv').onclick = exportCsv;

  $('viewport-title').innerHTML = 'Dinamica delle fasi<small>Trascina per ruotare · il colore è l\'intorno locale (q₆)</small>';
  $('viewport-legend').innerHTML = Object.entries({ crystal: 'cristallino', liquid: 'condensato disordinato', vapor: 'isolato (vapore)' })
    .map(([k, t]) => `<span><i class="swatch" style="background:${COLORS[k]}"></i>${t}</span>`).join('');
  $('viewport-note').textContent = 'LJ con forza traslata, rc = 2,5σ. Gas nobile modello; nessun legame chimico. Le particelle attraversano le facce periodiche della scatola.';
  $('analysis').innerHTML = `<h2>Modello e misure</h2><div class="analysis-grid">
    <article><h3>Forze e traiettorie</h3><p class="desc">U = 4ε[(σ/r)¹² − (σ/r)⁶]. Energia e forza sono portate a zero a 2,5σ con una correzione lineare. In NVE si usa velocity Verlet; in NVT Langevin BAOAB. Il calore del bagno e il lavoro di variazione del volume sono misurati separatamente.</p></article>
    <article><h3>Come si riconosce la fase</h3><p class="desc">Per ogni atomo si calcolano i parametri di Steinhardt q₆ₘ sui vicini entro 1,5σ (primo minimo di g(r)). Due vicini con orientazioni locali simili (prodotto scalare normalizzato > 0,7) formano un "legame solido"; con almeno 7 legami solidi l'atomo è cristallino (ten Wolde, Ruiz-Montero, Frenkel, 1996). Con meno di 3 vicini è isolato, cioè vapore. Per separare un liquido da un solido disordinato (vetro) si conta quanti primi vicini si perdono in 3τ: oltre il 20% le gabbie si rinnovano e il fluido scorre. Il criterio non dipende dall'orientazione del reticolo, quindi riconosce anche un cristallo nato dal liquido. Q₆ globale vale 0,575 per l'FCC perfetto e tende a 0 nel liquido. Sopra il punto critico liquido e vapore non si distinguono: la lettura "liquido + vapore" indica solo un fluido non uniforme.</p></article>
    <article><h3>Unità e dominio</h3><p class="desc">T*, ρ*, P*, E* e t* sono unità ridotte. Mappatura approssimata all'argon: σ = 3,405 Å, ε/k<sub>B</sub> = 119,8 K, τ ≈ 2,156 ps. Il troncamento e il numero finito di atomi spostano le transizioni: non sono temperature di fusione sperimentali. Nessuna estrapolazione a rame, acqua o sali.</p><a href="https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-properties" target="_blank" rel="noopener">Riferimenti NIST per il fluido LJ</a></article>
  </div>`;
  metrics();
}

function metrics() {
  const s = model.stats(), phase = phaseVerdict(order.fraction, mobile), f = order.fraction;
  const rows = [
    ['Tempo t*', s.time], ['Temperatura T*', s.T], ['Pressione P*', s.pressure], ['Densità ρ*', s.density],
    ['Energia per atomo E*/N', s.E / model.n], ['Calore dal bagno Q*/N', s.heat / model.n], ['Lavoro sul sistema W*/N', s.work / model.n],
    ['Errore (E−E₀−Q−W)/N', s.drift / model.n], ['Atomi cristallini (q₆)', f.crystal], ['Atomi isolati', f.vapor], ['Q₆ globale', order.Q6], ['Vicini persi in 3τ', mobility.lostFraction ?? NaN],
    ['Ordine FCC (200) iniziale', s.order], ['MSD / σ²', s.msd],
  ];
  $('panel-body').innerHTML = `<h3>Fase letta dalla struttura: ${phase.name}</h3><p class="hint">${phase.why}.</p>
    <h3>Misure istantanee</h3><dl class="info-list">${rows.map(([label, v]) => `<dt>${label}</dt><dd>${fmt(v, 4)}</dd>`).join('')}</dl>
    <p class="hint">Equivalente Ar approssimato: ${fmt(s.T * ARGON.epsilonK, 1)} K · ${fmt(s.time * ARGON.tauPs, 2)} ps. Le misure istantanee fluttuano.</p>`;
}

function chart(id, series, xlabel, ylabel, range = null) {
  const vals = series.flatMap(s => s.ys), xs = series[0].xs;
  const lo = range?.[0] ?? Math.min(...vals), hi = range?.[1] ?? Math.max(...vals), pad = range ? 0 : Math.max(0.001, (hi - lo) * 0.1);
  drawXY($(id), { series, xmin: Math.min(...xs), xmax: Math.max(0.01, ...xs), ymin: lo - pad, ymax: hi + pad, xlabel, ylabel });
}

function charts() {
  if (!model || !active) return;
  const x = history.map(s => s.time);
  $('radial-title').textContent = 'Conservazione dell\'energia';
  chart('chart-radial', [{ xs: x, ys: history.map(s => s.drift / model.n), color: '#ce8850', label: '(E−E₀−Q−W)/N' }], 't*', 'errore / ε');
  $('radial-note').textContent = 'Deve restare piccolo rispetto all\'energia per atomo. Non viene azzerato quando cambia il bagno.';
  $('levels-title').textContent = 'Ordine cristallino';
  chart('chart-levels', [
    { xs: x, ys: history.map(s => s.crystal), color: '#8b93ff', label: 'atomi cristallini (q₆)' },
    { xs: x, ys: history.map(s => s.order), color: '#67b7b0', label: 'Bragg (200) del reticolo iniziale' },
  ], 't*', 'frazione', [0, 1]);
  $('levels-note').textContent = 'La frazione q₆ riconosce qualunque cristallo; la coerenza di Bragg solo quello iniziale, con la sua orientazione.';
  const r = model.rdf();
  $('slice-title').textContent = 'Distribuzione radiale g(r)';
  chart('chart-slice', [{ xs: r.r, ys: r.g, color: '#67b7b0', label: 'g(r) istantanea' }], 'r / σ', 'g(r)');
  $('slice-note').textContent = 'Distribuzione delle distanze a coppie, normalizzata con il volume dei gusci e N(N−1).';
}

function exportCsv() {
  const keys = ['time', 'T', 'density', 'pressure', 'K', 'U', 'E', 'heat', 'work', 'drift', 'crystal', 'Q6', 'order', 'msd'];
  const text = `# LJ force-shift rc=2.5; reduced units; dt=0.002; N=256; seed=2026; crystal = q6 solid-bond fraction (d>0.7, >=7 bonds, rcut 1.5)\n${keys.join(',')}\n${history.map(s => keys.map(k => s[k]).join(',')).join('\n')}`;
  const url = URL.createObjectURL(new Blob([text], { type: 'text/csv' }));
  const a = document.createElement('a');
  a.href = url; a.download = 'fasi-md.csv'; a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
