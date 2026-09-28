// Modalità "Becher": chimica in soluzione. Si versano reagenti da uno scaffale e un solutore di equilibrio
// calcola pH, precipitati, complessi, reazioni redox, gas, calore e colore della soluzione.

import { Beaker, absorbance } from '../chem/aqueous.js';
import { REAGENTS, REAGENT_GROUPS, BEAKER_PRESETS, METALS } from '../chem/aqueousData.js';
import { drawXY } from './chemCharts.js';

const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const nf = (v, d = 2) => (v === null || v === undefined || !Number.isFinite(v) ? '—' : Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d, useGrouping: Math.abs(v) >= 1e4 }));
const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const sci = (v, d = 2) => {
  if (!Number.isFinite(v)) return '—';
  if (v === 0) return '0';
  const e = Math.floor(Math.log10(Math.abs(v)));
  if (e >= -2 && e <= 3) return nf(v, d + Math.max(0, -e));
  return `${nf(v / Math.pow(10, e), d)}·10${String(e).split('').map(c => SUP[c]).join('')}`;
};
const CAPACITY = 250; // mL

let active = false;
let raf = 0;
const B = {
  beaker: new Beaker(),
  preset: BEAKER_PRESETS[0],
  reagent: 'NaOH',
  amount: 5,
  drip: false,
  dripTimer: 0,
  particles: [],
  bubbles: [],
  anim: { level: 0, color: [235, 240, 245] },
  message: '',
};

export function activateBeaker() {
  active = true;
  document.body.dataset.mode = 'beaker';
  $('busy').hidden = true;
  if (!B.beaker.log.length && !B.started) { B.started = true; loadPreset(B.preset); }
  renderAll();
  loop();
}

export function deactivateBeaker() {
  active = false;
  cancelAnimationFrame(raf);
  stopDrip();
}

export function beakerRedraw() {
  if (active) { drawCharts(); }
}

function loadPreset(p) {
  B.preset = p;
  B.beaker = new Beaker();
  B.particles = [];
  B.bubbles = [];
  B.message = '';
  for (const [r, v] of p.steps) addReagent(r, v, false);
  const suggestion = { neutral: 'NaOH', silver: 'NH3', copper: 'NH3', iron: 'NaOH', lead: 'KI', zinc: 'HCl', 'cu-ag': 'AgNO3', buffer: 'HCl', marble: 'HCl', amphoteric: 'NaOH', chromate: 'HCl' }[p.id];
  if (suggestion) B.reagent = suggestion;
  B.amount = { neutral: 1, buffer: 1, amphoteric: 2, copper: 1, chromate: 2 }[p.id] ?? 5;
}

function reagentUnit(r) { return r.kind === 'metal' || r.kind === 'solid' ? 'g' : 'mL'; }

function addReagent(id, amount, render = true) {
  const r = REAGENTS.find(x => x.id === id);
  if (r.kind !== 'metal' && r.kind !== 'solid' && B.beaker.volumeML + amount > CAPACITY) {
    B.message = `Il becher (${CAPACITY} mL) traboccherebbe: svuotalo o aggiungi meno.`;
    stopDrip();
    if (render) renderAll();
    return;
  }
  const before = B.beaker.summary();
  const entry = B.beaker.add(id, amount);
  if (!entry) return;
  B.message = '';
  const after = B.beaker.summary();
  // particelle di precipitato: nuove particelle per ogni solido che è aumentato
  for (const s of after.solids) {
    const old = before.solids.find(x => x.id === s.id)?.mass ?? 0;
    const add = Math.min(220, Math.round((s.mass - old) * 900));
    for (let k = 0; k < add; k++) B.particles.push({ id: s.id, color: s.color, x: Math.random(), y: Math.random() * 0.8, vy: 0.002 + Math.random() * 0.004, settled: false, r: 1.2 + Math.random() * 1.8 });
  }
  for (const s of before.solids) {
    const now = after.solids.find(x => x.id === s.id)?.mass ?? 0;
    if (now < s.mass) {
      const frac = s.mass > 0 ? now / s.mass : 0;
      B.particles = B.particles.filter(p => p.id !== s.id || Math.random() < frac);
    }
  }
  if (B.particles.length > 1500) B.particles.splice(0, B.particles.length - 1500);
  const gas = (entry.h2 ?? 0) + (entry.co2 ?? 0);
  if (gas > 0) {
    const n = Math.min(160, 20 + Math.round(gas * 4000));
    for (let k = 0; k < n; k++) B.bubbles.push({ x: 0.2 + 0.6 * Math.random(), y: 1 + Math.random() * 1.5, r: 1.5 + Math.random() * 3, v: 0.004 + Math.random() * 0.006, gas: entry.h2 > 0 ? 'H₂' : 'CO₂' });
  }
  if (render) renderAll();
}

function stopDrip() {
  B.drip = false;
  clearInterval(B.dripTimer);
}

// ---------------------------------------------------------------------------
// Disegno del becher
// ---------------------------------------------------------------------------

function loop() {
  if (!active) return;
  drawBeaker();
  raf = requestAnimationFrame(loop);
}

function drawBeaker() {
  const canvas = $('lab-plot');
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  if (!w || !h) return;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
  const g = canvas.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  const s = B.beaker.summary();
  const ink = cssVar('--text'), muted = cssVar('--muted'), line = cssVar('--line');
  // geometria del becher
  const bh = Math.min(h * 0.72, w * 0.62), bw = bh * 0.72;
  const x0 = w * 0.40 - bw / 2, y1 = h * 0.9, y0 = y1 - bh;
  // animazione del livello e del colore
  const target = Math.min(1, s.V / CAPACITY);
  B.anim.level += (target - B.anim.level) * 0.12;
  B.anim.color = B.anim.color.map((c, k) => c + (s.color[k] - c) * 0.08);
  const liquidTop = y1 - B.anim.level * bh * 0.94;
  // banco
  g.fillStyle = cssVar('--panel-2');
  g.fillRect(0, y1 + 4, w, h - y1 - 4);
  g.strokeStyle = line; g.beginPath(); g.moveTo(0, y1 + 4); g.lineTo(w, y1 + 4); g.stroke();
  // liquido
  if (s.V > 0) {
    const [r, gg, b] = B.anim.color.map(Math.round);
    const grad = g.createLinearGradient(x0, 0, x0 + bw, 0);
    grad.addColorStop(0, `rgba(${r},${gg},${b},0.78)`);
    grad.addColorStop(0.5, `rgba(${r},${gg},${b},0.62)`);
    grad.addColorStop(1, `rgba(${r},${gg},${b},0.82)`);
    g.fillStyle = grad;
    g.beginPath();
    g.moveTo(x0 + 3, liquidTop);
    g.quadraticCurveTo(x0 + bw / 2, liquidTop + 4, x0 + bw - 3, liquidTop);
    g.lineTo(x0 + bw - 3, y1 - 8);
    g.quadraticCurveTo(x0 + bw - 3, y1 - 2, x0 + bw - 10, y1 - 2);
    g.lineTo(x0 + 10, y1 - 2);
    g.quadraticCurveTo(x0 + 3, y1 - 2, x0 + 3, y1 - 8);
    g.closePath();
    g.fill();
    // menisco
    g.strokeStyle = `rgba(${r},${gg},${b},0.95)`;
    g.lineWidth = 1.5;
    g.beginPath(); g.moveTo(x0 + 3, liquidTop - 2); g.quadraticCurveTo(x0 + bw / 2, liquidTop + 4, x0 + bw - 3, liquidTop - 2); g.stroke();
  }
  // precipitato: le particelle cadono e si depositano sul fondo
  const floorY = y1 - 5;
  const heapCount = new Map();
  for (const p of B.particles) {
    if (!p.settled) {
      p.y += p.vy;
      if (p.y >= 1) { p.y = 1; p.settled = true; }
    }
    const px = x0 + 8 + p.x * (bw - 16);
    let py = liquidTop + p.y * (floorY - liquidTop);
    if (p.settled) {
      const bucket = Math.floor(p.x * 40);
      const k = (heapCount.get(bucket) ?? 0) + 1;
      heapCount.set(bucket, k);
      py = floorY - k * 1.1 - p.r * 0.5;
    }
    g.fillStyle = p.color;
    g.globalAlpha = p.settled ? 0.95 : 0.8;
    g.beginPath(); g.arc(px, py, p.r, 0, 2 * Math.PI); g.fill();
  }
  g.globalAlpha = 1;
  // metalli sul fondo
  let mx = x0 + 18;
  for (const m of s.metals) {
    const n = Math.max(1, Math.min(14, Math.round(m.mass * 8)));
    for (let k = 0; k < n; k++) {
      g.fillStyle = m.color;
      g.strokeStyle = 'rgba(0,0,0,0.35)';
      const sz = 5 + (k % 3) * 2;
      g.beginPath();
      g.ellipse(mx + (k % 7) * 9, floorY - 4 - Math.floor(k / 7) * 7, sz, sz * 0.6, (k * 0.7) % 3, 0, 2 * Math.PI);
      g.fill(); g.stroke();
    }
    mx += Math.min(7, n) * 9 + 30;
  }
  // bolle di gas
  for (const bb of B.bubbles) {
    bb.y -= bb.v;
    const py = liquidTop + bb.y * (floorY - liquidTop);
    if (py < liquidTop) continue;
    g.strokeStyle = 'rgba(255,255,255,0.85)';
    g.lineWidth = 1;
    g.beginPath(); g.arc(x0 + 10 + bb.x * (bw - 20), py, bb.r, 0, 2 * Math.PI); g.stroke();
  }
  B.bubbles = B.bubbles.filter(bb => liquidTop + bb.y * (floorY - liquidTop) > liquidTop - 2);
  // vetro
  g.strokeStyle = cssVar('--scene-mode') === 'light' ? '#8b97a8' : '#c9d3e0';
  g.lineWidth = 2.2;
  g.beginPath();
  g.moveTo(x0 - 10, y0 - 6); g.quadraticCurveTo(x0 - 2, y0 - 4, x0, y0 + 6);
  g.lineTo(x0, y1 - 8); g.quadraticCurveTo(x0, y1, x0 + 8, y1);
  g.lineTo(x0 + bw - 8, y1); g.quadraticCurveTo(x0 + bw, y1, x0 + bw, y1 - 8);
  g.lineTo(x0 + bw, y0);
  g.stroke();
  // tacche graduate
  g.fillStyle = muted; g.strokeStyle = muted; g.lineWidth = 1;
  g.font = '11px "IBM Plex Mono", monospace'; g.textAlign = 'left';
  for (let v = 50; v <= 200; v += 50) {
    const y = y1 - (v / CAPACITY) * bh * 0.94;
    g.beginPath(); g.moveTo(x0 + bw - 22, y); g.lineTo(x0 + bw, y); g.stroke();
    g.fillText(`${v} mL`, x0 + bw - 60, y - 3);
  }
  // pHmetro
  const px0 = x0 + bw * 0.28;
  g.strokeStyle = muted; g.lineWidth = 3;
  g.beginPath(); g.moveTo(px0, y0 - 40); g.lineTo(px0, Math.max(liquidTop + 30, y0 + 20)); g.stroke();
  const mw = 132, mh = 58, mX = x0 + bw + 36, mY = y0 + 10;
  g.fillStyle = cssVar('--panel-2'); g.strokeStyle = line; g.lineWidth = 1;
  g.fillRect(mX, mY, mw, mh); g.strokeRect(mX, mY, mw, mh);
  g.fillStyle = muted; g.font = '11px "IBM Plex Sans", sans-serif'; g.fillText('pHmetro', mX + 8, mY + 15);
  g.fillStyle = ink; g.font = '600 24px "IBM Plex Mono", monospace';
  g.fillText(s.V > 0 ? nf(s.pH, 2) : '—', mX + 8, mY + 44);
  // termometro
  const tX = mX, tY = mY + mh + 16;
  g.fillStyle = cssVar('--panel-2'); g.fillRect(tX, tY, mw, mh); g.strokeRect(tX, tY, mw, mh);
  g.fillStyle = muted; g.font = '11px "IBM Plex Sans", sans-serif'; g.fillText('termometro', tX + 8, tY + 15);
  g.fillStyle = ink; g.font = '600 24px "IBM Plex Mono", monospace';
  g.fillText(`${nf(s.T - 273.15, 1)} °C`, tX + 8, tY + 44);
  // volume
  const vY = tY + mh + 16;
  g.fillStyle = cssVar('--panel-2'); g.fillRect(tX, vY, mw, mh); g.strokeRect(tX, vY, mw, mh);
  g.fillStyle = muted; g.font = '11px "IBM Plex Sans", sans-serif'; g.fillText('volume', tX + 8, vY + 15);
  g.fillStyle = ink; g.font = '600 24px "IBM Plex Mono", monospace';
  g.fillText(`${nf(s.V, 1)} mL`, tX + 8, vY + 44);
  // buretta durante la titolazione
  if (B.drip) {
    const bx = x0 + bw * 0.62;
    g.strokeStyle = muted; g.lineWidth = 2;
    g.strokeRect(bx - 5, y0 - 150, 10, 110);
    g.beginPath(); g.moveTo(bx, y0 - 40); g.lineTo(bx, y0 - 20); g.stroke();
    const t = (performance.now() / 380) % 1;
    g.fillStyle = cssVar('--phase-neg');
    g.beginPath(); g.arc(bx, y0 - 18 + t * (liquidTop - y0 + 18), 2.5, 0, 2 * Math.PI); g.fill();
  }
  // etichette dei precipitati
  g.textAlign = 'left';
  g.font = '12px "IBM Plex Sans", sans-serif';
  let ly = y1 + 22;
  for (const sd of s.solids) {
    g.fillStyle = sd.color;
    g.fillRect(x0, ly - 9, 10, 10);
    g.strokeStyle = line; g.strokeRect(x0, ly - 9, 10, 10);
    g.fillStyle = ink;
    g.fillText(`${sd.label}(s)↓ ${nf(sd.mass * 1000, 1)} mg`, x0 + 16, ly);
    ly += 16;
  }
  for (const m of s.metals) {
    g.fillStyle = m.color;
    g.fillRect(x0, ly - 9, 10, 10);
    g.strokeStyle = line; g.strokeRect(x0, ly - 9, 10, 10);
    g.fillStyle = ink;
    g.fillText(`${m.label}(s) ${nf(m.mass * 1000, 1)} mg`, x0 + 16, ly);
    ly += 16;
  }
  if (s.gasH2 > 0 || s.gasCO2 > 0) {
    g.fillStyle = muted;
    const Vm = 24.47; // L/mol a 25 °C e 1 atm
    g.fillText(`gas liberato: ${s.gasH2 > 0 ? `H₂ ${nf(s.gasH2 * Vm * 1000, 1)} mL ` : ''}${s.gasCO2 > 0 ? `CO₂ ${nf(s.gasCO2 * Vm * 1000, 1)} mL` : ''}`, x0 + bw + 36, vY + mh + 24);
  }
}

// ---------------------------------------------------------------------------
// Pannelli
// ---------------------------------------------------------------------------

function renderAll() {
  if (!active) return;
  renderSide();
  renderControls();
  renderPanelBody();
  drawCharts();
  renderAnalysis();
}

function renderSide() {
  const r = REAGENTS.find(x => x.id === B.reagent);
  const groups = REAGENT_GROUPS.map(([g, label]) => {
    const items = REAGENTS.filter(x => x.group === g);
    return `<label class="lbl">${label}</label><div class="mol-list">${items.map(x => `<button type="button" class="mol-chip ${x.id === B.reagent ? 'active' : ''}" data-r="${x.id}" title="${x.name}${x.c ? ` ${nf(x.c, x.c < 0.01 ? 4 : 2)} M` : ''}">${x.formula}</button>`).join('')}</div>`;
  }).join('');
  $('element-card').innerHTML = `
    <h3 class="side-h">Banco di chimica</h3>
    <label class="lbl" for="bk-preset">Esperienza</label>
    <select id="bk-preset">${BEAKER_PRESETS.map(p => `<option value="${p.id}" ${p.id === B.preset.id ? 'selected' : ''}>${p.name}</option>`).join('')}</select>
    <p class="mol-note">${B.preset.text}</p>
    <div class="btn-row"><button type="button" class="btn" id="bk-restart">Ricomincia</button><button type="button" class="btn" id="bk-empty">Svuota il becher</button></div>
    <h3 class="side-h" style="margin-top:14px">Scaffale dei reagenti</h3>
    ${groups}
    <p class="hint">Selezionato: <b>${r.name}</b>${r.c ? `, ${nf(r.c, r.c < 0.01 ? 4 : 2)} mol/L` : ''}.</p>`;
  $('bk-preset').addEventListener('change', (e) => { stopDrip(); loadPreset(BEAKER_PRESETS.find(p => p.id === e.target.value)); renderAll(); });
  $('bk-restart').addEventListener('click', () => { stopDrip(); loadPreset(B.preset); renderAll(); });
  $('bk-empty').addEventListener('click', () => { stopDrip(); B.beaker = new Beaker(); B.particles = []; B.bubbles = []; renderAll(); });
  $('element-card').querySelectorAll('[data-r]').forEach(b => b.addEventListener('click', () => {
    B.reagent = b.dataset.r;
    const rr = REAGENTS.find(x => x.id === B.reagent);
    if (rr.kind === 'ind') B.amount = 0.3;
    else if (rr.kind === 'metal' || rr.kind === 'solid') B.amount = 0.5;
    else if (B.amount < 0.5) B.amount = 5;
    renderSide(); renderControls();
  }));
}

function renderControls() {
  const r = REAGENTS.find(x => x.id === B.reagent);
  const unit = reagentUnit(r);
  const amounts = unit === 'g' ? [0.1, 0.5, 1, 2] : r.kind === 'ind' ? [0.1, 0.3, 1] : [0.5, 1, 2, 5, 10, 25, 50];
  $('controls').innerHTML = `
    <div class="ctl"><span class="lbl">Aggiungi ${r.formula}</span>
      <div class="seg" id="bk-amt">${amounts.map(a => `<button type="button" data-v="${a}" aria-pressed="${a === B.amount}">${nf(a, a < 1 ? 1 : 0)} ${unit}</button>`).join('')}</div>
      <div class="btn-row"><button type="button" class="btn" id="bk-add">Versa ${nf(B.amount, B.amount < 1 ? 1 : 0)} ${unit}</button>
      ${unit === 'mL' && r.kind === 'sol' ? `<button type="button" class="btn" id="bk-drip" aria-pressed="${B.drip}">${B.drip ? 'Ferma la buretta' : 'Titola goccia a goccia'}</button>` : ''}</div>
      ${B.message ? `<p class="hint warn">${B.message}</p>` : ''}
    </div>`;
  $('bk-amt').querySelectorAll('button').forEach(b => b.addEventListener('click', () => { B.amount = +b.dataset.v; renderControls(); }));
  $('bk-add').addEventListener('click', () => addReagent(B.reagent, B.amount));
  $('bk-drip')?.addEventListener('click', () => {
    if (B.drip) { stopDrip(); renderControls(); return; }
    B.drip = true;
    B.dripTimer = setInterval(() => {
      if (!active) { stopDrip(); return; }
      addReagent(B.reagent, 0.25);
      if (B.beaker.volumeML > CAPACITY - 1) stopDrip();
    }, 250);
    renderControls();
  });
}

function renderPanelBody() {
  const s = B.beaker.summary();
  const pOH = 14 - s.pH;
  const top = s.species.filter(x => x.c > 1e-9).slice(0, 14);
  $('viewport-title').innerHTML = `${B.preset.name}<small>${nf(s.V, 1)} mL · ${nf(s.T - 273.15, 1)} °C</small>`;
  $('viewport-legend').innerHTML = '';
  $('viewport-note').innerHTML = 'Il colore della soluzione è calcolato: luce bianca attraverso 1,5 cm di soluzione, legge di Lambert–Beer e funzioni colorimetriche CIE 1931.';
  $('panel-body').innerHTML = `
    <div><h3>Misure</h3>
    <dl class="info-list">
      <dt>pH = −log a(H⁺)</dt><dd>${nf(s.pH, 2)}</dd>
      <dt>pOH</dt><dd>${nf(pOH, 2)}</dd>
      <dt>[H⁺]</dt><dd>${sci(s.species.find(x => x.id === 'H+')?.c ?? 0)} M</dd>
      <dt>Temperatura</dt><dd>${nf(s.T - 273.15, 2)} °C</dd>
      <dt>Volume</dt><dd>${nf(s.V, 2)} mL</dd>
      <dt>Forza ionica I</dt><dd>${nf(s.I, 4)} M</dd>
      <dt>Calore liberato (totale)</dt><dd>${nf(s.heat, 3)} kJ</dd>
    </dl></div>
    ${s.solids.length ? `<div><h3>Precipitati</h3><dl class="info-list">${s.solids.map(x => `<dt><i class="swatch" style="background:${x.color}"></i> ${x.label} (K<sub>sp</sub> = 10<sup>${nf(x.Ksp, 2)}</sup>)</dt><dd>${nf(x.mass * 1000, 2)} mg</dd>`).join('')}</dl></div>` : ''}
    ${s.metals.length ? `<div><h3>Metalli</h3><dl class="info-list">${s.metals.map(x => `<dt><i class="swatch" style="background:${x.color}"></i> ${METALS[x.M].name}</dt><dd>${nf(x.mass * 1000, 1)} mg</dd>`).join('')}</dl></div>` : ''}
    <div><h3>Specie in soluzione</h3>
    <table class="data-table"><thead><tr><th>Specie</th><th class="num">concentrazione (M)</th></tr></thead><tbody>
      ${top.map(x => `<tr><td>${x.label}</td><td class="num">${sci(x.c)}</td></tr>`).join('')}
    </tbody></table></div>`;
}

function drawCharts() {
  if (!active) return;
  const s = B.beaker.summary();
  const eq = B.beaker.st.eq;
  // 1. spettro di assorbimento
  $('radial-title').textContent = 'Spettro di assorbimento (1 cm)';
  if (eq) {
    const xs = [], ys = [];
    for (let l = 350; l <= 800; l += 2) { xs.push(l); ys.push(absorbance(eq.species, l, 1)); }
    const ymax = Math.max(0.1, ...ys) * 1.15;
    drawXY($('chart-radial'), {
      series: [{ xs, ys, color: cssVar('--accent'), width: 2, label: 'A(λ) = Σ ε(λ) c ℓ', fill: `color-mix(in srgb, ${cssVar('--accent')} 18%, transparent)` }],
      xmin: 350, xmax: 800, ymin: 0, ymax: Math.min(ymax, 6), xlabel: 'λ (nm)', ylabel: 'assorbanza',
      hbands: [],
    });
    const [r, g, b] = s.color;
    $('radial-note').innerHTML = `Colore trasmesso <i class="swatch" style="background: rgb(${r},${g},${b}); width: 36px"></i>: la soluzione assorbe le lunghezze d'onda dei picchi e mostra il colore complementare.`;
  }
  // 2. pH in funzione del volume
  const H = B.beaker.history;
  $('levels-title').textContent = 'Curva di titolazione';
  if (H.length > 1) {
    const xs = H.map(p => p.V), ys = H.map(p => p.pH);
    drawXY($('chart-levels'), {
      series: [{ xs, ys, color: cssVar('--phase-neg'), width: 2.2, label: 'pH' }],
      xmin: Math.min(...xs), xmax: Math.max(...xs, Math.min(...xs) + 1), ymin: 0, ymax: 14, xlabel: 'volume totale (mL)', ylabel: 'pH',
      hlines: [{ y: 7, label: 'neutro', color: cssVar('--muted') }],
      legend: false,
    });
    $('levels-note').innerHTML = 'Ogni punto è un\'aggiunta. Con "Titola goccia a goccia" si costruisce la curva di titolazione: il salto di pH segna il punto equivalente.';
  }
  // 3. speciazione
  $('slice-title').textContent = 'Concentrazioni (scala logaritmica)';
  const top = s.species.filter(x => x.c > 1e-12 && x.id !== 'H+' && x.id !== 'OH-').slice(0, 10);
  if (top.length) {
    const c = $('chart-slice');
    drawBarsLog(c, [...top, { label: 'H⁺', c: s.species.find(x => x.id === 'H+')?.c ?? 1e-7 }, { label: 'OH⁻', c: s.species.find(x => x.id === 'OH-')?.c ?? 1e-7 }]);
    $('slice-note').innerHTML = 'Tutte le specie sono calcolate insieme: bilanci di massa, costanti di equilibrio, coefficienti di attività di Davies.';
  }
}

function drawBarsLog(canvas, items) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
  const g = canvas.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  const ink = cssVar('--text'), muted = cssVar('--muted'), grid = cssVar('--chart-grid');
  const pad = { l: 118, r: 14, t: 10, b: 26 };
  const pw = w - pad.l - pad.r;
  const rowH = Math.min(22, (h - pad.t - pad.b) / items.length);
  const lo = -14, hi = 1;
  const X = (v) => pad.l + (Math.max(lo, Math.min(hi, Math.log10(v))) - lo) / (hi - lo) * pw;
  g.font = '11px "IBM Plex Mono", monospace';
  g.strokeStyle = grid; g.fillStyle = muted; g.textAlign = 'center'; g.textBaseline = 'top';
  for (let e = lo; e <= hi; e += 2) {
    const x = X(Math.pow(10, e));
    g.beginPath(); g.moveTo(x, pad.t); g.lineTo(x, h - pad.b); g.stroke();
    g.fillText(`10${String(e).split('').map(c => SUP[c]).join('')}`, x, h - pad.b + 5);
  }
  items.forEach((it, k) => {
    const y = pad.t + k * rowH;
    g.fillStyle = k >= items.length - 2 ? cssVar('--phase-neg') : cssVar('--accent');
    g.fillRect(pad.l, y + 3, X(it.c) - pad.l, rowH - 6);
    g.fillStyle = ink; g.textAlign = 'right'; g.textBaseline = 'middle';
    g.font = '11.5px "IBM Plex Sans", sans-serif';
    g.fillText(it.label, pad.l - 6, y + rowH / 2);
  });
}

function renderAnalysis() {
  const log = B.beaker.log.slice().reverse().slice(0, 25);
  const s = B.beaker.summary();
  const present = new Set(s.solids.map(x => x.id));
  const near = Object.entries(s.SI).filter(([id, v]) => !present.has(id) && v > -3 && v <= 0.001).sort((a, b) => b[1] - a[1]).slice(0, 6);
  $('analysis').innerHTML = `
    <h2>Che cosa è successo nel becher</h2>
    <div class="analysis-grid">
      <div class="wide"><h3>Registro delle aggiunte</h3>
      ${log.length ? `<table class="data-table"><thead><tr><th>Aggiunta</th><th>Equazione netta osservata</th><th class="num">pH</th><th class="num">calore (J)</th></tr></thead><tbody>
      ${log.map(e => `<tr><td>${nf(e.amount, e.amount < 1 ? 2 : 1)} ${reagentUnit(e.reagent)} ${e.reagent.formula}</td><td class="eqn">${e.equation?.text ?? '<span class="desc-muted">nessuna reazione (diluizione)</span>'}</td><td class="num">${nf(e.pH, 2)}</td><td class="num">${e.q >= 0 ? '+' : '−'}${nf(Math.abs(e.q * 1000), 0)}</td></tr>`).join('')}
      </tbody></table>` : '<p class="desc-muted">Il becher è vuoto.</p>'}
      <p class="desc-muted">L'equazione netta confronta le quantità prima e dopo l'aggiunta (ioni spettatori esclusi). Quando avvengono più equilibri insieme compare la loro somma, con coefficienti frazionari. Calore positivo: reazione esotermica.</p></div>
      ${near.length ? `<div><h3>Vicino alla precipitazione</h3><p class="desc">Indice di saturazione SI = log(Q/K<sub>sp</sub>): precipita quando SI > 0.</p>
        <dl class="info-list">${near.map(([id, v]) => `<dt>${id.replace('(s)', '')}</dt><dd>${nf(v, 2)}</dd>`).join('')}</dl></div>` : ''}
      <article><h3>Come viene calcolato</h3>
        <p class="desc">Per ogni componente (ioni liberi, H⁺) vale il bilancio di massa; per ogni specie la legge di azione di massa con le attività:</p>
        <p class="eq">c<sub>s</sub> = β<sub>s</sub> Π (γ<sub>i</sub> c<sub>i</sub>)<sup>a<sub>si</sub></sup> / γ<sub>s</sub>,&nbsp; log γ = −0,509 z² [√I/(1+√I) − 0,3 I]</p>
        <p class="desc">Il sistema non lineare si risolve con il metodo di Newton. Un solido precipita quando il prodotto ionico supera K<sub>sp</sub>. I metalli reagiscono finché la f.e.m. di Nernst si annulla:</p>
        <p class="eq">E = E° − (0,05916/n) log Q,&nbsp; q = −ΔH = −Σ ν ΔfH°,&nbsp; ΔT = q/(m c<sub>p</sub>)</p>
      </article>
      <article><h3>Limiti</h3>
        <p class="desc-muted">Costanti a 25 °C (non corrette per la temperatura). Coefficienti di attività di Davies, affidabili fino a I ≈ 0,5 M. Le reazioni sono all'equilibrio: non c'è cinetica (lo zinco reagisce "subito"). Per alcuni complessi manca ΔfH° e il loro calore di formazione è trascurato. Colori da bande di assorbimento tipiche; i colori dei solidi sono quelli osservati.</p>
      </article>
    </div>`;
}

