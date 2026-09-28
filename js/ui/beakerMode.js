// Modalità "Becher": chimica in soluzione con reagenti liberi. Si compone qualunque soluzione (catione +
// anione alla concentrazione voluta), si versano metalli e solidi, si scalda con una piastra; il solutore di
// equilibrio (database MINTEQ v4) calcola pH, potenziale redox, precipitati, complessi, gas, calore e colore.

import { Beaker, absorbance, prettyName, phaseLabel, phaseName } from '../chem/aqueous.js';
import {
  CATIONS, ANIONS, SPECIAL, WATER, saltRecipe, solutionOf, solidOf, specialRecipe, quickShelf, beakerPresets,
} from '../chem/beakerReagents.js';
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
const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const CAPACITY = 250; // mL

let active = false;
let raf = 0;
const PRESETS = beakerPresets();
const SHELF = quickShelf();
const B = {
  beaker: new Beaker(),
  preset: PRESETS[0],
  recipe: null,
  amount: 1,
  composer: { kind: 'salt', cat: 2, an: 2, sp: 0, form: 'solution', c: 0.1 },
  drip: false, dripTimer: 0,
  particles: [], bubbles: [],
  anim: { level: 0, color: [235, 240, 245] },
  message: '',
  started: false,
};

export function activateBeaker() {
  active = true;
  document.body.dataset.mode = 'beaker';
  $('busy').hidden = true;
  if (!B.started) { B.started = true; loadPreset(B.preset); }
  renderAll();
  loop();
}
export function deactivateBeaker() { active = false; cancelAnimationFrame(raf); stopDrip(); }
export function beakerRedraw() { if (active) drawCharts(); }

function loadPreset(p) {
  stopDrip();
  B.preset = p;
  B.beaker = new Beaker();
  B.particles = []; B.bubbles = [];
  B.message = '';
  for (const [r, v] of p.steps) pour(r, v, false);
  if (p.next) { B.recipe = p.next; B.amount = p.amount ?? 1; }
}

const unitOf = (r) => (r.kind === 'solution' || r.kind === 'water' ? 'mL' : 'g');

function pour(recipe, amount, render = true) {
  if (unitOf(recipe) === 'mL' && B.beaker.volumeML + amount > CAPACITY) {
    B.message = `Il becher (${CAPACITY} mL) traboccherebbe: svuotalo o versa meno.`;
    stopDrip();
    if (render) renderAll();
    return;
  }
  if (unitOf(recipe) === 'g' && B.beaker.volumeML <= 0) {
    B.message = 'Prima versa un liquido (anche solo acqua distillata).';
    if (render) renderAll();
    return;
  }
  const before = B.beaker.summary();
  let entry;
  try { entry = B.beaker.add(recipe, amount); } catch (e) { B.message = `Errore nel calcolo: ${e.message}`; if (render) renderAll(); return; }
  if (!entry) return;
  B.message = '';
  const after = B.beaker.summary();
  for (const s of after.solids) {
    if (s.metal) continue;
    const old = before.solids.find(x => x.n === s.n)?.mass ?? 0;
    const add = Math.min(220, Math.round((s.mass - old) * 900));
    for (let k = 0; k < add; k++) B.particles.push({ n: s.n, color: s.color, x: Math.random(), y: Math.random() * 0.8, vy: 0.002 + Math.random() * 0.004, settled: false, r: 1.2 + Math.random() * 1.8 });
  }
  for (const s of before.solids) {
    const now = after.solids.find(x => x.n === s.n)?.mass ?? 0;
    if (now < s.mass) { const frac = s.mass > 0 ? now / s.mass : 0; B.particles = B.particles.filter(p => p.n !== s.n || Math.random() < frac); }
  }
  if (B.particles.length > 1500) B.particles.splice(0, B.particles.length - 1500);
  const gas = (entry.h2 ?? 0) + (entry.co2 ?? 0);
  if (gas > 0) {
    const n = Math.min(160, 20 + Math.round(gas * 4000));
    for (let k = 0; k < n; k++) B.bubbles.push({ x: 0.2 + 0.6 * Math.random(), y: 1 + Math.random() * 1.5, r: 1.5 + Math.random() * 3, v: 0.004 + Math.random() * 0.006 });
  }
  if (!after.converged) B.message = 'Attenzione: il calcolo dell\'equilibrio non è arrivato a piena convergenza.';
  if (render) renderAll();
}

function stopDrip() { B.drip = false; clearInterval(B.dripTimer); }

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
  const narrow = w < 560;
  const mw = narrow ? Math.max(92, w * 0.28) : 132, mh = narrow ? 50 : 58;
  const bh = Math.min(h * 0.68, (w - mw - 60) * 0.95 / 0.72, w * 0.62), bw = bh * 0.72;
  const x0 = narrow ? 28 : w * 0.40 - bw / 2, y1 = h * 0.86, y0 = y1 - bh;
  const bigFont = narrow ? '600 19px "IBM Plex Mono", monospace' : '600 24px "IBM Plex Mono", monospace';
  const target = Math.min(1, s.V / CAPACITY);
  B.anim.level += (target - B.anim.level) * 0.12;
  B.anim.color = B.anim.color.map((c, k) => c + (s.color[k] - c) * 0.08);
  const liquidTop = y1 - B.anim.level * bh * 0.94;
  // piastra
  const hot = s.thermostat && s.T > 303;
  g.fillStyle = hot ? `rgba(230, 80, 40, ${Math.min(0.9, (s.T - 298) / 80)})` : cssVar('--panel-2');
  g.fillRect(x0 - 24, y1 + 4, bw + 48, 10);
  g.fillStyle = cssVar('--panel-2');
  g.fillRect(0, y1 + 14, w, h - y1 - 14);
  g.strokeStyle = line; g.beginPath(); g.moveTo(0, y1 + 14); g.lineTo(w, y1 + 14); g.stroke();
  if (s.V > 0) {
    const [r, gg, b] = B.anim.color.map((c, k) => Math.round(0.88 * c + 0.12 * [150, 185, 225][k]));
    const grad = g.createLinearGradient(x0, 0, x0 + bw, 0);
    grad.addColorStop(0, `rgba(${r},${gg},${b},0.8)`);
    grad.addColorStop(0.5, `rgba(${r},${gg},${b},0.64)`);
    grad.addColorStop(1, `rgba(${r},${gg},${b},0.84)`);
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
    g.strokeStyle = `rgba(${Math.round(r * 0.7)},${Math.round(gg * 0.7)},${Math.round(b * 0.75)},0.95)`;
    g.lineWidth = 1.5;
    g.beginPath(); g.moveTo(x0 + 3, liquidTop - 2); g.quadraticCurveTo(x0 + bw / 2, liquidTop + 4, x0 + bw - 3, liquidTop - 2); g.stroke();
  }
  const floorY = y1 - 5;
  const heap = new Map();
  for (const p of B.particles) {
    if (!p.settled) { p.y += p.vy; if (p.y >= 1) { p.y = 1; p.settled = true; } }
    const px = x0 + 8 + p.x * (bw - 16);
    let py = liquidTop + p.y * (floorY - liquidTop);
    if (p.settled) {
      const bucket = Math.floor(p.x * 40);
      const k = (heap.get(bucket) ?? 0) + 1;
      heap.set(bucket, k);
      py = floorY - k * 1.1 - p.r * 0.5;
    }
    g.fillStyle = p.color;
    g.globalAlpha = p.settled ? 0.95 : 0.8;
    g.beginPath(); g.arc(px, py, p.r, 0, 2 * Math.PI); g.fill();
  }
  g.globalAlpha = 1;
  let mx = x0 + 18;
  for (const m of s.solids.filter(x => x.metal)) {
    const n = Math.max(1, Math.min(14, Math.round(m.mass * 8)));
    for (let k = 0; k < n; k++) {
      g.fillStyle = m.color; g.strokeStyle = 'rgba(0,0,0,0.35)';
      const sz = 5 + (k % 3) * 2;
      g.beginPath(); g.ellipse(mx + (k % 7) * 9, floorY - 4 - Math.floor(k / 7) * 7, sz, sz * 0.6, (k * 0.7) % 3, 0, 2 * Math.PI); g.fill(); g.stroke();
    }
    mx += Math.min(7, n) * 9 + 30;
  }
  for (const bb of B.bubbles) {
    bb.y -= bb.v;
    const py = liquidTop + bb.y * (floorY - liquidTop);
    if (py < liquidTop) continue;
    g.strokeStyle = 'rgba(255,255,255,0.85)'; g.lineWidth = 1;
    g.beginPath(); g.arc(x0 + 10 + bb.x * (bw - 20), py, bb.r, 0, 2 * Math.PI); g.stroke();
  }
  B.bubbles = B.bubbles.filter(bb => liquidTop + bb.y * (floorY - liquidTop) > liquidTop - 2);
  g.strokeStyle = cssVar('--scene-mode') === 'light' ? '#8b97a8' : '#c9d3e0';
  g.lineWidth = 2.2;
  g.beginPath();
  g.moveTo(x0 - 10, y0 - 6); g.quadraticCurveTo(x0 - 2, y0 - 4, x0, y0 + 6);
  g.lineTo(x0, y1 - 8); g.quadraticCurveTo(x0, y1, x0 + 8, y1);
  g.lineTo(x0 + bw - 8, y1); g.quadraticCurveTo(x0 + bw, y1, x0 + bw, y1 - 8);
  g.lineTo(x0 + bw, y0);
  g.stroke();
  g.fillStyle = muted; g.strokeStyle = muted; g.lineWidth = 1;
  g.font = '11px "IBM Plex Mono", monospace'; g.textAlign = 'left';
  for (let v = 50; v <= 200; v += 50) {
    const y = y1 - (v / CAPACITY) * bh * 0.94;
    g.beginPath(); g.moveTo(x0 + bw - 22, y); g.lineTo(x0 + bw, y); g.stroke();
    g.fillText(`${v} mL`, x0 + bw - 60, y - 3);
  }
  const px0 = x0 + bw * 0.28;
  g.strokeStyle = muted; g.lineWidth = 3;
  g.beginPath(); g.moveTo(px0, y0 - 40); g.lineTo(px0, Math.max(liquidTop + 30, y0 + 20)); g.stroke();
  const mX = narrow ? w - mw - 10 : x0 + bw + 36;
  const box = (y, title, value) => {
    g.fillStyle = cssVar('--panel-2'); g.strokeStyle = line; g.lineWidth = 1;
    g.fillRect(mX, y, mw, mh); g.strokeRect(mX, y, mw, mh);
    g.fillStyle = muted; g.font = '11px "IBM Plex Sans", sans-serif'; g.textAlign = 'left'; g.fillText(title, mX + 8, y + 15);
    g.fillStyle = ink; g.font = bigFont; g.fillText(value, mX + 8, y + (narrow ? 40 : 44));
  };
  let yy = y0 - 10;
  box(yy, 'pHmetro', s.V > 0 ? nf(s.pH, 2) : '—'); yy += mh + 12;
  box(yy, 'termometro', `${nf(s.T - 273.15, 1)} °C`); yy += mh + 12;
  box(yy, 'volume', `${nf(s.V, 1)} mL`); yy += mh + 12;
  if (s.Eh !== null) box(yy, 'potenziale redox Eh', `${s.Eh >= 0 ? '+' : '−'}${nf(Math.abs(s.Eh), 3)} V`);
  if (B.drip) {
    const bx = x0 + bw * 0.62;
    g.strokeStyle = muted; g.lineWidth = 2;
    g.strokeRect(bx - 5, y0 - 150, 10, 110);
    g.beginPath(); g.moveTo(bx, y0 - 40); g.lineTo(bx, y0 - 20); g.stroke();
    const t = (performance.now() / 380) % 1;
    g.fillStyle = cssVar('--phase-neg');
    g.beginPath(); g.arc(bx, y0 - 18 + t * (liquidTop - y0 + 18), 2.5, 0, 2 * Math.PI); g.fill();
  }
  g.textAlign = 'left';
  g.font = '12px "IBM Plex Sans", sans-serif';
  let ly = y1 + 32;
  for (const sd of s.solids) {
    g.fillStyle = sd.color; g.fillRect(x0, ly - 9, 10, 10);
    g.strokeStyle = line; g.strokeRect(x0, ly - 9, 10, 10);
    g.fillStyle = ink;
    g.fillText(`${sd.label}(s)${sd.metal ? '' : '↓'} ${nf(sd.mass * 1000, 1)} mg`, x0 + 16, ly);
    ly += 16;
  }
  if (s.gasH2 > 0 || s.gasCO2 > 0) {
    g.fillStyle = muted;
    const Vm = 24.47;
    g.fillText(`gas liberati: ${s.gasH2 > 0 ? `H₂ ${nf(s.gasH2 * Vm * 1000, 1)} mL ` : ''}${s.gasCO2 > 0 ? `CO₂ ${nf(s.gasCO2 * Vm * 1000, 1)} mL` : ''} (25 °C, 1 atm)`, x0, ly + 4);
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
  const c = B.composer;
  const salt = saltRecipe(c.cat, c.an);
  const opt = (list, cur) => list.map(([sp, f, z, name], k) => `<option value="${k}" ${k === cur ? 'selected' : ''}>${name} (${prettyName(sp)})</option>`).join('');
  let shelfHtml = '';
  let offset = 0;
  B.shelfIndex = [];
  for (const [group, list] of SHELF) {
    shelfHtml += `<label class="lbl">${group}</label><div class="mol-list">${list.map((r) => { const k = offset++; B.shelfIndex[k] = r; return `<button type="button" class="mol-chip ${B.recipe && B.recipe.label === r.label && B.recipe.kind === r.kind ? 'active' : ''}" data-shelf="${k}" title="${esc(r.name ?? r.label)}">${esc(r.label)}</button>`; }).join('')}</div>`;
  }
  $('element-card').innerHTML = `
    <h3 class="side-h">Banco di chimica</h3>
    <label class="lbl" for="bk-preset">Esperienza guidata</label>
    <select id="bk-preset">${PRESETS.map(p => `<option value="${p.id}" ${p.id === B.preset.id ? 'selected' : ''}>${p.name}</option>`).join('')}</select>
    <p class="mol-note">${B.preset.text}</p>
    <div class="btn-row"><button type="button" class="btn" id="bk-restart">Ricomincia</button><button type="button" class="btn" id="bk-empty">Svuota il becher</button></div>
    <h3 class="side-h" style="margin-top:14px">Prepara un reagente</h3>
    <div class="seg" id="bk-kind">
      <button type="button" data-v="salt" aria-pressed="${c.kind === 'salt'}">sale, acido, base</button>
      <button type="button" data-v="other" aria-pressed="${c.kind === 'other'}">molecola</button>
    </div>
    ${c.kind === 'salt' ? `<label class="lbl" for="bk-cat">Catione</label><select id="bk-cat">${opt(CATIONS, c.cat)}</select>
    <label class="lbl" for="bk-an">Anione</label><select id="bk-an">${opt(ANIONS, c.an)}</select>
    <p class="hint">Composto: <b>${salt.label}</b> · ${esc(salt.name)} · M = ${nf(salt.M, 2)} g/mol</p>
    <div class="seg" id="bk-form" style="margin-top:6px">
      <button type="button" data-v="solution" aria-pressed="${c.form === 'solution'}">soluzione</button>
      <button type="button" data-v="solid" aria-pressed="${c.form === 'solid'}">solido</button>
    </div>` : `<label class="lbl" for="bk-sp">Sostanza</label><select id="bk-sp">${SPECIAL.map((x, k) => `<option value="${k}" ${k === c.sp ? 'selected' : ''}>${esc(x.name)} (${esc(x.label)})</option>`).join('')}</select>
    <p class="hint">M = ${nf(SPECIAL[c.sp].M, 2)} g/mol</p>`}
    ${c.form === 'solution' || c.kind === 'other' ? `<label class="lbl" for="bk-c">Concentrazione (mol/L)</label>
      <div class="smiles-row"><input id="bk-c" type="number" min="0.000001" max="20" step="any" value="${c.c}"><button type="button" class="btn" id="bk-use">Usa</button></div>`
      : `<div class="btn-row"><button type="button" class="btn" id="bk-use">Usa il solido</button></div>
      <p class="hint">${salt.dsol !== null ? `Entalpia di soluzione ${salt.dsol > 0 ? '+' : ''}${nf(salt.dsol, 2)} kJ/mol (CRC).` : 'Per questo sale il calore di dissoluzione non è nei dati: l\'effetto termico della dissoluzione del solido non viene contato.'}</p>`}
    <h3 class="side-h" style="margin-top:14px">Scaffale rapido</h3>
    ${shelfHtml}`;
  $('bk-preset').addEventListener('change', (e) => { loadPreset(PRESETS.find(p => p.id === e.target.value)); renderAll(); });
  $('bk-restart').addEventListener('click', () => { loadPreset(B.preset); renderAll(); });
  $('bk-empty').addEventListener('click', () => { stopDrip(); B.beaker = new Beaker(); B.particles = []; B.bubbles = []; B.message = ''; renderAll(); });
  $('bk-kind').querySelectorAll('button').forEach(b => b.addEventListener('click', () => { c.kind = b.dataset.v; renderSide(); }));
  $('bk-cat')?.addEventListener('change', (e) => { c.cat = +e.target.value; renderSide(); });
  $('bk-an')?.addEventListener('change', (e) => { c.an = +e.target.value; renderSide(); });
  $('bk-sp')?.addEventListener('change', (e) => { c.sp = +e.target.value; c.c = SPECIAL[c.sp].c; renderSide(); });
  $('bk-form')?.querySelectorAll('button').forEach(b => b.addEventListener('click', () => { c.form = b.dataset.v; renderSide(); }));
  $('bk-use').addEventListener('click', () => {
    if (c.form === 'solution' || c.kind === 'other') {
      const v = parseFloat(String($('bk-c').value).replace(',', '.'));
      if (!(v > 0 && v <= 20)) { B.message = 'Concentrazione non valida (tra 0 e 20 mol/L).'; renderControls(); return; }
      c.c = v;
      B.recipe = c.kind === 'other' ? specialRecipe(SPECIAL[c.sp], v) : solutionOf(saltRecipe(c.cat, c.an), v);
      if (B.amount > 50) B.amount = 5;
    } else {
      B.recipe = solidOf(saltRecipe(c.cat, c.an));
      B.amount = 0.5;
    }
    renderSide(); renderControls();
  });
  $('element-card').querySelectorAll('[data-shelf]').forEach(b => b.addEventListener('click', () => {
    B.recipe = B.shelfIndex[+b.dataset.shelf];
    if (unitOf(B.recipe) === 'g') B.amount = 0.5;
    else if (B.recipe.species?.some(([s]) => /^H(PP|BTB|MO)$/.test(s))) B.amount = 0.3;
    else if (B.amount < 0.5 || B.amount > 50) B.amount = 5;
    renderSide(); renderControls();
  }));
}

function renderControls() {
  const r = B.recipe ?? WATER;
  const unit = unitOf(r);
  const st = B.beaker.st;
  $('controls').innerHTML = `
    <div class="ctl"><span class="lbl">Da versare: ${esc(r.label)}${r.name ? ` · ${esc(r.name)}` : ''}</span>
      <div class="range-row"><input id="bk-amt" type="number" min="0.01" step="any" value="${B.amount}" style="width:90px"> <span>${unit}</span>
      <button type="button" class="btn" id="bk-add">Versa</button>
      ${unit === 'mL' ? `<button type="button" class="btn" id="bk-drip" aria-pressed="${B.drip}">${B.drip ? 'Ferma la buretta' : 'Buretta (0,25 mL a goccia)'}</button>` : ''}</div>
      <div class="btn-row"><button type="button" class="btn" id="bk-water">+ 10 mL di acqua</button></div>
      ${B.message ? `<p class="hint warn">${B.message}</p>` : ''}
    </div>
    <div class="ctl"><span class="lbl">Temperatura</span>
      <div class="seg" id="bk-thermo">
        <button type="button" data-v="0" aria-pressed="${!st.thermostat}" title="Il calore di reazione resta nella soluzione (adiabatico)">isolato</button>
        <button type="button" data-v="1" aria-pressed="${st.thermostat}" title="Piastra/bagno termostatico: la soluzione resta alla temperatura impostata">piastra termostatata</button>
      </div>
      <div class="range-row"><input type="range" id="bk-T" min="0" max="95" step="1" value="${Math.round(st.Tset - 273.15)}" ${st.thermostat ? '' : 'disabled'}><output id="bk-T-out">${Math.round(st.Tset - 273.15)} °C</output></div>
      <p class="hint">Le costanti di equilibrio cambiano con T secondo van 't Hoff, con le entalpie di reazione del database.</p>
    </div>`;
  $('bk-amt').addEventListener('change', (e) => { const v = parseFloat(String(e.target.value).replace(',', '.')); if (v > 0) B.amount = v; });
  $('bk-add').addEventListener('click', () => { const v = parseFloat(String($('bk-amt').value).replace(',', '.')); if (v > 0) { B.amount = v; pour(r, v); } });
  $('bk-water').addEventListener('click', () => pour(WATER, 10));
  $('bk-drip')?.addEventListener('click', () => {
    if (B.drip) { stopDrip(); renderControls(); return; }
    B.drip = true;
    B.dripTimer = setInterval(() => {
      if (!active) { stopDrip(); return; }
      pour(r, 0.25);
      if (B.beaker.volumeML > CAPACITY - 1) stopDrip();
    }, 250);
    renderControls();
  });
  $('bk-thermo').querySelectorAll('button').forEach(b => b.addEventListener('click', () => {
    B.beaker.setThermostat(b.dataset.v === '1');
    renderAll();
  }));
  $('bk-T').addEventListener('change', (e) => { B.beaker.setThermostat(true, 273.15 + +e.target.value); renderAll(); });
  $('bk-T').addEventListener('input', (e) => { $('bk-T-out').textContent = `${e.target.value} °C`; });
}

function renderPanelBody() {
  const s = B.beaker.summary();
  const top = s.species.filter(x => x.c > 1e-10).slice(0, 16);
  $('viewport-title').innerHTML = `${esc(B.preset.name)}<small>${nf(s.V, 1)} mL · ${nf(s.T - 273.15, 1)} °C</small>`;
  $('viewport-legend').innerHTML = '';
  $('viewport-note').innerHTML = 'Equilibrio completo con il database termodinamico MINTEQ v4 (EPA/USGS). Il colore è calcolato dallo spettro: luce bianca attraverso 1,5 cm di soluzione, legge di Lambert–Beer, funzioni colorimetriche CIE 1931.';
  $('panel-body').innerHTML = `
    <div><h3>Misure</h3>
    <dl class="info-list">
      <dt>pH = −log a(H⁺)</dt><dd>${s.V > 0 ? nf(s.pH, 2) : '—'}</dd>
      ${s.Eh !== null ? `<dt>Potenziale redox Eh</dt><dd>${s.Eh >= 0 ? '+' : '−'}${nf(Math.abs(s.Eh), 3)} V</dd><dt>pe = −log a(e⁻)</dt><dd>${nf(s.pe, 2)}</dd>` : '<dt>Potenziale redox Eh</dt><dd title="Nessuna coppia redox con entrambe le forme presenti (≥ 1 µmol/L): un elettrodo di platino non misurerebbe un potenziale stabile.">non definito</dd>'}
      <dt>Temperatura</dt><dd>${nf(s.T - 273.15, 2)} °C</dd>
      <dt>Volume</dt><dd>${nf(s.V, 2)} mL</dd>
      <dt>Forza ionica I</dt><dd>${nf(s.I, 4)} mol/kg</dd>
      <dt>Calore ceduto (totale)</dt><dd>${nf(s.heat, 3)} kJ</dd>
    </dl></div>
    ${s.solids.length ? `<div><h3>Solidi nel becher</h3><dl class="info-list">${s.solids.map(x => `<dt><i class="swatch" style="background:${x.color}"></i> ${x.label} <span class="desc-muted">${esc(x.name)}${x.metal ? '' : `, log K<sub>sp</sub> = ${nf(x.logKsp, 2)}`}</span></dt><dd>${nf(x.mass * 1000, 2)} mg</dd>`).join('')}</dl></div>` : ''}
    <div><h3>Specie in soluzione</h3>
    <table class="data-table"><thead><tr><th>Specie</th><th class="num">c (mol/L)</th><th class="num">γ</th></tr></thead><tbody>
      ${top.map(x => `<tr title="${x.src ? `fonte dei dati: ${esc(x.src)}` : ''}"><td>${x.label}</td><td class="num">${sci(x.c)}</td><td class="num">${nf(x.c > 0 ? x.a / x.c : 1, 3)}</td></tr>`).join('')}
    </tbody></table>
    <p class="hint">Passa il mouse su una specie per vedere la fonte della sua costante (codici del database MINTEQ, per esempio NIST46 = NIST Critical Stability Constants).</p></div>`;
}

function drawCharts() {
  if (!active) return;
  const s = B.beaker.summary();
  const eq = B.beaker.st.eq;
  $('radial-title').textContent = 'Spettro di assorbimento (1 cm)';
  if (eq) {
    const xs = [], ys = [];
    for (let l = 350; l <= 800; l += 2) { xs.push(l); ys.push(absorbance(eq.species, l, 1)); }
    const ymax = Math.max(0.1, ...ys) * 1.15;
    drawXY($('chart-radial'), {
      series: [{ xs, ys, color: cssVar('--accent'), width: 2, label: 'A(λ) = Σ ε(λ) c ℓ', fill: `color-mix(in srgb, ${cssVar('--accent')} 18%, transparent)` }],
      xmin: 350, xmax: 800, ymin: 0, ymax: Math.min(ymax, 6), xlabel: 'λ (nm)', ylabel: 'assorbanza',
    });
    const [r, g, b] = s.color;
    $('radial-note').innerHTML = `Colore trasmesso <i class="swatch" style="background: rgb(${r},${g},${b}); width: 36px"></i>: la soluzione assorbe le lunghezze d'onda dei picchi e mostra il colore complementare.`;
  }
  const H = B.beaker.history;
  const redox = H.some(p => p.Eh !== undefined && p.Eh !== null);
  $('levels-title').textContent = 'Curva di titolazione';
  if (H.length > 1) {
    const xs = H.map(p => p.V), ys = H.map(p => p.pH);
    drawXY($('chart-levels'), {
      series: [{ xs, ys, color: cssVar('--phase-neg'), width: 2.2, label: 'pH' }],
      xmin: Math.min(...xs), xmax: Math.max(...xs, Math.min(...xs) + 1), ymin: 0, ymax: 14, xlabel: 'volume totale (mL)', ylabel: 'pH',
      hlines: [{ y: 7, label: 'neutro a 25 °C', color: cssVar('--muted') }], legend: false,
    });
    $('levels-note').innerHTML = `Ogni punto è un'aggiunta. Con la buretta si costruisce la curva di titolazione: il salto di pH segna il punto equivalente.${redox ? '' : ''}`;
  }
  $('slice-title').textContent = 'Concentrazioni (scala logaritmica)';
  const top = s.species.filter(x => x.c > 1e-12 && x.n !== 'H+' && x.n !== 'OH-').slice(0, 10);
  if (top.length) {
    drawBarsLog($('chart-slice'), [...top, { label: 'H⁺', c: s.species.find(x => x.n === 'H+')?.c ?? 1e-7 }, { label: 'OH⁻', c: s.species.find(x => x.n === 'OH-')?.c ?? 1e-7 }]);
    $('slice-note').innerHTML = 'Tutte le specie sono calcolate insieme: bilanci di massa, costanti di equilibrio (van \'t Hoff), coefficienti di attività (Debye–Hückel esteso / Davies).';
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
  const present = new Set(s.solids.map(x => x.n));
  const near = Object.entries(s.SI).filter(([n, v]) => !present.has(n) && n !== 'CO2(g)' && !/metal$/.test(n) && v > -3 && v <= 0.001).sort((a, b) => b[1] - a[1]).slice(0, 6);
  $('analysis').innerHTML = `
    <h2>Che cosa è successo nel becher</h2>
    <div class="analysis-grid">
      <div class="wide"><h3>Registro delle aggiunte</h3>
      ${log.length ? `<table class="data-table"><thead><tr><th>Aggiunta</th><th>Equazione netta osservata</th><th class="num">pH</th><th class="num">calore (J)</th></tr></thead><tbody>
      ${log.map(e => `<tr><td>${nf(e.amount, e.amount < 1 ? 2 : 1)} ${unitOf(e.recipe)} ${esc(e.recipe.label)}</td><td class="eqn">${e.equation?.text ?? '<span class="desc-muted">nessuna trasformazione apprezzabile</span>'}</td><td class="num">${nf(e.pH, 2)}</td><td class="num">${e.q >= 0 ? '+' : '−'}${nf(Math.abs(e.q * 1000), 0)}${e.heatKnown ? '' : '*'}</td></tr>`).join('')}
      </tbody></table>` : '<p class="desc-muted">Il becher è vuoto.</p>'}
      <p class="desc-muted">L'equazione netta confronta le quantità prima e dopo l'aggiunta (ioni spettatori esclusi): quando avvengono più equilibri insieme compare la loro somma. Calore positivo = reazione esotermica. * = per qualche specie coinvolta il database non riporta ΔrH e il suo contributo è trascurato.</p></div>
      ${near.length ? `<div><h3>Vicino alla precipitazione</h3><p class="desc">Indice di saturazione SI = log(Q/K<sub>sp</sub>): un solido precipita quando SI > 0.</p>
        <dl class="info-list">${near.map(([n, v]) => `<dt>${phaseLabel(n)} ${phaseName(n) !== n ? `<span class="desc-muted">${esc(phaseName(n))}</span>` : ''}</dt><dd>${nf(v, 2)}</dd>`).join('')}</dl></div>` : ''}
      <article><h3>Come viene calcolato</h3>
        <p class="desc">Per ogni componente (ioni principali, H⁺ ed elettroni) vale il bilancio di massa; per ogni specie la legge di azione di massa con le attività:</p>
        <p class="eq">c<sub>s</sub> = K<sub>s</sub>(T) Π (γ<sub>i</sub> c<sub>i</sub>)<sup>a<sub>si</sub></sup> / γ<sub>s</sub>,&nbsp; log γ = −A z² √I/(1 + B å √I) + b I</p>
        <p class="desc">Il sistema non lineare si risolve con il metodo di Newton; i solidi entrano ed escono secondo la regola delle fasi. Gli elettroni come componente danno l'equilibrio redox: la stessa costante di MINTEQ per Fe³⁺ + e⁻ = Fe²⁺ (log K = 13,03) descrive la titolazione con il permanganato e il potenziale di Nernst.</p>
        <p class="eq">log K(T) = log K(298) − ΔrH/(R ln10)·(1/T − 1/298,15),&nbsp; q = −Σ Δn ΔrH,&nbsp; ΔT = q/(m c<sub>p</sub>)</p>
      </article>
      <article><h3>Fonti e limiti</h3>
        <p class="desc-muted">Costanti ed entalpie: MINTEQ v4 (U.S. EPA), distribuito con PHREEQC (USGS); ammino-complessi mancanti dal database LLNL; FeSCN²⁺ da Inorg. Chim. Acta 2018; indicatori da Harris, Quantitative Chemical Analysis. Solidi ammessi: solo quelli che precipitano davvero in laboratorio (idrossidi amorfi, calcite, gesso…), non i minerali che si formano in tempi geologici. L'ossidazione e la riduzione dell'acqua sono escluse (sono lente): solo i metalli sviluppano H₂, con una sovratensione di 0,40 V. Nessuna cinetica: ogni aggiunta arriva subito all'equilibrio. La capacità termica è quella dell'acqua pura.</p>
      </article>
    </div>`;
}

