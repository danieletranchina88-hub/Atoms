// Modalità "Come nasce un legame": due atomi si avvicinano e a ogni distanza l'equazione di Schrödinger viene risolta
// (Hartree–Fock ristretto e non ristretto, interazione di configurazioni completa per H₂, GFN2-xTB). Si vedono gli
// orbitali molecolari che nascono dagli orbitali atomici, la densità elettronica che si sposta nella regione di legame,
// gli spin che si accoppiano, la curva di energia confrontata con le costanti spettroscopiche sperimentali e la
// scomposizione dell'energia in cinetica e potenziale (Ruedenberg) con il teorema del viriale.

import { DIATOMICS, experimentalWell, morseCurve } from '../chem/diatomicData.js';
import { geometry, methodsFor, METHOD_INFO, curveSummary, derivative, HARTREE_EV, ANGSTROM_TO_BOHR } from '../chem/bondFormation.js';
import { buildBasis, basisValues } from '../chem/integrals.js';
import { drawXY } from './chemCharts.js';

const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const nf = (v, d = 2) => (v === null || v === undefined || !Number.isFinite(v) ? '—' : Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d }));
const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const BOHR_A = 1 / ANGSTROM_TO_BOHR;
const FONT = '"IBM Plex Sans", system-ui, sans-serif';
const MONO = '"IBM Plex Mono", ui-monospace, monospace';

const METHOD_COLOR = { exact: '#f0b429', fci: '#f0b429', rhf: '#4aa3df', uhf: '#e0605e', gfn2: '#5cb87a' };
const VIEWS = [
  { id: 'mo', label: 'orbitale ψ', note: 'Orbitale molecolare scelto nella tabella o nel diagramma: rosso e blu sono i due segni di ψ, la linea chiara è il nodo (ψ = 0).' },
  { id: 'deform', label: 'densità di legame Δρ', note: 'Δρ = ρ(molecola) − ρ(atomi liberi e sferici, alla stessa distanza): rosso dove gli elettroni si accumulano, blu da dove se ne vanno.' },
  { id: 'density', label: 'densità ρ', note: 'Densità elettronica totale ρ(r) nel piano che contiene i nuclei, in scala logaritmica (curve di livello ogni fattore 2).' },
  { id: 'spin', label: 'spin ρα − ρβ', note: 'Densità di spin: rosso dove prevalgono gli elettroni α (↑), blu dove prevalgono i β (↓). Solo per le soluzioni non ristrette.' },
];

let active = false;
let worker = null;
const F = {
  mol: DIATOMICS.find(d => d.id === 'H2'),
  basis: '6-31G**',
  token: 0,
  grid: [],
  ctx: null,
  data: {},          // metodo → { points: [], Eatoms, done, ms }
  method: null,      // metodo mostrato nei dettagli
  idx: 0,            // indice della distanza nella griglia
  view: 'deform',
  sel: null,         // { label, spin: 'a' | 'b' }
  show: { gfn2: false },
  chart3: 'tv',
  playing: false,
  timer: null,
  prevC: null,
  error: null,
};

export function activateFormation() {
  active = true;
  document.body.dataset.mode = 'formation';
  $('busy').hidden = true;
  if (!F.ctx && !F.running) run();
  renderAll();
}

export function deactivateFormation() {
  active = false;
  stop();
  // un calcolo non finito si interrompe (consuma CPU): riparte quando si torna in questa modalità
  if (F.running && worker) { worker.terminate(); worker = null; F.running = false; F.ctx = null; }
}

export function formationRedraw() {
  if (!active) return;
  drawMap();
  drawCharts();
}

// ---------------------------------------------------------------------------
// Calcolo nel worker
// ---------------------------------------------------------------------------

function run() {
  stop();
  if (worker) worker.terminate();
  worker = new Worker(new URL('../chem/bondFormationWorker.js', import.meta.url), { type: 'module' });
  worker.onmessage = onMessage;
  const methods = methodsFor(F.mol);
  F.token++;
  F.ctx = null;
  F.grid = [];
  F.data = Object.fromEntries(methods.map(m => [m, { points: [], Eatoms: null, done: false, ms: 0 }]));
  F.method = methods.includes('fci') ? 'fci' : methods.includes('exact') ? 'exact' : methods.includes('rhf') ? 'rhf' : methods[0];
  F.sel = null;
  F.prevC = null;
  F.error = null;
  F.running = true;
  F.started = performance.now();
  worker.postMessage({ type: 'run', token: F.token, id: F.mol.id, basis: F.basis, methods });
}

let redrawPending = false;
function scheduleRedraw() {
  if (redrawPending) return;
  redrawPending = true;
  requestAnimationFrame(() => { redrawPending = false; if (active) renderAll(); });
}

function onMessage(ev) {
  const m = ev.data;
  if (m.token !== F.token) return;
  if (m.type === 'context') {
    F.ctx = m;
    F.grid = m.grid;
    // si parte con gli atomi già vicini all'equilibrio sperimentale
    F.idx = nearestIndex(F.mol.exp.re);
  } else if (m.type === 'point') {
    const d = F.data[m.method];
    d.points[m.i] = m.point;
    if (m.method === 'gfn2') d.Eatoms ??= null;
  } else if (m.type === 'done') {
    const d = F.data[m.method];
    d.done = true; d.Eatoms = m.Eatoms; d.ms = m.ms;
  } else if (m.type === 'finished') {
    F.running = false;
  } else if (m.type === 'error') {
    F.error = m.message; F.running = false;
  }
  scheduleRedraw();
}

const nearestIndex = (R) => F.grid.reduce((b, r, i) => (Math.abs(r - R) < Math.abs(F.grid[b] - R) ? i : b), 0);
const currentPoint = (method = F.method) => F.data[method]?.points[F.idx] ?? null;
const methodEatoms = (method) => (method === 'gfn2' ? F.data.gfn2?.Eatoms : F.ctx?.Eatoms);

/** Curva di un metodo: solo i punti già calcolati, ordinati per R crescente. */
function curve(method) {
  const d = F.data[method];
  if (!d) return { R: [], E: [], pts: [] };
  const pts = d.points.filter(Boolean).sort((a, b) => a.R - b.R);
  return { R: pts.map(p => p.R), E: pts.map(p => p.E), pts };
}

/**
 * Energia cinetica e potenziale rispetto agli atomi separati dalla curva E(R), con il teorema del viriale molecolare
 * (Slater 1933): ΔT = −ΔE − R dE/dR, ΔV = 2ΔE + R dE/dR (eV). Con una base finita il viriale della funzione d'onda non
 * è esatto (−V/T ≠ 2 di qualche millesimo, cioè alcuni eV sull'energia cinetica totale di un atomo come N): per questo
 * le variazioni dovute al legame si ricavano dalla curva, che il viriale lega esattamente a T e V.
 */
function virialCurve(method = F.method) {
  const c = curve(method), Ea = methodEatoms(method);
  if (c.R.length < 3 || Ea == null) return null;
  const d = derivative(c.R, c.E);
  const dE = c.E.map(e => e - Ea);
  return {
    R: c.R, pts: c.pts, dE: dE.map(v => v * HARTREE_EV),
    dT: dE.map((e, i) => (-e - c.R[i] * ANGSTROM_TO_BOHR * d[i]) * HARTREE_EV),
    dV: dE.map((e, i) => (2 * e + c.R[i] * ANGSTROM_TO_BOHR * d[i]) * HARTREE_EV),
  };
}

/** ΔT e ΔV al punto corrente (null se la curva non è ancora pronta attorno a esso). */
function virialAt(p) {
  const v = virialCurve();
  if (!v) return null;
  const i = v.pts.indexOf(p);
  if (i < 0) return null;
  // la derivata richiede i vicini sulla griglia già calcolati
  const gi = F.grid.indexOf(p.R), pts = F.data[F.method].points;
  const near = F.data[F.method].done || (pts[gi - 1] && pts[gi + 1]);
  return near ? { dT: v.dT[i], dV: v.dV[i] } : null;
}

// ---------------------------------------------------------------------------
// Animazione
// ---------------------------------------------------------------------------

function play() {
  if (!F.grid.length) return;
  F.playing = true;
  // dagli atomi separati verso l'interno; dopo una pausa si riprende da dove ci si era fermati
  if (!F.paused || F.idx >= F.grid.length - 1) F.idx = 0;
  F.paused = false;
  const step = () => {
    if (!F.playing || !active) return;
    if (F.idx < F.grid.length - 1) {
      F.idx++;
      renderDynamic();
      // rallenta vicino all'equilibrio, dove succede di più
      const R = F.grid[F.idx], re = F.mol.exp.re;
      F.timer = setTimeout(step, Math.abs(R - re) < 0.25 * re ? 260 : 150);
    } else {
      F.playing = false;
      renderControls();
    }
  };
  renderControls();
  F.timer = setTimeout(step, 120);
}

function stop() {
  F.playing = false;
  clearTimeout(F.timer);
}

function setIndex(i) {
  F.idx = Math.max(0, Math.min(F.grid.length - 1, i));
  renderDynamic();
}

// ---------------------------------------------------------------------------
// Pannelli
// ---------------------------------------------------------------------------

function renderAll() {
  if (!active) return;
  renderSide();
  renderControls();
  renderDynamic();
  renderAnalysis();
}

/** Ciò che cambia con la distanza. */
function renderDynamic() {
  if (!active) return;
  const s = $('fm-R');
  if (s && +s.value !== F.idx) s.value = F.idx;
  const lab = $('fm-Rlabel');
  if (lab) lab.textContent = F.grid.length ? `R = ${nf(F.grid[F.idx], 3)} Å` : '';
  renderPanel();
  drawMap();
  drawCharts();
}

function renderSide() {
  const mol = F.mol, ex = experimentalWell(mol);
  const methods = methodsFor(mol);
  const rows = [...methods.map(m => {
    const d = F.data[m];
    const c = curve(m);
    const sum = d?.done && c.R.length > 5 ? curveSummary(mol, c.R, c.E, methodEatoms(m)) : null;
    const pending = !d?.done;
    const De = sum ? (sum.re ? sum.De : null) : null;
    const unrel = m === 'gfn2';
    return `<tr><td><i class="swatch" style="background:${METHOD_COLOR[m]}"></i> ${METHOD_INFO[m].short}</td>
      <td class="num">${pending ? '…' : sum?.bound ? nf(sum.re, 3) : '—'}</td>
      <td class="num">${pending ? '…' : sum?.bound ? `${nf(De, De < 0.01 ? 4 : 2)}${unrel ? '<sup>†</sup>' : ''}` : 'non legato'}</td>
      <td class="num">${pending ? '…' : sum?.bound ? nf(sum.we, 0) : '—'}</td></tr>`;
  }), `<tr class="exp-row"><td><i class="swatch dash"></i> esperimento</td><td class="num">${nf(ex.re, 4)}</td><td class="num">${nf(ex.De, ex.De < 0.01 ? 5 : 3)}</td><td class="num">${nf(ex.we, 0)}</td></tr>`];
  $('element-card').innerHTML = `
    <h3 class="side-h">Come nasce un legame</h3>
    <div class="fm-mols" role="group" aria-label="Molecola">${DIATOMICS.map(d => `<button type="button" class="chip" data-mol="${d.id}" aria-pressed="${d.id === mol.id}">${d.name}</button>`).join('')}</div>
    <p class="mol-note"><b>${mol.name}</b>, stato ${mol.state}. ${mol.blurb}</p>
    <label class="lbl" for="fm-basis">Base di funzioni gaussiane</label>
    <select id="fm-basis">${['STO-3G', '6-31G', '6-31G**'].map(b => `<option ${b === F.basis ? 'selected' : ''}>${b}</option>`).join('')}</select>
    <h4 class="fm-h">Calcolo contro esperimento</h4>
    <table class="data-table fm-table"><thead><tr><th>Metodo</th><th class="num">r<sub>e</sub> (Å)</th><th class="num">D<sub>e</sub> (eV)</th><th class="num">ω<sub>e</sub> (cm⁻¹)</th></tr></thead>
    <tbody>${rows.join('')}</tbody></table>
    <p class="hint">r<sub>e</sub>: distanza di equilibrio; D<sub>e</sub>: profondità della buca rispetto agli atomi separati; ω<sub>e</sub>: frequenza di vibrazione armonica dalla curvatura nel minimo.${methods.includes('gfn2') ? ' <sup>†</sup>GFN2-xTB è parametrizzato vicino all\'equilibrio: r<sub>e</sub> e ω<sub>e</sub> sono buoni, l\'energia di dissociazione no.' : ''}</p>
    <p class="hint">Esperimento: ${esc(mol.exp.D0src)}; r<sub>e</sub>, ω<sub>e</sub>: Huber e Herzberg (1979), NIST WebBook.${ex.zpe ? ` D<sub>e</sub> = D<sub>0</sub> + G(0), con G(0) = ω<sub>e</sub>/2 − ω<sub>e</sub>x<sub>e</sub>/4 = ${nf(ex.zpe, 3)} eV (energia di punto zero).` : ''}${mol.exp.note ? ` ${esc(mol.exp.note)}` : ''}</p>
    ${F.error ? `<p class="hint warn">Errore nel calcolo: ${esc(F.error)}</p>` : ''}`;
  $('element-card').querySelectorAll('[data-mol]').forEach(b => b.addEventListener('click', () => {
    const m = DIATOMICS.find(d => d.id === b.dataset.mol);
    if (m === F.mol) return;
    F.mol = m;
    run();
    renderAll();
  }));
  $('fm-basis').addEventListener('change', (e) => { F.basis = e.target.value; run(); renderAll(); });
}

function renderControls() {
  const methods = methodsFor(F.mol);
  const done = Object.values(F.data).filter(d => d.done).length;
  const progress = F.running ? `Calcolo in corso: ${done} di ${methods.length} metodi completati, ${Object.values(F.data).reduce((a, d) => a + d.points.filter(Boolean).length, 0)} punti…` : `${F.grid.length} distanze per metodo · ${Object.entries(F.data).map(([m, d]) => `${METHOD_INFO[m].short} ${nf(d.ms / 1000, 1)} s`).join(' · ')}`;
  $('controls').innerHTML = `
    <div class="ctl"><span class="lbl">Distanza fra i nuclei</span>
      <div class="range-row">
        <button type="button" class="btn" id="fm-play">${F.playing ? '❚❚ ferma' : '▶ avvicina gli atomi'}</button>
        <span id="fm-Rlabel" class="fm-rlabel">${F.grid.length ? `R = ${nf(F.grid[F.idx], 3)} Å` : ''}</span>
      </div>
      <input type="range" id="fm-R" min="0" max="${Math.max(0, F.grid.length - 1)}" step="1" value="${F.idx}" aria-label="Distanza fra i nuclei (da lontano a vicino)">
      <div class="fm-range-ends"><span>atomi separati</span><span>nuclei vicini</span></div>
      <p class="hint">${progress}</p>
    </div>
    <div class="ctl"><span class="lbl">Funzione d'onda mostrata</span>
      <div class="seg">${methods.filter(m => m !== 'gfn2').map(m => `<button type="button" data-method="${m}" aria-pressed="${m === F.method}">${METHOD_INFO[m].short}</button>`).join('')}</div>
    </div>
    <div class="ctl"><span class="lbl">Che cosa vedere</span>
      <div class="seg">${VIEWS.map(v => `<button type="button" data-view="${v.id}" aria-pressed="${v.id === F.view}" ${v.id === 'spin' && !hasSpin() ? 'disabled' : ''}>${v.label}</button>`).join('')}</div>
    </div>`;
  $('fm-play').addEventListener('click', () => { if (F.playing) { stop(); F.paused = true; renderControls(); } else play(); });
  $('fm-R').addEventListener('input', (e) => { stop(); F.paused = false; setIndex(+e.target.value); $('fm-play').textContent = '▶ avvicina gli atomi'; });
  $('controls').querySelectorAll('[data-method]').forEach(b => b.addEventListener('click', () => { F.method = b.dataset.method; F.sel = null; F.prevC = null; if (F.view === 'spin' && !hasSpin()) F.view = 'deform'; renderControls(); renderDynamic(); }));
  $('controls').querySelectorAll('[data-view]').forEach(b => b.addEventListener('click', () => { F.view = b.dataset.view; renderControls(); drawMap(); }));
}

const hasSpin = () => F.method === 'uhf' || F.method === 'exact' || F.mol.mult !== 1;

/** Orbitale scelto per la mappa nel punto corrente (seguito per nome al variare di R). */
function selectedMO(p) {
  if (!p) return null;
  const list = F.sel?.spin === 'b' && p.mosB ? p.mosB : p.mos;
  if (F.sel) {
    const same = list.filter(o => o.label === F.sel.label && (o.axis || '') === (F.sel.axis || ''));
    if (same.length) return { mo: same[0], spin: F.sel.spin };
  }
  // per difetto: l'orbitale occupato più alto con carattere di legame, altrimenti l'HOMO
  const occ = list.filter(o => o.occ > 0.5);
  const homo = occ[occ.length - 1];
  const bond = [...occ].reverse().find(o => o.character === 'legante' && (o.axis === '' || o.axis === 'x'));
  return { mo: bond ?? homo ?? list[0], spin: 'a' };
}

function renderPanel() {
  const p = currentPoint();
  const body = $('panel-body');
  if (!F.ctx || !p) {
    body.innerHTML = `<p class="hint">${F.error ? esc(F.error) : 'Calcolo della funzione d\'onda a questa distanza…'}</p>`;
    return;
  }
  const Ea = F.ctx.Eatoms;
  const dE = (p.E - Ea) * HARTREE_EV;
  const vir = virialAt(p);
  const sel = selectedMO(p);
  const ex = experimentalWell(F.mol);
  const sym = (s) => s.replace(/([gu])$/, '<sub>$1</sub>');
  const moRows = (list, spin) => list.filter(o => o.e * HARTREE_EV > -60 || o.occ > 0).map(o => {
    const isSel = sel && sel.mo === o && sel.spin === spin;
    const bars = o.comp.filter(c => Math.abs(c.w) > 0.04).slice(0, 4).map(c => `<span class="fm-comp ${c.atom ? 'b' : 'a'}${c.extra ? ' x' : ''}" style="flex:${Math.max(0.05, Math.abs(c.w))}" title="${esc(atomSym(c.atom))} ${c.label}: ${nf(100 * c.w, 0)} %">${atomSym(c.atom)} ${c.label}</span>`).join('');
    const occTxt = F.method === 'fci' ? nf(o.occ, 3) : o.occ > 1.5 ? '↑↓' : o.occ > 0.5 ? (spin === 'b' ? '↓' : '↑') : '';
    return `<tr class="${isSel ? 'sel' : ''}" data-mo="${esc(o.label)}" data-axis="${o.axis}" data-spin="${spin}" tabindex="0">
      <td>${sym(o.label)}${o.star ? '*' : ''}</td><td class="num">${nf(o.e * HARTREE_EV, 2)}</td><td class="fm-occ">${occTxt}</td>
      <td><span class="fm-char ${o.character === 'legante' ? 'bond' : o.character === 'antilegante' ? 'anti' : ''}">${o.character}</span></td>
      <td><div class="fm-bars">${bars}</div></td></tr>`;
  }).join('');
  const moTable = (list, spin, title) => `<table class="data-table fm-mo"><thead><tr><th>${title}</th><th class="num">ε (eV)</th><th>e⁻</th><th>carattere</th><th>da quali orbitali atomici</th></tr></thead><tbody>${moRows(list, spin)}</tbody></table>`;
  const ov = p.overlap?.slice(0, 3).map(o => `S(${atomSym(0)} ${o.a}, ${atomSym(1)} ${o.b}) = ${nf(o.S, 3)}`).join(' · ');
  const fciTxt = p.natOcc ? `<dt>Occupazioni naturali</dt><dd>${p.natOcc.slice(0, 3).map(v => nf(v, 3)).join(' · ')}</dd>
    <dt>Configurazioni principali</dt><dd>${p.configs.slice(0, 3).map(c => `${nf(100 * c.w, 1)} % ${configLabel(p, c)}`).join(' · ')}</dd>` : '';
  const spinTxt = p.S2 != null ? `<dt>⟨S²⟩ (esatto ${nf(F.mol.mult === 1 ? 0 : (F.mol.mult - 1) / 2 * ((F.mol.mult - 1) / 2 + 1), 2)})</dt><dd>${nf(p.S2, 3)}${Math.abs(p.spinA) > 0.05 ? ` · spin su ${atomSym(0)}: ${nf(p.spinA, 2)}` : ''}</dd>` : '';
  body.innerHTML = `
    <div><h3>A ${nf(p.R, 3)} Å (${nf(p.R / ex.re, 2)} r<sub>e</sub>)</h3>
    <p class="fm-story">${story(p)}</p>
    <dl class="info-list">
      <dt>Energia rispetto agli atomi separati</dt><dd>${nf(dE, 3)} eV</dd>
      <dt>Cinetica ΔT · potenziale ΔV (viriale)</dt><dd>${vir ? `${nf(vir.dT, 2)} · ${nf(vir.dV, 2)} eV` : '…'}</dd>
      <dt>−⟨V⟩/⟨T⟩ della funzione d'onda (esatto: 2 nel minimo)</dt><dd>${nf(-p.V / p.T, 4)}</dd>
      <dt>Ordine di legame (Mayer)</dt><dd>${nf(p.bondOrder, 2)}</dd>
      <dt>Cariche di Mulliken</dt><dd>${atomSym(0)} ${signed(p.charges[0])} · ${atomSym(1)} ${signed(p.charges[1])}</dd>
      ${spinTxt}${fciTxt}
      ${ov ? `<dt>Sovrapposizione fra orbitali atomici</dt><dd>${ov}</dd>` : ''}
    </dl></div>
    <div><h3>Orbitali molecolari</h3>
    ${moTable(p.mos, 'a', p.mosB ? 'spin α' : 'OM')}
    ${p.mosB ? moTable(p.mosB, 'b', 'spin β') : ''}
    <p class="hint">Carattere dalla popolazione di sovrapposizione di Mulliken fra i due atomi (positiva: legante). Composizione: peso degli orbitali degli atomi liberi, calcolati nella stessa base${F.mol.Z[0] !== F.mol.Z[1] ? '' : ''}. Clic su una riga per vederlo nella mappa.</p></div>`;
  body.querySelectorAll('[data-mo]').forEach(tr => {
    const pick = () => { F.sel = { label: tr.dataset.mo, axis: tr.dataset.axis, spin: tr.dataset.spin }; F.prevC = null; F.view = 'mo'; renderControls(); renderDynamic(); };
    tr.addEventListener('click', pick);
    tr.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pick(); } });
  });
}

const atomSym = (k) => (F.ctx ? (k ? F.ctx.B.symbol : F.ctx.A.symbol) : '');
const signed = (q) => { const r = Math.round(q * 100) / 100; return `${r > 0 ? '+' : r < 0 ? '−' : ''}${nf(Math.abs(r), 2)}`; };

function configLabel(p, c) {
  const name = (k) => { const o = p.mos[k]; return o ? `${o.label.replace(/([gu])$/, '<sub>$1</sub>')}` : `φ${k + 1}`; };
  return c.p === c.q ? `(${name(c.p)})²` : `${name(c.q)}·${name(c.p)}`;
}

/** Descrizione di ciò che accade alla distanza corrente, ricavata dai numeri calcolati. */
function story(p) {
  const ex = experimentalWell(F.mol), re = ex.re, R = p.R;
  const dE = (p.E - F.ctx.Eatoms) * HARTREE_EV;
  const dT = virialAt(p)?.dT ?? NaN;
  const c = curve(F.method);
  const sum = c.R.length > 6 ? curveSummary(F.mol, c.R, c.E, F.ctx.Eatoms) : null;
  const reCalc = sum?.bound ? sum.re : re;
  const s0 = p.overlap?.[0]?.S ?? 0;
  const parts = [];
  if (F.mol.id === 'He2') {
    parts.push(R > 2.2 ? 'Le nuvole 1s quasi non si toccano.' : `Le nuvole 1s si sovrappongono (S = ${nf(s0, 3)}): σg² e σu² sono entrambi pieni, l'effetto legante dell'uno è annullato (anzi superato) dall'antilegante dell'altro.`);
    parts.push(dE > 0.01 ? `L'energia sale di ${nf(dE, 2)} eV: è la repulsione di Pauli. L'energia cinetica cresce (${nf(dT, 2)} eV) perché gli elettroni con lo stesso spin non possono occupare la stessa regione e i loro orbitali devono deformarsi per restare ortogonali.` : 'Nessun legame chimico: resta solo la debolissima attrazione di dispersione (van der Waals), che Hartree–Fock non contiene.');
    return parts.join(' ');
  }
  if (sum && !sum.bound && R > 0.75 * re && R < 1.6 * re) {
    parts.push(`Con ${METHOD_INFO[F.method].short} la curva non scende sotto l'energia degli atomi separati: questo legame esiste solo grazie alla correlazione elettronica, cioè al modo in cui ogni elettrone evita gli altri istante per istante, che un singolo determinante di Slater non descrive (${nf(dE, 2)} eV rispetto agli atomi, contro −${nf(ex.De, 2)} eV sperimentali).`);
    return parts.join(' ');
  }
  if (R > 1.9 * reCalc) {
    parts.push(`Gli atomi sono quasi indipendenti: la sovrapposizione fra i loro orbitali è piccola (S = ${nf(s0, 3)}).`);
    if (p.method === 'rhf') parts.push(`Hartree–Fock ristretto costringe però i due elettroni di legame nello stesso orbitale anche qui: metà della funzione d'onda mette entrambi su un atomo (A⁻B⁺), e l'energia resta ${nf(dE, 1)} eV sopra quella degli atomi separati. È il fallimento noto di RHF nella dissociazione.`);
    else if (p.method === 'uhf' && F.mol.mult !== 1) parts.push(`Ogni atomo di ossigeno ha due elettroni spaiati; nella molecola ne restano due, con spin paralleli (tripletto). Un solo determinante con M<sub>S</sub> = 1 non può però descrivere due atomi O(³P) indipendenti: lontano dall'equilibrio l'energia UHF resta ${nf(dE, 2)} eV sopra quella degli atomi separati.`);
    else if (p.method === 'uhf' && p.S2 > 0.3) parts.push(`UHF lascia invece ogni elettrone spaiato sul suo atomo, con spin opposti (popolazione di spin ${nf(p.spinA, 2)} su ${atomSym(0)}): l'energia è giusta, ma il determinante mescola singoletto e tripletto (⟨S²⟩ = ${nf(p.S2, 2)} invece di 0).`);
    else if (p.method === 'fci') parts.push(`La soluzione esatta è metà (σg)² e metà (σu)² con segno opposto: è proprio la funzione di Heitler e London, un elettrone su ciascun atomo e spin accoppiati.`);
    else if (p.method === 'exact') parts.push('L\'elettrone è condiviso a metà fra i due protoni, ma l\'energia è quasi quella di H + H⁺: lo scambio fra i due nuclei è ancora lento.');
  } else if (R > 1.12 * reCalc) {
    parts.push(`Le nuvole si sovrappongono (S = ${nf(s0, 2)}) e gli orbitali atomici si combinano: in fase (orbitale legante, energia più bassa) e in opposizione di fase (antilegante).`);
    parts.push(dT < 0 && dE < 0 ? `L'energia scende (${nf(dE, 2)} eV) e scende anche l'energia cinetica (${nf(dT, 2)} eV): gli elettroni possono muoversi su entrambi i nuclei e la loro funzione d'onda si "distende" lungo il legame. Secondo Ruedenberg è questo il primo motore del legame covalente.` : `L'energia scende (${nf(dE, 2)} eV).`);
  } else if (R > 0.92 * reCalc) {
    parts.push(`Vicino al minimo (r<sub>e</sub> calcolata ${nf(reCalc, 3)} Å, sperimentale ${nf(re, 3)} Å). Qui gli orbitali si contraggono attorno ai nuclei: l'energia potenziale scende molto e quella cinetica risale, fino a ΔT ≈ −ΔE e ΔV ≈ 2ΔE come vuole il teorema del viriale nel minimo (ΔT = ${nf(dT, 2)}, ΔE = ${nf(dE, 2)} eV).`);
  } else {
    parts.push(`Più vicini dell'equilibrio: la repulsione fra i nuclei e la crescita dell'energia cinetica (${nf(dT, 2)} eV, gli elettroni compressi in poco spazio) prevalgono, e la curva risale ripida.`);
  }
  return parts.join(' ');
}

// ---------------------------------------------------------------------------
// Mappa nel piano del legame
// ---------------------------------------------------------------------------

const hexToRgb = (h) => { const s = h.replace('#', ''); const v = parseInt(s.length === 3 ? s.split('').map(c => c + c).join('') : s, 16); return [(v >> 16) & 255, (v >> 8) & 255, v & 255]; };

function setupCanvas(canvas) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

/** Campo nel piano xz (asse del legame orizzontale). */
function computeField(p, view, nx, ny, zhalf, xhalf) {
  const Rb = p.R * ANGSTROM_TO_BOHR;
  const atoms = geometry(F.mol.Z, Rb);
  const basis = buildBasis(atoms, F.basis);
  const n = basis.nbf;
  const chi = new Float64Array(n);
  const out = new Float32Array(nx * ny);
  let coef = null, P = null;
  if (view === 'mo') {
    const s = selectedMO(p);
    if (!s) return null;
    const C = s.spin === 'b' && p.Cb ? p.Cb : p.Ca;
    coef = new Float64Array(n);
    for (let i = 0; i < n; i++) coef[i] = C[i * n + s.mo.index];
    // segno coerente con il fotogramma precedente (la fase di un orbitale è arbitraria)
    if (F.prevC && F.prevC.length === n) { let d = 0; for (let i = 0; i < n; i++) d += coef[i] * F.prevC[i]; if (d < 0) for (let i = 0; i < n; i++) coef[i] = -coef[i]; }
    else { let big = 0; for (let i = 0; i < n; i++) if (Math.abs(coef[i]) > Math.abs(big)) big = coef[i]; if (big < 0) for (let i = 0; i < n; i++) coef[i] = -coef[i]; }
    F.prevC = coef;
  } else {
    P = new Float64Array(n * n);
    for (let i = 0; i < n * n; i++) {
      if (view === 'spin') P[i] = p.Pa[i] - p.Pb[i];
      else P[i] = p.Pa[i] + p.Pb[i] - (view === 'deform' ? F.ctx.Ppro[i] : 0);
    }
  }
  const tmp = new Float64Array(n);
  for (let j = 0; j < ny; j++) {
    const x = (xhalf - (j + 0.5) * 2 * xhalf / ny) * ANGSTROM_TO_BOHR;
    for (let i = 0; i < nx; i++) {
      const z = (-zhalf + (i + 0.5) * 2 * zhalf / nx) * ANGSTROM_TO_BOHR;
      basisValues(basis, x, 0, z, chi);
      let v = 0;
      if (coef) for (let a = 0; a < n; a++) v += coef[a] * chi[a];
      else {
        for (let a = 0; a < n; a++) { let s = 0; const ra = a * n; for (let b = 0; b < n; b++) s += P[ra + b] * chi[b]; tmp[a] = s; }
        for (let a = 0; a < n; a++) v += chi[a] * tmp[a];
      }
      out[j * nx + i] = v;
    }
  }
  return out;
}

/** Curve di livello con i quadrati in marcia (marching squares), in coordinate della griglia. */
function contour(field, nx, ny, level) {
  const seg = [];
  const at = (i, j) => field[j * nx + i] - level;
  for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
    const a = at(i, j), b = at(i + 1, j), c = at(i + 1, j + 1), d = at(i, j + 1);
    const pts = [];
    const edge = (v1, v2, x1, y1, x2, y2) => { if ((v1 > 0) !== (v2 > 0)) { const t = v1 / (v1 - v2); pts.push([x1 + (x2 - x1) * t, y1 + (y2 - y1) * t]); } };
    edge(a, b, i, j, i + 1, j); edge(b, c, i + 1, j, i + 1, j + 1); edge(c, d, i + 1, j + 1, i, j + 1); edge(d, a, i, j + 1, i, j);
    if (pts.length === 2) seg.push(pts);
    else if (pts.length === 4) { seg.push([pts[0], pts[1]]); seg.push([pts[2], pts[3]]); }
  }
  return seg;
}

function drawMap() {
  const canvas = $('lab-plot');
  if (!canvas || !active) return;
  const { ctx, w, h } = setupCanvas(canvas);
  const p = currentPoint();
  const ex = experimentalWell(F.mol);
  const view = F.view === 'spin' && !hasSpin() ? 'deform' : F.view;
  $('viewport-title').innerHTML = `${F.mol.name}<small>${F.grid.length ? `R = ${nf(F.grid[F.idx], 3)} Å · ` : ''}${METHOD_INFO[F.method]?.name ?? ''} / ${F.basis}</small>`;
  $('viewport-note').innerHTML = VIEWS.find(v => v.id === view).note + (view === 'mo' && p ? ` Mostrato: ${moName(p)}.` : '');
  if (!p || !F.ctx) {
    ctx.fillStyle = cssVar('--muted'); ctx.font = `13px ${FONT}`; ctx.textAlign = 'center';
    ctx.fillText(F.error ? `Errore: ${F.error}` : 'Calcolo in corso…', w / 2, h / 2);
    return;
  }
  // finestra fissa per la molecola: gli atomi si muovono davvero quando cambia R
  const Rmax = F.grid[0];
  const zhalf = Math.min(Rmax / 2 + 1.7, Math.max(ex.re / 2 + 2.6, 3.2));
  const xhalf = zhalf * h / w;
  const scale = Math.max(2, Math.round(w / 230));
  const nx = Math.round(w / scale), ny = Math.round(h / scale);
  const field = computeField(p, view, nx, ny, zhalf, xhalf);
  if (!field) return;
  const bg = hexToRgb(cssVar('--slice-bg')), pos = hexToRgb(cssVar('--phase-pos')), neg = hexToRgb(cssVar('--phase-neg')), acc = hexToRgb(cssVar('--accent'));
  let max = 0;
  for (let i = 0; i < field.length; i++) max = Math.max(max, Math.abs(field[i]));
  // scala fissa per Δρ e spin (si confrontano le distanze), adattiva per ψ
  const signedView = view !== 'density';
  const vmax = view === 'deform' ? 0.12 : view === 'spin' ? 0.12 : max || 1;
  const img = ctx.createImageData(nx, ny);
  const lmax = Math.log(Math.max(max, 1e-6)), lmin = lmax - Math.log(1e4);
  for (let k = 0; k < field.length; k++) {
    const v = field[k];
    let t, col;
    if (signedView) { t = Math.min(1, Math.pow(Math.abs(v) / vmax, 0.5)); col = v >= 0 ? pos : neg; }
    else { t = v > 0 ? Math.max(0, Math.min(1, (Math.log(v) - lmin) / (lmax - lmin))) : 0; col = acc; }
    img.data[4 * k] = bg[0] + (col[0] - bg[0]) * t;
    img.data[4 * k + 1] = bg[1] + (col[1] - bg[1]) * t;
    img.data[4 * k + 2] = bg[2] + (col[2] - bg[2]) * t;
    img.data[4 * k + 3] = 255;
  }
  const tmp = document.createElement('canvas');
  tmp.width = nx; tmp.height = ny;
  tmp.getContext('2d').putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(tmp, 0, 0, w, h);
  // curve di livello
  const sx = w / nx, sy = h / ny;
  const drawLevels = (levels, color, width, dash = []) => {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash);
    ctx.beginPath();
    for (const lv of levels) for (const [[x1, y1], [x2, y2]] of contour(field, nx, ny, lv)) { ctx.moveTo((x1 + 0.5) * sx, (y1 + 0.5) * sy); ctx.lineTo((x2 + 0.5) * sx, (y2 + 0.5) * sy); }
    ctx.stroke(); ctx.setLineDash([]);
  };
  const ink = cssVar('--slice-node');
  ctx.globalAlpha = 0.55;
  if (view === 'density') drawLevels(Array.from({ length: 12 }, (_, k) => 0.001 * 2 ** k).filter(v => v < max), ink, 0.8);
  else if (view === 'mo') { const lv = [0.02, 0.05, 0.1, 0.2, 0.4].filter(v => v < max); drawLevels(lv, ink, 0.8); drawLevels(lv.map(v => -v), ink, 0.8, [3, 3]); }
  else { const lv = [0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32]; drawLevels(lv, ink, 0.8); drawLevels(lv.map(v => -v), ink, 0.8, [3, 3]); }
  ctx.globalAlpha = 0.9;
  if (signedView) drawLevels([0], ink, 1.3);
  ctx.globalAlpha = 1;
  // nuclei
  const Xz = (zA) => (zA + zhalf) / (2 * zhalf) * w;
  const y0 = h / 2;
  ctx.font = `600 13px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
  [[-p.R / 2, 0], [p.R / 2, 1]].forEach(([z, k]) => {
    ctx.fillStyle = cssVar('--nucleus');
    ctx.beginPath(); ctx.arc(Xz(z), y0, 3.5, 0, 2 * Math.PI); ctx.fill();
    ctx.fillStyle = cssVar('--text');
    ctx.fillText(atomSym(k), Xz(z), y0 - 7);
  });
  // righello della distanza
  ctx.strokeStyle = cssVar('--muted'); ctx.fillStyle = cssVar('--muted'); ctx.lineWidth = 1;
  const yr = h - 22;
  ctx.beginPath(); ctx.moveTo(Xz(-p.R / 2), yr); ctx.lineTo(Xz(p.R / 2), yr); ctx.moveTo(Xz(-p.R / 2), yr - 4); ctx.lineTo(Xz(-p.R / 2), yr + 4); ctx.moveTo(Xz(p.R / 2), yr - 4); ctx.lineTo(Xz(p.R / 2), yr + 4); ctx.stroke();
  ctx.font = `11px ${MONO}`; ctx.textBaseline = 'top';
  ctx.fillText(`${nf(p.R, 2)} Å`, (Xz(-p.R / 2) + Xz(p.R / 2)) / 2, yr + 4);
  ctx.textAlign = 'right';
  ctx.fillText(`finestra ${nf(2 * zhalf, 1)} × ${nf(2 * xhalf, 1)} Å`, w - 10, h - 16);
  // legenda della scala
  ctx.textAlign = 'left';
  const unit = view === 'mo' ? `|ψ| max ${nf(max, 3)} bohr⁻³ᐟ²` : view === 'density' ? `ρ max ${nf(max, 1)} e/bohr³ (logaritmica)` : `scala ±${nf(vmax, 2)} e/bohr³ (max |Δ| ${nf(max, 3)})`;
  ctx.fillText(unit, 10, h - 16);
}

function moName(p) {
  const s = selectedMO(p);
  if (!s) return '';
  const o = s.mo;
  return `${o.label}${o.star ? '*' : ''}${p.mosB ? ` (spin ${s.spin === 'b' ? 'β' : 'α'})` : ''}, ε = ${nf(o.e * HARTREE_EV, 2)} eV, ${o.character}`;
}

// ---------------------------------------------------------------------------
// Grafici
// ---------------------------------------------------------------------------

function drawCharts() {
  if (!active) return;
  drawEnergyChart();
  drawDiagram();
  drawChart3();
}

let energyAxis = null;

function drawEnergyChart() {
  const canvas = $('chart-radial');
  $('radial-title').textContent = 'Energia mentre gli atomi si avvicinano';
  if (!F.grid.length) { drawXY(canvas, { series: [], xmin: 0, xmax: 1, ymin: 0, ymax: 1 }); return; }
  const ex = experimentalWell(F.mol);
  const morse = morseCurve(F.mol);
  const series = [];
  const Rmin = Math.min(...F.grid), Rmax = Math.max(...F.grid);
  const deep = Math.max(ex.De, 0.002);
  for (const m of Object.keys(F.data)) {
    if (m === 'gfn2' && !F.show.gfn2) continue;
    const c = curve(m), Ea = methodEatoms(m);
    if (!c.R.length || Ea == null) continue;
    series.push({ xs: c.R, ys: c.E.map(e => (e - Ea) * HARTREE_EV), color: METHOD_COLOR[m], width: m === F.method ? 2.6 : 1.6, label: METHOD_INFO[m].short });
  }
  if (morse) {
    const xs = [], ys = [];
    for (let k = 0; k <= 160; k++) { const R = Rmin + (Rmax - Rmin) * k / 160; xs.push(R); ys.push(morse.V(R)); }
    series.push({ xs, ys, color: cssVar('--text'), width: 1.4, dash: [6, 4], label: 'esperimento (Morse)' });
  }
  const isHe = F.mol.id === 'He2';
  const ymin = -1.35 * deep, ymax = isHe ? 0.004 : Math.max(0.9 * deep, 1.5);
  const R = F.grid[F.idx];
  const points = [{ x: ex.re, y: -ex.De, color: cssVar('--text'), r: 3.5 }];
  const p = currentPoint();
  if (p) points.push({ x: p.R, y: (p.E - F.ctx.Eatoms) * HARTREE_EV, color: METHOD_COLOR[F.method], r: 5 });
  energyAxis = { xmin: Rmin, xmax: Rmax, padLeft: 52 };
  drawXY(canvas, {
    series, xmin: Rmin, xmax: Rmax, ymin, ymax, xlabel: 'R (Å)', ylabel: 'E − E(atomi) (eV)',
    hlines: [{ y: 0, label: 'atomi separati', dash: [2, 3] }], vlines: [{ x: R, color: METHOD_COLOR[F.method], dash: [] }], points,
  });
  const c = curve(F.method);
  const sum = c.R.length > 6 ? curveSummary(F.mol, c.R, c.E, F.ctx?.Eatoms) : null;
  const cf = coulsonFischer();
  $('radial-note').innerHTML = `Punti calcolati a ${F.grid.length} distanze; la curva tratteggiata è il potenziale di Morse costruito con D<sub>e</sub>, r<sub>e</sub> e ω<sub>e</sub> sperimentali. ${sum?.bound ? `${METHOD_INFO[F.method].short}: r<sub>e</sub> = ${nf(sum.re, 3)} Å, D<sub>e</sub> = ${nf(sum.De, 2)} eV (esperimento ${nf(ex.re, 3)} Å, ${nf(ex.De, 2)} eV).` : ''}${cf ? ` Punto di Coulson–Fischer: sotto ${nf(cf, 2)} Å la soluzione UHF coincide con RHF.` : ''} Clic sul grafico per scegliere la distanza.`;
}

/** Distanza oltre la quale la soluzione UHF di singoletto rompe la simmetria di spin. */
function coulsonFischer() {
  if (!F.data.uhf || !F.data.rhf) return null;
  const pts = curve('uhf').pts;
  let last = null;
  for (const p of pts) { if (p.S2 != null && p.S2 < 0.01) last = p.R; else if (last != null) break; }
  return last != null && pts.some(p => p.S2 > 0.05) ? last : null;
}

/** Diagramma degli orbitali molecolari alla distanza corrente, con gli orbitali atomici ai lati. */
function drawDiagram() {
  const canvas = $('chart-levels');
  $('levels-title').textContent = 'Dagli orbitali atomici agli orbitali molecolari';
  const { ctx, w, h } = setupCanvas(canvas);
  const p = currentPoint();
  if (!p || !F.ctx) return;
  const ink = cssVar('--text'), muted = cssVar('--muted'), acc = cssVar('--accent');
  // orbitali atomici di valenza (niente core profondo) e orbitali molecolari nella stessa finestra di energia
  const atomLevels = (at) => at.orbitals.filter(o => !o.extra);
  const lists = [p.mos, ...(p.mosB ? [p.mosB] : [])];
  // si lasciano fuori i gusci interni (1s di Li…Ne, 1s–2p di Cl), che non partecipano al legame
  const coreTop = Math.max(-Infinity, ...[F.ctx.A, F.ctx.B].flatMap(at => at.orbitals.filter(o => o.core).map(o => o.e)));
  const valE = [...atomLevels(F.ctx.A), ...atomLevels(F.ctx.B)].filter(o => !o.core).map(o => o.e);
  const homoE = Math.max(...lists.flat().filter(o => o.occ > 0.5).map(o => o.e));
  const lowMO = lists.flat().map(o => o.e).filter(e => e > coreTop + 0.3);
  const coreCut = Math.min(...valE, ...lowMO) - 0.08;
  const virt = lists.flat().filter(o => o.occ < 0.5 && o.e > homoE).map(o => o.e).sort((a, b) => a - b);
  const top = Math.max(homoE + 0.3, Math.min(virt[1] ?? virt[0] ?? homoE, 0.8), ...valE.filter(e => e < 0.6)) + 0.08;
  const inWin = (e) => e >= coreCut && e <= top;
  const pad = { t: 16, b: 22, l: 40, r: 10 };
  const Y = (e) => pad.t + (top - e) / (top - coreCut) * (h - pad.t - pad.b);
  const colA = pad.l + (w - pad.l - pad.r) * 0.12, colB = pad.l + (w - pad.l - pad.r) * 0.88, colM = pad.l + (w - pad.l - pad.r) * 0.5;
  const lw = Math.min(46, (w - pad.l - pad.r) * 0.1);
  // asse
  ctx.font = `10.5px ${MONO}`; ctx.fillStyle = muted; ctx.strokeStyle = cssVar('--chart-grid'); ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  const stepEv = (top - coreCut) * HARTREE_EV > 30 ? 10 : 5;
  for (let v = Math.ceil(coreCut * HARTREE_EV / stepEv) * stepEv; v <= top * HARTREE_EV; v += stepEv) {
    const y = Y(v / HARTREE_EV);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
    ctx.fillText(`${v}`, pad.l - 4, y);
  }
  ctx.save(); ctx.translate(10, h / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.fillText('ε (eV)', 0, 0); ctx.restore();
  // livelli: gruppi degeneri disposti in orizzontale
  const groups = (items, keyE) => {
    const out = [];
    for (const it of items) {
      const g = out.find(x => Math.abs(x.e - keyE(it)) < 2e-4);
      if (g) g.items.push(it); else out.push({ e: keyE(it), items: [it] });
    }
    return out;
  };
  const placeLevel = (cx, width, g, draw) => {
    const k = g.items.length, gap = 4, each = (width - gap * (k - 1)) / k;
    g.items.forEach((it, i) => draw(it, cx - width / 2 + i * (each + gap), each));
  };
  const arrows = (x, y, wid, occ, spin) => {
    ctx.fillStyle = ink; ctx.font = `12px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    const t = spin === 'a' ? '↑' : spin === 'b' ? '↓' : occ > 1.5 ? '↑↓' : occ > 0.5 ? '↑' : '';
    if (t) ctx.fillText(t, x + wid / 2, y - 1);
  };
  const aoPos = new Map();
  [[F.ctx.A, colA, 0], [F.ctx.B, colB, 1]].forEach(([at, cx, k]) => {
    for (const g of groups(atomLevels(at).filter(o => inWin(o.e)), o => o.e)) {
      placeLevel(cx, lw * 1.5, g, (o, x, wid) => {
        const y = Y(o.e);
        ctx.strokeStyle = acc; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + wid, y); ctx.stroke();
        aoPos.set(`${k}:${o.label}`, { x: x + wid / 2, y, side: k });
        // elettroni dell'atomo libero (media sferica: frazioni distribuite sui p)
        const occ = Math.round(o.occ * 2) / 2;
        ctx.fillStyle = ink; ctx.font = `11px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
        ctx.fillText(occ >= 1.75 ? '↑↓' : occ >= 0.75 ? '↑' : occ > 0.2 ? '½' : '', x + wid / 2, y - 1);
      });
      ctx.fillStyle = muted; ctx.font = `10.5px ${MONO}`; ctx.textBaseline = 'middle';
      ctx.textAlign = k ? 'left' : 'right';
      ctx.fillText(g.items[0].label.replace(/[xyz]$/, ''), k ? cx + lw * 0.8 : cx - lw * 0.8, Y(g.e));
    }
    ctx.fillStyle = ink; ctx.font = `600 12px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    ctx.fillText(at.symbol, cx, h - pad.b + 4);
  });
  ctx.fillStyle = ink; ctx.font = `600 12px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  ctx.fillText(`${F.mol.name} a ${nf(p.R, 2)} Å`, colM, h - pad.b + 4);
  // orbitali molecolari (per UHF: α a sinistra, β a destra della colonna centrale)
  const sel = selectedMO(p);
  F.diagramHits = [];
  lists.forEach((list, li) => {
    const spin = p.mosB ? (li ? 'b' : 'a') : null;
    const cx = p.mosB ? colM + (li ? 1 : -1) * lw * 0.95 : colM;
    for (const g of groups(list.filter(o => inWin(o.e)), o => o.e)) {
      placeLevel(cx, lw * (p.mosB ? 1.6 : 2), g, (o, x, wid) => {
        const y = Y(o.e);
        // linee verso gli orbitali atomici che contribuiscono
        for (const c of o.comp) {
          if (Math.abs(c.w) < 0.08) continue;
          const a = aoPos.get(`${c.atom}:${c.label}`);
          if (!a) continue;
          ctx.strokeStyle = muted; ctx.globalAlpha = Math.min(0.85, 0.15 + Math.abs(c.w)); ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
          ctx.beginPath(); ctx.moveTo(a.side ? a.x - lw * 0.75 : a.x + lw * 0.75, a.y); ctx.lineTo(a.side ? x + wid : x, y); ctx.stroke();
          ctx.setLineDash([]); ctx.globalAlpha = 1;
        }
        const isSel = sel && sel.mo.label === o.label && sel.mo.axis === o.axis && (sel.spin === (spin ?? 'a'));
        ctx.strokeStyle = o.character === 'legante' ? cssVar('--phase-pos') : o.character === 'antilegante' ? cssVar('--phase-neg') : ink;
        ctx.lineWidth = isSel ? 4 : 2.2;
        ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + wid, y); ctx.stroke();
        if (F.method === 'fci') { ctx.fillStyle = ink; ctx.font = `10px ${MONO}`; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; if (o.occ > 0.002) ctx.fillText(nf(o.occ, 2), x + wid / 2, y - 1); }
        else arrows(x, y, wid, o.occ, o.occ > 0.5 ? spin : null);
        F.diagramHits.push({ x0: x - 2, x1: x + wid + 2, y, label: o.label, axis: o.axis, spin: spin ?? 'a' });
      });
      ctx.fillStyle = muted; ctx.font = `10.5px ${MONO}`; ctx.textBaseline = 'middle';
      const o = g.items[0];
      if (!p.mosB || li === 1) { ctx.textAlign = 'left'; ctx.fillText(`${o.label}${o.star ? '*' : ''}`, cx + lw * (p.mosB ? 0.85 : 1.05), Y(g.e)); }
    }
  });
  $('levels-note').innerHTML = `Ai lati gli orbitali degli atomi liberi (calcolati, media sferica), al centro quelli della molecola a questa distanza: rosso legante, blu antilegante, nero non legante. Le linee tratteggiate collegano ogni orbitale molecolare agli orbitali atomici che lo compongono (più scure se pesano di più).${p.mosB ? ' Con UHF gli elettroni α (↑) e β (↓) hanno orbitali diversi.' : ''}${F.method === 'fci' ? ' Per FCI i numeri sono le occupazioni medie: la correlazione porta una frazione di elettrone negli orbitali antileganti.' : ''} Clic su un livello per vederlo nella mappa.`;
}

function drawChart3() {
  const canvas = $('chart-slice');
  const sel = $('fm-chart3');
  if (sel && sel.value !== F.chart3) sel.value = F.chart3;
  if (F.chart3 === 'eps') return drawEpsChart(canvas);
  $('slice-title').firstChild.textContent = 'Energia cinetica e potenziale (Ruedenberg)';
  const v = virialCurve();
  if (!v || v.R.length < 4 || !F.ctx) { drawXY(canvas, { series: [], xmin: 0, xmax: 1, ymin: 0, ymax: 1 }); return; }
  const series = [
    { xs: v.R, ys: v.dT, color: cssVar('--phase-neg'), width: 2.2, label: 'ΔT cinetica' },
    { xs: v.R, ys: v.dV, color: cssVar('--phase-pos'), width: 2.2, label: 'ΔV potenziale' },
    { xs: v.R, ys: v.dE, color: cssVar('--text'), width: 2.2, label: 'ΔE = ΔT + ΔV' },
  ];
  // stesso calcolo dalla curva sperimentale con il teorema del viriale: T = −E − R dE/dR, V = 2E + R dE/dR
  const morse = morseCurve(F.mol);
  const ex = experimentalWell(F.mol);
  if (morse) {
    const xs = [], ts = [], vs = [];
    const Rmin = Math.min(...F.grid), Rmax = Math.max(...F.grid);
    for (let k = 0; k <= 120; k++) {
      const R = Rmin + (Rmax - Rmin) * k / 120, h = 1e-4;
      const E = morse.V(R), dEdR = (morse.V(R + h) - morse.V(R - h)) / (2 * h);
      xs.push(R); ts.push(-E - R * dEdR); vs.push(2 * E + R * dEdR);
    }
    series.push({ xs, ys: ts, color: cssVar('--phase-neg'), width: 1.2, dash: [5, 4], label: '' }, { xs, ys: vs, color: cssVar('--phase-pos'), width: 1.2, dash: [5, 4], label: '' });
  }
  const lim = 3.2 * Math.max(ex.De, 0.5);
  drawXY(canvas, {
    series, xmin: Math.min(...F.grid), xmax: Math.max(...F.grid), ymin: -lim, ymax: lim, xlabel: 'R (Å)', ylabel: 'rispetto agli atomi (eV)',
    hlines: [{ y: 0, dash: [2, 3] }], vlines: [{ x: F.grid[F.idx], color: METHOD_COLOR[F.method], dash: [] }],
  });
  const p = currentPoint();
  const wf = p ? `−⟨V⟩/⟨T⟩ = ${nf(-p.V / p.T, 4)}` : '';
  $('slice-note').innerHTML = `Linee piene: dal calcolo ${METHOD_INFO[F.method].short}; tratteggiate: dalla curva sperimentale. In entrambi i casi T e V si ottengono dalla curva E(R) con il teorema del viriale molecolare, ΔT = −ΔE − R dE/dR e ΔV = 2ΔE + R dE/dR (Slater, J. Chem. Phys. 1, 687, 1933), esatto per la soluzione esatta: con una base finita il viriale della funzione d'onda sbaglia di qualche millesimo (qui ${wf}), un errore sull'energia cinetica totale più grande delle variazioni dovute al legame. Avvicinandosi l'energia cinetica prima cala (gli elettroni si delocalizzano su due nuclei) e poi cresce molto, mentre la potenziale crolla (gli orbitali si contraggono): nel minimo ΔT = −ΔE > 0 e ΔV = 2ΔE (Ruedenberg, Rev. Mod. Phys. 34, 326, 1962; Bacskay e Nordholm, J. Phys. Chem. A 117, 7946, 2013).`;
}

function drawEpsChart(canvas) {
  $('slice-title').firstChild.textContent = 'Energie degli orbitali al variare di R';
  const c = curve(F.method);
  if (c.R.length < 2) { drawXY(canvas, { series: [], xmin: 0, xmax: 1, ymin: 0, ymax: 1 }); return; }
  const p = currentPoint() ?? c.pts[c.pts.length - 1];
  const lists = p.mosB ? ['a', 'b'] : ['a'];
  const series = [];
  let lo = Infinity, hi = -Infinity;
  for (const spin of lists) {
    const ref = (spin === 'b' ? p.mosB : p.mos);
    const homo = ref.filter(o => o.occ > 0.5).length;
    const keep = ref.filter((o, k) => k >= Math.max(0, homo - 5) && k < homo + 2 && (o.axis === '' || o.axis === 'x'));
    for (const o of keep) {
      const ys = c.pts.map(q => { const l = spin === 'b' ? q.mosB : q.mos; const m = l?.find(x => x.label === o.label && x.axis === o.axis); return m ? m.e * HARTREE_EV : NaN; });
      ys.forEach(v => { if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } });
      series.push({ xs: c.R, ys, color: o.character === 'legante' ? cssVar('--phase-pos') : o.character === 'antilegante' ? cssVar('--phase-neg') : cssVar('--muted'), width: o.occ > 0.5 ? 2.2 : 1.2, dash: spin === 'b' ? [4, 3] : [], label: '' });
    }
  }
  lo = Math.max(lo, -45); hi = Math.min(hi, 15);
  const padY = (hi - lo) * 0.06 || 1;
  drawXY(canvas, { series, xmin: Math.min(...F.grid), xmax: Math.max(...F.grid), ymin: lo - padY, ymax: hi + padY, xlabel: 'R (Å)', ylabel: 'ε (eV)', legend: false, vlines: [{ x: F.grid[F.idx], color: METHOD_COLOR[F.method], dash: [] }] });
  $('slice-note').innerHTML = 'Diagramma di correlazione: come ogni orbitale molecolare (linee spesse se occupato) cambia energia mentre gli atomi si avvicinano. A grande distanza i livelli coincidono con quelli degli atomi; vicino all\'equilibrio i leganti (rossi) scendono e gli antileganti (blu) salgono. I nomi sono quelli del punto scelto.';
}

// ---------------------------------------------------------------------------
// Testo di approfondimento con le fonti
// ---------------------------------------------------------------------------

function renderAnalysis() {
  $('analysis').innerHTML = `
    <h2>Che cosa succede quando nasce un legame</h2>
    <div class="analysis-grid">
      <article><h3>Il calcolo, distanza per distanza</h3>
        <p>Per ogni distanza R fra i nuclei, tenuti fermi (approssimazione di Born e Oppenheimer, 1927), si risolve l'equazione di Schrödinger per gli elettroni. La funzione d'onda è scritta come combinazione di funzioni gaussiane centrate sugli atomi (basi di Pople: Hehre, Ditchfield e Pople, J. Chem. Phys. 56, 2257, 1972; Hariharan e Pople, Theor. Chim. Acta 28, 213, 1973). Gli integrali sono calcolati esattamente e il campo autoconsistente di Hartree–Fock (Roothaan, Rev. Mod. Phys. 23, 69, 1951) è accelerato con DIIS (Pulay, 1980). Ogni distanza parte dalla soluzione della distanza vicina, per seguire con continuità lo stesso stato elettronico.</p>
        <p>Per H₂ si calcola anche l'interazione di configurazioni completa: tutte le configurazioni dei due elettroni negli orbitali della base, cioè la soluzione esatta dell'equazione di Schrödinger in quella base (Szabo e Ostlund, <i>Modern Quantum Chemistry</i>, cap. 4). Per H₂⁺, con un solo elettrone, Hartree–Fock è già esatto.</p>
      </article>
      <article><h3>Orbitali che nascono da orbitali</h3>
        <p>Gli orbitali molecolari sono combinazioni lineari degli orbitali atomici (LCAO: Lennard-Jones, Trans. Faraday Soc. 25, 668, 1929; Mulliken, Phys. Rev. 32, 186, 1928). Qui la composizione non è disegnata a mano: ogni orbitale molecolare calcolato è riscritto esattamente nella base degli orbitali dei due atomi liberi, e il peso di ciascuno è la sua popolazione di Mulliken (J. Chem. Phys. 23, 1833, 1955). Il carattere legante o antilegante viene dal segno della popolazione di sovrapposizione fra i due atomi; l'ordine di legame è quello di Mayer (Chem. Phys. Lett. 97, 270, 1983).</p>
        <p>Si vede così il mescolamento s–p: in N₂ e CO l'orbitale 3σ<sub>g</sub> (o 5σ) ha un forte contributo 2s, che lo spinge sopra i π, mentre in O₂ e F₂ i livelli tornano nell'ordine "da manuale".</p>
      </article>
      <article><h3>Spin che si accoppiano</h3>
        <p>Heitler e London (Z. Phys. 44, 455, 1927) spiegarono H₂ con due elettroni, uno per atomo, a spin opposti. Hartree–Fock ristretto mette invece entrambi gli elettroni nello stesso orbitale anche quando gli atomi sono lontani: metà della funzione d'onda è ionica (H⁺ H⁻) e la dissociazione è sbagliata di molti eV. La soluzione non ristretta (Pople e Nesbet, J. Chem. Phys. 22, 571, 1954) lascia separare gli spin: sotto una distanza critica, il punto di Coulson–Fischer (Phil. Mag. 40, 386, 1949), coincide con quella ristretta; sopra, gli elettroni α e β si localizzano su atomi diversi. La vista "spin" mostra questo disaccoppiamento; FCI lo descrive senza rompere la simmetria.</p>
      </article>
      <article><h3>Perché l'energia scende: cinetica e potenziale</h3>
        <p>Nel minimo il teorema del viriale impone ΔT = −ΔE > 0 e ΔV = 2ΔE: l'energia cinetica degli elettroni è più alta nella molecola che negli atomi. Ruedenberg (Rev. Mod. Phys. 34, 326, 1962) mostrò che il legame nasce comunque da un abbassamento dell'energia cinetica, dovuto alla delocalizzazione degli elettroni su due nuclei, a cui segue la contrazione degli orbitali che abbassa il potenziale (vedi anche Bacskay e Nordholm, J. Phys. Chem. A 117, 7946, 2013). Il terzo grafico ricava ΔT e ΔV con il viriale molecolare (Slater, J. Chem. Phys. 1, 687, 1933) sia dalla curva calcolata sia da quella sperimentale.</p>
      </article>
      <article><h3>Dati sperimentali</h3>
        <p class="desc-muted">r<sub>e</sub>, ω<sub>e</sub>, ω<sub>e</sub>x<sub>e</sub>: K. P. Huber e G. Herzberg, <i>Molecular Spectra and Molecular Structure IV. Constants of Diatomic Molecules</i> (1979), dal NIST Chemistry WebBook (SRD 69). Energie di dissociazione: Huber e Herzberg; H₂ da Liu et al., J. Chem. Phys. 130, 174306 (2009); Li₂ da Le Roy et al., J. Chem. Phys. 131, 204309 (2009); HF da di Lonardo e Douglas (1973); Cl₂ da Le Roy (1973); He₂ da Aziz, Janzen e Moldover, Phys. Rev. Lett. 74, 1586 (1995). La curva "sperimentale" è il potenziale di Morse (Phys. Rev. 34, 57, 1929) con quelle costanti: esatto vicino al minimo e nel limite di dissociazione, approssimato in mezzo.</p>
      </article>
      <article><h3>Limiti</h3>
        <p class="desc-muted">Nuclei fermi (niente moto di punto zero nel calcolo: per questo si confronta con D<sub>e</sub>, non con D<sub>0</sub>); nessun effetto relativistico; basi finite, che danno anche un piccolo errore di sovrapposizione (in He₂ un minimo spurio di pochi μeV). Hartree–Fock non contiene la correlazione elettronica: sottostima D<sub>e</sub> di circa 1–4 eV e non lega F₂. UHF dissocia bene ma non è un autostato dello spin. Per O₂ un solo determinante con M<sub>S</sub> = 1 non può descrivere due atomi O(³P) accoppiati: la curva UHF termina sopra gli atomi separati. GFN2-xTB (Bannwarth, Ehlert e Grimme, J. Chem. Theory Comput. 15, 1652, 2019) è un metodo semiempirico parametrizzato sulle geometrie di equilibrio: buono per r<sub>e</sub> e ω<sub>e</sub>, non per la dissociazione; in He₂ ritrova invece la buca di van der Waals grazie alla dispersione D4.</p>
      </article>
    </div>`;
}

// ---------------------------------------------------------------------------
// Interazione con i grafici
// ---------------------------------------------------------------------------

export function initFormationMode() {
  $('chart-radial').addEventListener('click', (e) => {
    if (!active || !energyAxis || !F.grid.length) return;
    const r = e.currentTarget.getBoundingClientRect();
    const pw = r.width - energyAxis.padLeft - 14;
    const R = energyAxis.xmin + (e.clientX - r.left - energyAxis.padLeft) / pw * (energyAxis.xmax - energyAxis.xmin);
    stop();
    setIndex(nearestIndex(R));
    renderControls();
  });
  $('chart-levels').addEventListener('click', (e) => {
    if (!active || !F.diagramHits) return;
    const r = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - r.left, y = e.clientY - r.top;
    let best = null, bd = 9;
    for (const hgt of F.diagramHits) if (x >= hgt.x0 && x <= hgt.x1 && Math.abs(y - hgt.y) < bd) { bd = Math.abs(y - hgt.y); best = hgt; }
    if (!best) return;
    F.sel = { label: best.label, axis: best.axis, spin: best.spin };
    F.prevC = null;
    F.view = 'mo';
    renderControls();
    renderDynamic();
  });
  $('fm-chart3')?.addEventListener('change', (e) => { F.chart3 = e.target.value; drawChart3(); });
  $('fm-gfn2')?.addEventListener('change', (e) => { F.show.gfn2 = e.target.checked; drawEnergyChart(); });
}

