// Modalità "Reazioni": termochimica calcolata dai principi primi.
//   ΔE = Σ ν E(prodotti) − Σ ν E(reagenti)        (energia elettronica, HF o MP2)
//   ΔH(T) = ΔE + Δ[ZPE + E_termica + RT]           (termodinamica statistica)
//   ΔG(T) = ΔH − TΔS,   K = exp(−ΔG°/RT)

import { MOLECULES, REACTIONS, BARRIERS, moleculeById } from '../chem/library.js';
import { LIBRARY_DATA } from '../chem/libraryData.js';
import { parseSmiles, hillFormula } from '../chem/smiles.js';
import { analyzeStructure, connectivityFromGeometry } from '../chem/structure.js';
import { thermochemistry } from '../chem/vibrations.js';
import { buildMolecule } from '../render/moleculeView.js';
import { textSprite } from '../render/viewer.js';
import { drawEnthalpy, drawBars, drawProfile } from './chemCharts.js';
import { drawLineChart } from './charts.js';
import { FREQ_SCALE } from './moleculeMode.js';

const BOHR_ANG = 0.52917721090;
const HARTREE_KJ = 2625.4996394799;
const R_GAS = 8.314462618e-3; // kJ/(mol K)
const $ = (id) => document.getElementById(id);
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const nf = (v, d = 1) => (v === null || v === undefined || !Number.isFinite(v) ? '—' : Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d }));

let viewer = null;
let active = false;
const R = {
  reaction: REACTIONS[0],
  custom: { left: [[1, 'H2'], [1, 'Cl2']], right: [[2, 'HCl']] },
  isCustom: false,
  kind: 'thermo',          // 'thermo' | 'barrier'
  barrier: BARRIERS[0],
  method: 'mp2',
  scaled: true,
  T: 298.15,
};

export function initReactionMode(v) {
  viewer = v;
}

export function activateReaction() {
  active = true;
  document.body.dataset.mode = 'reaction';
  renderAll();
}

export function deactivateReaction() {
  active = false;
}

export function reactionRedraw() {
  if (active) { render3D(); drawCharts(); }
}

const formulaCache = new Map();
function formula(id) {
  if (id.startsWith('ts-')) return '‡';
  if (!formulaCache.has(id)) {
    const m = moleculeById(id);
    const g = parseSmiles(m.smiles);
    formulaCache.set(id, m.formula ?? hillFormula(g.atoms, g.atoms.reduce((s, a) => s + a.charge, 0)));
  }
  return formulaCache.get(id);
}

function sideText(side) {
  return side.map(([c, id]) => `${c > 1 ? `${c} ` : ''}${formula(id)}`).join(' + ');
}

function current() {
  return R.isCustom ? { name: 'Reazione personalizzata', ...R.custom, dH: null } : R.reaction;
}

/** Proprietà termodinamiche di una specie alla temperatura T (kJ/mol, J/(mol K)). */
function species(id, T) {
  const d = LIBRARY_DATA[id];
  if (!d) return null;
  const atoms = d.geometry.map(([Z, x, y, z]) => ({ Z, xyz: [x / BOHR_ANG, y / BOHR_ANG, z / BOHR_ANG] }));
  const scale = R.scaled ? (FREQ_SCALE[d.basis] ?? 0.9) : 1;
  const modes = d.freqs.map(f => ({ freq: f.freq * scale }));
  const th = thermochemistry(atoms, modes, d.energy, { T, multiplicity: d.multiplicity });
  const Eel = R.method === 'mp2' ? d.energy + (d.mp2 ?? 0) : d.energy;
  return {
    id, basis: d.basis, E: Eel * HARTREE_KJ, Ehf: d.energy, Emp2: d.energy + (d.mp2 ?? 0),
    zpe: th.zpe, Hcorr: th.Hcorr * HARTREE_KJ, Gcorr: th.Gcorr * HARTREE_KJ, S: th.S,
    pointGroup: d.pointGroup, sigma: th.sigma, freqs: d.freqs.length,
  };
}

/** Bilancio della reazione: atomi e carica devono conservarsi. */
function balance(rx) {
  const count = new Map();
  let charge = 0;
  const add = (side, sign) => side.forEach(([c, id]) => {
    const g = parseSmiles(moleculeById(id).smiles);
    g.atoms.forEach(a => count.set(a.Z, (count.get(a.Z) ?? 0) + sign * c));
    charge += sign * c * g.atoms.reduce((s, a) => s + a.charge, 0);
  });
  add(rx.left, -1);
  add(rx.right, 1);
  const bad = [...count.entries()].filter(([, v]) => v !== 0);
  return { ok: bad.length === 0 && charge === 0, bad, charge };
}

/** Grandezze di reazione alla temperatura T. */
function thermo(rx, T) {
  const sum = (side) => {
    const out = { E: 0, H: 0, G: 0, S: 0, zpe: 0, Ehf: 0, Emp2: 0, missing: [] };
    for (const [c, id] of side) {
      const s = species(id, T);
      if (!s) { out.missing.push(id); continue; }
      out.E += c * s.E; out.H += c * (s.E + s.Hcorr); out.G += c * (s.E + s.Gcorr); out.S += c * s.S; out.zpe += c * s.zpe;
      out.Ehf += c * s.Ehf * HARTREE_KJ; out.Emp2 += c * s.Emp2 * HARTREE_KJ;
    }
    return out;
  };
  const L = sum(rx.left), Rt = sum(rx.right);
  const missing = [...L.missing, ...Rt.missing];
  const dG = Rt.G - L.G;
  // variazione del numero di moli di gas (per Kc)
  const dn = rx.right.reduce((s, [c]) => s + c, 0) - rx.left.reduce((s, [c]) => s + c, 0);
  return {
    missing,
    dE: Rt.E - L.E, dEhf: Rt.Ehf - L.Ehf, dEmp2: Rt.Emp2 - L.Emp2,
    dZPE: Rt.zpe - L.zpe,
    dH: Rt.H - L.H, dS: Rt.S - L.S, dG,
    lnK: -dG / (R_GAS * T), dn,
  };
}

function renderAll() {
  if (!active) return;
  renderSide();
  renderControls();
  renderResults();
  render3D();
  drawCharts();
  renderAnalysis();
}

function speciesOptions(selected) {
  return MOLECULES.filter(m => LIBRARY_DATA[m.id]).map(m => `<option value="${m.id}" ${m.id === selected ? 'selected' : ''}>${formula(m.id)} · ${m.name}</option>`).join('');
}

function renderSide() {
  const list = REACTIONS.map(r => `<button type="button" class="rx-item ${!R.isCustom && R.reaction.id === r.id ? 'active' : ''}" data-id="${r.id}">
      <span class="rx-name">${r.name}</span><span class="rx-eq">${sideText(r.left)} → ${sideText(r.right)}</span></button>`).join('');
  const sideEditor = (side, key) => side.map(([c, id], i) => `
    <div class="rx-row">
      <input type="number" min="1" max="9" value="${c}" data-side="${key}" data-i="${i}" class="rx-coef" aria-label="Coefficiente">
      <select data-side="${key}" data-i="${i}" class="rx-sp" aria-label="Specie">${speciesOptions(id)}</select>
      <button type="button" class="linkish rx-del" data-side="${key}" data-i="${i}" aria-label="Rimuovi">✕</button>
    </div>`).join('');
  const bal = balance(R.custom);
  const blist = BARRIERS.map(b => `<button type="button" class="rx-item ${R.kind === 'barrier' && R.barrier.id === b.id ? 'active' : ''}" data-barrier="${b.id}">
      <span class="rx-name">${b.name}</span><span class="rx-eq">${sideText(b.left)} → [‡] → ${sideText(b.right)}</span></button>`).join('');
  $('element-card').innerHTML = `
    <h3 class="side-h">Termochimica: reazioni in fase gassosa</h3>
    <div class="rx-list">${list}</div>
    <h3 class="side-h">Cinetica: stati di transizione</h3>
    <div class="rx-list short">${blist}</div>
    <h3 class="side-h">Costruisci una reazione</h3>
    <div class="rx-editor">
      <p class="lbl">Reagenti</p>${sideEditor(R.custom.left, 'left')}
      <button type="button" class="linkish" id="rx-add-left">+ aggiungi reagente</button>
      <p class="lbl">Prodotti</p>${sideEditor(R.custom.right, 'right')}
      <button type="button" class="linkish" id="rx-add-right">+ aggiungi prodotto</button>
      <p class="hint">${bal.ok ? 'Reazione bilanciata.' : `Non bilanciata: ${bal.bad.map(([Z, v]) => `${ELEMENT_SYMBOL(Z)} ${v > 0 ? '+' : ''}${v}`).join(', ')}${bal.charge ? `, carica ${bal.charge}` : ''}.`}</p>
      <button type="button" class="btn" id="rx-use" ${bal.ok ? '' : 'disabled'}>Calcola questa reazione</button>
    </div>`;
  $('element-card').querySelectorAll('.rx-item[data-id]').forEach(b => b.addEventListener('click', () => {
    R.reaction = REACTIONS.find(r => r.id === b.dataset.id);
    R.isCustom = false;
    R.kind = 'thermo';
    renderAll();
  }));
  $('element-card').querySelectorAll('.rx-item[data-barrier]').forEach(b => b.addEventListener('click', () => {
    R.barrier = BARRIERS.find(x => x.id === b.dataset.barrier);
    R.kind = 'barrier';
    renderAll();
  }));
  const card = $('element-card');
  card.querySelectorAll('.rx-coef').forEach(inp => inp.addEventListener('change', () => {
    R.custom[inp.dataset.side][+inp.dataset.i][0] = Math.max(1, Math.min(9, +inp.value || 1));
    renderSide();
  }));
  card.querySelectorAll('.rx-sp').forEach(sel => sel.addEventListener('change', () => {
    R.custom[sel.dataset.side][+sel.dataset.i][1] = sel.value;
    renderSide();
  }));
  card.querySelectorAll('.rx-del').forEach(b => b.addEventListener('click', () => {
    const side = R.custom[b.dataset.side];
    if (side.length > 1) side.splice(+b.dataset.i, 1);
    renderSide();
  }));
  $('rx-add-left').addEventListener('click', () => { R.custom.left.push([1, 'H2']); renderSide(); });
  $('rx-add-right').addEventListener('click', () => { R.custom.right.push([1, 'H2O']); renderSide(); });
  $('rx-use').addEventListener('click', () => { R.isCustom = true; R.kind = 'thermo'; renderAll(); });
}

const ELEMENT_SYMBOL = (Z) => ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'][Z] ?? `Z${Z}`;

function renderControls() {
  $('controls').innerHTML = `
    <div class="ctl"><span class="lbl">Energia elettronica</span>
      <div class="seg" id="rx-method">
        <button type="button" data-v="hf" aria-pressed="${R.method === 'hf'}">Hartree–Fock</button>
        <button type="button" data-v="mp2" aria-pressed="${R.method === 'mp2'}">MP2 (con correlazione)</button>
      </div>
    </div>
    <div class="ctl"><span class="lbl">Temperatura</span>
      <div class="range-row"><input type="range" id="rx-T" min="100" max="2000" step="5" value="${R.T}"><output id="rx-T-out">${nf(R.T, 0)} K</output></div>
    </div>
    <label class="mini-toggle"><input type="checkbox" id="rx-scaled" ${R.scaled ? 'checked' : ''}> frequenze scalate (correzione empirica di Scott e Radom)</label>`;
  $('rx-method').querySelectorAll('button').forEach(b => b.addEventListener('click', () => { R.method = b.dataset.v; renderAll(); }));
  const T = $('rx-T');
  T.addEventListener('input', () => { $('rx-T-out').textContent = `${T.value} K`; });
  T.addEventListener('change', () => { R.T = +T.value; renderResults(); drawCharts(); renderAnalysis(); });
  $('rx-scaled').addEventListener('change', (e) => { R.scaled = e.target.checked; renderAll(); });
}

function renderResults() {
  if (R.kind === 'barrier') { renderBarrierResults(); return; }
  const rx = current();
  const t = thermo(rx, R.T);
  const body = $('panel-body');
  if (t.missing.length) {
    body.innerHTML = `<p class="desc">Dati non ancora calcolati per: ${t.missing.map(formula).join(', ')}.</p>`;
    return;
  }
  const K = t.lnK;
  const Ktxt = Math.abs(K) < 700 ? formatK(Math.exp(K)) : `e^${nf(K, 0)}`;
  const t298 = thermo(rx, 298.15);
  const err = rx.dH !== null && rx.dH !== undefined ? t298.dH - rx.dH : null;
  const Tinv = t.dS !== 0 ? t.dH / (t.dS / 1000) : null;
  body.innerHTML = `
    <div>
      <h3>${rx.name}</h3>
      <p class="rx-big">${sideText(rx.left)} → ${sideText(rx.right)}</p>
      <dl class="info-list">
        <dt>ΔE elettronica (${R.method === 'mp2' ? 'MP2' : 'HF'})</dt><dd>${nf(t.dE)} kJ/mol</dd>
        <dt>ΔZPE</dt><dd>${nf(t.dZPE)} kJ/mol</dd>
        <dt>ΔH° (${nf(R.T, 0)} K)</dt><dd><b>${nf(t.dH)}</b> kJ/mol</dd>
        <dt>ΔS°</dt><dd>${nf(t.dS)} J/(mol·K)</dd>
        <dt>ΔG° = ΔH° − TΔS°</dt><dd><b>${nf(t.dG)}</b> kJ/mol</dd>
        <dt>K<sub>p</sub> = e<sup>−ΔG°/RT</sup></dt><dd>${Ktxt}</dd>
        <dt>Spontanea (ΔG° &lt; 0)</dt><dd>${t.dG < 0 ? 'sì' : 'no'}</dd>
        ${Tinv && Tinv > 0 && Math.sign(t.dH) === Math.sign(t.dS) ? `<dt>Inversione (ΔG° = 0)</dt><dd>≈ ${nf(Tinv, 0)} K</dd>` : ''}
      </dl>
    </div>
    ${rx.dH !== null && rx.dH !== undefined ? `
    <div>
      <h3>Confronto con l'esperimento (298 K)</h3>
      <dl class="info-list">
        <dt>ΔH° calcolata</dt><dd>${nf(t298.dH)} kJ/mol</dd>
        <dt>ΔH° sperimentale</dt><dd>${nf(rx.dH)} kJ/mol</dd>
        <dt>Errore</dt><dd>${err >= 0 ? '+' : ''}${nf(err)} kJ/mol</dd>
      </dl>
      <p class="desc-muted" style="margin-top:6px">${R.method === 'hf' ? 'Hartree–Fock trascura la correlazione elettronica: gli errori su rotture e formazioni di legami possono arrivare a 100 kJ/mol. Prova MP2.' : 'MP2 recupera gran parte della correlazione; con basi piccole restano errori di qualche decina di kJ/mol, soprattutto per O₂ e i radicali.'}</p>
    </div>` : ''}`;
}

function formatK(K) {
  if (K === 0) return '0';
  const e = Math.floor(Math.log10(K));
  if (e >= -3 && e <= 4) return nf(K, 3);
  const m = K / Math.pow(10, e);
  const sup = String(e).split('').map(c => ({ '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' }[c])).join('');
  return `${nf(m, 2)}·10${sup}`;
}

function render3D() {
  if (R.kind === 'barrier') { renderBarrier3D(); return; }
  const rx = current();
  viewer.clear();
  viewer.setAxesVisible(false);
  // molecole in fila: reagenti → prodotti
  const items = [];
  const push = (side) => side.forEach(([c, id]) => items.push({ c, id }));
  push(rx.left);
  items.push({ arrow: true });
  push(rx.right);
  const geoms = items.map(it => {
    if (it.arrow) return { width: 3 };
    const d = LIBRARY_DATA[it.id];
    const g = parseSmiles(moleculeById(it.id).smiles);
    const atoms = d ? d.geometry.map(([Z, x, y, z]) => ({ Z, xyz: [x / BOHR_ANG, y / BOHR_ANG, z / BOHR_ANG] })) : [{ Z: 1, xyz: [0, 0, 0] }];
    const c = [0, 1, 2].map(k => atoms.reduce((s, a) => s + a.xyz[k], 0) / atoms.length);
    const centered = atoms.map(a => ({ Z: a.Z, xyz: a.xyz.map((v, k) => v - c[k]) }));
    const w = Math.max(2.5, ...centered.map(a => Math.abs(a.xyz[0]))) * 2 + 2;
    const st = analyzeStructure({ atoms: g.atoms, bonds: g.bonds });
    return { atoms: centered, bonds: g.bonds.map((b, k) => ({ a: b.a, b: b.b, order: st.bonds[k].averageOrder })), width: w, coef: it.c };
  });
  const total = geoms.reduce((s, g) => s + g.width, 0);
  let x = -total / 2;
  const color = cssVar('--text');
  geoms.forEach((g, i) => {
    const cx = x + g.width / 2;
    if (items[i].arrow) {
      const s = textSprite('→', color, 48);
      s.position.set(cx, 0, 0);
      s.scale.set(2 * s.userData.aspect, 2, 1);
      viewer.add(s);
    } else {
      const shifted = g.atoms.map(a => ({ Z: a.Z, xyz: [a.xyz[0] + cx, a.xyz[1], a.xyz[2]] }));
      const mol = buildMolecule(shifted, g.bonds, {});
      viewer.add(mol.group);
      if (g.coef > 1) {
        const s = textSprite(`${g.coef}×`, color, 40);
        s.position.set(cx - g.width / 2 + 0.6, 0, 2.6);
        s.scale.set(1.3 * s.userData.aspect, 1.3, 1);
        viewer.add(s);
      }
    }
    x += g.width;
  });
  viewer.frame(Math.max(6, total * 0.42), false);
  viewer.camera.position.set(0, -total * 1.05, total * 0.35);
  viewer.controls.target.set(0, 0, 0);
  $('viewport-title').innerHTML = `${rx.name}<small>${sideText(rx.left)} → ${sideText(rx.right)}</small>`;
  $('viewport-legend').innerHTML = '';
  $('viewport-note').textContent = 'Geometrie ottimizzate con Hartree–Fock. Ogni specie è calcolata isolata, in fase gassosa.';
}

function drawCharts() {
  if (R.kind === 'barrier') { drawBarrierCharts(); return; }
  const rx = current();
  const t = thermo(rx, R.T);
  if (t.missing.length) return;
  const t298 = thermo(rx, 298.15);
  const hfOnly = { ...R };
  // entalpia con entrambi i metodi a 298 K
  const saved = R.method;
  R.method = 'hf';
  const hf298 = thermo(rx, 298.15);
  R.method = 'mp2';
  const mp298 = thermo(rx, 298.15);
  R.method = saved;
  void hfOnly;
  drawEnthalpy($('chart-radial'), {
    left: sideText(rx.left), right: sideText(rx.right),
    values: [
      { label: 'HF', value: hf298.dH, color: cssVar('--phase-neg') },
      { label: 'MP2', value: mp298.dH, color: cssVar('--accent') },
      ...(rx.dH !== null && rx.dH !== undefined ? [{ label: 'sperim.', value: rx.dH, color: cssVar('--text'), dashed: true }] : []),
    ],
  });
  $('radial-title').textContent = 'Diagramma entalpico (298 K)';
  $('radial-note').textContent = `ΔH° ${t298.dH < 0 ? '< 0: reazione esotermica' : '> 0: reazione endotermica'}. Legge di Hess: l'entalpia è una funzione di stato, conta solo la differenza tra prodotti e reagenti.`;

  // ΔG(T) e ΔH(T)
  const Ts = [];
  for (let T = 100; T <= 2000; T += 25) Ts.push(T);
  const dG = Ts.map(T => thermo(rx, T).dG);
  const dH = Ts.map(T => thermo(rx, T).dH);
  const TdS = Ts.map((T, i) => dH[i] - dG[i]);
  drawLineChart($('chart-levels'), {
    series: [
      { xs: Ts, ys: dG, color: cssVar('--accent'), label: 'ΔG°(T)', width: 2.2, keepZero: true },
      { xs: Ts, ys: dH, color: cssVar('--phase-neg'), label: 'ΔH°(T)', dash: [5, 4], keepZero: true },
      { xs: Ts, ys: TdS.map(v => -v), color: cssVar('--muted'), label: '−TΔS°', dash: [2, 3], keepZero: true },
    ],
    log: false, signed: true, xlabel: 'T (K)', ylabel: 'kJ/mol',
    markers: [{ x: R.T, label: `${nf(R.T, 0)} K`, color: cssVar('--accent') }],
  });
  $('levels-title').textContent = 'Energia libera in funzione di T';
  $('levels-note').textContent = `Dove ΔG° < 0 la reazione è favorita (K > 1). ${t.dS < 0 ? 'ΔS° < 0: l\'aumento di T sfavorisce i prodotti (principio di Le Châtelier).' : 'ΔS° > 0: l\'aumento di T favorisce i prodotti.'}`;

  // contributi
  drawBars($('chart-slice'), {
    unit: 'kJ/mol',
    items: [
      { label: 'ΔE Hartree–Fock', value: t.dEhf },
      { label: 'Δ correlazione MP2', value: t.dEmp2 - t.dEhf },
      { label: 'ΔZPE', value: t.dZPE },
      { label: 'Δ(H − E − ZPE)', value: t.dH - t.dE - t.dZPE },
      { label: '−TΔS°', value: -R.T * t.dS / 1000 },
      { label: 'ΔG°', value: t.dG, color: cssVar('--accent') },
    ],
  });
  $('slice-title').textContent = `Contributi a ΔG° (${nf(R.T, 0)} K)`;
  $('slice-note').textContent = 'ΔG° = ΔE_el + ΔZPE + Δ(energia termica + RT) − TΔS°. Ogni termine viene da un pezzo diverso della fisica: struttura elettronica, vibrazioni, rotazioni e traslazioni.';
}

function renderAnalysis() {
  const rx = R.kind === 'barrier' ? { left: R.barrier.left, right: [...R.barrier.right, [1, `ts-${R.barrier.id}`]] } : current();
  const ids = [...new Set([...rx.left, ...rx.right].map(([, id]) => id))];
  const rows = ids.map(id => {
    const s = species(id, R.T);
    const name = id.startsWith('ts-') ? 'stato di transizione ‡' : formula(id);
    if (!s) return `<tr><td>${name}</td><td colspan="9">non calcolata</td></tr>`;
    return `<tr><td><b>${name}</b></td><td>${s.basis}</td><td class="num">${nf(s.Ehf, 6)}</td><td class="num">${nf(s.Emp2, 6)}</td>
      <td class="num">${nf(s.zpe, 2)}</td><td class="num">${nf(s.Hcorr, 2)}</td><td class="num">${nf(s.S, 2)}</td><td class="num">${nf(s.Gcorr, 2)}</td>
      <td>${s.pointGroup}</td><td class="num">${s.sigma}</td></tr>`;
  }).join('');
  $('analysis').innerHTML = `
    <h2>Dati delle specie (${nf(R.T, 0)} K, 1 atm)</h2>
    <div class="table-wrap"><table class="data-table">
      <thead><tr><th>Specie</th><th>Base</th><th class="num">E HF (Ha)</th><th class="num">E MP2 (Ha)</th><th class="num">ZPE (kJ/mol)</th><th class="num">H − E (kJ/mol)</th><th class="num">S° (J/mol·K)</th><th class="num">G − E (kJ/mol)</th><th>Gruppo</th><th class="num">σ</th></tr></thead>
      <tbody>${rows}</tbody>
    </table></div>
    <p class="desc-muted">Per ogni specie: geometria ottimizzata, frequenze armoniche dall'hessiana analitica-numerica, funzioni di partizione traslazionale q<sub>t</sub> = (2πmkT/h²)<sup>3/2</sup>kT/p, rotazionale q<sub>r</sub> = √π/σ · √(T³/Θ<sub>A</sub>Θ<sub>B</sub>Θ<sub>C</sub>), vibrazionale q<sub>v</sub> = Π 1/(1 − e<sup>−hν/kT</sup>). L'energia MP2 aggiunge la correlazione elettronica alla geometria HF.</p>`;
}

// ---------------------------------------------------------------------------
// Cinetica: teoria dello stato di transizione (Eyring, Evans e Polanyi, 1935)
//   k(T) = κ(T) · (k_B T / h) · (RT/p°)^(1−m) · exp(−ΔG‡°/RT)      (m = molecolarità)
//   κ(T) = 1 + (1/24)(hν‡/k_B T)²        correzione per l'effetto tunnel di Wigner
// ---------------------------------------------------------------------------

const KB = 1.380649e-23, HP = 6.62607015e-34, CL = 2.99792458e10;

function barrierData(T) {
  const b = R.barrier;
  const tsId = `ts-${b.id}`;
  const ts = species(tsId, T);
  const d = LIBRARY_DATA[tsId];
  if (!ts || !d) return null;
  const sum = (side) => side.reduce((acc, [c, id]) => {
    const s = species(id, T);
    if (!s) { acc.missing = true; return acc; }
    acc.E += c * s.E; acc.H += c * (s.E + s.Hcorr); acc.G += c * (s.E + s.Gcorr); acc.S += c * s.S; acc.zpe += c * s.zpe;
    acc.Ehf += c * s.Ehf * HARTREE_KJ; acc.Emp2 += c * s.Emp2 * HARTREE_KJ;
    return acc;
  }, { E: 0, H: 0, G: 0, S: 0, zpe: 0, Ehf: 0, Emp2: 0, missing: false });
  const L = sum(b.left), P = sum(b.right);
  if (L.missing || P.missing) return null;
  const molecularity = b.left.reduce((s, [c]) => s + c, 0);
  const nuImag = Math.abs(Math.min(...d.freqs.map(f => f.freq))) * (R.scaled ? (FREQ_SCALE[d.basis] ?? 0.9) : 1);
  const u = HP * CL * nuImag / (KB * T);
  const kappa = 1 + u * u / 24;
  const dG = ts.E + ts.Gcorr - L.G;
  let k = kappa * KB * T / HP * Math.exp(-dG * 1000 / (8.314462618 * T));
  if (molecularity === 2) k *= KB * T / 1e5 * 1e6; // → cm³ molecola⁻¹ s⁻¹ (stato standard 1 bar)
  return {
    dE: ts.E - L.E, dEhf: ts.Ehf * HARTREE_KJ - L.Ehf, dEmp2: ts.Emp2 * HARTREE_KJ - L.Emp2,
    dZPE: ts.zpe - L.zpe, dH: ts.E + ts.Hcorr - L.H, dS: ts.S - L.S, dG,
    dHr: P.H - L.H, dEr: P.E - L.E,
    k, kappa, nuImag, molecularity,
  };
}

function renderBarrierResults() {
  const b = R.barrier;
  const t = barrierData(R.T);
  const body = $('panel-body');
  if (!t) { body.innerHTML = '<p class="desc">Dati dello stato di transizione non ancora calcolati.</p>'; return; }
  const unit = t.molecularity === 2 ? 'cm³ molecola⁻¹ s⁻¹' : 's⁻¹';
  const half = t.molecularity === 1 ? Math.LN2 / t.k : null;
  const t298 = barrierData(298.15);
  body.innerHTML = `
    <div>
      <h3>${b.name}</h3>
      <p class="rx-big">${sideText(b.left)} → [‡] → ${sideText(b.right)}</p>
      <dl class="info-list">
        <dt>Frequenza immaginaria ν‡</dt><dd>${nf(t.nuImag, 0)}i cm⁻¹</dd>
        <dt>ΔE‡ elettronica (${R.method === 'mp2' ? 'MP2' : 'HF'})</dt><dd>${nf(t.dE)} kJ/mol</dd>
        <dt>ΔZPE‡</dt><dd>${nf(t.dZPE)} kJ/mol</dd>
        <dt>ΔH‡ (${nf(R.T, 0)} K)</dt><dd><b>${nf(t.dH)}</b> kJ/mol</dd>
        <dt>ΔS‡</dt><dd>${nf(t.dS)} J/(mol·K)</dd>
        <dt>ΔG‡</dt><dd><b>${nf(t.dG)}</b> kJ/mol</dd>
        <dt>Correzione tunnel κ (Wigner)</dt><dd>${nf(t.kappa, 3)}</dd>
        <dt>Costante di velocità k</dt><dd>${formatK(t.k)} ${unit}</dd>
        ${half !== null ? `<dt>Tempo di dimezzamento</dt><dd>${formatTime(half)}</dd>` : ''}
        <dt>ΔH di reazione</dt><dd>${nf(t.dHr)} kJ/mol</dd>
      </dl>
    </div>
    <div>
      <h3>Confronto con il riferimento</h3>
      <dl class="info-list">
        <dt>Barriera calcolata (ΔE‡ + ΔZPE‡)</dt><dd>${nf(t298.dE + t298.dZPE)} kJ/mol</dd>
        <dt>Barriera classica ΔE‡</dt><dd>${nf(t298.dE)} kJ/mol</dd>
        <dt>Riferimento</dt><dd>${nf(b.expBarrier)} kJ/mol</dd>
      </dl>
      <p class="mol-note" style="margin-top:10px">${b.note}</p>
    </div>`;
}

function formatTime(s) {
  if (s < 1e-9) return `${nf(s * 1e12, 2)} ps`;
  if (s < 1e-6) return `${nf(s * 1e9, 2)} ns`;
  if (s < 1e-3) return `${nf(s * 1e6, 2)} µs`;
  if (s < 1) return `${nf(s * 1e3, 2)} ms`;
  if (s < 3600) return `${nf(s, 2)} s`;
  if (s < 86400 * 365) return `${nf(s / 3600, 2)} h`;
  return `${formatK(s / (86400 * 365.25))} anni`;
}

function renderBarrier3D() {
  const b = R.barrier;
  const d = LIBRARY_DATA[`ts-${b.id}`];
  viewer.clear();
  viewer.setAxesVisible(false);
  if (!d) { $('viewport-title').innerHTML = `${b.name}<small>dati non disponibili</small>`; return; }
  const atoms = d.geometry.map(([Z, x, y, z]) => ({ Z, xyz: [x / BOHR_ANG, y / BOHR_ANG, z / BOHR_ANG] }));
  const c = [0, 1, 2].map(k => atoms.reduce((s, a) => s + a.xyz[k], 0) / atoms.length);
  atoms.forEach(a => { a.xyz = a.xyz.map((v, k) => v - c[k]); });
  // legami: pieni se corti, tratteggiati se allungati (legami che si formano o si rompono)
  const bonds = connectivityFromGeometry(atoms).map(bd => {
    const A = atoms[bd.a].xyz, B = atoms[bd.b].xyz;
    const r = Math.hypot(A[0] - B[0], A[1] - B[1], A[2] - B[2]) * BOHR_ANG * 100;
    const ref = (atoms[bd.a].Z === 1 ? 32 : 72) + (atoms[bd.b].Z === 1 ? 32 : 72);
    return { ...bd, order: r > 1.12 * ref ? 0.5 : 1 };
  });
  const mol = buildMolecule(atoms, bonds, {});
  viewer.add(mol.group);
  const imag = d.freqs.find(f => f.freq < 0);
  if (imag) {
    const x0 = atoms.map(a => a.xyz);
    viewer.frameCallbacks.add((t) => {
      const s = 0.45 * Math.sin(2 * Math.PI * t / 1.6);
      mol.update(x0.map((p, i) => [p[0] + s * imag.d[3 * i], p[1] + s * imag.d[3 * i + 1], p[2] + s * imag.d[3 * i + 2]]));
    });
  }
  let ext = 0;
  for (const a of atoms) ext = Math.max(ext, Math.hypot(...a.xyz));
  viewer.frame(ext + 3, false);
  viewer.resetView();
  $('viewport-title').innerHTML = `Stato di transizione ‡<small>${b.name} · ${d.pointGroup} · HF/${d.basis}</small>`;
  $('viewport-legend').innerHTML = '';
  $('viewport-note').textContent = `L'animazione segue il modo normale con frequenza immaginaria (${nf(Math.abs(imag?.freq ?? 0), 0)}i cm⁻¹): è la coordinata di reazione, la direzione in cui la struttura scivola verso reagenti o prodotti. Tratteggiati i legami che si stanno formando o rompendo.`;
}

function drawBarrierCharts() {
  const b = R.barrier;
  const t = barrierData(R.T);
  if (!t) return;
  const saved = R.method;
  R.method = 'hf';
  const hf = barrierData(298.15);
  R.method = 'mp2';
  const mp = barrierData(298.15);
  R.method = saved;
  drawProfile($('chart-radial'), {
    left: sideText(b.left), right: sideText(b.right),
    curves: [
      { label: 'HF', ts: hf.dH, product: hf.dHr, color: cssVar('--phase-neg') },
      { label: 'MP2', ts: mp.dH, product: mp.dHr, color: cssVar('--accent') },
    ],
    reference: b.expBarrier,
  });
  $('radial-title').textContent = 'Profilo energetico (ΔH, 298 K)';
  $('radial-note').textContent = 'Il massimo lungo la coordinata di reazione è lo stato di transizione. La linea tratteggiata indica la barriera di riferimento.';

  // grafico di Arrhenius: ln k contro 1000/T
  const Ts = [];
  for (let T = 200; T <= 2000; T += 20) Ts.push(T);
  const xs = Ts.map(T => 1000 / T);
  const lnk = Ts.map(T => Math.log(barrierData(T).k));
  // energia di attivazione di Arrhenius dalla pendenza a T corrente: Ea = −R d(ln k)/d(1/T)
  const h = 1;
  const Ea = -8.314462618e-3 * (Math.log(barrierData(R.T + h).k) - Math.log(barrierData(R.T - h).k)) / (1 / (R.T + h) - 1 / (R.T - h));
  drawLineChart($('chart-levels'), {
    series: [{ xs, ys: lnk, color: cssVar('--accent'), label: 'ln k (Eyring + tunnel)', width: 2.2, keepZero: true }],
    log: false, signed: true, xlabel: '1000/T (K⁻¹)', ylabel: `ln k (${t.molecularity === 2 ? 'cm³ s⁻¹' : 's⁻¹'})`,
    markers: [{ x: 1000 / R.T, label: `${nf(R.T, 0)} K`, color: cssVar('--accent') }],
  });
  $('levels-title').textContent = 'Grafico di Arrhenius';
  $('levels-note').textContent = `ln k = ln A − Ea/RT: la pendenza è −Ea/R. A ${nf(R.T, 0)} K, Ea = ${nf(Ea, 1)} kJ/mol (≈ ΔH‡ + ${t.molecularity === 2 ? '2' : ''}RT).`;

  drawBars($('chart-slice'), {
    unit: 'kJ/mol',
    items: [
      { label: 'ΔE‡ Hartree–Fock', value: t.dEhf },
      { label: 'Δ correlazione MP2', value: t.dEmp2 - t.dEhf },
      { label: 'ΔZPE‡', value: t.dZPE },
      { label: 'Δ(H − E − ZPE)‡', value: t.dH - t.dE - t.dZPE },
      { label: '−TΔS‡', value: -R.T * t.dS / 1000 },
      { label: 'ΔG‡', value: t.dG, color: cssVar('--accent') },
    ],
  });
  $('slice-title').textContent = `Contributi a ΔG‡ (${nf(R.T, 0)} K)`;
  $('slice-note').textContent = 'La barriera di energia libera determina la velocità: ogni 5,7 kJ/mol in più a 298 K rallentano la reazione di un fattore 10.';
}
