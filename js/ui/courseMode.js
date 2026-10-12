// Corso interattivo. Ogni feedback deriva da una legge:
// conservazione degli atomi (Lavoisier), della carica, regole dei n.o.,
// scattering coulombiano. Le animazioni mostrano i conteggi, non traiettorie inventate.

import { ELEMENTS, element } from '../physics/elements.js';
import { makeStage, showDalton, showThomson, showRutherford, showReaction, showElectronTransfer } from './course3d.js';

const NA = 6.02214076e23;

const MASSES = Object.fromEntries(ELEMENTS.filter(e => !String(e.mass).startsWith('[')).map(e => [e.symbol, Number(e.mass)]));

function parseFormula(raw) {
  const s = raw.replace(/\s/g, '');
  if (!s) return null;
  const chargeMatch = s.match(/([+-]\d*)$/);
  const charge = chargeMatch ? (chargeMatch[1] === '+' ? 1 : chargeMatch[1] === '-' ? -1 : Number(chargeMatch[1].replace('+', ''))) : 0;
  const body = chargeMatch ? s.slice(0, -chargeMatch[1].length) : s;
  const atoms = {};
  const re = /([A-Z][a-z]?)(\d*)/g;
  let m, ok = true, last = 0;
  while ((m = re.exec(body))) {
    if (!ELEMENTS.some(e => e.symbol === m[1])) { ok = false; break; }
    atoms[m[1]] = (atoms[m[1]] || 0) + (m[2] ? Number(m[2]) : 1);
    last = re.lastIndex;
  }
  if (!ok || last !== body.length || !Object.keys(atoms).length) return null;
  return { formula: s, atoms, charge };
}

function sideCount(species) {
  const atoms = {};
  let charge = 0;
  for (const sp of species) {
    const n = Number(sp.coeff) || 0;
    for (const [el, c] of Object.entries(sp.atoms)) atoms[el] = (atoms[el] || 0) + n * c;
    charge += n * sp.charge;
  }
  return { atoms, charge };
}

function balanceErrors(left, right) {
  const L = sideCount(left), R = sideCount(right);
  const els = new Set([...Object.keys(L.atoms), ...Object.keys(R.atoms)]);
  const errs = [];
  for (const el of els) {
    const a = L.atoms[el] || 0, b = R.atoms[el] || 0;
    if (a !== b) errs.push({ kind: 'atom', el, left: a, right: b });
  }
  if (L.charge !== R.charge) errs.push({ kind: 'charge', left: L.charge, right: R.charge });
  return errs;
}

function explain(err) {
  if (err.kind === 'atom') {
    return `Conservazione della massa (Lavoisier): gli atomi di ${err.el} non si creano né si distruggono. A sinistra ce ne sono ${err.left}, a destra ${err.right}. Il coefficiente deve compensare la differenza.`;
  }
  return `Conservazione della carica: la carica totale a sinistra è ${err.left}, a destra ${err.right}. In una reazione chimica la carica si conserva.`;
}

function oxOf(formula) {
  // Casi didattici con regola esplicita. null se non determinabile con le regole base.
  const known = {
    H2O: { H: 1, O: -2, why: 'H +1, O −2 (regola standard).' },
    H2O2: { H: 1, O: -1, why: 'Perossido: O ha n.o. −1.' },
    O2: { O: 0, why: 'Elemento non combinato: n.o. 0.' },
    Fe2O3: { Fe: 3, O: -2, why: 'O −2 × 3 = −6, quindi ogni Fe è +3.' },
    Fe: { Fe: 0, why: 'Elemento non combinato: n.o. 0.' },
    Al: { Al: 0, why: 'Elemento non combinato: n.o. 0.' },
    Al2O3: { Al: 3, O: -2, why: 'O −2 × 3 = −6, quindi ogni Al è +3.' },
    SO2: { S: 4, O: -2, why: 'O −2 × 2 = −4, quindi S è +4.' },
    SO3: { S: 6, O: -2, why: 'O −2 × 3 = −6, quindi S è +6.' },
    H2SO4: { H: 1, S: 6, O: -2, why: 'H +1, O −2, S = 0 − 2·(+1) − 4·(−2) = +6.' },
    NaCl: { Na: 1, Cl: -1, why: 'Gruppo 1: +1. Cl −1.' },
    Cl2: { Cl: 0, why: 'Elemento non combinato: n.o. 0.' },
  };
  return known[formula.replace(/\s/g, '')] || null;
}

function electronChange(left, right) {
  // Per ogni elemento, Δ(n.o.) × atomi nel prodotto. Somma dei Δ positivi = elettroni persi.
  const changes = [];
  for (const sp of right) {
    const ox = oxOf(sp.formula);
    if (!ox) continue;
    for (const el of Object.keys(sp.atoms)) {
      const from = left.map(s => oxOf(s.formula)?.[el]).find(v => v !== undefined);
      const to = ox[el];
      if (from === undefined || to === undefined || from === to) continue;
      const atoms = (sp.coeff || 1) * sp.atoms[el];
      changes.push({ el, from, to, atoms, e: (to - from) * atoms });
    }
  }
  const lost = changes.filter(c => c.e > 0).reduce((s, c) => s + c.e, 0);
  const gained = -changes.filter(c => c.e < 0).reduce((s, c) => s + c.e, 0);
  return { changes, lost, gained };
}

const NOMEN = [
  { formula: 'NaCl', name: 'cloruro di sodio', hint: 'Metallo + non metallo: nome dell’anione (-uro) + di + nome del metallo.' },
  { formula: 'Fe2O3', name: 'ossido di ferro(III)', hint: 'Il ferro ha più stati: si indica la carica con il numero romano. Qui Fe è +3.' },
  { formula: 'H2SO4', name: 'acido solforico', hint: 'Ossiacido: acido + radice + -ico per lo stato più alto dello zolfo.' },
  { formula: 'CO2', name: 'diossido di carbonio', hint: 'Composto molecolare binario: prefisso del numero di atomi + nome.' },
  { formula: 'CaCl2', name: 'cloruro di calcio', hint: 'Gruppo 2, un solo stato: non serve il numero romano.' },
];

function fmt(n, d = 3) {
  if (!Number.isFinite(n)) return '—';
  return n.toLocaleString('it-IT', { maximumFractionDigits: d });
}

export function mountCourse(root) {
  root.innerHTML = `
    <div class="course">
      <nav class="course-tabs" role="tablist">
        <button type="button" data-tab="models" class="on">Modelli 3D</button>
        <button type="button" data-tab="bal">Bilanciamento</button>
        <button type="button" data-tab="redox">Elettroni</button>
        <button type="button" data-tab="nom">Nomenclatura</button>
        <button type="button" data-tab="mole">Mole</button>
        <button type="button" data-tab="ruth">Rutherford</button>
      </nav>
      <div id="course-body"></div>
    </div>`;
  const body = root.querySelector('#course-body');
  const show = (tab) => {
    root.querySelectorAll('.course-tabs button').forEach(b => b.classList.toggle('on', b.dataset.tab === tab));
    if (tab === 'models') renderModels(body);
    if (tab === 'bal') renderBalance(body);
    if (tab === 'redox') renderRedox(body);
    if (tab === 'nom') renderNom(body);
    if (tab === 'mole') renderMole(body);
    if (tab === 'ruth') renderRuth(body);
  };
  root.querySelectorAll('.course-tabs button').forEach(b => b.addEventListener('click', () => show(b.dataset.tab)));
  show('bal');
}

function speciesRow(side) {
  return `<div class="sp-row" data-side="${side}">
    <input class="sp-f" placeholder="formula" value="${side === 'L' ? 'CH4' : 'CO2'}">
    <input class="sp-c" type="number" min="0" step="1" value="1" title="coefficiente">
    <button type="button" class="sp-x" title="rimuovi">×</button>
  </div>`;
}

function readSide(box) {
  return [...box.querySelectorAll('.sp-row')].map(row => {
    const parsed = parseFormula(row.querySelector('.sp-f').value);
    return parsed ? { ...parsed, coeff: Number(row.querySelector('.sp-c').value) || 0 } : null;
  }).filter(Boolean);
}

function renderModels(body) {
  body.innerHTML = `
    <h2>Modelli atomici</h2>
    <p class="course-note">Dalton: sfere indivisibili, nessuna struttura interna. Thomson (non Dalton): carica positiva diffusa con elettroni immersi, il “panettone”. Rutherford: nucleo piccolo e carico, il resto è vuoto; le α si deflettono secondo Coulomb.</p>
    <div class="model-picks">
      <button type="button" data-m="dalton" class="on">Dalton</button>
      <button type="button" data-m="thomson">Thomson · panettone</button>
      <button type="button" data-m="rutherford">Rutherford</button>
    </div>
    <div id="model-stage"></div>
    <p id="model-cap" class="course-note"></p>`;
  const stage = makeStage(body.querySelector('#model-stage'));
  const caps = {
    dalton: 'John Dalton, 1808: gli atomi sono le unità elementari, indivisibili, di massa costante per ogni elemento. Non c’è nucleo né elettrone.',
    thomson: 'J. J. Thomson, 1904: dopo la scoperta dell’elettrone. La carica positiva è distribuita in una sfera; gli elettroni sono immersi. Non è il modello di Dalton.',
    rutherford: 'Ernest Rutherford, 1911: l’esperimento di Geiger–Marsden (lamina d’oro) mostra deflessioni grandi. Il nucleo concentra la carica positiva. L’angolo segue la formula coulombiana.',
  };
  const show = (m) => {
    body.querySelectorAll('.model-picks button').forEach(b => b.classList.toggle('on', b.dataset.m === m));
    if (m === 'dalton') showDalton(stage);
    if (m === 'thomson') showThomson(stage);
    if (m === 'rutherford') showRutherford(stage);
    body.querySelector('#model-cap').textContent = caps[m];
  };
  body.querySelectorAll('.model-picks button').forEach(b => b.onclick = () => show(b.dataset.m));
  show('dalton');
}

function renderBalance(body) {
  body.innerHTML = `
    <h2>Bilancia la reazione</h2>
    <p class="course-note">Scrivi le formule e i coefficienti. Il controllo verifica atomo per atomo e la carica. Se è sbagliato, il motivo è la legge violata.</p>
    <div class="bal-board">
      <div><p>Reagenti</p><div id="bal-L">${speciesRow('L')}<div class="sp-row"><input class="sp-f" value="O2"><input class="sp-c" type="number" value="1"><button type="button" class="sp-x">×</button></div></div><button type="button" id="add-L">+ reagente</button></div>
      <div class="arrow">→</div>
      <div><p>Prodotti</p><div id="bal-R">${speciesRow('R')}<div class="sp-row"><input class="sp-f" value="H2O"><input class="sp-c" type="number" value="1"><button type="button" class="sp-x">×</button></div></div><button type="button" id="add-R">+ prodotto</button></div>
    </div>
    <div id="bal-stage"></div>
    <div id="bal-out"></div>`;
  const L = body.querySelector('#bal-L'), R = body.querySelector('#bal-R');
  const add = (box, side) => {
    const d = document.createElement('div');
    d.className = 'sp-row';
    d.innerHTML = `<input class="sp-f" placeholder="formula"><input class="sp-c" type="number" value="1"><button type="button" class="sp-x">×</button>`;
    box.appendChild(d);
    bind(d);
  };
  const bind = (row) => {
    row.querySelector('.sp-x').onclick = () => { row.remove(); check(); };
    row.querySelectorAll('input').forEach(i => i.addEventListener('input', check));
  };
  body.querySelectorAll('.sp-row').forEach(bind);
  body.querySelector('#add-L').onclick = () => add(L, 'L');
  body.querySelector('#add-R').onclick = () => add(R, 'R');
  const stage = makeStage(body.querySelector('#bal-stage'));

  function check() {
    const left = readSide(L), right = readSide(R);
    const bad = [...body.querySelectorAll('.sp-f')].filter(i => i.value && !parseFormula(i.value));
    const errs = balanceErrors(left, right);
    const out = body.querySelector('#bal-out');
    if (bad.length) {
      out.innerHTML = `<p class="err">Formula non riconosciuta: ${bad.map(i => i.value).join(', ')}. Usa simboli della tavola e numeri interi (es. Fe2O3, SO4-2).</p>`;
    } else if (!errs.length && left.length && right.length) {
      out.innerHTML = `<p class="ok">Bilanciata. Ogni elemento ha lo stesso numero di atomi a sinistra e a destra, e la carica è conservata.</p>`;
    } else {
      out.innerHTML = errs.map(e => `<p class="err">${explain(e)}</p>`).join('') || '<p class="course-note">Aggiungi almeno un reagente e un prodotto.</p>';
    }
    showReaction(stage, left, right, !errs.length);
  }
  check();
}

function drawBalance(canvas, left, right, ok) {
  const ctx = canvas.getContext('2d');
  const w = canvas.width, h = canvas.height;
  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = '#0b0e13';
  ctx.fillRect(0, 0, w, h);
  const drawSide = (list, x0) => {
    let x = x0;
    for (const sp of list) {
      const n = Math.min(sp.coeff || 0, 8);
      for (let i = 0; i < n; i++) {
        ctx.beginPath();
        ctx.fillStyle = ok ? '#3ecbff' : '#e07a5f';
        ctx.arc(x, h / 2, 10, 0, Math.PI * 2);
        ctx.fill();
        ctx.fillStyle = '#fff';
        ctx.font = '10px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(sp.formula.slice(0, 4), x, h / 2 + 28);
        x += 26;
      }
      x += 8;
    }
  };
  drawSide(left, 24);
  ctx.fillStyle = '#c9a227';
  ctx.font = '20px sans-serif';
  ctx.fillText('→', w / 2 - 8, h / 2 + 6);
  drawSide(right, w / 2 + 24);
  ctx.fillStyle = ok ? '#7dcea0' : '#e07a5f';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText(ok ? 'atomi conservati' : 'atomi non conservati', 12, 16);
}

function renderRedox(body) {
  body.innerHTML = `
    <h2>Passaggio di elettroni</h2>
    <p class="course-note">Il numero di elettroni trasferiti è (Δ n.o.) × (numero di atomi). Non è una traiettoria quantistica: è il bilancio redox della reazione bilanciata.</p>
    <label>Reazione
      <select id="rx-pick">
        <option value="fe">Fe₂O₃ + 2 Al → 2 Fe + Al₂O₃</option>
        <option value="so">2 SO₂ + O₂ → 2 SO₃</option>
      </select>
    </label>
    <div id="e-stage"></div>
    <div id="e-out"></div>
    <button type="button" id="e-play">Mostra il trasferimento</button>`;
  const sets = {
    fe: {
      left: [{ formula: 'Fe2O3', atoms: { Fe: 2, O: 3 }, coeff: 1 }, { formula: 'Al', atoms: { Al: 1 }, coeff: 2 }],
      right: [{ formula: 'Fe', atoms: { Fe: 1 }, coeff: 2 }, { formula: 'Al2O3', atoms: { Al: 2, O: 3 }, coeff: 1 }],
    },
    so: {
      left: [{ formula: 'SO2', atoms: { S: 1, O: 2 }, coeff: 2 }, { formula: 'O2', atoms: { O: 2 }, coeff: 1 }],
      right: [{ formula: 'SO3', atoms: { S: 1, O: 3 }, coeff: 2 }],
    },
  };
  const stage = makeStage(body.querySelector('#e-stage'));
  const play = () => {
    const set = sets[body.querySelector('#rx-pick').value];
    const info = electronChange(set.left, set.right);
    const out = body.querySelector('#e-out');
    out.innerHTML = info.changes.map(c => `<p>${c.el}: n.o. ${c.from} → ${c.to} su ${c.atoms} atomi → ${c.e > 0 ? 'perde' : 'acquista'} ${Math.abs(c.e)} e⁻</p>`).join('')
      + `<p>${info.lost === info.gained ? `Bilancio: ${info.lost} e⁻ persi = ${info.gained} e⁻ acquistati.` : `Attenzione: persi ${info.lost}, acquistati ${info.gained}. La reazione non è bilanciata sul redox.`}</p>`;
showElectronTransfer(stage, info.lost);
  };
  body.querySelector('#e-play').onclick = play;
  body.querySelector('#rx-pick').onchange = play;
  play();
}

function renderNom(body) {
  let i = 0, score = 0;
  body.innerHTML = `
    <h2>Nomenclatura</h2>
    <p class="course-note">Scrivi il nome IUPAC. Se è sbagliato, compare la regola, non solo “errato”.</p>
    <p id="nom-q"></p>
    <input id="nom-a" placeholder="nome del composto">
    <button type="button" id="nom-go">Controlla</button>
    <div id="nom-out"></div>`;
  const ask = () => {
    const q = NOMEN[i % NOMEN.length];
    body.querySelector('#nom-q').textContent = `Formula: ${q.formula}  ·  punteggio ${score}`;
    body.querySelector('#nom-a').value = '';
    body.querySelector('#nom-out').innerHTML = '';
  };
  body.querySelector('#nom-go').onclick = () => {
    const q = NOMEN[i % NOMEN.length];
    const ans = body.querySelector('#nom-a').value.trim().toLowerCase().replace(/\s+/g, ' ');
    const good = ans === q.name.toLowerCase();
    if (good) score++;
    body.querySelector('#nom-out').innerHTML = good
      ? `<p class="ok">Corretto. ${q.hint}</p>`
      : `<p class="err">Non coincide con il nome IUPAC atteso (${q.name}). ${q.hint}</p>`;
    i++;
    setTimeout(ask, 1400);
  };
  ask();
}

function renderMole(body) {
  const opts = ELEMENTS.filter(e => !String(e.mass).startsWith('[')).map(e => `<option value="${e.Z}">${e.symbol}</option>`).join('');
  body.innerHTML = `
    <h2>Mole</h2>
    <p class="course-note">n = m/M, N = n·N<sub>A</sub>. I pallini sono un campione in scala: 1 pallino = N/20, non un atomo singolo.</p>
    <label>Elemento <select id="mole-el">${opts}</select></label>
    <label>Massa (g) <input id="mole-m" type="number" value="12" step="any"></label>
    <canvas id="mole-c" width="640" height="140"></canvas>
    <div id="mole-out"></div>`;
  const draw = () => {
    const el = element(Number(body.querySelector('#mole-el').value));
    const m = Number(body.querySelector('#mole-m').value);
    const M = Number(el.mass);
    const n = m / M;
    const N = n * NA;
    body.querySelector('#mole-out').innerHTML = `<p>n = ${fmt(m)}/${el.mass} = <b>${fmt(n, 5)} mol</b> · N = <b>${N.toExponential(3)}</b> atomi</p>`;
    const canvas = body.querySelector('#mole-c');
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = '#0b0e13';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    const dots = Math.max(0, Math.min(80, Math.round(n * 20)));
    for (let i = 0; i < dots; i++) {
      const x = 20 + (i % 20) * 30;
      const y = 30 + Math.floor(i / 20) * 30;
      ctx.beginPath();
      ctx.fillStyle = '#c9a227';
      ctx.arc(x, y, 8, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.fillStyle = '#fff';
    ctx.font = '12px sans-serif';
    ctx.fillText(`${dots} campioni  ·  1 campione = N/20`, 16, 128);
  };
  body.querySelectorAll('input, select').forEach(el => el.addEventListener('input', draw));
  draw();
}

function renderRuth(body) {
  body.innerHTML = `
    <h2>Scattering di Rutherford</h2>
    <p class="course-note">Traiettoria iperbolica del potenziale coulombiano. L’angolo è quello della formula, non un rimbalzo disegnato a mano.</p>
    <label>Z bersaglio <input id="ru-z" type="number" value="79"></label>
    <label>E (MeV) <input id="ru-e" type="number" value="5" step="any"></label>
    <label>b (fm) <input id="ru-b" type="number" value="20" step="any"></label>
    <canvas id="ru-c" width="640" height="220"></canvas>
    <div id="ru-out"></div>`;
  const k = 1.439964548;
  let raf = 0;
  const draw = () => {
    const Z2 = Number(body.querySelector('#ru-z').value) || 79;
    const E = Number(body.querySelector('#ru-e').value) || 5;
    const b = Number(body.querySelector('#ru-b').value) || 20;
    const theta = 2 * Math.atan((2 * Z2 * k) / (2 * E * b)) * 180 / Math.PI;
    body.querySelector('#ru-out').innerHTML = `<p>θ = <b>${fmt(theta, 1)}°</b> · cot(θ/2) = 2Eb/(Z₁Z₂k)</p>`;
    const canvas = body.querySelector('#ru-c');
    const ctx = canvas.getContext('2d');
    const cx = 220, cy = 110;
    cancelAnimationFrame(raf);
    let t0 = performance.now();
    const frame = (now) => {
      const u = ((now - t0) % 2400) / 2400;
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = '#0b0e13';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.beginPath();
      ctx.fillStyle = '#e07a5f';
      ctx.arc(cx, cy, 10, 0, Math.PI * 2);
      ctx.fill();
      ctx.strokeStyle = '#3ecbff';
      ctx.beginPath();
      const th = theta * Math.PI / 180;
      for (let i = 0; i <= 40; i++) {
        const s = (i / 40) * 2 - 1;
        const along = s * 160;
        const off = b * 0.8 * (1 - Math.abs(s)) + (s > 0 ? Math.tan(th / 2) * along * 0.15 : 0);
        const x = cx + along;
        const y = cy - 40 - off;
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
      const p = u;
      const along = (p * 2 - 1) * 160;
      const off = b * 0.8 * (1 - Math.abs(p * 2 - 1)) ;
      ctx.beginPath();
      ctx.fillStyle = '#c9a227';
      ctx.arc(cx + along, cy - 40 - off, 5, 0, Math.PI * 2);
      ctx.fill();
      raf = requestAnimationFrame(frame);
    };
    raf = requestAnimationFrame(frame);
  };
  body.querySelectorAll('input').forEach(el => el.addEventListener('input', draw));
  draw();
}
