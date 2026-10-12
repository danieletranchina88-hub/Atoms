import { ELEMENTS } from '../physics/elements.js';
import {
  number,
  parseFormula,
  parseEquation,
  molarMass,
  checkBalance,
  solveBalance,
  equationText,
  stoichiometry,
  oxidation,
  redoxLedger,
  ionicFormula,
  dilution,
  isotopeMean,
  significantFigures,
  empiricalFormula,
  rutherford,
  chromatography,
  measurementStatistics,
  waterPhase,
  measuredCalculation,
  NA,
} from './chemistry.js';
import {
  REACTIONS,
  COMPOUNDS,
  CATIONS,
  ANIONS,
  MIX_REAGENTS,
  MIX_RULES,
  QUESTIONS,
  normalizeName,
} from './content.js';
import {
  canvasAnimation,
  reactionDiagram,
  dilutionDiagram,
  scatterDiagram,
  isotopeDiagram,
} from './visuals.js';
import {
  makeStage,
  showDalton,
  showThomson,
  showRutherford,
  showIsomers,
} from '../ui/course3d.js';

export const esc = (value) =>
  String(value).replace(
    /[&<>"']/g,
    (c) =>
      ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[
        c
      ],
  );
export const fmt = (n, d = 5) =>
  Number.isFinite(n)
    ? n.toLocaleString('it-IT', { maximumSignificantDigits: d })
    : '—';
export function formulaHTML(raw) {
  const s = typeof raw === 'string' ? parseFormula(raw) : raw;
  const body = s.body
    .split(/[·:]/)
    .map((p, i) => {
      const m = i ? p.match(/^\d+/) : null;
      return (
        (m ? esc(m[0]) : '') +
        esc(m ? p.slice(m[0].length) : p).replace(
          /\d+/g,
          (n) => `<sub>${n}</sub>`,
        )
      );
    })
    .join('·');
  return (
    body +
    (s.charge
      ? `<sup>${Math.abs(s.charge) === 1 ? '' : Math.abs(s.charge)}${s.charge > 0 ? '+' : '−'}</sup>`
      : '') +
    (s.phase ? `<small>(${esc(s.phase)})</small>` : '')
  );
}
export function equationHTML(eq) {
  const side = (arr) =>
    arr
      .map((s) => `${s.coeff === 1 ? '' : s.coeff + ' '}${formulaHTML(s)}`)
      .join(' <span class="cl-plus">+</span> ');
  return side(eq.left) + ' <span class="cl-arrow">→</span> ' + side(eq.right);
}
const input = (id, label, value, options = '') =>
  `<label class="cl-field" for="${id}"><span>${label}</span><input id="${id}" value="${esc(value)}" ${options}></label>`;
const select = (id, label, options) =>
  `<label class="cl-field" for="${id}"><span>${label}</span><select id="${id}">${options}</select></label>`;
const options = (arr) =>
  arr.map((x) => `<option value="${esc(x)}">${esc(x)}</option>`).join('');
const btn = (id, text, primary = false) =>
  `<button type="button" id="${id}" class="cl-btn${primary ? ' primary' : ''}">${text}</button>`;
const table = (headers, rows) =>
  `<div class="cl-table-wrap"><table class="cl-table"><thead><tr>${headers.map((x) => `<th scope="col">${x}</th>`).join('')}</tr></thead><tbody>${rows.map((r) => `<tr>${r.map((x) => `<td>${x}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;
const metrics = (arr) =>
  `<div class="cl-metrics">${arr.map(([label, val]) => `<div><small>${label}</small><strong>${val}</strong></div>`).join('')}</div>`;
const notice = (text, ok = false) =>
  `<div class="cl-feedback ${ok ? 'success' : ''}" role="status">${text}</div>`;
const visual = (id, caption) =>
  `<figure class="cl-visual"><canvas id="${id}" aria-label="${esc(caption)}"></canvas><figcaption>${caption}</figcaption></figure>`;
const shell = (root, title, lead, html) => {
  root.innerHTML = `<div class="cl-tool-head"><span class="cl-eyebrow">B A N C O &nbsp; D I &nbsp; L A V O R O</span><h3>${title}</h3><p>${lead}</p></div>${html}`;
};
const lookup = (id) => REACTIONS.find((r) => r.id === id) || REACTIONS[0];
const reactionOptions = (excludeUnbalanced = false) =>
  REACTIONS.filter((r) => !excludeUnbalanced || r.id !== 'emati')
    .map(
      (r) => `<option value="${r.id}">${esc(r.name)} · ${esc(r.type)}</option>`,
    )
    .join('');
function guarded(out, fn, onFail = () => {}) {
  try {
    fn();
  } catch (e) {
    out.innerHTML = notice(esc(e.message));
    onFail();
  }
}
function listen(root, selector, fn, event = 'input') {
  root
    .querySelectorAll(selector)
    .forEach((el) => el.addEventListener(event, fn));
}
const val = (root, id) => root.querySelector('#' + id).value;
const getEq = (root, id) => parseEquation(val(root, id));
function atomsKey(raw) {
  const s = parseFormula(raw);
  return (
    JSON.stringify(Object.entries(s.atoms).sort()) +
    ':' +
    s.charge +
    ':' +
    s.phase
  );
}
function sameFormula(raw, expected) {
  const a = parseFormula(raw),
    b = parseFormula(expected);
  return (
    JSON.stringify(Object.entries(a.atoms).sort()) ===
      JSON.stringify(Object.entries(b.atoms).sort()) &&
    a.charge === b.charge &&
    (!a.phase || a.phase === b.phase)
  );
}

function renderBalance(root, ctx) {
  shell(
    root,
    'Un bilancio che puoi costruire',
    'Scegli uno schema o scrivi le tue specie. Cambia i coefficienti, mai gli indici: il conteggio si aggiorna a ogni modifica.',
    `
    <div class="cl-form-row">${select('cl-rx', 'Schema di partenza', reactionOptions())}${btn('cl-reset', 'Tutti i coefficienti a 1')}</div>
    <details class="cl-custom"><summary>Scrivi una reazione libera</summary><p>Usa spazio + spazio tra le specie. Parentesi, idrati e cariche sono supportati: <code>Fe^3+</code>, <code>SO4^2-</code>, <code>Ca(OH)2</code>, <code>CuSO4·5H2O</code>.</p>${input('cl-custom-equation', 'Schema, con una freccia ->', 'Fe + O2 -> Fe2O3', 'maxlength="1000"')}${btn('cl-load', 'Carica schema')}</details>
    <div id="cl-editor" class="cl-equation-editor"></div>
    <div class="cl-actions">${btn('cl-solve', 'Mostra i coefficienti minimi')}${btn('cl-to-stoich', 'Calcola le quantità', true)}</div>
    <div id="cl-balance-out" aria-live="polite"></div><p class="cl-model-note">Conservare atomi e carica non dimostra spontaneità o velocità. Il bilanciatore restituisce solo soluzioni univoche, positive e intere.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id),
    out = $('cl-balance-out');
  let eq;
  const bound = new WeakSet();
  const row = (s, side, i) =>
    `<div class="cl-species" data-side="${side}"><label><span class="cl-sr-only">Coefficiente della specie ${i + 1} ${side === 'L' ? 'reagente' : 'prodotto'}</span><input class="cl-coeff" value="${s.coeff}" type="number" min="1" max="1000000" step="1"></label><label><span class="cl-sr-only">Formula della specie ${i + 1} ${side === 'L' ? 'reagente' : 'prodotto'}</span><input class="cl-formula" value="${esc(s.formula + (s.phase ? '(' + s.phase + ')' : ''))}" maxlength="200"></label><button type="button" class="cl-remove" aria-label="Rimuovi specie">×</button></div>`;
  const read = () => {
    const side = (s) =>
      [...$('cl-editor').querySelectorAll(`[data-side="${s}"]`)].map((el) => ({
        ...parseFormula(el.querySelector('.cl-formula').value),
        coeff: number(el.querySelector('.cl-coeff').value, {
          min: 1,
          max: 1000000,
          integer: true,
        }),
      }));
    return { left: side('L'), right: side('R') };
  };
  const check = () =>
    guarded(
      out,
      () => {
        eq = read();
        const info = checkBalance(eq.left, eq.right);
        ctx.shared.equation = equationText(eq.left, eq.right);
        out.innerHTML =
          `<div class="cl-equation">${equationHTML(eq)}</div>` +
          notice(
            info.balanced
              ? 'Atomi e carica conservati. Il bilancio è corretto.'
              : 'Confronta le colonne: ogni differenza deve diventare zero.',
            info.balanced,
          ) +
          table(
            ['Grandezza', 'Reagenti', 'Prodotti', 'Differenza'],
            info.rows.map((r) => [
              esc(r.el),
              r.left,
              r.right,
              `<span class="${r.left === r.right ? 'cl-ok' : 'cl-error'}">${r.right - r.left}</span>`,
            ]),
          );
        ctx.snapshot = {
          tipo: 'bilanciamento',
          equazione: ctx.shared.equation,
          bilanciata: info.balanced,
        };
      },
      () => {
        eq = null;
        ctx.snapshot = null;
      },
    );
  const bind = () => {
    root.querySelectorAll('.cl-species').forEach((row) => {
      if (bound.has(row)) return;
      bound.add(row);
      row
        .querySelectorAll('input')
        .forEach((i) => i.addEventListener('input', check));
    });
    root.querySelectorAll('.cl-remove').forEach(
      (b) =>
        (b.onclick = () => {
          b.closest('.cl-species').remove();
          check();
        }),
    );
    ['L', 'R'].forEach(
      (side) =>
        ($('cl-add-' + side).onclick = () => {
          const box = $('cl-side-' + side);
          if (box.children.length >= 8) {
            out.innerHTML = notice('Al massimo 8 specie per lato.');
            return;
          }
          box.insertAdjacentHTML(
            'beforeend',
            row(
              { ...parseFormula(side === 'L' ? 'O2' : 'H2O'), coeff: 1 },
              side,
              box.children.length,
            ),
          );
          bind();
          check();
        }),
    );
  };
  const show = (data) => {
    $('cl-editor').innerHTML =
      `<section><h4>Reagenti</h4><div id="cl-side-L">${data.left.map((s, i) => row(s, 'L', i)).join('')}</div><button type="button" id="cl-add-L" class="cl-btn">+ Specie</button></section><span class="cl-big-arrow">→</span><section><h4>Prodotti</h4><div id="cl-side-R">${data.right.map((s, i) => row(s, 'R', i)).join('')}</div><button type="button" id="cl-add-R" class="cl-btn">+ Specie</button></section>`;
    bind();
    check();
  };
  $('cl-rx').value = ctx.shared.rx || 'metano';
  const loadPreset = () => {
    const r = lookup(val(root, 'cl-rx'));
    ctx.shared.rx = r.id;
    const data = parseEquation(r.equation);
    [...data.left, ...data.right].forEach((s) => (s.coeff = 1));
    show(data);
  };
  $('cl-rx').onchange = loadPreset;
  $('cl-load').onclick = () =>
    guarded(out, () => show(parseEquation(val(root, 'cl-custom-equation'))));
  $('cl-reset').onclick = () => {
    root.querySelectorAll('.cl-coeff').forEach((i) => (i.value = 1));
    check();
  };
  $('cl-solve').onclick = () =>
    guarded(out, () => {
      const data = read();
      show(solveBalance(data.left, data.right));
      out.insertAdjacentHTML(
        'beforeend',
        '<p class="cl-hint">Metodo: un vincolo per ogni elemento e per la carica, soluzione intera ridotta al minimo comune divisore.</p>',
      );
    });
  $('cl-to-stoich').onclick = () =>
    guarded(out, () => {
      const data = read();
      if (!checkBalance(data.left, data.right).balanced)
        throw Error('Bilancia prima di passare alle quantità.');
      ctx.shared.equation = equationText(data.left, data.right);
      ctx.showTool('stechiometria');
    });
  loadPreset();
  return () => {};
}

function renderStoichiometry(root, ctx, experiment = false) {
  shell(
    root,
    experiment
      ? 'Il tuo banco delle reazioni'
      : 'Dai coefficienti alle quantità',
    experiment
      ? 'Scegli una trasformazione, oppure prova a combinare i reagenti. Modifica le quantità e fai avanzare il bilancio.'
      : 'Lavora in moli o grammi, trova il limitante e controlla dove resta la materia.',
    `
    ${experiment ? `<details class="cl-custom"><summary>Mescola liberamente due reagenti</summary><div class="cl-form-row">${select('cl-mix-a', 'Primo reagente', options(MIX_REAGENTS))}${select('cl-mix-b', 'Secondo reagente', options(MIX_REAGENTS))}${btn('cl-mix', 'Prova la combinazione', true)}</div><p class="cl-hint">Le combinazioni riconosciute usano un repertorio di reazioni verificato. Una combinazione non prevista non viene trasformata in una reazione inventata.</p><div id="cl-mix-result" aria-live="polite"></div></details>` : ''}
    <div class="cl-form-row">${select('cl-st-rx', 'Reazione', reactionOptions(true))}${btn('cl-st-load', 'Carica esempio')}</div>
    ${input('cl-st-equation', 'Equazione libera, già bilanciata', ctx.shared.equation || lookup(ctx.shared.rx || 'metano').equation, 'maxlength="1000"')}
    <div class="cl-actions">${btn('cl-st-use', 'Usa questa equazione')}${btn('cl-st-balance', 'Bilancia qui')}</div>
    <div id="cl-st-inputs" class="cl-amounts"></div>
    <div class="cl-form-row">${input('cl-yield', 'Resa di isolamento (%)', 100, 'type="number" min="0" max="100" step="any"')}
      <label class="cl-field cl-grow" for="cl-extent"><span>Avanzamento: <output id="cl-extent-label">100%</output></span><input id="cl-extent" type="range" min="0" max="100" value="100" step="0.1"></label>${btn('cl-play', '▶ Anima')}</div>
    ${visual('cl-rx-canvas', 'Rappresentazione schematica: marcatori e colori didattici. Il cursore controlla l’avanzamento, non il tempo fisico o la cinetica.')}
    <div id="cl-st-out" aria-live="polite"></div><div id="cl-rx-description"></div>`,
  );
  const $ = (id) => root.querySelector('#' + id),
    out = $('cl-st-out');
  let eq = null,
    model = null,
    playing = false,
    last = 0,
    frame = 0,
    rx = lookup(ctx.shared.rx || 'metano'),
    custom = false;
  $('cl-st-rx').value = rx.id;
  const update = () =>
    guarded(
      out,
      () => {
        if (!eq) return;
        const fraction =
            number(val(root, 'cl-extent'), { min: 0, max: 100 }) / 100,
          yieldPercent = number(val(root, 'cl-yield'), { min: 0, max: 100 });
        const amounts = [
          ...$('cl-st-inputs').querySelectorAll('.cl-amount'),
        ].map((el) => ({
          value: el.querySelector('input').value,
          unit: el.querySelector('select').value,
        }));
        const result = stoichiometry(
          eq.left,
          eq.right,
          amounts,
          fraction,
          yieldPercent,
        );
        $('cl-extent-label').textContent = fmt(fraction * 100, 4) + '%';
        model = { rx, left: eq.left, right: eq.right, result, fraction };
        const limiting = result.limiting
          .map((i) => formulaHTML(eq.left[i]))
          .join(' e ');
        out.innerHTML =
          `<div class="cl-equation">${equationHTML(eq)}</div>` +
          metrics([
            ['Limitante a conversione completa', limiting],
            ['Avanzamento massimo', fmt(result.maximum) + ' mol'],
            ['Avanzamento attuale', fmt(result.extent) + ' mol'],
          ]) +
          table(
            ['Reagente', 'Iniziali (mol)', 'Consumate (mol)', 'Residue (mol)'],
            eq.left.map((s, i) => [
              formulaHTML(s),
              fmt(result.initial[i]),
              fmt(result.initial[i] - result.remaining[i]),
              fmt(result.remaining[i]),
            ]),
          ) +
          table(
            [
              'Prodotto',
              'Formate ora (mol)',
              'Massa ora (g)',
              'Teorica finale (g)',
              'Isolabile finale (g)',
            ],
            result.products.map((p, i) => [
              formulaHTML(eq.right[i]),
              fmt(p.mol),
              fmt(p.grams),
              fmt(p.theoreticalGrams),
              fmt(p.isolatedGrams),
            ]),
          );
        out.insertAdjacentHTML(
          'beforeend',
          '<p class="cl-hint">Passaggi: n = m/M se usi grammi → ξ massimo = min(nᵢ/νᵢ) → n consumate = νᵢξ. La resa di isolamento riguarda il prodotto recuperato alla fine: non cambia i bilanci della reazione.</p>',
        );
        ctx.snapshot = {
          tipo: experiment ? 'esperimento' : 'stechiometria',
          equazione: equationText(eq.left, eq.right),
          quantita: amounts,
          avanzamento: result.extent,
          moliResidue: result.remaining,
          prodotti: result.products,
          resaPercentuale: yieldPercent,
          modello: 'Bilancio stechiometrico; nessuna cinetica calcolata',
        };
      },
      () => {
        model = null;
        playing = false;
        ctx.snapshot = null;
        $('cl-play').textContent = '▶ Anima';
      },
    );
  const load = (raw, asCustom = false) =>
    guarded(
      out,
      () => {
        const next = parseEquation(raw);
        if (!checkBalance(next.left, next.right).balanced)
          throw Error(
            'L’equazione non è bilanciata. Premi “Bilancia qui” per calcolare i coefficienti.',
          );
        // Do not leave an old model active behind an invalid new equation.
        eq = next;
        custom = asCustom;
        ctx.shared.equation = equationText(eq.left, eq.right);
        if (custom)
          rx = {
            name: 'Reazione libera',
            id: 'custom',
            type: 'Bilancio impostato dall’utente',
            visual: 'mix',
            source: 'Schema libero',
            note: 'Il bilancio verifica atomi e carica. Prodotti, spontaneità, meccanismo e condizioni devono essere stabiliti da te: non sono dedotti dal bilanciamento.',
          };
        $('cl-st-equation').value = ctx.shared.equation;
        $('cl-st-inputs').innerHTML = eq.left
          .map(
            (s, i) =>
              `<div class="cl-amount">${input('cl-amount-' + i, formulaHTML(s) + ' · quantità iniziale', s.coeff, 'type="number" min="0" max="1000000000" step="any"')}<label class="cl-field" for="cl-unit-${i}"><span>Unità</span><select id="cl-unit-${i}"><option value="mol">mol</option>${s.electron ? '' : '<option value="g">g</option>'}</select></label></div>`,
          )
          .join('');
        $('cl-rx-description').innerHTML =
          `<p class="cl-source-inline">${esc(rx.source)} · ${esc(rx.type)}</p><p class="cl-model-note">${esc(rx.note)}</p>${rx.net ? `<p class="cl-hint">Equazione ionica netta</p><div class="cl-equation">${equationHTML(parseEquation(rx.net))}</div>` : ''}`;
        listen(root, '#cl-st-inputs input,#cl-st-inputs select', update);
        $('cl-extent').value = 100;
        playing = false;
        $('cl-play').textContent = '▶ Anima';
        update();
      },
      () => {
        eq = null;
        model = null;
        playing = false;
        ctx.snapshot = null;
        $('cl-st-inputs').innerHTML = '';
        $('cl-rx-description').innerHTML = '';
      },
    );
  $('cl-st-load').onclick = () => {
    rx = lookup(val(root, 'cl-st-rx'));
    ctx.shared.rx = rx.id;
    load(rx.equation);
  };
  $('cl-st-rx').onchange = $('cl-st-load').onclick;
  $('cl-st-use').onclick = () => load(val(root, 'cl-st-equation'), true);
  $('cl-st-equation').oninput = () => {
    eq = null;
    model = null;
    playing = false;
    ctx.snapshot = null;
    $('cl-play').textContent = '▶ Anima';
    $('cl-st-inputs')
      .querySelectorAll('input,select')
      .forEach((el) => (el.disabled = true));
    out.innerHTML = notice(
      'Equazione modificata: premi “Usa questa equazione” per validarla.',
    );
  };
  $('cl-st-balance').onclick = () =>
    guarded(out, () => {
      const data = getEq(root, 'cl-st-equation');
      const solved = solveBalance(data.left, data.right);
      load(equationText(solved.left, solved.right), true);
    });
  listen(root, '#cl-extent,#cl-yield', () => {
    playing = false;
    $('cl-play').textContent = '▶ Anima';
    update();
  });
  $('cl-play').onclick = () => {
    if (!eq || !model) return;
    playing = !playing;
    if (playing && Number(val(root, 'cl-extent')) >= 100)
      $('cl-extent').value = 0;
    $('cl-play').textContent = playing ? 'Ⅱ Pausa' : '▶ Anima';
    update();
  };
  const reduced = matchMedia('(prefers-reduced-motion: reduce)');
  const loop = (now) => {
    frame = requestAnimationFrame(loop);
    const dt = last ? Math.min((now - last) / 1000, 0.05) : 0;
    last = now;
    if (playing && !ctx.paused && !document.hidden) {
      const progress = Math.min(
        100,
        Number(val(root, 'cl-extent')) + dt * (reduced.matches ? 150 : 12),
      );
      $('cl-extent').value = progress;
      update();
      if (progress >= 100) {
        playing = false;
        $('cl-play').textContent = '↺ Ripeti';
      }
    }
  };
  frame = requestAnimationFrame(loop);
  if (experiment) {
    $('cl-mix-b').value = 'NaOH';
    $('cl-mix').onclick = () => {
      const a = val(root, 'cl-mix-a'),
        b = val(root, 'cl-mix-b'),
        rule = MIX_RULES.find(
          (r) => (r[0] === a && r[1] === b) || (r[0] === b && r[1] === a),
        );
      if (rule) {
        rx = lookup(rule[2]);
        ctx.shared.rx = rx.id;
        $('cl-st-rx').value = rx.id;
        load(rx.equation);
        $('cl-mix-result').innerHTML = notice(
          'Combinazione riconosciuta: ' +
            esc(rx.name) +
            '. Modifica ora le quantità sotto.',
          true,
        );
      } else {
        playing = false;
        model = null;
        eq = null;
        ctx.snapshot = null;
        $('cl-st-inputs').innerHTML = '';
        $('cl-rx-description').innerHTML = '';
        out.innerHTML = '';
        $('cl-play').textContent = '▶ Anima';
        const copper = [a, b].includes('Cu') && [a, b].includes('HCl');
        $('cl-st-equation').value = '';
        $('cl-mix-result').innerHTML = notice(
          copper
            ? 'Cu + HCl diluito, senza ossidanti: non si prevede dissoluzione spontanea del rame con produzione di H₂. Il cloruro rameico non si ottiene da questo semplice bilancio.'
            : 'Nessuna trasformazione è prevista dal repertorio per questa coppia. Questo non dimostra inerzia in ogni condizione: puoi studiare un’equazione nota nello spazio libero.',
        );
      }
    };
  }
  const stopCanvas = canvasAnimation(
    $('cl-rx-canvas'),
    (c, t, p) => reactionDiagram(c, t, p, model),
    { getPaused: () => ctx.paused },
  );
  const shared = ctx.shared.equation,
    known = REACTIONS.find((r) => r.equation === shared);
  if (known) rx = known;
  load(shared || rx.equation, !!shared && !known);
  return () => {
    cancelAnimationFrame(frame);
    stopCanvas();
  };
}

function renderNomenclature(root, ctx) {
  shell(
    root,
    'Un nome, una formula, una regola',
    'Allenati in entrambe le direzioni. Stock, sistematica e tradizionale sono richieste distinte, con varianti accettate esplicitamente.',
    `
    <div class="cl-form-row">${select('cl-nom-scheme', 'Convenzione', '<option value="stock">Stock</option><option value="systematic">Sistematica con prefissi</option><option value="traditional">Tradizionale</option>')}${select('cl-nom-direction', 'Direzione', '<option value="name">Formula → nome</option><option value="formula">Nome → formula</option>')}</div>
    <div class="cl-question"><span class="cl-eyebrow">S C R I V I &nbsp; L A &nbsp; R I S P O S T A</span><h4 id="cl-nom-question"></h4>${input('cl-nom-answer', 'La tua risposta', '', 'autocomplete="off" maxlength="150"')}<div class="cl-actions">${btn('cl-nom-check', 'Controlla', true)}${btn('cl-nom-next', 'Prossimo composto')}${btn('cl-nom-hint', 'Spiega la regola')}</div><div id="cl-nom-out" aria-live="polite"></div></div>
    <details class="cl-custom"><summary>Costruisci un sale dagli ioni</summary><div class="cl-form-row">${select('cl-cation', 'Catione', options(CATIONS))}${select('cl-anion', 'Anione', options(ANIONS))}</div><div id="cl-ion-out" aria-live="polite"></div></details>
    <details class="cl-custom"><summary>Consulta il repertorio dei nomi</summary>${input('cl-nom-search', 'Cerca formula o nome', '', 'type="search" placeholder="es. FeCl3, solfato, acido"')}<div id="cl-nom-list"></div></details>
    <p class="cl-source-inline">L3 · pp. 18–26. Le voci aggiuntive applicano le regole della lezione.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id);
  let index = ctx.shared.nomIndex || 0,
    q,
    lastAnswer = '';
  const ask = () => {
    const scheme = val(root, 'cl-nom-scheme'),
      list = COMPOUNDS.filter((c) => c[scheme].length);
    q = list[index % list.length];
    lastAnswer = '';
    $('cl-nom-question').innerHTML =
      val(root, 'cl-nom-direction') === 'name'
        ? formulaHTML(q.formula)
        : esc(q[scheme][0]);
    $('cl-nom-answer').value = '';
    $('cl-nom-out').innerHTML = '';
    ctx.snapshot = null;
  };
  $('cl-nom-check').onclick = () => {
    const answer = val(root, 'cl-nom-answer').trim();
    if (!answer) {
      $('cl-nom-out').innerHTML = notice(
        'Scrivi una risposta prima di controllare.',
      );
      return;
    }
    if (answer === lastAnswer) return;
    lastAnswer = answer;
    const scheme = val(root, 'cl-nom-scheme'),
      dir = val(root, 'cl-nom-direction');
    let good = false;
    try {
      good =
        dir === 'name'
          ? q[scheme].some((n) => normalizeName(n) === normalizeName(answer))
          : sameFormula(answer, q.formula);
    } catch {}
    const expected =
      dir === 'name' ? q[scheme].map(esc).join(' / ') : formulaHTML(q.formula);
    $('cl-nom-out').innerHTML = notice(
      (good
        ? 'Corretto. '
        : 'Da rivedere. Risposta attesa: ' + expected + '. ') + esc(q.rule),
      good,
    );
    ctx.recordAttempt(
      'nom:' + scheme + ':' + dir + ':' + q.formula,
      good,
      'composti',
    );
    ctx.snapshot = {
      tipo: 'nomenclatura',
      formula: q.formula,
      convenzione: scheme,
      risposta: answer,
      corretta: good,
      regola: q.rule,
    };
  };
  $('cl-nom-answer').onkeydown = (e) => {
    if (e.key === 'Enter') $('cl-nom-check').click();
  };
  $('cl-nom-next').onclick = () => {
    index++;
    ctx.shared.nomIndex = index;
    ask();
  };
  $('cl-nom-hint').onclick = () => {
    $('cl-nom-out').innerHTML = notice(esc(q.rule));
  };
  $('cl-nom-scheme').onchange = ask;
  $('cl-nom-direction').onchange = ask;
  const build = () =>
    guarded($('cl-ion-out'), () => {
      const data = ionicFormula(val(root, 'cl-cation'), val(root, 'cl-anion'));
      $('cl-ion-out').innerHTML =
        metrics([
          ['Formula neutra', formulaHTML(data.formula)],
          ['Cationi', data.cationCount],
          ['Anioni', data.anionCount],
        ]) +
        `<p>${data.cationCount} × ${parseFormula(val(root, 'cl-cation')).charge} + ${data.anionCount} × (${parseFormula(val(root, 'cl-anion')).charge}) = 0</p>`;
    });
  listen(root, '#cl-cation,#cl-anion', build);
  const search = () => {
    const term = normalizeName(val(root, 'cl-nom-search'));
    const rows = COMPOUNDS.filter((c) =>
      normalizeName(
        c.formula +
          ' ' +
          [...c.stock, ...c.systematic, ...c.traditional].join(' '),
      ).includes(term),
    );
    $('cl-nom-list').innerHTML = table(
      ['Formula', 'Stock', 'Sistematica', 'Tradizionale'],
      rows.map((c) => [
        formulaHTML(c.formula),
        esc(c.stock.join(' / ') || '—'),
        esc(c.systematic.join(' / ') || '—'),
        esc(c.traditional.join(' / ') || '—'),
      ]),
    );
  };
  $('cl-nom-search').oninput = search;
  ask();
  build();
  search();
  return () => {};
}

function renderRedox(root, ctx) {
  shell(
    root,
    'Segui gli elettroni nel bilancio',
    'Calcola i numeri di ossidazione, confronta le specie e verifica che elettroni persi e acquistati coincidano.',
    `
    <div class="cl-form-row">${select('cl-red-rx', 'Reazione redox', reactionOptions(true))}${btn('cl-red-load', 'Carica esempio')}</div>${input('cl-red-equation', 'Equazione completa bilanciata', lookup('termite').equation, 'maxlength="1000"')}
    <div class="cl-actions">${btn('cl-red-check', 'Analizza il bilancio', true)}${btn('cl-red-solve', 'Bilancia atomi e carica')}</div><div id="cl-red-out" aria-live="polite"></div>
    <details class="cl-custom" open><summary>Calcola il numero di ossidazione di una specie</summary>${input('cl-ox-formula', 'Formula o ione', 'H2SO4', 'maxlength="200" placeholder="es. H2O2, MnO4^-, Fe3O4"')}<div id="cl-ox-out" aria-live="polite"></div></details>
    <details class="cl-custom"><summary>Semireazioni: esercizio in ambiente acido</summary><p><b>1. Riduzione:</b> MnO₄⁻ → Mn²⁺. Aggiungi 4 H₂O a destra, 8 H⁺ a sinistra e 5 e⁻ a sinistra.</p><p><b>2. Ossidazione:</b> Fe²⁺ → Fe³⁺ + e⁻. Moltiplica questa semireazione per 5.</p><p><b>3. Somma:</b> cancella 5 e⁻ su entrambi i lati. Nell’equazione completa verifica Mn, Fe, O, H e la carica.</p><p>Puoi inserire anche una singola semireazione con <code>e^-</code>: il motore controlla atomi e carica, e distingue gli elettroni esterni dal bilancio completo.</p>${btn('cl-half', 'Prova la semireazione')}</details>
    <p class="cl-model-note">Il n.o. è una convenzione, non una carica atomica calcolata. Nessuna animazione mostra traiettorie elettroniche: la densità quantistica resta nella tavola periodica e nel sandbox.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id),
    out = $('cl-red-out');
  $('cl-red-rx').value = 'termite';
  const check = () =>
    guarded(
      out,
      () => {
        const eq = getEq(root, 'cl-red-equation'),
          data = redoxLedger(eq.left, eq.right);
        out.innerHTML =
          `<div class="cl-equation">${equationHTML(eq)}</div>` +
          metrics([
            ['Elettroni ceduti', fmt(data.lost) + ' e⁻'],
            ['Elettroni acquistati', fmt(data.gained) + ' e⁻'],
            [
              'Bilancio',
              data.external
                ? 'Semireazione'
                : data.balanced
                  ? 'Conservato'
                  : 'Da verificare',
            ],
          ]) +
          table(
            [
              'Elemento',
              'n.o. iniziale',
              'n.o. finale',
              'Atomi coinvolti',
              'Variazione',
            ],
            data.changes.map((c) => [
              esc(c.el),
              fmt(c.from),
              fmt(c.to),
              fmt(c.atoms),
              `<span class="${c.electrons > 0 ? 'cl-error' : 'cl-ok'}">${c.electrons > 0 ? 'Ossidazione: cede' : 'Riduzione: acquista'} ${fmt(Math.abs(c.electrons))} e⁻</span>`,
            ]),
          );
        if (!data.changes.length)
          out.insertAdjacentHTML(
            'beforeend',
            notice(
              'Nessun numero di ossidazione cambia: non è una redox.',
              true,
            ),
          );
        if (data.external)
          out.insertAdjacentHTML(
            'beforeend',
            notice(
              'Semireazione con e⁻ espliciti: l’altra semireazione fornirà il bilancio elettronico complessivo.',
            ),
          );
        out.insertAdjacentHTML(
          'beforeend',
          '<p class="cl-hint">Passaggi: assegna i n.o. → elimina gli atomi che non cambiano stato → moltiplica Δn.o. per gli atomi coinvolti. Non attribuire lo stato dell’ossigeno elementare agli O già presenti negli altri reagenti.</p>',
        );
        ctx.snapshot = {
          tipo: 'redox',
          equazione: equationText(eq.left, eq.right),
          variazioni: data.changes,
          elettroniCeduti: data.lost,
          elettroniAcquistati: data.gained,
        };
      },
      () => (ctx.snapshot = null),
    );
  $('cl-red-check').onclick = check;
  $('cl-red-load').onclick = () => {
    $('cl-red-equation').value = lookup(val(root, 'cl-red-rx')).equation;
    check();
  };
  $('cl-red-rx').onchange = $('cl-red-load').onclick;
  $('cl-red-equation').oninput = () => {
    ctx.snapshot = null;
    out.innerHTML = notice(
      'Premi “Analizza il bilancio” per verificare la nuova equazione.',
    );
  };
  $('cl-red-solve').onclick = () =>
    guarded(out, () => {
      const eq = getEq(root, 'cl-red-equation'),
        sol = solveBalance(eq.left, eq.right);
      $('cl-red-equation').value = equationText(sol.left, sol.right);
      check();
    });
  $('cl-half').onclick = () => {
    $('cl-red-equation').value = 'MnO4^- + H^+ + e^- -> Mn^2+ + H2O';
    $('cl-red-solve').click();
  };
  const ox = () =>
    guarded($('cl-ox-out'), () => {
      const data = oxidation(val(root, 'cl-ox-formula'));
      $('cl-ox-out').innerHTML =
        table(
          ['Elemento', 'Atomi', 'n.o.', 'Contributo'],
          Object.entries(data.sp.atoms).map(([el, n]) => [
            esc(el),
            n,
            data.states[el] === undefined ? 'incognito' : fmt(data.states[el]),
            data.states[el] === undefined ? '—' : fmt(n * data.states[el]),
          ]),
        ) + notice(esc(data.note), !data.unresolved);
      if (data.unresolved) ctx.snapshot = null;
      else
        ctx.snapshot = {
          tipo: 'numero di ossidazione',
          formula: data.sp.formula,
          stati: data.states,
        };
    });
  $('cl-ox-formula').oninput = ox;
  ox();
  check();
  return () => {};
}

function renderSolutions(root, ctx) {
  shell(
    root,
    'Diluire senza perdere il soluto',
    'Modifica concentrazione e volumi: le moli restano costanti. Il solvente aggiunto è la differenza fra i volumi del modello.',
    `
    <div class="cl-form-row">${input('cl-sol-C', 'Concentrazione iniziale (mol/L)', 2, 'type="number" min="0" step="any"')}${input('cl-sol-Vi', 'Volume iniziale (mL)', 50, 'type="number" min="0.001" step="any"')}${input('cl-sol-Vf', 'Volume finale (mL)', 200, 'type="number" min="0.001" step="any"')}</div>
    ${visual('cl-sol-canvas', 'Stesso numero di marcatori del soluto nei due recipienti. Il modello assume soluzione omogenea e assenza di reazione o perdita di soluto.')}
    <div id="cl-sol-out" aria-live="polite"></div>
    <details class="cl-custom" open><summary>Preparare una soluzione da un solido</summary><div class="cl-form-row">${input('cl-prep-formula', 'Formula del soluto', 'NaCl', 'maxlength="200"')}${input('cl-prep-C', 'Concentrazione desiderata (mol/L)', 0.1, 'type="number" min="0" step="any"')}${input('cl-prep-V', 'Volume finale (mL)', 250, 'type="number" min="0.001" step="any"')}</div><div id="cl-prep-out" aria-live="polite"></div></details>
    <p class="cl-source-inline">L4 · pp. 18–23. Il volume è quello finale della soluzione, non il solo volume del solvente.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id);
  let model = null;
  const update = () =>
    guarded(
      $('cl-sol-out'),
      () => {
        const C = number(val(root, 'cl-sol-C'), { min: 0 }),
          Vi = number(val(root, 'cl-sol-Vi'), { min: 0.001 }),
          Vf = number(val(root, 'cl-sol-Vf'), { min: 0.001 }),
          out = dilution(C, Vi, Vf);
        model = { C, Vi, Vf, out };
        $('cl-sol-out').innerHTML =
          metrics([
            ['Concentrazione finale', fmt(out.finalConcentration) + ' mol/L'],
            ['Soluto conservato', fmt(out.mol) + ' mol'],
            ['Incremento di volume', fmt(out.addedMl) + ' mL'],
          ]) +
          `<div class="cl-calculation">c𝒇 = cᵢ × Vᵢ / V𝒇 = ${fmt(C)} × ${fmt(Vi)} / ${fmt(Vf)} = <b>${fmt(out.finalConcentration)} mol/L</b></div>`;
        ctx.snapshot = {
          tipo: 'diluizione',
          concentrazioneIniziale: C,
          volumeInizialeMl: Vi,
          volumeFinaleMl: Vf,
          ...out,
        };
      },
      () => {
        model = null;
        ctx.snapshot = null;
      },
    );
  const prep = () =>
    guarded(
      $('cl-prep-out'),
      () => {
        const M = molarMass(val(root, 'cl-prep-formula')),
          C = number(val(root, 'cl-prep-C'), { min: 0, max: 1000 }),
          V = number(val(root, 'cl-prep-V'), { min: 0.001, max: 1e9 });
        $('cl-prep-out').innerHTML =
          `<div class="cl-calculation">n = cV = ${fmt(C)} × ${fmt(V / 1000)} = ${fmt((C * V) / 1000)} mol<br>m = nM = ${fmt((C * V) / 1000)} × ${fmt(M)} = <b>${fmt(((C * V) / 1000) * M)} g</b></div><p class="cl-hint">Calcolo formale: verifica la solubilità reale del soluto e le condizioni prima di interpretarlo come una soluzione ottenibile.</p>`;
        ctx.snapshot = {
          tipo: 'preparazione soluzione',
          formula: val(root, 'cl-prep-formula'),
          concentrazione: C,
          volumeMl: V,
          massaG: ((C * V) / 1000) * M,
        };
      },
      () => (ctx.snapshot = null),
    );
  listen(root, '#cl-sol-C,#cl-sol-Vi,#cl-sol-Vf', update);
  listen(root, '#cl-prep-formula,#cl-prep-C,#cl-prep-V', prep);
  prep();
  update();
  const stop = canvasAnimation(
    $('cl-sol-canvas'),
    (c, t, p) => dilutionDiagram(c, t, p, model),
    { getPaused: () => ctx.paused },
  );
  return stop;
}

function renderMeasures(root, ctx) {
  shell(
    root,
    'Misure con le unità in vista',
    'Converti, misura la densità e verifica la precisione dichiarata. I risultati si aggiornano mentre modifichi i dati.',
    `
    <div class="cl-grid-two"><section class="cl-card"><h4>Temperatura</h4><div class="cl-form-row">${input('cl-temp', 'Temperatura', 25, 'type="number" step="any"')}${select('cl-temp-unit', 'Unità', '<option value="C">°C</option><option value="K">K</option><option value="F">°F</option>')}</div><div id="cl-temp-out" aria-live="polite"></div></section>
    <section class="cl-card"><h4>Densità</h4><div class="cl-form-row">${input('cl-density-m', 'Massa (g)', 54, 'type="number" min="0" step="any"')}${input('cl-density-v', 'Volume (mL)', 20, 'type="number" min="0.001" step="any"')}</div><div id="cl-density-out" aria-live="polite"></div></section>
    <section class="cl-card"><h4>Volumi</h4><div class="cl-form-row">${input('cl-vol', 'Valore', 250, 'type="number" min="0" step="any"')}${select('cl-vol-unit', 'Da', '<option value="mL">mL = cm³</option><option value="L">L</option><option value="m3">m³</option>')}</div><div id="cl-vol-out" aria-live="polite"></div></section>
    <section class="cl-card"><h4>Cifre significative</h4>${input('cl-significant', 'Numero scritto con la sua precisione', '0.00450', 'inputmode="decimal" maxlength="50"')}<div id="cl-significant-out" aria-live="polite"></div><p class="cl-hint">Gli zeri finali di un intero senza decimali sono ambigui: usa la notazione scientifica per dichiararli.</p></section></div>
    <details class="cl-custom"><summary>Precisione e accuratezza: una serie di misure</summary><div class="cl-form-row">${input('cl-stat-reference', 'Valore di riferimento (g)', 100, 'type="number" step="any"')}${input('cl-stat-bias', 'Errore sistematico (g)', 2, 'type="number" step="any"')}${input('cl-stat-spread', 'Dispersione del modello (g)', 0.3, 'type="number" min="0" step="any"')}</div><div id="cl-stat-out" aria-live="polite"></div><p class="cl-hint">Serie sintetica controllata, non dati sperimentali reali. Una piccola dispersione non elimina lo scostamento sistematico dal riferimento.</p></details>
    <details class="cl-custom"><summary>Calcola e arrotonda con la precisione dei dati</summary><div class="cl-form-row">${input('cl-calc-a', 'Primo dato', '1.20', 'maxlength="50"')}${select('cl-calc-op', 'Operazione', '<option value="multiply">×</option><option value="divide">÷</option><option value="add">+</option><option value="subtract">−</option>')}${input('cl-calc-b', 'Secondo dato', '2.3', 'maxlength="50"')}</div><div id="cl-calc-out" aria-live="polite"></div></details>
    <p class="cl-source-inline">L1 · pp. 26–36.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id);
  const temperature = () =>
    guarded(
      $('cl-temp-out'),
      () => {
        const x = number(val(root, 'cl-temp')),
          unit = val(root, 'cl-temp-unit'),
          K =
            unit === 'K'
              ? x
              : unit === 'C'
                ? x + 273.15
                : ((x - 32) * 5) / 9 + 273.15;
        if (K < 0)
          throw Error(
            'Temperatura inferiore allo zero assoluto: controlla il valore.',
          );
        $('cl-temp-out').innerHTML = metrics([
          ['Kelvin', fmt(K, 7) + ' K'],
          ['Celsius', fmt(K - 273.15, 7) + ' °C'],
          ['Fahrenheit', fmt(((K - 273.15) * 9) / 5 + 32, 7) + ' °F'],
        ]);
        ctx.snapshot = { tipo: 'conversione temperatura', kelvin: K };
      },
      () => (ctx.snapshot = null),
    );
  const density = () =>
    guarded(
      $('cl-density-out'),
      () => {
        const m = number(val(root, 'cl-density-m'), { min: 0 }),
          v = number(val(root, 'cl-density-v'), { min: 0.001 });
        $('cl-density-out').innerHTML =
          `<div class="cl-calculation">ρ = ${fmt(m)} / ${fmt(v)} = <b>${fmt(m / v)} g/mL</b><br>= ${fmt((1000 * m) / v)} kg/m³</div>`;
        ctx.snapshot = {
          tipo: 'densita',
          massaG: m,
          volumeMl: v,
          densitaGml: m / v,
        };
      },
      () => (ctx.snapshot = null),
    );
  const volume = () =>
    guarded(
      $('cl-vol-out'),
      () => {
        const x = number(val(root, 'cl-vol'), { min: 0 }),
          litres = x * { mL: 0.001, L: 1, m3: 1000 }[val(root, 'cl-vol-unit')];
        $('cl-vol-out').innerHTML = metrics([
          ['Litri', fmt(litres) + ' L'],
          ['Millilitri', fmt(litres * 1000) + ' mL'],
          ['Metri cubi', fmt(litres / 1000) + ' m³'],
        ]);
        ctx.snapshot = { tipo: 'conversione volume', litri: litres };
      },
      () => (ctx.snapshot = null),
    );
  const significant = () =>
    guarded(
      $('cl-significant-out'),
      () => {
        $('cl-significant-out').innerHTML = notice(
          `<b>${significantFigures(val(root, 'cl-significant'))} cifre significative.</b> La precisione dipende da come scrivi il dato.`,
          true,
        );
        ctx.snapshot = {
          tipo: 'cifre significative',
          dato: val(root, 'cl-significant'),
          cifre: significantFigures(val(root, 'cl-significant')),
        };
      },
      () => (ctx.snapshot = null),
    );
  const statistics = () =>
    guarded(
      $('cl-stat-out'),
      () => {
        const reference = number(val(root, 'cl-stat-reference')),
          bias = number(val(root, 'cl-stat-bias')),
          spread = number(val(root, 'cl-stat-spread'), { min: 0 }),
          values = [-1.5, -1.1, -0.8, -0.4, -0.1, 0.1, 0.4, 0.8, 1.1, 1.5].map(
            (z) => reference + bias + z * spread,
          ),
          data = measurementStatistics(values);
        $('cl-stat-out').innerHTML =
          metrics([
            ['Media', fmt(data.mean, 7) + ' g'],
            ['Scostamento dal riferimento', fmt(data.mean - reference) + ' g'],
            ['Deviazione standard campionaria', fmt(data.sd) + ' g'],
          ]) +
          table(
            ['Misura', 'Valore (g)'],
            values.map((v, i) => [i + 1, fmt(v, 7)]),
          );
        ctx.snapshot = {
          tipo: 'misure sintetiche',
          riferimento: reference,
          erroreSistematico: bias,
          misure: values,
          ...data,
        };
      },
      () => (ctx.snapshot = null),
    );
  const calculation = () =>
    guarded(
      $('cl-calc-out'),
      () => {
        const data = measuredCalculation(
          val(root, 'cl-calc-a'),
          val(root, 'cl-calc-b'),
          val(root, 'cl-calc-op'),
        );
        $('cl-calc-out').innerHTML =
          metrics([
            ['Risultato non arrotondato', fmt(data.value, 12)],
            ['Risultato da riportare', esc(data.reported.replace('.', ','))],
          ]) + notice(esc(data.rule), true);
        ctx.snapshot = {
          tipo: 'calcolo con precisione',
          datoA: val(root, 'cl-calc-a'),
          datoB: val(root, 'cl-calc-b'),
          operazione: val(root, 'cl-calc-op'),
          ...data,
        };
      },
      () => (ctx.snapshot = null),
    );
  listen(root, '#cl-stat-reference,#cl-stat-bias,#cl-stat-spread', statistics);
  listen(root, '#cl-calc-a,#cl-calc-op,#cl-calc-b', calculation);
  statistics();
  calculation();
  listen(root, '#cl-temp,#cl-temp-unit', temperature);
  listen(root, '#cl-density-m,#cl-density-v', density);
  listen(root, '#cl-vol,#cl-vol-unit', volume);
  listen(root, '#cl-significant', significant);
  temperature();
  density();
  volume();
  significant();
  return () => {};
}

function renderMole(root, ctx) {
  shell(
    root,
    'Dalla formula alle entità',
    'Calcola qualsiasi formula riconosciuta: composizione, massa molare, moli e numero di entità. Per i sali le entità sono unità formula.',
    `
    <div class="cl-form-row">${input('cl-mole-formula', 'Formula', 'Ca(OH)2', 'maxlength="200"')}${input('cl-mole-mass', 'Massa (g)', 7.4092, 'type="number" min="0" max="1000000000" step="any"')}</div><div id="cl-mole-out" aria-live="polite"></div>
    <details class="cl-custom"><summary>Trova una formula empirica</summary><p>Inserisci masse o percentuali su una stessa base. Il controllo cerca piccoli rapporti interi con scarto inferiore al 2%: verifica sempre la compatibilità con la precisione dei dati.</p><div class="cl-form-row">${input('cl-emp-C', 'Massa C (g)', 40, 'type="number" min="0" step="any"')}${input('cl-emp-H', 'Massa H (g)', 6.71, 'type="number" min="0" step="any"')}${input('cl-emp-O', 'Massa O (g)', 53.29, 'type="number" min="0" step="any"')}</div><div id="cl-emp-out" aria-live="polite"></div></details>
    <details class="cl-custom"><summary>Analisi di combustione di un composto C/H/O</summary><div class="cl-form-row">${input('cl-comb-sample', 'Massa del campione (g)', 0.3, 'type="number" min="0" step="any"')}${input('cl-comb-co2', 'CO₂ ottenuta (g)', 0.44, 'type="number" min="0" step="any"')}${input('cl-comb-h2o', 'H₂O ottenuta (g)', 0.18, 'type="number" min="0" step="any"')}</div><div id="cl-comb-out" aria-live="polite"></div></details>
    <p class="cl-source-inline">L2 · pp. 31–33; L3 · pp. 28–33. <a href="https://www.bipm.org/en/si-base-units/mole" target="_blank" rel="noopener">Definizione SI attuale della mole · BIPM</a>.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id);
  const update = () =>
    guarded(
      $('cl-mole-out'),
      () => {
        const sp = parseFormula(val(root, 'cl-mole-formula')),
          M = molarMass(sp),
          m = number(val(root, 'cl-mole-mass'), { min: 0, max: 1e9 }),
          n = m / M;
        $('cl-mole-out').innerHTML =
          `<div class="cl-equation">${formulaHTML(sp)}</div>` +
          metrics([
            ['Massa molare', fmt(M, 7) + ' g/mol'],
            ['Quantità di sostanza', fmt(n, 7) + ' mol'],
            ['Entità specificate', (n * NA).toExponential(5)],
          ]) +
          table(
            [
              'Elemento',
              'Atomi per entità',
              'Contributo (g/mol)',
              'Frazione in massa',
            ],
            Object.entries(sp.atoms).map(([el, c]) => [
              esc(el),
              c,
              fmt(molarMass(el) * c, 7),
              fmt(((molarMass(el) * c) / M) * 100, 5) + '%',
            ]),
          ) +
          `<div class="cl-calculation">n = m/M = ${fmt(m, 7)} / ${fmt(M, 7)} = ${fmt(n, 7)} mol<br>N = n × 6,02214076 × 10²³ = <b>${(n * NA).toExponential(5)} entità</b></div>`;
        ctx.snapshot = {
          tipo: 'mole',
          formula: sp.formula,
          massaG: m,
          massaMolare: M,
          moli: n,
          entita: n * NA,
        };
      },
      () => (ctx.snapshot = null),
    );
  const empirical = () =>
    guarded(
      $('cl-emp-out'),
      () => {
        const entries = ['C', 'H', 'O'].map((el) => ({
            el,
            mass: val(root, 'cl-emp-' + el),
          })),
          data = empiricalFormula(entries);
        $('cl-emp-out').innerHTML =
          metrics([
            ['Formula empirica', formulaHTML(data.formula)],
            ['Scarto dei rapporti', fmt(data.error * 100) + '%'],
          ]) +
          table(
            ['Elemento', 'Moli', 'Rapporto n/n minimo'],
            entries.map((e, i) => [
              e.el,
              fmt(data.mol[i]),
              fmt(data.ratios[i]),
            ]),
          );
        ctx.snapshot = { tipo: 'formula empirica', dati: entries, ...data };
      },
      () => (ctx.snapshot = null),
    );
  const combustion = () =>
    guarded(
      $('cl-comb-out'),
      () => {
        const sample = number(val(root, 'cl-comb-sample'), { min: 1e-12 }),
          co2 = number(val(root, 'cl-comb-co2'), { min: 1e-12 }),
          water = number(val(root, 'cl-comb-h2o'), { min: 1e-12 });
        const mC = (co2 / molarMass('CO2')) * molarMass('C'),
          mH = (water / molarMass('H2O')) * 2 * molarMass('H'),
          mO = sample - mC - mH;
        if (mO <= 0)
          throw Error(
            'La massa residua di O non è positiva. Questo calcolo vale per un composto contenente C, H e O: verifica i dati.',
          );
        const data = empiricalFormula([
          { el: 'C', mass: mC },
          { el: 'H', mass: mH },
          { el: 'O', mass: mO },
        ]);
        $('cl-comb-out').innerHTML =
          metrics([
            ['C nel campione', fmt(mC) + ' g'],
            ['H nel campione', fmt(mH) + ' g'],
            ['O per differenza', fmt(mO) + ' g'],
            ['Formula empirica', formulaHTML(data.formula)],
          ]) +
          `<p class="cl-hint">Ipotesi: combustione completa, composto contenente soltanto C/H/O, prodotti raccolti quantitativamente. Una formula molecolare richiede anche la massa molare del composto.</p>`;
        ctx.snapshot = {
          tipo: 'analisi combustione',
          massaCampioneG: sample,
          massaCO2G: co2,
          massaH2OG: water,
          formulaEmpirica: data.formula,
        };
      },
      () => (ctx.snapshot = null),
    );
  listen(root, '#cl-mole-formula,#cl-mole-mass', update);
  listen(root, '#cl-emp-C,#cl-emp-H,#cl-emp-O', empirical);
  listen(root, '#cl-comb-sample,#cl-comb-co2,#cl-comb-h2o', combustion);
  empirical();
  combustion();
  update();
  return () => {};
}

function renderAtoms(root, ctx) {
  shell(
    root,
    'Atomi: composizione, storia e modello',
    'Costruisci il conteggio di un nuclide, cambia una miscela isotopica e confronta i modelli storici con la simulazione quantistica già presente.',
    `
    <div class="cl-grid-two"><section class="cl-card"><h4>Protoni, neutroni, elettroni</h4><div class="cl-form-row">${input('cl-nuc-Z', 'Z', 17, 'type="number" min="1" max="118" step="1"')}${input('cl-nuc-A', 'A', 37, 'type="number" min="1" max="350" step="1"')}${input('cl-nuc-q', 'Carica q', -1, 'type="number" step="1" min="-20" max="118"')}</div><div id="cl-nuc-out" aria-live="polite"></div><p class="cl-hint">Il conteggio non verifica che il nuclide esista o sia stabile.</p></section>
    <section class="cl-card"><h4>Il ponte verso la quantistica</h4><p>Gli orbitali sono densità di probabilità. Esplorali nel motore DFT-LDA già presente, con nuvole, superfici e dinamica coerente dell’idrogeno.</p><div class="cl-actions"><button type="button" class="cl-btn primary" data-atom="H">Apri l’idrogeno</button><button type="button" class="cl-btn" data-atom="C">Apri il carbonio</button><button type="button" class="cl-btn" data-open-sandbox>Apri il sandbox</button></div><p class="cl-hint">Le orbite classiche dei modelli storici non descrivono gli elettroni reali.</p></section></div>
    <details class="cl-custom" open><summary>Miscela isotopica del cloro</summary><div class="cl-form-row">${input('cl-iso-m1', 'Massa ³⁵Cl usata (u)', 35, 'type="number" min="0.001" step="any"')}${input('cl-iso-m2', 'Massa ³⁷Cl usata (u)', 37, 'type="number" min="0.001" step="any"')}<label class="cl-field cl-grow" for="cl-iso-x"><span>Abbondanza di ³⁵Cl (%)</span><input id="cl-iso-x" type="range" min="0" max="100" step="0.1" value="77.5"></label></div>${visual('cl-iso-canvas', 'I marcatori visualizzano l’abbondanza; i numeri 35 e 37 sono masse approssimate, come nell’esempio della lezione.')}<div id="cl-iso-out" aria-live="polite"></div></details>
    <details class="cl-custom"><summary>Esplora i modelli storici in 3D</summary><div class="cl-actions">${btn('cl-dalton', 'Dalton')}${btn('cl-thomson', 'Thomson')}${btn('cl-rutherford', 'Rutherford')}${btn('cl-isomeri', 'Isomeri C₂H₆O')}${btn('cl-model-pause', 'Ⅱ Pausa')}</div><div id="cl-history-stage"></div><div id="cl-history-out"></div></details>
    <details class="cl-custom"><summary>Esperimento di Rutherford: parametri liberi</summary><div class="cl-form-row">${input('cl-ru-Z', 'Z del bersaglio', 79, 'type="number" min="1" max="118" step="1"')}${input('cl-ru-E', 'Energia α (MeV)', 5, 'type="number" min="0.01" max="1000" step="any"')}${input('cl-ru-b', 'Parametro d’urto b (fm)', 20, 'type="number" min="0" max="10000" step="any"')}</div>${visual('cl-ru-canvas', 'Traiettoria α ottenuta dalla soluzione iperbolica del potenziale coulombiano repulsivo; animazione rallentata con tempo parametrico fisico.')}<div id="cl-ru-out" aria-live="polite"></div></details>
    <p class="cl-source-inline">L2 · pp. 4–24. I modelli storici sono didattici; lo scattering usa un nucleo puntiforme fisso e trascura interazioni nucleari e schermaggio.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id);
  let isotope = null,
    scattering = null,
    stage = null,
    paused = false,
    dead = false;
  const nucleus = () =>
    guarded(
      $('cl-nuc-out'),
      () => {
        const Z = number(val(root, 'cl-nuc-Z'), {
            min: 1,
            max: 118,
            integer: true,
          }),
          A = number(val(root, 'cl-nuc-A'), {
            min: Z,
            max: 350,
            integer: true,
          }),
          q = number(val(root, 'cl-nuc-q'), {
            min: -20,
            max: Z,
            integer: true,
          }),
          el = ELEMENTS.find((e) => e.Z === Z);
        $('cl-nuc-out').innerHTML = metrics([
          [
            'Nuclide formale',
            `<sup>${A}</sup><sub>${Z}</sub>${esc(el.symbol)}`,
          ],
          ['Protoni', Z],
          ['Neutroni', A - Z],
          ['Elettroni', Z - q],
        ]);
        ctx.snapshot = {
          tipo: 'conteggio nuclide',
          Z,
          A,
          carica: q,
          neutroni: A - Z,
          elettroni: Z - q,
        };
      },
      () => (ctx.snapshot = null),
    );
  const isotopes = () =>
    guarded(
      $('cl-iso-out'),
      () => {
        const x = number(val(root, 'cl-iso-x'), { min: 0, max: 100 }),
          m1 = number(val(root, 'cl-iso-m1'), { min: 0.001 }),
          m2 = number(val(root, 'cl-iso-m2'), { min: 0.001 }),
          mean = isotopeMean(m1, m2, x);
        isotope = x;
        $('cl-iso-out').innerHTML =
          metrics([
            ['Abbondanza ³⁵Cl', fmt(x) + '%'],
            ['Abbondanza ³⁷Cl', fmt(100 - x) + '%'],
            ['Massa media', fmt(mean, 7) + ' u'],
          ]) +
          `<div class="cl-calculation">${fmt(m1)} × ${fmt(x / 100)} + ${fmt(m2)} × ${fmt(1 - x / 100)} = <b>${fmt(mean, 7)} u</b></div>`;
        ctx.snapshot = {
          tipo: 'miscela isotopica',
          masseU: [m1, m2],
          percentuale35: x,
          mediaU: mean,
        };
      },
      () => {
        isotope = null;
        ctx.snapshot = null;
      },
    );
  const scatter = () =>
    guarded(
      $('cl-ru-out'),
      () => {
        scattering = rutherford({
          Z: val(root, 'cl-ru-Z'),
          E: val(root, 'cl-ru-E'),
          b: val(root, 'cl-ru-b'),
        });
        $('cl-ru-out').innerHTML =
          metrics([
            [
              'Angolo asintotico',
              fmt((scattering.theta * 180) / Math.PI) + '°',
            ],
            ['Minima distanza', fmt(scattering.closest) + ' fm'],
          ]) +
          `<div class="cl-calculation">θ = 2 arctan[Z₁Z₂ k / (2Eb)], con Z₁ = 2 e k = 1,439964548 MeV·fm.</div><p class="cl-hint">La formula è classica. A distanze comparabili al raggio nucleare il modello puntiforme non descrive le interazioni nucleari.</p>`;
        ctx.snapshot = {
          tipo: 'scattering Rutherford',
          Z: Number(val(root, 'cl-ru-Z')),
          energiaMeV: Number(val(root, 'cl-ru-E')),
          bFm: Number(val(root, 'cl-ru-b')),
          angoloGradi: (scattering.theta * 180) / Math.PI,
        };
      },
      () => {
        scattering = null;
        ctx.snapshot = null;
      },
    );
  const historical = (m) =>
    guarded($('cl-history-out'), () => {
      if (!stage)
        stage = makeStage($('cl-history-stage'), {
          getPaused: () => ctx.paused,
        });
      const captions = {
        isomeri:
          'A sinistra etanolo (C−C−O−H), a destra dimetiletere (C−O−C): stessa formula C₂H₆O, connettività diversa. Modelli a sfere e bastoncini schematici: colori convenzionali e geometrie non ottimizzate, non densità quantistiche.',
        dalton:
          'Dalton: sfere senza struttura interna. Diversi elementi sono distinguibili nel modello; isotopi e particelle subatomiche ne superano i limiti.',
        thomson:
          'Thomson: elettroni in un fondo positivo diffuso. È un modello storico, con posizioni schematiche statiche; non è una densità quantistica.',
        rutherford:
          'Rutherford: nucleo e scattering α coulombiano. Il nucleo è ingrandito; la traiettoria è calcolata, non una curva disegnata a mano.',
      };
      if (m === 'dalton') showDalton(stage);
      if (m === 'thomson') showThomson(stage);
      if (m === 'rutherford') showRutherford(stage);
      if (m === 'isomeri') showIsomers(stage);
      $('cl-history-out').innerHTML = notice(captions[m], true);
      ctx.snapshot = { tipo: 'modello storico', modello: m };
    });
  ['dalton', 'thomson', 'rutherford', 'isomeri'].forEach(
    (m) => ($('cl-' + m).onclick = () => historical(m)),
  );
  $('cl-model-pause').onclick = () => {
    paused = !paused;
    stage?.setPaused(paused);
    $('cl-model-pause').textContent = paused ? '▶ Riprendi' : 'Ⅱ Pausa';
  };
  const historyDetails = $('cl-history-stage').closest('details');
  historyDetails.addEventListener('toggle', () => {
    if (dead) return;
    if (historyDetails.open && !stage) historical('dalton');
    else if (!historyDetails.open) {
      stage?.dispose();
      stage = null;
    }
  });
  root
    .querySelectorAll('[data-atom]')
    .forEach((b) => (b.onclick = () => ctx.openAtom(b.dataset.atom)));
  root.querySelector('[data-open-sandbox]').onclick = () => ctx.openSandbox();
  listen(root, '#cl-nuc-Z,#cl-nuc-A,#cl-nuc-q', nucleus);
  listen(root, '#cl-iso-x,#cl-iso-m1,#cl-iso-m2', isotopes);
  listen(root, '#cl-ru-Z,#cl-ru-E,#cl-ru-b', scatter);
  nucleus();
  scatter();
  isotopes();
  const stop1 = canvasAnimation(
      $('cl-iso-canvas'),
      (c, t, p) => {
        if (isotope !== null) isotopeDiagram(c, t, p, isotope);
      },
      { getPaused: () => ctx.paused },
    ),
    stop2 = canvasAnimation(
      $('cl-ru-canvas'),
      (c, t, p) => scatterDiagram(c, t, p, scattering),
      { getPaused: () => ctx.paused },
    );
  return () => {
    dead = true;
    stage?.dispose();
    stage = null;
    stop1();
    stop2();
  };
}

function renderExercises(root, ctx) {
  shell(
    root,
    'Allenati, poi spiega il risultato',
    'Risposte aperte, calcoli e formule. Il feedback resta visibile finché decidi tu di cambiare esercizio.',
    `
    <div class="cl-form-row">${select('cl-question-unit', 'Argomento', `<option value="all">Tutti gli argomenti</option>${ctx.units.map((u) => `<option value="${u.id}">${esc(u.title)}</option>`).join('')}`)}${btn('cl-question-new', 'Nuovo esercizio')}</div>
    <div class="cl-question"><span id="cl-question-ref" class="cl-source-inline"></span><h4 id="cl-question-text"></h4>${input('cl-question-answer', 'Risposta (solo valore numerico quando richiesto)', '', 'autocomplete="off" maxlength="200"')}<div class="cl-actions">${btn('cl-question-check', 'Verifica', true)}${btn('cl-question-reveal', 'Mostra spiegazione')}</div><div id="cl-question-feedback" aria-live="polite"></div></div><div id="cl-question-stats"></div>`,
  );
  const $ = (id) => root.querySelector('#' + id);
  $('cl-question-unit').value = ctx.unit;
  let index = ctx.shared.questionIndex || 0,
    q,
    key,
    lastAnswer = '',
    revealed = false;
  const generated = (unit) => {
    const seed = index + 1;
    if (unit === 'soluzioni') {
      const C = 1 + (seed % 4),
        Vi = 25 * (1 + (seed % 5)),
        factor = 2 + (seed % 5);
      return {
        unit,
        q: `Diluisci ${Vi} mL a ${C} mol/L fino a ${Vi * factor} mL. Concentrazione finale in mol/L?`,
        answer: C / factor,
        tolerance: 0.0005,
        why: `c𝒇 = cᵢVᵢ/V𝒇 = ${C} × ${Vi} / ${Vi * factor} = ${fmt(C / factor, 7)} mol/L.`,
        ref: 'Esercizio generato · L4 p. 21',
      };
    }
    if (unit === 'misure') {
      const m = 15 + (seed % 41),
        v = 5 + (seed % 17);
      return {
        unit,
        q: `Massa ${m} g, volume ${v} mL. Calcola la densità in g/mL (almeno 4 cifre significative).`,
        answer: m / v,
        tolerance: 0.001,
        why: `ρ = m/V = ${m}/${v} = ${fmt(m / v, 7)} g/mL.`,
        ref: 'Esercizio generato · L1 pp. 29–30',
      };
    }
    if (unit === 'mole') {
      const mol = 1 + (seed % 9);
      return {
        unit,
        q: `${mol * 18} g di una sostanza con M = 18 g/mol corrispondono a quante moli?`,
        answer: mol,
        tolerance: 0,
        why: `n = m/M = ${mol * 18}/18 = ${mol} mol.`,
        ref: 'Esercizio generato · L2 pp. 31–33',
      };
    }
    if (unit === 'reazioni') {
      const h = 1 + (seed % 12),
        o = 1 + ((seed * 3) % 9),
        product = Math.min(h, 2 * o);
      return {
        unit,
        q: `2 H₂ + O₂ → 2 H₂O. Con ${h} mol H₂ e ${o} mol O₂, quante moli di H₂O puoi ottenere teoricamente?`,
        answer: product,
        tolerance: 0,
        why: `ξ massimo = min(${h}/2, ${o}/1); n(H₂O) = 2ξ = ${product} mol.`,
        ref: 'Esercizio generato · L4 p. 14',
      };
    }
    return null;
  };
  const ask = () => {
    const unit = val(root, 'cl-question-unit'),
      list = QUESTIONS.filter((x) => unit === 'all' || x.unit === unit);
    q =
      index % 3 === 2
        ? generated(
            unit === 'all'
              ? ['soluzioni', 'misure', 'mole', 'reazioni'][index % 4]
              : unit,
          )
        : null;
    q = q || list[index % list.length];
    key = 'exercise:' + q.unit + ':' + q.q;
    lastAnswer = '';
    revealed = false;
    $('cl-question-text').textContent = q.q;
    $('cl-question-ref').textContent = q.ref;
    $('cl-question-answer').value = '';
    $('cl-question-feedback').innerHTML = '';
    stats();
    ctx.snapshot = null;
  };
  const stats = () => {
    $('cl-question-stats').innerHTML = metrics([
      ['Tentativi', ctx.progress.attempts || 0],
      [
        'Esercizi risolti distinti',
        Object.keys(ctx.progress.solved || {}).length,
      ],
      [
        'Lezioni segnate come studiate',
        Object.keys(ctx.progress.studied || {}).length,
      ],
    ]);
  };
  $('cl-question-check').onclick = () => {
    const answer = val(root, 'cl-question-answer').trim();
    if (!answer) {
      $('cl-question-feedback').innerHTML = notice(
        'Inserisci prima la tua risposta.',
      );
      return;
    }
    if (answer === lastAnswer) return;
    lastAnswer = answer;
    let good = false;
    try {
      good =
        typeof q.answer === 'number'
          ? Math.abs(number(answer) - q.answer) <= (q.tolerance || 0)
          : q.formula
            ? atomsKey(answer) === atomsKey(q.answer)
            : [q.answer, ...(q.aliases || [])].some(
                (x) => normalizeName(x) === normalizeName(answer),
              );
    } catch {}
    const expected =
      typeof q.answer === 'number' ? fmt(q.answer, 9) : esc(q.answer);
    $('cl-question-feedback').innerHTML = notice(
      (good ? 'Corretto. ' : 'Da rivedere. Risposta: ' + expected + '. ') +
        esc(q.why) +
        (revealed
          ? ' (Spiegazione già consultata: questo tentativo non aumenta i progressi.)'
          : ''),
      good,
    );
    if (!revealed) ctx.recordAttempt(key, good, q.unit);
    stats();
    ctx.snapshot = {
      tipo: 'esercizio',
      domanda: q.q,
      risposta: answer,
      corretta: good,
      spiegazione: q.why,
    };
  };
  $('cl-question-answer').onkeydown = (e) => {
    if (e.key === 'Enter') $('cl-question-check').click();
  };
  $('cl-question-reveal').onclick = () => {
    revealed = true;
    $('cl-question-feedback').innerHTML = notice(esc(q.why));
  };
  $('cl-question-new').onclick = () => {
    index++;
    ctx.shared.questionIndex = index;
    ask();
  };
  $('cl-question-unit').onchange = () => {
    index = 0;
    ask();
  };
  ask();
  return () => {};
}

function renderMatter(root, ctx) {
  shell(
    root,
    'Materia: osserva prima di classificare',
    'Esplora le fasi dell’acqua e una separazione cromatografica. Gli schemi delle particelle rappresentano molecole, non elettroni.',
    `
    <section class="cl-card"><h4>Acqua a pressione atmosferica</h4><label class="cl-field" for="cl-phase-T"><span>Temperatura (°C)</span><input id="cl-phase-T" type="range" min="-20" max="120" step="1" value="20"></label>${visual('cl-phase-canvas', 'Schema cinetico classico: posizioni e moti illustrativi delle molecole. Non è una dinamica molecolare calcolata. Condizioni di riferimento: acqua pura a 1 atm; soglie 0 e 100 °C arrotondate, proporzioni di coesistenza illustrative.')}<div id="cl-phase-out" aria-live="polite"></div></section>
    <details class="cl-custom" open><summary>Cromatografia: due soluti, due affinità</summary><div class="cl-form-row">${input('cl-chrom-K1', 'K del soluto A', 1, 'type="number" min="0" max="10000" step="any"')}${input('cl-chrom-K2', 'K del soluto B', 4, 'type="number" min="0" max="10000" step="any"')}${input('cl-chrom-ratio', 'Rapporto delle fasi Vₛ/Vₘ', 1, 'type="number" min="0.001" max="10000" step="any"')}<label class="cl-field" for="cl-chrom-front"><span>Avanzamento del fronte (%)</span><input id="cl-chrom-front" type="range" min="0" max="100" step="1" value="75"></label></div>${visual('cl-chrom-canvas', 'Modello di ripartizione lineare all’equilibrio: K = cₛ/cₘ e Rf = 1/(1 + K·Vₛ/Vₘ). Non include allargamento delle bande, adsorbimento non lineare o cinetica.')}<div id="cl-chrom-out" aria-live="polite"></div></details>
    <section class="cl-question"><h4>Classifica un campione</h4><div class="cl-form-row">${select('cl-matter-sample', 'Campione', options(['Ossigeno O2', 'Acqua distillata', 'Acqua salata limpida', 'Acqua e sabbia', 'Latte', 'Aria filtrata']))}${select('cl-matter-class', 'La tua classificazione', options(['Sostanza elementare', 'Sostanza composta', 'Miscela omogenea', 'Miscela eterogenea']))}</div>${btn('cl-matter-check', 'Verifica la classificazione', true)}<div id="cl-matter-out" aria-live="polite"></div></section><p class="cl-source-inline">L1 · pp. 11–25; la ripartizione cromatografica è un modello aggiuntivo dell’esempio di p. 22.</p>`,
  );
  const $ = (id) => root.querySelector('#' + id);
  let phase = null,
    chrom = null,
    lastAnswer = '';
  const phases = () =>
    guarded($('cl-phase-out'), () => {
      const T = number(val(root, 'cl-phase-T'), { min: -20, max: 120 }),
        state = waterPhase(T);
      phase = { T, state };
      $('cl-phase-out').innerHTML =
        metrics([
          ['Temperatura', fmt(T) + ' °C'],
          ['Fase di riferimento', esc(state)],
          ['Trasformazione', 'Fisica'],
        ]) +
        `<p class="cl-hint">La specie resta H₂O. Alle temperature di transizione possono coesistere due fasi; per determinare la loro proporzione serve anche il bilancio energetico. Le particelle del diagramma non sono in scala di densità o velocità.</p>`;
      ctx.snapshot = {
        tipo: 'fase acqua',
        temperaturaC: T,
        fase: state,
        pressioneAtm: 1,
      };
    });
  const chromatographyUpdate = () =>
    guarded(
      $('cl-chrom-out'),
      () => {
        const K1 = number(val(root, 'cl-chrom-K1'), { min: 0, max: 10000 }),
          K2 = number(val(root, 'cl-chrom-K2'), { min: 0, max: 10000 }),
          ratio = number(val(root, 'cl-chrom-ratio'), {
            min: 0.001,
            max: 10000,
          }),
          front =
            number(val(root, 'cl-chrom-front'), { min: 0, max: 100 }) / 100,
          A = chromatography(K1, ratio, front),
          B = chromatography(K2, ratio, front);
        chrom = { A, B, front };
        $('cl-chrom-out').innerHTML =
          metrics([
            ['Rf soluto A', fmt(A.Rf)],
            ['Rf soluto B', fmt(B.Rf)],
            [
              'Separazione fra i centri',
              fmt(Math.abs(A.fractionAlongPaper - B.fractionAlongPaper) * 100) +
                '% della corsa',
            ],
          ]) +
          `<p class="cl-hint">Maggiore K significa maggiore affinità per la fase stazionaria e minore migrazione. A K uguali il modello non separa i due soluti.</p>`;
        ctx.snapshot = {
          tipo: 'cromatografia',
          K1,
          K2,
          rapportoFasi: ratio,
          fronte: front,
          RfA: A.Rf,
          RfB: B.Rf,
        };
      },
      () => {
        chrom = null;
        ctx.snapshot = null;
      },
    );
  const classifications = [
    ['Sostanza elementare', 'Una sola specie elementare: O₂.'],
    ['Sostanza composta', 'La sostanza è H₂O, formata da elementi diversi.'],
    [
      'Miscela omogenea',
      'Soluto e solvente sono uniformi alla scala macroscopica.',
    ],
    [
      'Miscela eterogenea',
      'Solido e liquido costituiscono fasi distinguibili.',
    ],
    [
      'Miscela eterogenea',
      'Il latte è un colloide: l’apparente uniformità non lo rende una soluzione molecolare.',
    ],
    ['Miscela omogenea', 'Miscela gassosa omogenea di più sostanze.'],
  ];
  $('cl-matter-check').onclick = () => {
    const sample = val(root, 'cl-matter-sample'),
      answer = val(root, 'cl-matter-class'),
      entry = classifications[$('cl-matter-sample').selectedIndex],
      good = answer === entry[0];
    $('cl-matter-out').innerHTML = notice(
      (good ? 'Corretto. ' : 'La categoria attesa è ' + entry[0] + '. ') +
        entry[1],
      good,
    );
    const key = sample + answer;
    if (key !== lastAnswer)
      ctx.recordAttempt('materia:' + sample, good, 'materia');
    lastAnswer = key;
  };
  listen(root, '#cl-phase-T', phases);
  listen(
    root,
    '#cl-chrom-K1,#cl-chrom-K2,#cl-chrom-ratio,#cl-chrom-front',
    chromatographyUpdate,
  );
  phases();
  chromatographyUpdate();
  const stop1 = canvasAnimation(
    $('cl-phase-canvas'),
    (c, t, p) => {
      if (!phase) return;
      c.strokeStyle = p.line;
      c.lineWidth = 2;
      c.strokeRect(220, 55, 520, 230);
      for (let i = 0; i < 48; i++) {
        const col = i % 8,
          row = Math.floor(i / 8),
          gas =
            phase.state === 'gas' ||
            (phase.state === 'liquido + gas' && i % 2 === 0),
          solid =
            phase.state === 'solido' ||
            (phase.state === 'solido + liquido' && i % 2 === 0),
          moving = !solid,
          x = gas
            ? 240 + ((i * 0.618 + t * 0.024) % 1) * 480
            : 330 +
              col * 42 +
              (moving ? Math.sin(t * 0.6 + i) * 9 : Math.sin(t * 1.5 + i) * 2),
          y = gas
            ? 75 + ((i * 0.414 + t * 0.017) % 1) * 190
            : solid
              ? 70 + row * 34
              : 90 +
                row * 30 +
                (moving
                  ? Math.cos(t * 0.5 + i) * 8
                  : Math.cos(t * 1.5 + i) * 2);
        c.fillStyle = '#5ab8e7';
        c.beginPath();
        c.arc(x, y, 5, 0, Math.PI * 2);
        c.fill();
      }
      c.fillStyle = p.text;
      c.font = '22px sans-serif';
      c.textAlign = 'center';
      c.fillText(phase.state + ' · ' + phase.T + ' °C', 480, 330);
    },
    { getPaused: () => ctx.paused },
  );
  const stop2 = canvasAnimation(
    $('cl-chrom-canvas'),
    (c, t, p) => {
      if (!chrom) return;
      c.strokeStyle = p.line;
      c.lineWidth = 2;
      c.strokeRect(90, 55, 780, 230);
      c.strokeStyle = p.muted;
      c.setLineDash([5, 5]);
      c.beginPath();
      c.moveTo(120 + chrom.front * 710, 55);
      c.lineTo(120 + chrom.front * 710, 285);
      c.stroke();
      c.setLineDash([]);
      for (const [entry, color, y, label] of [
        [chrom.A, '#f0b429', 125, 'A'],
        [chrom.B, '#b4a0ed', 225, 'B'],
      ]) {
        const x = 120 + entry.fractionAlongPaper * 710;
        c.fillStyle = color;
        c.beginPath();
        c.ellipse(x, y, 10, 28, 0, 0, 2 * Math.PI);
        c.fill();
        c.font = '18px sans-serif';
        c.fillText(label, x + 20, y + 7);
      }
      c.fillStyle = p.muted;
      c.font = '16px sans-serif';
      c.fillText('origine', 95, 320);
      c.fillText('fronte solvente', 600, 320);
    },
    { getPaused: () => ctx.paused },
  );
  return () => {
    stop1();
    stop2();
  };
}

export function mountTool(root, id, ctx) {
  const handlers = {
    materia: renderMatter,
    bilancio: renderBalance,
    stechiometria: renderStoichiometry,
    esperimenti: (r, c) => renderStoichiometry(r, c, true),
    nomenclatura: renderNomenclature,
    redox: renderRedox,
    soluzioni: renderSolutions,
    misure: renderMeasures,
    mole: renderMole,
    atomi: renderAtoms,
    esercizi: renderExercises,
  };
  ctx.snapshot = null;
  return (handlers[id] || renderBalance)(root, ctx);
}
