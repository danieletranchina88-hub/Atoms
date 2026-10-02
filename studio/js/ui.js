// Componenti riusati da tutti gli argomenti: progressi, classificatori, quiz, palestra di esercizi, soluzioni passo passo.

import { grade } from './numeri.js';

const KEY = 'studio-chimica-generale-v1';

export function loadProgress() {
  try { return JSON.parse(localStorage.getItem(KEY)) ?? {}; } catch { return {}; }
}
export function saveProgress(update) {
  const s = loadProgress();
  update(s);
  try { localStorage.setItem(KEY, JSON.stringify(s)); } catch { /* memoria del browser non disponibile */ }
  return s;
}
export const topicProgress = (id) => loadProgress()[id] ?? {};

export const rand = (a, b) => a + Math.random() * (b - a);
export const randInt = (a, b) => Math.floor(rand(a, b + 1));
export const pick = (arr) => arr[Math.floor(Math.random() * arr.length)];
export function shuffle(arr) {
  const a = [...arr];
  for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(Math.random() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
  return a;
}
const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// Classificatore: per ogni voce si sceglie una categoria e si legge subito il perché.
export function sorter(root, { items, options, label }) {
  let order = shuffle(items), score = 0, done = 0;
  const draw = () => {
    score = 0; done = 0;
    root.innerHTML = `
      <div class="sorter-head"><span class="sorter-score" aria-live="polite"></span>
        <button type="button" class="btn-quiet" data-reset>Ricomincia</button></div>
      <ol class="sorter-list">${order.map((it, i) => `
        <li class="sorter-item" data-i="${i}">
          <span class="sorter-name">${it.name}</span>
          <span class="choices" role="group" aria-label="${esc(label ?? 'Scegli la categoria')}: ${esc(it.name.replace(/<[^>]+>/g, ''))}">${options.map(o => `<button type="button" class="choice" data-v="${o.id}">${o.label}</button>`).join('')}</span>
          <p class="why" hidden></p>
        </li>`).join('')}</ol>`;
    update();
  };
  const update = () => {
    root.querySelector('.sorter-score').textContent = done ? `${score} giuste su ${done} · ${order.length - done} da fare` : `${order.length} voci da classificare`;
  };
  root.addEventListener('click', (ev) => {
    if (ev.target.closest('[data-reset]')) { order = shuffle(items); draw(); return; }
    const b = ev.target.closest('.choice');
    if (!b) return;
    const li = b.closest('.sorter-item');
    if (li.dataset.answered) return;
    const it = order[+li.dataset.i];
    const right = b.dataset.v === it.answer;
    li.dataset.answered = right ? 'ok' : 'ko';
    li.querySelectorAll('.choice').forEach(c => {
      c.disabled = true;
      if (c.dataset.v === it.answer) c.classList.add('is-right');
    });
    if (!right) b.classList.add('is-wrong');
    const why = li.querySelector('.why');
    why.innerHTML = `<b>${right ? 'Giusto.' : `No: ${options.find(o => o.id === it.answer).label.toLowerCase()}.`}</b> ${it.why}`;
    why.hidden = false;
    done++; if (right) score++;
    update();
  });
  draw();
}

// Quiz a risposta multipla con spiegazioni; salva il miglior punteggio dell'argomento.
export function quiz(root, questions, { topic, onScore }) {
  const draw = () => {
    const qs = shuffle(questions).map(q => ({ ...q, opts: shuffle(q.options.map((t, i) => ({ t, ok: i === q.correct }))) }));
    const best = topicProgress(topic).quiz;
    root.innerHTML = `
      <form class="quiz" novalidate>
        <ol class="quiz-list">${qs.map((q, i) => `
          <li class="quiz-q" data-i="${i}">
            <fieldset><legend><span class="quiz-n">${i + 1}</span>${q.q}</legend>
              ${q.opts.map((o, j) => `<label class="quiz-opt"><input type="radio" name="q${i}" id="q${i}-${j}" value="${j}"><span>${o.t}</span></label>`).join('')}
            </fieldset>
            <p class="why" hidden></p>
          </li>`).join('')}</ol>
        <div class="quiz-foot">
          <button type="submit" class="btn">Correggi</button>
          <span class="quiz-result" aria-live="polite">${best ? `Miglior risultato finora: ${best.score}/${best.total}` : `${qs.length} domande`}</span>
        </div>
      </form>`;
    const form = root.querySelector('form');
    form.addEventListener('submit', (ev) => {
      ev.preventDefault();
      if (form.dataset.graded) { draw(); return; }
      let score = 0;
      qs.forEach((q, i) => {
        const li = form.querySelector(`.quiz-q[data-i="${i}"]`);
        const sel = form.querySelector(`input[name="q${i}"]:checked`);
        const okIndex = q.opts.findIndex(o => o.ok);
        const right = sel && +sel.value === okIndex;
        if (right) score++;
        li.dataset.answered = right ? 'ok' : 'ko';
        li.querySelectorAll('.quiz-opt').forEach((l, j) => { if (j === okIndex) l.classList.add('is-right'); else if (sel && +sel.value === j) l.classList.add('is-wrong'); });
        li.querySelectorAll('input').forEach(x => { x.disabled = true; });
        const why = li.querySelector('.why');
        why.innerHTML = `<b>${right ? 'Giusto.' : sel ? 'Sbagliato.' : 'Senza risposta.'}</b> ${q.why}`;
        why.hidden = false;
      });
      form.dataset.graded = '1';
      const p = saveProgress(s => {
        const t = s[topic] ?? (s[topic] = {});
        if (!t.quiz || score > t.quiz.score) t.quiz = { score, total: qs.length };
      });
      const verdict = score / qs.length >= 0.9 ? 'Argomento pronto per l\'esame.' : score / qs.length >= 0.7 ? 'Buona base: rileggi le spiegazioni delle risposte sbagliate.' : 'Ripassa le sezioni collegate alle domande sbagliate e riprova.';
      form.querySelector('.quiz-result').innerHTML = `<b>${score}/${qs.length}</b> · ${verdict} Miglior risultato: ${p[topic].quiz.score}/${p[topic].quiz.total}`;
      form.querySelector('button[type=submit]').textContent = 'Rifai il quiz con l\'ordine cambiato';
      onScore?.(score, qs.length);
    });
  };
  draw();
}

// Palestra: esercizi generati a caso, corretti con le regole delle cifre significative.
export function practice(root, { kinds }) {
  let kind = 'misto', current = null, tried = 0, right = 0;
  root.innerHTML = `
    <div class="chips" role="group" aria-label="Tipo di esercizio">
      <button type="button" class="chip" data-k="misto" aria-pressed="true">Misti</button>
      ${kinds.map(k => `<button type="button" class="chip" data-k="${k.id}" aria-pressed="false">${k.label}</button>`).join('')}
    </div>
    <div class="ex-card">
      <p class="ex-kind"></p>
      <div class="ex-prompt"></div>
      <form class="ex-form" novalidate>
        <label class="sr-only" for="ex-answer">Risposta</label>
        <input id="ex-answer" class="ex-input" inputmode="decimal" autocomplete="off" spellcheck="false" placeholder="Risposta">
        <span class="ex-unit"></span>
        <button type="submit" class="btn">Verifica</button>
      </form>
      <p class="ex-feedback" aria-live="polite"></p>
      <div class="ex-actions">
        <button type="button" class="btn-quiet" data-sol>Mostra la soluzione</button>
        <button type="button" class="btn-quiet" data-new>Nuovo esercizio</button>
        <span class="ex-count"></span>
      </div>
      <div class="ex-solution" hidden></div>
    </div>`;
  const $ = (s) => root.querySelector(s);
  const next = () => {
    const k = kind === 'misto' ? pick(kinds) : kinds.find(x => x.id === kind);
    current = { ...k.gen(), kindLabel: k.label, counted: false };
    $('.ex-kind').textContent = k.label;
    $('.ex-prompt').innerHTML = current.prompt;
    $('.ex-unit').innerHTML = current.unit ?? '';
    $('.ex-input').value = '';
    $('.ex-feedback').textContent = '';
    $('.ex-feedback').dataset.status = '';
    $('.ex-solution').hidden = true;
    $('.ex-solution').innerHTML = `<h4>Soluzione</h4>${current.steps}`;
  };
  const count = () => { $('.ex-count').textContent = tried ? `${right} giusti su ${tried}` : ''; };
  root.addEventListener('click', (ev) => {
    const chip = ev.target.closest('.chip');
    if (chip) {
      kind = chip.dataset.k;
      root.querySelectorAll('.chip').forEach(c => c.setAttribute('aria-pressed', String(c === chip)));
      next();
    }
    if (ev.target.closest('[data-new]')) next();
    if (ev.target.closest('[data-sol]')) {
      $('.ex-solution').hidden = false;
      if (!current.counted) { current.counted = true; tried++; count(); }
    }
  });
  $('.ex-form').addEventListener('submit', (ev) => {
    ev.preventDefault();
    const g = grade($('.ex-input').value, current.expected);
    const fb = $('.ex-feedback');
    fb.dataset.status = g.status;
    fb.innerHTML = g.status === 'ok' ? `<b>${g.msg}</b>` : g.msg;
    if (g.status !== 'invalid' && !current.counted) { current.counted = true; tried++; if (g.status === 'ok') right++; count(); }
    if (g.status === 'ok' || g.status === 'almost') $('.ex-solution').hidden = false;
  });
  next();
}

// Esercizi svolti: i passaggi si scoprono uno alla volta, per provare prima da soli.
export function steppers(root) {
  root.querySelectorAll('.svolto').forEach(box => {
    const steps = [...box.querySelectorAll('.passi > li')];
    let shown = 0;
    const bar = document.createElement('div');
    bar.className = 'svolto-bar';
    bar.innerHTML = '<button type="button" class="btn-quiet" data-next>Mostra il primo passaggio</button><button type="button" class="btn-quiet" data-all>Mostra tutto</button>';
    box.append(bar);
    const sync = () => {
      steps.forEach((s, i) => { s.hidden = i >= shown; });
      const nb = bar.querySelector('[data-next]');
      nb.hidden = shown >= steps.length;
      bar.querySelector('[data-all]').hidden = shown >= steps.length;
      nb.textContent = shown ? `Passaggio successivo (${shown}/${steps.length})` : 'Mostra il primo passaggio';
      box.querySelector('.passi').hidden = shown === 0;
    };
    bar.addEventListener('click', (ev) => {
      if (ev.target.closest('[data-next]')) shown++;
      if (ev.target.closest('[data-all]')) shown = steps.length;
      sync();
    });
    sync();
  });
}

// Obiettivi d'esame spuntabili, salvati per argomento.
export function checklist(root, topic) {
  const st = topicProgress(topic).obiettivi ?? {};
  root.querySelectorAll('input[type=checkbox]').forEach(cb => {
    cb.checked = !!st[cb.id];
    cb.addEventListener('change', () => saveProgress(s => { const t = s[topic] ?? (s[topic] = {}); (t.obiettivi ??= {})[cb.id] = cb.checked; }));
  });
}
