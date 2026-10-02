// Navigazione: #programma (o nessun indirizzo) mostra il programma, #<id-argomento> apre l'argomento.

import { PROGRAMMA, PARTI, argomento } from './programma.js';
import { loadProgress } from './ui.js';

const main = document.getElementById('main');
const MODULI = { materia: () => import('./argomenti/materia.js') };

function tile(a, { link = true } = {}) {
  const pr = loadProgress()[a.id]?.quiz;
  const pct = pr ? Math.round(100 * pr.score / pr.total) : 0;
  const inner = `
    <span class="tile-top"><span>${a.n}</span><span>${a.ore} h</span></span>
    <span class="tile-s">${a.simbolo}</span>
    <span class="tile-t">${a.titolo}</span>
    ${a.pronto ? `<span class="tile-p" title="${pr ? `Quiz: ${pr.score}/${pr.total}` : 'Quiz non ancora fatto'}"><i style="width:${pct}%"></i></span>` : '<span class="tile-badge">in arrivo</span>'}`;
  if (a.pronto && link) return `<a class="tile is-ready" href="#${a.id}" aria-label="${a.n}. ${a.titolo}, ${a.ore} ore${pr ? `, quiz ${pr.score} su ${pr.total}` : ''}">${inner}</a>`;
  return `<div class="tile ${a.pronto ? 'is-ready' : 'is-locked'}" ${a.pronto ? '' : 'aria-disabled="true"'}>${inner}</div>`;
}

function home() {
  document.title = 'Quaderno di Chimica Generale';
  const prog = loadProgress();
  const ready = PROGRAMMA.filter(a => a.pronto);
  const passed = ready.filter(a => prog[a.id]?.quiz && prog[a.id].quiz.score / prog[a.id].quiz.total >= 0.9).length;
  const ore = PROGRAMMA.reduce((s, a) => s + a.ore, 0);
  main.innerHTML = `
    <section class="intro">
      <p class="eyebrow">Chimica generale e inorganica · ${PROGRAMMA.length} argomenti · ${ore} ore di lezione</p>
      <h1>Quaderno di Chimica Generale</h1>
      <p class="lead">Il programma del corso, un argomento alla volta. Ogni argomento ha la teoria, gli strumenti per vedere i concetti, esercizi svolti passo passo, una palestra di esercizi sempre nuovi, le domande d'orale e un quiz finale.</p>
      <div class="progress-line">
        <span><b>${ready.length}</b> di ${PROGRAMMA.length} argomenti pronti</span>
        <span><b>${passed}</b> quiz superati con almeno il 90%</span>
      </div>
    </section>
    ${PARTI.map(p => {
      const items = PROGRAMMA.filter(a => a.parte === p.id);
      return `<section class="part" aria-labelledby="parte-${p.id}">
        <div class="part-head"><h2 id="parte-${p.id}">${p.titolo}</h2><span>${p.sottotitolo} · ${items.reduce((s, a) => s + a.ore, 0)} ore</span></div>
        <ul class="tiles">${items.map(a => `<li>${tile(a)}</li>`).join('')}</ul>
      </section>`;
    }).join('')}
    <section class="how">
      <h2>Come studiare un argomento</h2>
      <ol>
        <li>Leggi la teoria e usa gli strumenti interattivi finché il concetto è chiaro.</li>
        <li>Rifai gli esercizi svolti coprendo la soluzione: scopri un passaggio solo quando sei bloccato.</li>
        <li>Allenati in palestra finché sbagli raramente, anche sulle cifre significative.</li>
        <li>Rispondi a voce alle domande d'orale, poi confronta con la risposta modello.</li>
        <li>Fai il quiz: punta ad almeno il 90% prima di passare all'argomento successivo.</li>
      </ol>
      <p class="note">I progressi (obiettivi spuntati e punteggi dei quiz) restano salvati solo in questo browser.</p>
    </section>
    <footer class="foot">Fonti: R. H. Petrucci, F. G. Herring, J. D. Madura, C. Bissonnette, <i>Chimica generale. Principi ed applicazioni moderne</i>, Piccin; BIPM, <i>Il Sistema Internazionale di Unità</i>, 9ª edizione (2019); CRC Handbook of Chemistry and Physics per i dati numerici.</footer>`;
}

async function topic(a) {
  document.title = `${a.titolo} · Quaderno di Chimica Generale`;
  main.innerHTML = '<p class="note" style="padding-top:40px">Caricamento…</p>';
  const mod = await MODULI[a.id]();
  main.innerHTML = `
    <header class="topic-head">
      ${tile(a, { link: false })}
      <div>
        <p class="eyebrow">Argomento ${a.n} di ${PROGRAMMA.length} · Parte ${a.parte} · ${a.ore} ore di lezione</p>
        <h1>${a.titolo}</h1>
        ${a.fonte ? `<p class="topic-meta">Riferimento: ${a.fonte}.</p>` : ''}
      </div>
    </header>
    <nav class="toc" aria-label="Sezioni dell'argomento">${mod.SEZIONI.map(([id, t]) => `<button type="button" data-s="${id}">${t}</button>`).join('')}</nav>
    <div id="topic-body"></div>
    ${nextPrev(a)}`;
  mod.render(document.getElementById('topic-body'));
  const toc = main.querySelector('.toc');
  toc.addEventListener('click', (ev) => {
    const b = ev.target.closest('[data-s]');
    if (b) document.getElementById(b.dataset.s)?.scrollIntoView({ behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth' });
  });
  const obs = new IntersectionObserver((entries) => {
    entries.filter(e => e.isIntersecting).forEach(e => {
      toc.querySelectorAll('button').forEach(b => {
        const on = b.dataset.s === e.target.id;
        b.setAttribute('aria-current', String(on));
        if (on) b.scrollIntoView({ block: 'nearest', inline: 'nearest' });
      });
    });
  }, { rootMargin: '-130px 0px -65% 0px' });
  main.querySelectorAll('.sec').forEach(s => obs.observe(s));
  window.scrollTo(0, 0);
}

function nextPrev(a) {
  const next = PROGRAMMA.find(x => x.n === a.n + 1);
  return `<footer class="foot">
    <p><a href="#programma">← Torna al programma</a>${next ? ` · Prossimo argomento: <b>${next.titolo}</b>${next.pronto ? ` (<a href="#${next.id}">apri</a>)` : ' (in preparazione)'}` : ''}</p>
    <p>Fonti: Petrucci et al., <i>Chimica generale</i>, Piccin; BIPM, SI 9ª ed. (2019); CRC Handbook of Chemistry and Physics.</p>
  </footer>`;
}

function route() {
  const id = location.hash.slice(1);
  const a = argomento(id);
  if (a?.pronto && MODULI[a.id]) topic(a);
  else { home(); window.scrollTo(0, 0); }
}

window.addEventListener('hashchange', route);
route();
