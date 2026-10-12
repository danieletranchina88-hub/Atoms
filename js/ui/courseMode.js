// The course is mounted and disposed independently of the existing scientific views.
import { UNITS, TOOLS, SOURCES } from '../course/content.js';
import { mountTool, esc } from '../course/widgets.js';

const STORAGE = 'atoms-general-chemistry-course-v1';
function readState() {
  try {
    const raw = JSON.parse(localStorage.getItem(STORAGE) || '{}');
    const object = (v) =>
      v && typeof v === 'object' && !Array.isArray(v) ? v : {};
    return {
      unit: UNITS.some((u) => u.id === raw.unit) ? raw.unit : 'introduzione',
      tool: TOOLS.some((t) => t[0] === raw.tool) ? raw.tool : null,
      attempts:
        Number.isSafeInteger(raw.attempts) && raw.attempts >= 0
          ? raw.attempts
          : 0,
      solved: object(raw.solved),
      studied: object(raw.studied),
      notes: Array.isArray(raw.notes)
        ? raw.notes
            .filter(
              (n) => n && typeof n.text === 'string' && n.text.length <= 4000,
            )
            .slice(0, 40)
        : [],
    };
  } catch {
    return {
      unit: 'introduzione',
      tool: null,
      attempts: 0,
      solved: {},
      studied: {},
      notes: [],
    };
  }
}
function route(hash) {
  const parts = String(hash).replace(/^#/, '').split('/');
  return parts[0] === 'corso'
    ? {
        unit: UNITS.some((u) => u.id === parts[1]) ? parts[1] : null,
        tool: TOOLS.some((t) => t[0] === parts[2]) ? parts[2] : null,
      }
    : {};
}
export function mountCourse(
  root,
  {
    openAtom = () => {},
    openSandbox = () => {},
    initialRoute = location.hash,
  } = {},
) {
  const progress = readState(),
    start = route(initialRoute),
    shared = { rx: 'argento' };
  let current =
      UNITS.find((u) => u.id === (start.unit || progress.unit)) || UNITS[0],
    tool = start.tool || progress.tool || current.tool,
    cleanup = () => {},
    disposed = false,
    storageOK = true;
  root.innerHTML = `<div class="course-app">
    <header class="cl-course-hero"><div><p class="cl-eyebrow">A T O M S &nbsp; / &nbsp; C H I M I C A &nbsp; G E N E R A L E</p><h2>Capire. Provare.<br><em>Fare chimica.</em></h2><p>Un percorso dalle tue lezioni al banco di lavoro.<br>Le formule diventano esperimenti, gli errori diventano spiegazioni.</p><div class="cl-hero-links"><button type="button" id="cl-pause-animations" class="cl-btn" aria-pressed="false">Ⅱ Pausa animazioni</button><button type="button" class="cl-btn" data-course-sources>Fonti del corso ↗</button><button type="button" class="cl-btn" data-course-sandbox>Sandbox quantistico ↗</button></div></div><div class="cl-hero-symbol" aria-hidden="true"><span>n</span><div class="cl-orbit cl-orbit-a"></div><div class="cl-orbit cl-orbit-b"></div><div class="cl-orbit cl-orbit-c"></div><small>osserva · calcola · verifica</small></div></header>
    <div class="cl-layout"><aside class="cl-syllabus"><div class="cl-syllabus-title"><span class="cl-eyebrow">I L &nbsp; P E R C O R S O</span><span id="cl-progress-count"></span></div><nav aria-label="Argomenti del corso">${UNITS.map((u, i) => `<button type="button" data-course-unit="${u.id}" class="cl-unit"><span class="cl-unit-number">${String(i + 1).padStart(2, '0')}</span><span><strong>${esc(u.title)}</strong><small>${esc(u.subtitle)}</small></span><span class="cl-unit-done" aria-label="Studiato">✓</span></button>`).join('')}</nav><div class="cl-syllabus-foot"><span class="cl-progress-line"><i id="cl-progress-fill"></i></span><p>Progressi e quaderno salvati su questo dispositivo.</p><p id="cl-storage-warning" hidden>Il browser non consente il salvataggio: i progressi restano disponibili solo durante questa visita.</p><button type="button" id="cl-open-notebook" class="cl-btn">Il mio quaderno <span id="cl-note-count">0</span></button></div></aside>
    <main class="cl-main"><section class="cl-lesson" id="cl-lesson"></section><section class="cl-workbench"><div class="cl-tool-navigation"><span class="cl-eyebrow">S C E G L I &nbsp; U N ’ A T T I V I T À</span><nav aria-label="Strumenti del corso">${TOOLS.map(([id, label, n]) => `<button type="button" data-course-tool="${id}"><small>${n}</small>${esc(label)}</button>`).join('')}</nav></div><div class="cl-tool" id="cl-tool"></div><footer class="cl-notebook-bar"><p>Annota cosa hai provato e conserva i dati del calcolo.</p><button type="button" id="cl-note-current" class="cl-btn">+ Nel quaderno</button><span id="cl-note-status" role="status"></span></footer></section></main></div>
    <dialog class="cl-dialog" id="cl-source-dialog" aria-labelledby="cl-source-title"><div class="cl-dialog-top"><h3 id="cl-source-title">Materiale e modelli del corso</h3><button type="button" data-dialog-close class="cl-btn" aria-label="Chiudi fonti">×</button></div><p>Spiegazioni originali basate sulle lezioni di <b>Dario Duca</b>, Dipartimento di Fisica e Chimica, Università di Palermo, Chimica Generale ed Inorganica I, A.A. 2026–2027. I riferimenti usano il numero di pagina del PDF.</p>${SOURCES.map((s) => `<article class="cl-source-card"><span>${s.id}</span><div><h4>${esc(s.title)}</h4><p>${esc(s.file)} · ${s.pages} pagine · ${esc(s.date)}</p></div></article>`).join('')}<p>Le slide originali non sono ripubblicate qui. Ogni lezione indica dove ritrovare l’argomento nel materiale fornito.</p><h4>Aggiornamenti e approfondimenti dichiarati</h4><ul><li>Definizione attuale della mole: <a href="https://www.bipm.org/en/si-base-units/mole" target="_blank" rel="noopener">BIPM, SI base unit: mole</a>. Il riferimento ai 12 g di ¹²C nelle slide è storico.</li><li>Semireazioni e reazioni aggiuntive applicano le regole delle lezioni e sono segnate come approfondimenti.</li><li>Lo scattering α usa la soluzione classica del potenziale coulombiano repulsivo: nucleo puntiforme fisso, nessuna interazione nucleare o schermaggio.</li><li>Gli esperimenti del banco usano bilanci stechiometrici con quantità calcolate; l’avanzamento animato non è una cinetica simulata. I colori e i marcatori sono rappresentazioni didattiche.</li><li>La simulazione DFT-LDA della tavola periodica e i modelli quantistici del sandbox restano indipendenti dal corso.</li></ul></dialog>
    <dialog class="cl-dialog" id="cl-notebook-dialog" aria-labelledby="cl-notebook-title"><div class="cl-dialog-top"><h3 id="cl-notebook-title">Il mio quaderno di chimica</h3><button type="button" data-dialog-close class="cl-btn" aria-label="Chiudi quaderno">×</button></div><label class="cl-field" for="cl-note-text"><span>Che cosa hai osservato o vuoi ricordare?</span><textarea id="cl-note-text" rows="4" maxlength="4000" placeholder="Es. Ho dimezzato le moli di NaCl: cambia il limitante e resta AgNO₃ in eccesso."></textarea></label><div class="cl-actions"><button type="button" id="cl-save-note" class="cl-btn primary">Salva con i dati attuali</button><button type="button" id="cl-export-notes" class="cl-btn">Esporta JSON</button></div><p class="cl-hint">Al massimo 40 osservazioni su questo browser. I dati personali non vengono inviati a un server.</p><div id="cl-notebook-list"></div><p id="cl-notebook-status" role="status"></p></dialog>
  </div>`;
  const $ = (id) => root.querySelector('#' + id);
  const save = () => {
    progress.unit = current.id;
    progress.tool = tool;
    try {
      localStorage.setItem(STORAGE, JSON.stringify(progress));
    } catch {
      storageOK = false;
      $('cl-storage-warning').hidden = false;
    }
    updateProgress();
  };
  const updateProgress = () => {
    const studied = UNITS.filter((u) => progress.studied[u.id]).length;
    $('cl-progress-count').textContent = `${studied}/${UNITS.length}`;
    $('cl-progress-fill').style.width = (100 * studied) / UNITS.length + '%';
    $('cl-note-count').textContent = progress.notes.length;
    root.querySelectorAll('[data-course-unit]').forEach((b) => {
      b.classList.toggle('selected', b.dataset.courseUnit === current.id);
      b.classList.toggle('studied', !!progress.studied[b.dataset.courseUnit]);
      b.setAttribute(
        'aria-current',
        b.dataset.courseUnit === current.id ? 'step' : 'false',
      );
    });
    $('cl-storage-warning').hidden = storageOK;
  };
  const recordAttempt = (key, correct, unit) => {
    progress.attempts++;
    if (correct) {
      if (Object.keys(progress.solved).length >= 500 && !progress.solved[key])
        delete progress.solved[Object.keys(progress.solved)[0]];
      progress.solved[key] = { unit, date: new Date().toISOString() };
    }
    save();
  };
  const ctx = {
    shared,
    progress,
    units: UNITS,
    unit: current.id,
    snapshot: null,
    paused: false,
    recordAttempt,
    openAtom,
    openSandbox,
    showTool: (id) => showTool(id),
  };
  const writeHash = () => {
    try {
      history.replaceState(null, '', `#corso/${current.id}/${tool}`);
    } catch {}
  };
  const showTool = (id, persist = true) => {
    if (disposed) return;
    cleanup();
    ctx.snapshot = null;
    tool = TOOLS.some((t) => t[0] === id) ? id : current.tool;
    root.querySelectorAll('[data-course-tool]').forEach((b) => {
      const on = b.dataset.courseTool === tool;
      b.classList.toggle('selected', on);
      b.setAttribute('aria-pressed', String(on));
    });
    ctx.unit = current.id;
    cleanup = mountTool($('cl-tool'), tool, ctx) || (() => {});
    $('cl-note-status').textContent = '';
    if (persist) {
      save();
      writeHash();
    }
  };
  const showUnit = (id, nextTool) => {
    current = UNITS.find((u) => u.id === id) || UNITS[0];
    ctx.unit = current.id;
    const i = UNITS.indexOf(current);
    $('cl-lesson').innerHTML =
      `<div class="cl-lesson-header"><div><span class="cl-eyebrow">L E Z I O N E &nbsp; ${String(i + 1).padStart(2, '0')} / ${UNITS.length}</span><h3>${esc(current.title)}</h3><p>${esc(current.subtitle)}</p></div><button type="button" id="cl-mark-studied" class="cl-btn${progress.studied[current.id] ? ' primary' : ''}">${progress.studied[current.id] ? '✓ Studiata' : 'Segna come studiata'}</button></div><div class="cl-objectives">${current.goals.map((g) => `<span>↗ ${esc(g)}</span>`).join('')}</div><details class="cl-reading" open><summary>Leggi la spiegazione <span>${esc(current.source)} · pp. ${esc(current.pages)}</span></summary><div class="cl-reading-grid">${current.sections.map(([title, text]) => `<article><h4>${esc(title)}</h4><p>${esc(text)}</p></article>`).join('')}</div><button type="button" class="cl-source-link" data-course-sources>Riferimenti al materiale ↗</button></details><div class="cl-lesson-actions"><button type="button" id="cl-unit-practice" class="cl-btn">Esercitati su questo argomento</button>${i < UNITS.length - 1 ? `<button type="button" id="cl-unit-next" class="cl-btn">Lezione successiva →</button>` : ''}</div>`;
    $('cl-mark-studied').onclick = () => {
      if (progress.studied[current.id]) delete progress.studied[current.id];
      else progress.studied[current.id] = new Date().toISOString();
      save();
      $('cl-mark-studied').textContent = progress.studied[current.id]
        ? '✓ Studiata'
        : 'Segna come studiata';
      $('cl-mark-studied').classList.toggle(
        'primary',
        !!progress.studied[current.id],
      );
    };
    $('cl-unit-practice').onclick = () => showTool('esercizi');
    if ($('cl-unit-next'))
      $('cl-unit-next').onclick = () => showUnit(UNITS[i + 1].id);
    $('cl-lesson').querySelector('[data-course-sources]').onclick = () =>
      openDialog($('cl-source-dialog'));
    showTool(nextTool || current.tool);
    updateProgress();
  };
  const openDialog = (dialog) => {
    if (typeof dialog.showModal === 'function') dialog.showModal();
    else dialog.setAttribute('open', '');
  };
  const closeDialog = (dialog) => {
    if (typeof dialog.close === 'function') dialog.close();
    else dialog.removeAttribute('open');
  };
  const renderNotes = () => {
    $('cl-notebook-list').innerHTML = progress.notes.length
      ? progress.notes
          .map(
            (n, i) =>
              `<article class="cl-note-entry"><div><span>${esc(n.unit || 'Chimica')} · ${esc(n.date || '')}</span><button type="button" class="cl-btn" data-note-delete="${i}" aria-label="Elimina osservazione ${i + 1}">×</button></div><p>${esc(n.text)}</p>${n.snapshot ? `<details><summary>Dati del calcolo</summary><pre>${esc(JSON.stringify(n.snapshot, null, 2))}</pre></details>` : ''}</article>`,
          )
          .join('')
      : '<p class="cl-hint">Il quaderno è vuoto. Svolgi un calcolo, cambia un parametro e annota cosa succede.</p>';
    $('cl-notebook-list')
      .querySelectorAll('[data-note-delete]')
      .forEach(
        (b) =>
          (b.onclick = () => {
            progress.notes.splice(Number(b.dataset.noteDelete), 1);
            save();
            renderNotes();
          }),
      );
  };
  const openNotes = () => {
    renderNotes();
    $('cl-notebook-status').textContent = '';
    openDialog($('cl-notebook-dialog'));
  };
  $('cl-pause-animations').onclick = () => {
    ctx.paused = !ctx.paused;
    $('cl-pause-animations').textContent = ctx.paused
      ? '▶ Riprendi animazioni'
      : 'Ⅱ Pausa animazioni';
    $('cl-pause-animations').setAttribute('aria-pressed', String(ctx.paused));
  };
  root
    .querySelectorAll('[data-course-unit]')
    .forEach((b) => (b.onclick = () => showUnit(b.dataset.courseUnit)));
  root
    .querySelectorAll('[data-course-tool]')
    .forEach((b) => (b.onclick = () => showTool(b.dataset.courseTool)));
  root
    .querySelectorAll('[data-course-sources]')
    .forEach((b) => (b.onclick = () => openDialog($('cl-source-dialog'))));
  root.querySelector('[data-course-sandbox]').onclick = openSandbox;
  root
    .querySelectorAll('[data-dialog-close]')
    .forEach((b) => (b.onclick = () => closeDialog(b.closest('dialog'))));
  root.querySelectorAll('dialog').forEach((d) =>
    d.addEventListener('click', (e) => {
      if (e.target === d) {
        const r = d.getBoundingClientRect();
        if (
          e.clientX < r.left ||
          e.clientX > r.right ||
          e.clientY < r.top ||
          e.clientY > r.bottom
        )
          closeDialog(d);
      }
    }),
  );
  $('cl-open-notebook').onclick = openNotes;
  $('cl-note-current').onclick = openNotes;
  $('cl-save-note').onclick = () => {
    const text = $('cl-note-text').value.trim();
    if (!text) {
      $('cl-notebook-status').textContent = 'Scrivi almeno un’osservazione.';
      return;
    }
    if (progress.notes.length >= 40) {
      $('cl-notebook-status').textContent =
        'Quaderno pieno: esporta o elimina un’osservazione prima di aggiungerne un’altra.';
      return;
    }
    progress.notes.unshift({
      text,
      unit: current.title,
      date: new Date().toLocaleString('it-IT'),
      snapshot: ctx.snapshot ? JSON.parse(JSON.stringify(ctx.snapshot)) : null,
    });
    save();
    $('cl-note-text').value = '';
    $('cl-notebook-status').textContent = storageOK
      ? 'Osservazione salvata.'
      : 'Osservazione disponibile per questa visita; salvataggio locale non consentito.';
    renderNotes();
    $('cl-note-status').textContent = 'Osservazione aggiunta.';
  };
  $('cl-export-notes').onclick = () => {
    const blob = new Blob(
      [
        JSON.stringify(
          {
            corso: 'Atoms · Chimica generale',
            versione: 1,
            esportato: new Date().toISOString(),
            osservazioni: progress.notes,
            progressi: { studiati: progress.studied, risolti: progress.solved },
          },
          null,
          2,
        ),
      ],
      { type: 'application/json' },
    );
    const url = URL.createObjectURL(blob),
      a = document.createElement('a');
    a.href = url;
    a.download = 'atoms-quaderno-chimica.json';
    a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
  showUnit(current.id, tool);
  return {
    navigate(hash) {
      const next = route(hash);
      if (next.unit) showUnit(next.unit, next.tool);
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      cleanup();
      root.querySelectorAll('dialog[open]').forEach(closeDialog);
      root.replaceChildren();
    },
  };
}
