// Interfaccia verso il worker molecolare: richieste con promesse, avanzamento e annullamento.

let worker = null;
let nextId = 1;
const pending = new Map();

function start() {
  worker = new Worker(new URL('./moleculeWorker.js', import.meta.url), { type: 'module' });
  worker.onmessage = (ev) => {
    const { id, result, error, progress, ...extra } = ev.data;
    const p = pending.get(id);
    if (!p) return;
    if (progress !== undefined && result === undefined && error === undefined) {
      p.onProgress?.(progress, extra);
      return;
    }
    pending.delete(id);
    if (error) p.reject(new Error(error)); else p.resolve(result);
  };
  worker.onerror = (e) => {
    for (const p of pending.values()) p.reject(new Error(e.message || 'Errore nel worker'));
    pending.clear();
  };
}

export function request(type, payload = {}, onProgress = null) {
  if (!worker) start();
  const id = nextId++;
  return new Promise((resolve, reject) => {
    pending.set(id, { resolve, reject, onProgress });
    worker.postMessage({ id, type, ...payload });
  });
}

/** Interrompe tutti i calcoli in corso (il worker viene ricreato alla richiesta successiva). */
export function cancelAll() {
  if (!worker) return;
  worker.terminate();
  worker = null;
  for (const p of pending.values()) p.reject(new Error('annullato'));
  pending.clear();
}

export function busy() {
  return pending.size > 0;
}
