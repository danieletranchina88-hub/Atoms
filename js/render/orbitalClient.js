import { prepareOrbitalScene } from './orbitalScene.js';

let worker;
let nextId = 0;
const pending = new Map();

export function requestOrbitalScene(payload) {
  if (worker === undefined) {
    try {
      worker = new Worker(new URL('./orbitalWorker.js', import.meta.url), { type: 'module' });
      worker.onmessage = ({ data: { id, scene, error } }) => {
        const p = pending.get(id);
        if (!p) return;
        pending.delete(id);
        if (error) p.reject(new Error(error)); else p.resolve(scene);
      };
      worker.onerror = () => {
        worker.terminate();
        worker = null;
        for (const p of pending.values()) prepareOrbitalScene(p.payload).then(p.resolve, p.reject);
        pending.clear();
      };
    } catch { worker = null; }
  }
  if (!worker) return prepareOrbitalScene(payload);
  return new Promise((resolve, reject) => {
    const id = ++nextId;
    pending.set(id, { resolve, reject, payload });
    worker.postMessage({ id, payload });
  });
}
