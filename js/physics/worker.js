// Web Worker: esegue il calcolo autoconsistente senza bloccare l'interfaccia.
import { computeAtom } from './atom.js';

self.onmessage = (ev) => {
  const { id, Z } = ev.data;
  try {
    const atom = computeAtom(Z);
    self.postMessage({ id, atom });
  } catch (err) {
    self.postMessage({ id, error: String(err && err.message || err) });
  }
};
