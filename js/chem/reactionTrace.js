// Observables of an accepted MD trajectory. No reaction templates, inferred electron paths,
// oxidation states or transition-state assignments are used here.

export const TRACE_LIMITS = { frames: 600, bytes: 32 * 1024 * 1024 };
const pairKey = (i, j) => `${Math.min(i, j)}:${Math.max(i, j)}`;

export function bondOrder(bonds, i, j) {
  for (let k = 0; k < bonds.length; k += 4)
    if ((bonds[k] === i && bonds[k + 1] === j) || (bonds[k] === j && bonds[k + 1] === i)) return bonds[k + 3];
  return 0;
}

/** Distance Å, angle degrees, signed IUPAC-style torsion degrees. Degenerate geometry → null. */
export function geometry(pos, ids) {
  const valid = ids.filter(i => Number.isInteger(i) && i >= 0 && 3 * i + 2 < pos.length);
  if (valid.length !== ids.length || new Set(ids).size !== ids.length) return { distance: null, angle: null, dihedral: null };
  const vec = (i, j) => [0, 1, 2].map(c => pos[3 * j + c] - pos[3 * i + c]);
  const dot = (a, b) => a.reduce((s, v, c) => s + v * b[c], 0);
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = a => Math.hypot(...a);
  const out = { distance: null, angle: null, dihedral: null };
  if (ids.length < 2) return out;
  out.distance = norm(vec(ids[0], ids[1]));
  if (ids.length < 3) return out;
  const a = vec(ids[1], ids[0]), b = vec(ids[1], ids[2]);
  if (norm(a) * norm(b) > 1e-12) out.angle = Math.acos(Math.max(-1, Math.min(1, dot(a, b) / (norm(a) * norm(b))))) * 180 / Math.PI;
  if (ids.length < 4 || norm(b) < 1e-12) return out;
  const c = vec(ids[2], ids[3]), u = b.map(v => v / norm(b));
  const p = a.map((v, k) => v - dot(a, u) * u[k]), q = c.map((v, k) => v - dot(c, u) * u[k]);
  if (norm(p) * norm(q) > 1e-12) out.dihedral = Math.atan2(dot(cross(u, p), q), dot(p, q)) * 180 / Math.PI;
  return out;
}

/** Population changes at fixed atomic identities: Δn_A = −Δq_A. Not pairwise current. */
export function populationChange(q, reference) {
  if (!reference || q.length !== reference.length) return null;
  const delta = Array.from(q, (v, i) => reference[i] - v);
  return { delta, total: delta.reduce((s, v) => s + v, 0), gained: delta.reduce((s, v) => s + Math.max(v, 0), 0) };
}

export class BondMonitor {
  constructor({ form = 0.55, breakAt = 0.35 } = {}) {
    if (!(form > breakAt && breakAt >= 0)) throw new Error('Invalid bond hysteresis.');
    this.form = form; this.breakAt = breakAt; this.reset();
  }
  reset() { this.previous = null; this.connected = new Set(); this.serial = 0; }
  sample({ t, pos, q, bonds }) {
    const current = new Map();
    for (let k = 0; k < bonds.length; k += 4) current.set(pairKey(bonds[k], bonds[k + 1]), bonds[k + 3]);
    const events = [];
    const keys = new Set([...current.keys(), ...this.connected]);
    for (const key of keys) {
      const old = this.connected.has(key), order = current.get(key) ?? 0;
      const connected = order > this.form || (old && order > this.breakAt);
      if (connected) this.connected.add(key); else this.connected.delete(key);
      if (this.previous && old !== connected) {
        const [i, j] = key.split(':').map(Number);
        events.push({ serial: ++this.serial, t, tBefore: this.previous.t, i, j, kind: connected ? 'formed' : 'broken',
          before: this.previous.orders.get(key) ?? 0, after: order,
          distance: geometry(pos, [i, j]).distance,
          deltaPopulation: [this.previous.q[i] - q[i], this.previous.q[j] - q[j]],
          at: [0, 1, 2].map(c => (pos[3 * i + c] + pos[3 * j + c]) / 2) });
      }
    }
    this.previous = { t, q: Array.from(q), orders: current };
    return events;
  }
}

const bytesOf = x => {
  if (ArrayBuffer.isView(x)) return x.byteLength;
  if (Array.isArray(x)) return x.reduce((s, v) => s + bytesOf(v), 0);
  if (x && typeof x === 'object') return Object.values(x).reduce((s, v) => s + bytesOf(v), 0);
  return typeof x === 'number' ? 8 : typeof x === 'string' ? 2 * x.length : 4;
};

/** Ring buffer with real samples (no interpolation). Stored wavefunctions reproduce the SCF branch. */
export class TrajectoryRecorder {
  constructor({ maxFrames = TRACE_LIMITS.frames, maxBytes = TRACE_LIMITS.bytes } = {}) {
    this.maxFrames = maxFrames; this.maxBytes = maxBytes; this.epoch = 0; this.reset();
  }
  reset(reason = 'Nuovo esperimento') {
    this.epoch++; this.reason = reason; this.frames = []; this.bytes = 0; this.nextId = 0; this.evicted = 0;
  }
  add(frame, wave = null) {
    const sample = structuredClone({ ...frame, id: this.nextId++, epoch: this.epoch, wave });
    const bytes = bytesOf(sample);
    if (bytes > this.maxBytes) { sample.wave = null; }
    sample.bytes = bytesOf(sample);
    this.frames.push(sample); this.bytes += sample.bytes;
    while (this.frames.length > 1 && (this.frames.length > this.maxFrames || this.bytes > this.maxBytes)) {
      this.bytes -= this.frames.shift().bytes; this.evicted++;
    }
    return sample;
  }
  get(id) { return this.frames.find(f => f.id === id) ?? null; }
  nearest(t) {
    return this.frames.reduce((a, b) => !a || Math.abs(b.stats.t - t) < Math.abs(a.stats.t - t) ? b : a, null);
  }
  meta() {
    return { epoch: this.epoch, reason: this.reason, count: this.frames.length, bytes: this.bytes, evicted: this.evicted,
      first: this.frames[0]?.id ?? null, last: this.frames.at(-1)?.id ?? null };
  }
}
