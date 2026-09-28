// Dinamica molecolare della sandbox: equazioni di Newton integrate con l'algoritmo velocity Verlet,
// termostato stocastico di Bussi–Donadio–Parrinello (distribuzione canonica esatta), pareti morbide
// che misurano la pressione come flusso di quantità di moto, riconoscimento automatico delle specie
// chimiche e delle reazioni.
//
// Unità: Å, fs, amu, eV, K.

import { ReactiveFF } from './reactive.js';
import { ELEMENTS } from '../physics/elements.js';
import { hillFormula } from './smiles.js';
import { KB_EV } from './reactiveData.js';

export const ACC = 0.00964853321;        // (eV/Å)/amu → Å/fs²
export const MV2 = 1 / ACC;              // amu·Å²/fs² → eV
export const EV_A3_TO_BAR = 1.602176634e6;
const K_WALL = 10;                       // costante elastica delle pareti (eV/Å²)
const K_GRAB = 3;                        // pinzetta ottica (eV/Å²)

// numeri casuali: generatore di Marsaglia (xorshift) con distribuzioni normale e gamma
export function makeRng(seed = 12345) {
  let s = seed >>> 0 || 1;
  const uni = () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return (s + 0.5) / 4294967296; };
  let spare = null;
  const gauss = () => {
    if (spare !== null) { const v = spare; spare = null; return v; }
    let u, v, q;
    do { u = 2 * uni() - 1; v = 2 * uni() - 1; q = u * u + v * v; } while (q >= 1 || q === 0);
    const f = Math.sqrt(-2 * Math.log(q) / q);
    spare = v * f;
    return u * f;
  };
  // gamma(a, 1) di Marsaglia–Tsang
  const gamma = (a) => {
    if (a < 1) return gamma(a + 1) * Math.pow(uni(), 1 / a);
    const d = a - 1 / 3, c = 1 / Math.sqrt(9 * d);
    for (;;) {
      let x, v;
      do { x = gauss(); v = 1 + c * x; } while (v <= 0);
      v = v * v * v;
      const u = uni();
      if (u < 1 - 0.0331 * x * x * x * x) return d * v;
      if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
    }
  };
  return { uni, gauss, gamma };
}

export class Simulation {
  constructor({ box = 24, T = 300, dt = 0.4, seed = 2024 } = {}) {
    this.ff = new ReactiveFF();
    this.rng = makeRng(seed);
    this.Z = [];
    this.pos = new Float64Array(0);
    this.vel = new Float64Array(0);
    this.F = new Float64Array(0);
    this.mass = new Float64Array(0);
    this.box = box;
    this.T = T;
    this.thermostat = true;
    this.tau = 100;          // fs
    this.dt = dt;
    this.time = 0;
    this.res = null;
    this.Ewall = 0;
    this.wallForce = 0;
    this.pAccum = 0; this.pSamples = 0; this.pressure = 0;
    this.heatBath = 0;       // energia ceduta dal termostato (eV)
    this.work = 0;           // energia fornita da scintille, luce, pinzetta (eV)
    this.grab = null;        // { i, target: [x,y,z] }
    this.bondedPrev = new Set();
    this.fragOf = null;
    this.events = [];
    this.history = [];
    this.censusEvery = 50;
    this.stepCount = 0;
    this.species = new Map();
    this.clamped = 0;
    this.provider = null;    // forze alternative (per esempio Hartree–Fock ab initio)
  }

  get N() { return this.Z.length; }

  resize(N) {
    const grow = (a) => { const b = new Float64Array(3 * N); b.set(a.subarray(0, Math.min(a.length, 3 * N))); return b; };
    this.pos = grow(this.pos); this.vel = grow(this.vel); this.F = new Float64Array(3 * N);
    this.mass = Float64Array.from(this.Z, z => ELEMENTS[z - 1].mass);
    this.ff.splitPrev.clear();
    this.res = null;
    this.fragOf = null;
    this.bondedPrev = new Set();
  }

  clear() {
    this.Z = [];
    this.resize(0);
    this.time = 0; this.heatBath = 0; this.work = 0; this.events = []; this.history = []; this.species = new Map();
    this.grab = null;
  }

  /** Velocità di Maxwell–Boltzmann alla temperatura T per l'atomo i. */
  thermalize(i, T) {
    const s = Math.sqrt(KB_EV * T / (this.mass[i] * MV2));
    for (let c = 0; c < 3; c++) this.vel[3 * i + c] = s * this.rng.gauss();
  }

  randomRotation() {
    // quaternione uniforme (Shoemake)
    const { uni } = this.rng;
    const u1 = uni(), u2 = uni(), u3 = uni();
    const a = Math.sqrt(1 - u1), b = Math.sqrt(u1);
    const q = [a * Math.sin(2 * Math.PI * u2), a * Math.cos(2 * Math.PI * u2), b * Math.sin(2 * Math.PI * u3), b * Math.cos(2 * Math.PI * u3)];
    const [x, y, z, w] = q;
    return [
      [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
      [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
      [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ];
  }

  /**
   * Inserisce `count` copie di una molecola (coordinate in Å) in posizioni libere, orientate a caso.
   * Restituisce quante copie sono state inserite (la scatola può essere piena).
   */
  addMolecule(template, count = 1, T = this.T, at = null) {
    const n0 = template.Z.length;
    // centro della molecola
    const c = [0, 1, 2].map(k => template.pos.reduce((s, _, i) => (i % 3 === k ? s + template.pos[i] : s), 0) / n0);
    const local = template.Z.map((_, i) => [0, 1, 2].map(k => template.pos[3 * i + k] - c[k]));
    const radius = Math.max(0, ...local.map(v => Math.hypot(...v)));
    const half = this.box / 2 - 1.2 - radius;
    let placed = 0;
    for (let copy = 0; copy < count; copy++) {
      let ok = false;
      for (let attempt = 0; attempt < 400 && !ok; attempt++) {
        const Rm = this.randomRotation();
        const center = at && attempt === 0 ? at : [0, 1, 2].map(() => (2 * this.rng.uni() - 1) * Math.max(half, 0));
        const coords = local.map(v => [0, 1, 2].map(k => center[k] + Rm[k][0] * v[0] + Rm[k][1] * v[1] + Rm[k][2] * v[2]));
        const minD = at ? 1.6 : 2.4;
        ok = coords.every(p => {
          for (let i = 0; i < this.N; i++) {
            const d = Math.hypot(p[0] - this.pos[3 * i], p[1] - this.pos[3 * i + 1], p[2] - this.pos[3 * i + 2]);
            if (d < minD) return false;
          }
          return p.every(x => Math.abs(x) < this.box / 2 - 0.5);
        });
        if (!ok) continue;
        const start = this.N;
        this.Z = [...this.Z, ...template.Z];
        this.resize(this.Z.length);
        coords.forEach((p, k) => { this.pos.set(p, 3 * (start + k)); });
        // velocità del centro di massa di Maxwell–Boltzmann per l'intera molecola, più moto interno termico
        for (let k = 0; k < n0; k++) this.thermalize(start + k, T);
      }
      if (ok) placed++;
    }
    return placed;
  }

  /** Inserisce atomi con posizioni (Å) e velocità (Å/fs) assegnate. */
  addAtoms(list) {
    const start = this.N;
    this.Z = [...this.Z, ...list.map(a => a.Z)];
    this.resize(this.Z.length);
    list.forEach((a, k) => {
      this.pos.set(a.pos, 3 * (start + k));
      this.vel.set(a.vel ?? [0, 0, 0], 3 * (start + k));
    });
  }

  remove(indices) {
    const drop = new Set(indices);
    const keep = this.Z.map((_, i) => i).filter(i => !drop.has(i));
    const pos = new Float64Array(3 * keep.length), vel = new Float64Array(3 * keep.length);
    keep.forEach((i, k) => { pos.set(this.pos.subarray(3 * i, 3 * i + 3), 3 * k); vel.set(this.vel.subarray(3 * i, 3 * i + 3), 3 * k); });
    this.Z = keep.map(i => this.Z[i]);
    this.resize(this.Z.length);
    this.pos = pos; this.vel = vel;
    if (this.grab && drop.has(this.grab.i)) this.grab = null;
    else if (this.grab) this.grab.i = keep.indexOf(this.grab.i);
  }

  kinetic() {
    let K = 0;
    for (let i = 0; i < this.N; i++) {
      const vx = this.vel[3 * i], vy = this.vel[3 * i + 1], vz = this.vel[3 * i + 2];
      K += 0.5 * this.mass[i] * (vx * vx + vy * vy + vz * vz);
    }
    return K * MV2;
  }

  temperature() {
    return this.N ? 2 * this.kinetic() / (3 * this.N * KB_EV) : 0;
  }

  forces() {
    const { N, pos, F } = this;
    if (!N) { this.res = { E: 0, parts: { bond: 0, angle: 0, vdw: 0, es: 0 }, q: new Float64Array(0), bonds: [] }; this.Ewall = 0; return; }
    this.res = this.provider ? this.provider.compute(this.Z, pos, F) : this.ff.compute(this.Z, pos, F);
    // pareti morbide: E = ½ k d² per ogni atomo oltre il bordo della scatola
    const h = this.box / 2;
    let Ew = 0, Fw = 0;
    for (let k = 0; k < 3 * N; k++) {
      const x = pos[k];
      if (x > h) { const d = x - h; F[k] -= K_WALL * d; Ew += 0.5 * K_WALL * d * d; Fw += K_WALL * d; }
      else if (x < -h) { const d = -h - x; F[k] += K_WALL * d; Ew += 0.5 * K_WALL * d * d; Fw += K_WALL * d; }
    }
    this.Ewall = Ew;
    this.wallForce = Fw;
    if (this.grab && this.grab.i < N) {
      const i = this.grab.i;
      for (let c = 0; c < 3; c++) F[3 * i + c] += K_GRAB * (this.grab.target[c] - pos[3 * i + c]);
    }
  }

  /** Termostato CSVR (Bussi, Donadio, Parrinello, J. Chem. Phys. 2007). */
  bussi() {
    const N = this.N;
    if (!N || !this.thermostat) return;
    const nf = 3 * N;
    const K = this.kinetic();
    if (K <= 0) { for (let i = 0; i < N; i++) this.thermalize(i, this.T); return; }
    const Kt = 0.5 * nf * KB_EV * this.T;
    const c = Math.exp(-this.dt / this.tau);
    const r1 = this.rng.gauss();
    const s2 = nf > 1 ? 2 * this.rng.gamma((nf - 1) / 2) : 0;
    const Knew = K + (1 - c) * (Kt * (r1 * r1 + s2) / nf - K) + 2 * r1 * Math.sqrt(c * (1 - c) * Kt * K / nf);
    const alpha = Math.sqrt(Math.max(Knew, 0) / K);
    for (let k = 0; k < 3 * N; k++) this.vel[k] *= alpha;
    this.heatBath += Knew - K;
  }

  step() {
    const { N, dt, pos, vel, mass } = this;
    if (!this.res) this.forces();
    const F = this.F;
    // velocity Verlet: v(t+½dt) = v + ½ a dt;  x(t+dt) = x + v dt;  v(t+dt) = v(t+½dt) + ½ a(t+dt) dt
    const maxStep = 0.25;
    for (let i = 0; i < N; i++) {
      const a = 0.5 * dt * ACC / mass[i];
      for (let c = 0; c < 3; c++) vel[3 * i + c] += a * F[3 * i + c];
      // limite di sicurezza allo spostamento per passo (atomi sovrapposti, urti violentissimi)
      const vx = vel[3 * i], vy = vel[3 * i + 1], vz = vel[3 * i + 2];
      const disp = Math.sqrt(vx * vx + vy * vy + vz * vz) * dt;
      if (disp > maxStep) {
        const f = maxStep / disp;
        const before = 0.5 * mass[i] * (vx * vx + vy * vy + vz * vz) * MV2;
        vel[3 * i] *= f; vel[3 * i + 1] *= f; vel[3 * i + 2] *= f;
        this.work -= before * (1 - f * f);
        this.clamped++;
      }
      for (let c = 0; c < 3; c++) pos[3 * i + c] += dt * vel[3 * i + c];
    }
    let Wgrab = 0;
    if (this.grab && this.grab.i < N) {
      const i = this.grab.i;
      for (let c = 0; c < 3; c++) Wgrab += K_GRAB * (this.grab.target[c] - pos[3 * i + c]) * vel[3 * i + c] * dt;
    }
    this.forces();
    for (let i = 0; i < N; i++) {
      const a = 0.5 * dt * ACC / mass[i];
      for (let c = 0; c < 3; c++) vel[3 * i + c] += a * F[3 * i + c];
    }
    this.work += Wgrab;
    this.bussi();
    this.time += dt;
    this.stepCount++;
    // media mobile esponenziale della forza sulle pareti (memoria ≈ 2 ps)
    const decay = Math.exp(-dt / 2000);
    this.pAccum = this.pAccum * decay + this.wallForce;
    this.pSamples = this.pSamples * decay + 1;
    if (this.stepCount % this.censusEvery === 0) this.census();
  }

  /** Pressione sulle pareti: forza totale media / area totale (6 L²), in bar. */
  measurePressure() {
    if (this.pSamples) this.pressure = this.pAccum / this.pSamples / (6 * this.box * this.box) * EV_A3_TO_BAR;
    return this.pressure;
  }

  /** Scintilla: gli atomi entro `radius` da `center` vengono portati alla temperatura Tspark. */
  spark(center, radius = 4, Tspark = 6000) {
    const K0 = this.kinetic();
    for (let i = 0; i < this.N; i++) {
      const d = Math.hypot(this.pos[3 * i] - center[0], this.pos[3 * i + 1] - center[1], this.pos[3 * i + 2] - center[2]);
      if (d < radius) this.thermalize(i, Tspark);
    }
    this.work += this.kinetic() - K0;
  }

  /**
   * Assorbimento di un fotone di lunghezza d'onda λ (nm): l'energia E = hc/λ va in un legame scelto a caso,
   * come moto relativo dei due atomi lungo l'asse del legame. Se E supera l'energia di legame, il legame si rompe
   * (fotolisi, come Cl₂ + hν → 2 Cl·).
   */
  photon(lambdaNm, filter = null) {
    const E = 1239.841984 / lambdaNm; // eV
    const bonds = (this.res?.bonds ?? []).filter(b => b.n * b.w > 0.5 && (!filter || filter(this.Z[b.i], this.Z[b.j])));
    if (!bonds.length) return null;
    const b = bonds[Math.floor(this.rng.uni() * bonds.length)];
    const { i, j } = b;
    const u = [0, 1, 2].map(c => this.pos[3 * j + c] - this.pos[3 * i + c]);
    const r = Math.hypot(...u);
    u.forEach((_, c) => { u[c] /= r; });
    const mi = this.mass[i], mj = this.mass[j], mu = mi * mj / (mi + mj);
    const dv = Math.sqrt(2 * E / (mu * MV2));
    const K0 = this.kinetic();
    for (let c = 0; c < 3; c++) {
      this.vel[3 * i + c] -= mu / mi * dv * u[c];
      this.vel[3 * j + c] += mu / mj * dv * u[c];
    }
    this.work += this.kinetic() - K0;
    return { i, j, E };
  }

  /**
   * Specie chimiche: frammenti connessi dai legami (ordine > 0,55 per formarsi, < 0,35 per rompersi:
   * l'isteresi evita di contare come reazioni le vibrazioni). Le reazioni sono i cambi di connettività.
   */
  census() {
    const N = this.N;
    const bonds = this.res?.bonds ?? [];
    const now = new Set();
    const parent = Int32Array.from({ length: N }, (_, i) => i);
    const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    for (const b of bonds) {
      const key = b.i * N + b.j;
      const s = b.n * b.w;
      if (s > 0.55 || (s > 0.35 && this.bondedPrev.has(key))) {
        now.add(key);
        parent[find(b.i)] = find(b.j);
      }
    }
    this.bondedPrev = now;
    const groups = new Map();
    for (let i = 0; i < N; i++) {
      const r = find(i);
      if (!groups.has(r)) groups.set(r, []);
      groups.get(r).push(i);
    }
    const fragOf = new Array(N);
    const frags = [];
    for (const atoms of groups.values()) {
      const formula = hillFormula(atoms.map(i => ({ Z: this.Z[i] })));
      const id = atoms.join(',');
      const f = { atoms, formula, id };
      frags.push(f);
      for (const i of atoms) fragOf[i] = f;
    }
    // reazioni: frammenti scomparsi → frammenti nuovi, raggruppati per atomi in comune
    if (this.fragOf && this.fragOf.length === N) {
      const oldIds = new Set(this.fragOf.map(f => f.id));
      const newIds = new Set(frags.map(f => f.id));
      const gone = [...new Map(this.fragOf.filter(f => !newIds.has(f.id)).map(f => [f.id, f])).values()];
      const born = frags.filter(f => !oldIds.has(f.id));
      if (gone.length || born.length) {
        const up = new Map();
        const fnd = (x) => { while (up.get(x) !== x) x = up.get(x); return x; };
        const all = [...gone.map(f => 'o' + f.id), ...born.map(f => 'n' + f.id)];
        all.forEach(k => up.set(k, k));
        for (const f of born) {
          for (const i of f.atoms) {
            const o = this.fragOf[i];
            if (!newIds.has(o.id)) up.set(fnd('o' + o.id), fnd('n' + f.id));
          }
        }
        const ev = new Map();
        for (const f of gone) { const r = fnd('o' + f.id); if (!ev.has(r)) ev.set(r, { reactants: [], products: [] }); ev.get(r).reactants.push(f.formula); }
        for (const f of born) { const r = fnd('n' + f.id); if (!ev.has(r)) ev.set(r, { reactants: [], products: [] }); ev.get(r).products.push(f.formula); }
        for (const e of ev.values()) {
          e.reactants.sort(); e.products.sort();
          const same = e.reactants.join('+') === e.products.join('+');
          this.events.push({ t: this.time, reactants: e.reactants, products: e.products, exchange: same });
        }
        if (this.events.length > 400) this.events.splice(0, this.events.length - 400);
      }
    }
    this.fragOf = fragOf;
    const counts = new Map();
    for (const f of frags) counts.set(f.formula, (counts.get(f.formula) ?? 0) + 1);
    this.species = counts;
    this.frags = frags;
    const E = this.totalEnergy();
    this.history.push({ t: this.time, counts: Object.fromEntries(counts), T: this.temperature(), E, Ep: (this.res?.E ?? 0) + this.Ewall, Ec: E - this.heatBath - this.work, P: this.pressure });
    if (this.history.length > 2000) this.history.splice(0, this.history.length - 2000);
  }

  totalEnergy() {
    return this.kinetic() + (this.res?.E ?? 0) + this.Ewall;
  }
}
