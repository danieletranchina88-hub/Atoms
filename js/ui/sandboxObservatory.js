import * as THREE from 'three';
import { geometry, bondOrder, populationChange } from '../chem/reactionTrace.js';
import { ELEMENTS } from '../physics/elements.js';
import { drawXY } from './chemCharts.js';

const $ = id => document.getElementById(id);
const fmt = (v, d = 3) => Number.isFinite(v) ? v.toLocaleString('it-IT', { maximumFractionDigits: d, minimumFractionDigits: d }) : '—';
const sign = v => `${v >= 0 ? '+' : '−'}${fmt(Math.abs(v))}`;
const label = (Z, i) => `${ELEMENTS[Z[i] - 1]?.symbol ?? '?'}${i + 1}`;
const colors = ['#4b9dff', '#f0b429', '#e277a8', '#69c998'];
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

export class SandboxObservatory {
  constructor({ post, display, live, select, view, reference }) {
    Object.assign(this, { post, display, showLive: live, select, view, setReference: reference });
    this.samples = []; this.events = []; this.epoch = null; this.replay = false; this.cursor = null;
    this.ids = []; this.options = { labels: true, trails: false, vectors: 'none' };
    this.reqId = 0; this.lastDraw = 0; this.lastLabelTime = -Infinity; this.orbitals = []; this.slice = null;
  }
  receive(trace, frame) {
    if (!trace) return;
    const epochChanged = this.epoch !== trace.epoch;
    if (epochChanged) {
      this.epoch = trace.epoch; this.samples = []; this.events = []; this.cursor = null;
      this.replay = false; this.ids = frame.N === 3 ? [0, 1, 2] : frame.N > 1 ? [0, 1] : frame.N ? [0] : [];
      this.setReference(null); this.orbitals = []; this.slice = null; this.orbitalKey = '';
    }
    this.meta = trace; this.liveFrame = frame;
    const seen = new Set(this.samples.map(f => f.id));
    for (const f of trace.samples) if (!seen.has(f.id)) this.samples.push(f);
    this.samples = this.samples.filter(f => f.id >= trace.first);
    this.events = trace.events;
    if (this.replay && !this.samples.some(f => f.id === this.cursor)) this.live();
    if (!this.replay) this.current = frame;
    if (epochChanged && this.root?.isConnected) this.mount(this.root);
    else this.refresh();
  }
  inspect(id) {
    if (!Number.isInteger(id) || !this.samples.some(f => f.id === id)) return;
    this.replay = true; this.cursor = id; this.reqId++;
    this.post({ type: 'trace-inspect', id, epoch: this.epoch, reqId: this.reqId });
  }
  accept(m) {
    if (m.reqId !== this.reqId || m.frame.epoch !== this.epoch || !this.replay) return;
    this.current = m.frame; this.display(m.frame); this.refresh(true);
  }
  live() {
    this.reqId++; this.replay = false; this.cursor = null;
    this.current = this.liveFrame; this.showLive(this.liveFrame); this.refresh(true);
  }
  pick(i) {
    if (!this.current || i < 0 || i >= this.current.N) return;
    if (this.ids[0] === i) return;
    this.ids[0] = i;
    this.ids = [...new Set(this.ids)];
    this.mount(this.root);
  }
  mount(root) {
    this.root = root;
    if (!root) return;
    const f = this.current ?? this.liveFrame;
    const atoms = f ? Array.from(f.Z, (_, i) => `<option value="${i}">${label(f.Z, i)}</option>`).join('') : '';
    root.innerHTML = `<div class="sbo-heading"><h3>Microscopio di reazione</h3><span id="sbo-status" class="desc-muted"></span></div>
      <p class="hint">Segui gli stessi nuclei lungo i passi calcolati. I grafici e il replay usano campioni reali della traiettoria.</p>
      <div class="sbo-toolbar">
        <button class="btn" id="sbo-prev" title="Campione precedente">←</button>
        <input id="sbo-time" type="range" min="0" max="0" value="0" aria-label="Fotogramma della traiettoria">
        <button class="btn" id="sbo-next" title="Campione successivo">→</button>
        <button class="btn" id="sbo-live">Vista attuale</button>
        <button class="btn" id="sbo-reference">Fissa riferimento Δρ</button>
        <button class="btn" id="sbo-download">Traiettoria JSON</button>
      </div>
      <div class="sbo-toolbar">
        <label><input id="sbo-record" type="checkbox" ${this.meta?.config.enabled !== false ? 'checked' : ''}> Registra</label>
        <label><input id="sbo-waves" type="checkbox" ${this.meta?.config.waves !== false ? 'checked' : ''}> Salva densità (≤40 atomi)</label>
        <label><input id="sbo-stop" type="checkbox" ${this.meta?.config.pauseOnBond ? 'checked' : ''}> Pausa a un cambio di legame</label>
        <label for="sbo-stride">Campiona ogni</label><select id="sbo-stride">${[[0, 'automatico'], [1, '1 passo'], [5, '5 passi'], [20, '20 passi'], [100, '100 passi']].map(([n, t]) => `<option value="${n}" ${(this.meta?.config.stride ?? 0) === n ? 'selected' : ''}>${t}</option>`).join('')}</select>
        <button class="btn" id="sbo-reset">Nuova registrazione</button>
      </div>
      <div class="sbo-toolbar sbo-measure">${['A', 'B', 'C', 'D'].map((t, k) => `<label for="sbo-atom-${k}">${t}</label><select id="sbo-atom-${k}"><option value="">—</option>${atoms}</select>`).join('')}
        <label><input id="sbo-labels" type="checkbox" ${this.options.labels ? 'checked' : ''}> Etichette</label>
        <label><input id="sbo-trails" type="checkbox" ${this.options.trails ? 'checked' : ''}> Tracce dei nuclei</label>
        <select id="sbo-vectors" aria-label="Vettori"><option value="none">Nessun vettore</option><option value="forces">Forze F = −∇U</option><option value="vel">Velocità dei nuclei</option></select>
      </div>
      <div id="sbo-readout" class="sbo-readout"></div>
      <div class="sbo-charts">
        <div><h4>Ordini di legame fra i nuclei scelti</h4><canvas id="sbo-bonds"></canvas></div>
        <div><h4>Popolazioni elettroniche rispetto al primo campione nel buffer</h4><canvas id="sbo-population"></canvas></div>
        <div><h4>Energia lungo la traiettoria</h4><canvas id="sbo-energy"></canvas></div>
        <div><h4>Distanza A–B lungo la traiettoria</h4><canvas id="sbo-distance"></canvas></div>
      </div>
      <p class="hint">Δn<sub>A</sub> = −Δq<sub>A</sub>: valori positivi indicano aumento della popolazione attribuita all'atomo. Le cariche parziali dipendono dall'analisi del modello. Le curve di energia comprendono bagno, pareti e lavoro esterno; un massimo temporale non identifica uno stato di transizione.</p>
      <div class="sbo-toolbar"><label for="sbo-mo">Orbitale da mostrare</label><select id="sbo-mo"><option>Carica una vista elettronica…</option></select>
        <label for="sbo-slice-axis">Sezione della densità</label><select id="sbo-slice-axis"><option value="z">xy (z costante)</option><option value="y">xz (y costante)</option><option value="x">yz (x costante)</option></select>
        <input id="sbo-slice-offset" type="number" step="0.1" value="${this.sliceOffset ?? 0}" aria-label="Posizione della sezione in angstrom"><span>Å</span>
      </div>
      <div class="sbo-density"><canvas id="sbo-slice"></canvas><div><h4>Sezione elettronica</h4><p id="sbo-slice-note" class="hint">Scegli nuvola elettronica, Δρ, spin o orbitali nella vista 3D. La sezione usa la stessa funzione d'onda.</p><p class="hint">Δρ confronta due configurazioni sulla stessa griglia di laboratorio: comprende sia ridistribuzione elettronica sia spostamento dei nuclei. Non è una corrente, né una traiettoria di elettroni.</p></div></div>
      <h4>Passaggi di legame con identità atomiche</h4><div id="sbo-event-list" class="table-scroll short"></div>
      <p class="hint">Il rilevatore controlla ogni passo accettato, anche i legami che cambiano dentro una singola molecola. Soglie di visualizzazione: formazione >0,55; rottura ≤0,35. Non assegnano un meccanismo o un prodotto sperimentale. Buffer: massimo 600 campioni o 32 MiB; i più vecchi vengono rimossi.</p>`;
    $('sbo-time').oninput = e => this.inspect(this.samples[+e.target.value]?.id);
    $('sbo-prev').onclick = () => this.move(-1);
    $('sbo-next').onclick = () => this.move(1);
    $('sbo-live').onclick = () => this.live();
    $('sbo-reference').onclick = () => this.post({ type: 'wave', purpose: 'reference', reqId: -1, ...(this.replay ? { traceId: this.cursor } : {}) });
    $('sbo-download').onclick = () => this.post({ type: 'trace-export' });
    for (const [id, key] of [['sbo-record', 'enabled'], ['sbo-waves', 'waves'], ['sbo-stop', 'pauseOnBond']])
      $(id).onchange = e => this.post({ type: 'trace-config', [key]: e.target.checked });
    $('sbo-stride').onchange = e => this.post({ type: 'trace-config', stride: +e.target.value || null });
    $('sbo-reset').onclick = () => { this.live(); this.post({ type: 'trace-reset' }); };
    for (let k = 0; k < 4; k++) {
      $(`sbo-atom-${k}`).value = this.ids[k] ?? '';
      $(`sbo-atom-${k}`).onchange = () => {
        this.ids = [0, 1, 2, 3].map(k => $(`sbo-atom-${k}`).value).filter(v => v !== '').map(Number);
        this.ids = [...new Set(this.ids)]; this.select(this.ids[0] ?? -1); this.mount(this.root);
      };
    }
    for (const key of ['labels', 'trails']) $(`sbo-${key}`).onchange = e => { this.options[key] = e.target.checked; this.refresh(true); };
    $('sbo-vectors').value = this.options.vectors;
    $('sbo-vectors').onchange = e => { this.options.vectors = e.target.value; this.refresh(true); };
    $('sbo-mo').onchange = e => { this.orbitalKey = e.target.value; this.view('orbital', this.orbitalKey); };
    $('sbo-slice-axis').value = this.sliceAxis ?? 'z';
    $('sbo-slice-axis').onchange = e => { this.sliceAxis = e.target.value; this.view(null); };
    $('sbo-slice-offset').onchange = e => { if (e.target.checkValidity()) { this.sliceOffset = +e.target.value; this.view(null); } };
    this.refresh(true); this.updateOrbitals(); this.drawSlice();
  }
  move(d) {
    const k = this.replay ? this.samples.findIndex(f => f.id === this.cursor) : this.samples.length - 1;
    this.inspect(this.samples[Math.max(0, Math.min(this.samples.length - 1, k + d))]?.id);
  }
  refresh(force = false) {
    if (!this.root?.isConnected) return;
    const f = this.current ?? this.liveFrame;
    if (this.root.querySelectorAll('#sbo-atom-0 option').length !== (f?.N ?? 0) + 1) { this.mount(this.root); return; }
    const now = performance.now();
    if (!force && now - this.lastDraw < 200) return;
    this.lastDraw = now;
    for (const [id, key] of [['sbo-record', 'enabled'], ['sbo-waves', 'waves'], ['sbo-stop', 'pauseOnBond']]) $(id).checked = !!this.meta?.config[key];
    const k = this.replay ? this.samples.findIndex(s => s.id === this.cursor) : this.samples.length - 1;
    $('sbo-time').max = Math.max(0, this.samples.length - 1); $('sbo-time').value = Math.max(0, k);
    $('sbo-status').textContent = `${this.replay ? 'Replay' : 'Attuale'} · ${fmt(f?.stats.t ?? 0, 2)} fs · ${this.samples.length} campioni · ${fmt((this.meta?.bytes ?? 0) / 1048576, 1)} MiB${this.meta?.evicted ? ` · ${this.meta.evicted} rimossi` : ''}`;
    $('sbo-prev').disabled = this.samples.length < 2 || (this.replay && k <= 0);
    $('sbo-next').disabled = this.samples.length < 2 || k >= this.samples.length - 1;
    $('sbo-reference').disabled = !f?.N;
    $('sbo-download').disabled = !this.samples.length;
    const g = f ? geometry(f.pos, this.ids) : {};
    const delta = f ? populationChange(f.q, this.samples[0]?.q) : null;
    const fields = [`A–B: ${fmt(g.distance)} Å`, `∠ABC: ${fmt(g.angle, 1)}°`, `τABCD: ${fmt(g.dihedral, 1)}°`];
    if (f && this.ids.length >= 2) fields.push(`BO(A–B): ${fmt(bondOrder(f.bonds, this.ids[0], this.ids[1]))}`);
    if (delta) fields.push(`ΣΔn: ${sign(delta.total)} e`);
    if (this.options.vectors !== 'none') fields.push(`Frecce riscalate · ${this.options.vectors === 'forces' ? 'eV/Å' : 'Å/fs'} nelle etichette`);
    $('sbo-readout').textContent = fields.join(' · ');
    this.drawCharts(); this.drawEvents();
  }
  drawCharts() {
    if (!this.samples.length) return;
    const xs = this.samples.map(f => f.stats.t), xmin = xs[0], xmax = Math.max(xs.at(-1), xmin + 1);
    const chart = (id, series, ylabel) => {
      const values = series.flatMap(s => s.ys).filter(Number.isFinite);
      const lo = values.length ? Math.min(...values) : 0, hi = values.length ? Math.max(...values) : 1;
      const pad = Math.max(0.02, (hi - lo) * 0.1);
      drawXY($(id), { series, xmin, xmax, ymin: lo - pad, ymax: hi + pad, xlabel: 't (fs)', ylabel,
        vlines: this.replay ? [{ x: this.current?.stats.t, color: css('--text'), label: 'replay' }] : [] });
    };
    const pairs = [];
    for (let a = 0; a < this.ids.length; a++) for (let b = a + 1; b < this.ids.length; b++) pairs.push([this.ids[a], this.ids[b]]);
    chart('sbo-bonds', pairs.map(([i, j], k) => ({ xs, ys: this.samples.map(f => bondOrder(f.bonds, i, j)), color: colors[k % 4], label: `${label(this.current.Z, i)}–${label(this.current.Z, j)}` })), 'BO (modello)');
    const initial = this.samples[0];
    chart('sbo-population', this.ids.map((i, k) => ({ xs, ys: this.samples.map(f => initial.q[i] - f.q[i]), color: colors[k], label: label(this.current.Z, i) })), 'Δn (e)');
    chart('sbo-energy', [['Epot', 'potenziale', colors[0]], ['Ekin', 'cinetica', colors[1]], ['Etot', 'totale', colors[2]]].map(([key, name, color]) => ({ xs, ys: this.samples.map(f => f.stats[key] - initial.stats[key]), color, label: `Δ${name}` })), 'ΔE (eV/scatola)');
    chart('sbo-distance', this.ids.length >= 2 ? [{ xs, ys: this.samples.map(f => geometry(f.pos, this.ids.slice(0, 2)).distance), color: colors[0], label: 'A–B' }] : [], 'r (Å)');
  }
  drawEvents() {
    const f = this.current;
    const events = this.events.slice(-30).reverse();
    $('sbo-event-list').innerHTML = events.length && f ? `<table class="data-table"><thead><tr><th>t (fs)</th><th>Passaggio</th><th>Legame</th><th>BO prima → dopo</th><th>r (Å)</th><th></th></tr></thead><tbody>${events.map(e => `<tr><td>${fmt(e.t, 2)}</td><td>${e.kind === 'formed' ? 'formazione' : 'rottura'}</td><td>${label(f.Z, e.i)}–${label(f.Z, e.j)}</td><td>${fmt(e.before, 2)} → ${fmt(e.after, 2)}</td><td>${fmt(e.distance)}</td><td><button class="linkish" data-time="${e.t}" data-atoms="${e.i},${e.j}">Rivedi</button></td></tr>`).join('')}</tbody></table>` : '<p class="hint">Nessun attraversamento delle soglie di legame rilevato.</p>';
    $('sbo-event-list').querySelectorAll('[data-time]').forEach(b => b.onclick = () => {
      const t = +b.dataset.time;
      const sample = this.samples.reduce((a, c) => !a || Math.abs(c.stats.t - t) < Math.abs(a.stats.t - t) ? c : a, null);
      this.ids = b.dataset.atoms.split(',').map(Number); this.select(this.ids[0]); this.mount(this.root);
      if (sample) this.inspect(sample.id);
    });
  }
  wave(w, t) {
    this.orbitals = w?.orbitals ?? []; this.eUnit = w?.eUnit ?? 1; this.waveTime = t;
    this.updateOrbitals();
  }
  updateOrbitals() {
    const el = $('sbo-mo');
    if (!el) return;
    el.innerHTML = '<option value="homo">HOMO</option><option value="lumo">LUMO</option>' + this.orbitals.map(o => `<option value="${o.id}">MO ${o.index + 1} · ${o.spin} · ${fmt(o.e * this.eUnit, 2)} eV · occ ${fmt(o.occ, 2)}</option>`).join('');
    el.value = this.orbitalKey || 'homo';
  }
  densitySlice(slice, meta) {
    this.slice = slice; this.sliceMeta = meta; this.drawSlice();
  }
  drawSlice() {
    const c = $('sbo-slice'), g = this.slice;
    if (!c || !g) return;
    c.width = 360; c.height = 360;
    const ctx = c.getContext('2d'), temp = document.createElement('canvas'); temp.width = temp.height = g.res;
    const t = temp.getContext('2d'), image = t.createImageData(g.res, g.res);
    const signed = this.sliceMeta.what !== 'density';
    const max = Math.max(1e-12, ...g.values.map(Math.abs));
    for (let y = 0; y < g.res; y++) for (let x = 0; x < g.res; x++) {
      const v = g.values[x + g.res * y];
      const w = Math.log1p(Math.abs(v) / max * 100) / Math.log(101);
      const dst = 4 * (x + g.res * (g.res - 1 - y));
      const rgb = signed && v < 0 ? [235, 92, 108] : [73, 166, 248];
      for (let k = 0; k < 3; k++) image.data[dst + k] = Math.round(18 * (1 - w) + rgb[k] * w);
      image.data[dst + 3] = 255;
    }
    t.putImageData(image, 0, 0); ctx.drawImage(temp, 0, 0, 360, 360);
    const axes = g.axis === 'z' ? [0, 1, 2] : g.axis === 'y' ? [0, 2, 1] : [1, 2, 0];
    const p = this.sliceMeta.pos;
    ctx.font = '12px system-ui'; ctx.textAlign = 'center'; ctx.fillStyle = '#fff';
    for (let i = 0; i < this.sliceMeta.Z.length; i++) if (Math.abs(p[3 * i + axes[2]] - g.offset) < 0.5) {
      const x = 180 + p[3 * i + axes[0]] / (2 * g.half) * 360, y = 180 - p[3 * i + axes[1]] / (2 * g.half) * 360;
      ctx.beginPath(); ctx.arc(x, y, 2, 0, 2 * Math.PI); ctx.fill(); ctx.fillText(label(this.sliceMeta.Z, i), x, y - 6);
    }
    ctx.textAlign = 'left'; ctx.fillText(`−${fmt(g.half, 1)} Å`, 8, 350); ctx.textAlign = 'right'; ctx.fillText(`+${fmt(g.half, 1)} Å`, 352, 350);
    const names = { density: 'ρ', spin: 'ρα−ρβ', difference: 'Δρ', orbital: 'ψ' };
    $('sbo-slice-note').textContent = `${names[this.sliceMeta.what]} · ${g.axis} = ${fmt(g.offset, 2)} Å · t = ${fmt(this.sliceMeta.t, 2)} fs. Scala cromatica logaritmica simmetrica; max |valore| = ${fmt(max, 4)} ${this.sliceMeta.what === 'orbital' ? 'bohr⁻³ᐟ²' : 'e/bohr³'}. ${signed ? 'Blu: positivo; rosso: negativo.' : 'Blu: densità crescente.'} Nuclei marcati se entro 0,5 Å dal piano.`;
  }
  attach3D(parent) {
    this.labelKey = ''; this.lastLabelTime = -Infinity;
    this.group = new THREE.Group(); parent.add(this.group);
    this.labels = new THREE.Group(); this.group.add(this.labels);
    this.lines = new THREE.Group(); this.group.add(this.lines);
  }
  update3D(f, selected, viewer) {
    if (!this.group || !f) return;
    const dispose = group => { for (const c of [...group.children]) { c.traverse(o => { o.geometry?.dispose(); o.material?.map?.dispose(); o.material?.dispose(); }); group.remove(c); } };
    dispose(this.lines);
    const ids = [...new Set([...this.ids, selected])].filter(i => i >= 0 && i < f.N).slice(0, 8);
    const point = (p, i) => new THREE.Vector3(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
    for (const i of ids) {
      if (this.options.trails) {
        const samples = this.samples.filter(s => s.stats.t <= f.stats.t).slice(-100);
        if (samples.length > 1) this.lines.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(samples.map(s => point(s.pos, i))), new THREE.LineBasicMaterial({ color: colors[ids.indexOf(i) % 4], transparent: true, opacity: 0.65 })));
      }
    }
    const vector = f[this.options.vectors];
    if (vector) {
      const maximum = Math.max(1e-12, ...ids.map(i => point(vector, i).length()));
      for (const i of ids) {
        const v = point(vector, i), n = v.length();
        if (n > 1e-12) this.lines.add(new THREE.ArrowHelper(v.normalize(), point(f.pos, i), n / maximum * Math.min(2, f.stats.box / 4), this.options.vectors === 'forces' ? '#ffad66' : '#69c998', 0.15, 0.08));
      }
    }
    // Explicit marker for a recently detected bond crossing, including broken bonds now absent from the renderer.
    for (const e of this.events.filter(e => e.t <= f.stats.t && f.stats.t - e.t < 3).slice(-8)) {
      const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints([point(f.pos, e.i), point(f.pos, e.j)]), new THREE.LineDashedMaterial({ color: e.kind === 'formed' ? '#69ee9b' : '#ff6b57', dashSize: 0.12, gapSize: 0.08 }));
      line.computeLineDistances(); this.lines.add(line);
    }
    const key = `${ids}|${this.options.labels}|${this.options.vectors}|${f.stats.t}|${f.epoch}`;
    if (this.labelKey !== key && performance.now() - this.lastLabelTime > 100) {
      this.labelKey = key; this.lastLabelTime = performance.now(); dispose(this.labels);
      if (this.options.labels) for (const i of ids) {
        const c = document.createElement('canvas'); c.width = 512; c.height = 160;
        const ctx = c.getContext('2d'); ctx.fillStyle = 'rgba(15,22,33,0.86)'; ctx.fillRect(0, 0, 512, 160);
        ctx.fillStyle = '#fff'; ctx.font = '50px system-ui'; ctx.fillText(`${label(f.Z, i)}  q ${sign(f.q[i])} e`, 12, 60);
        ctx.font = '34px system-ui';
        const vec = f[this.options.vectors];
        ctx.fillText(vec ? `${this.options.vectors === 'forces' ? '|F|' : '|v|'} ${fmt(point(vec, i).length(), 4)} ${this.options.vectors === 'forces' ? 'eV/Å' : 'Å/fs'}` : f.spin ? `nα−nβ = ${sign(f.spin[i])}` : 'identità del nucleo', 12, 125);
        const texture = new THREE.CanvasTexture(c); texture.colorSpace = THREE.SRGBColorSpace;
        const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, depthTest: false, depthWrite: false, transparent: true }));
        sprite.renderOrder = 50;
        sprite.userData.i = i; sprite.scale.set(3.2, 1, 1); this.labels.add(sprite);
      }
    }
    const up = new THREE.Vector3(0, 1, 0).applyQuaternion(viewer.camera.quaternion);
    for (const [k, sprite] of this.labels.children.entries()) {
      const nucleus = point(f.pos, sprite.userData.i);
      sprite.position.copy(nucleus).addScaledVector(up, k % 2 ? -1.3 : 1.3);
      const leader = new THREE.Line(new THREE.BufferGeometry().setFromPoints([nucleus, sprite.position]), new THREE.LineDashedMaterial({ color: '#778899', dashSize: 0.06, gapSize: 0.06, transparent: true, opacity: 0.65 }));
      leader.computeLineDistances(); this.lines.add(leader);
    }
  }
}

export function downloadTrajectory(data) {
  const body = JSON.stringify(data, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v);
  const url = URL.createObjectURL(new Blob([body], { type: 'application/json' }));
  const a = document.createElement('a'); a.href = url; a.download = 'atoms-sandbox-traiettoria.json'; a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
