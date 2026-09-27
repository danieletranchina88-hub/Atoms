// Grafici per molecole e reazioni: struttura di Lewis, diagramma degli OM, spettro IR, entalpia.

import { ELEMENTS } from '../physics/elements.js';

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

function setup(canvas) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

const FONT = '"IBM Plex Sans", system-ui, sans-serif';
const MONO = '"IBM Plex Mono", ui-monospace, monospace';
const fmt = (v, d = 2) => Number(v).toLocaleString('it-IT', { minimumFractionDigits: d, maximumFractionDigits: d });

// ---------------------------------------------------------------------------
// Struttura di Lewis
// ---------------------------------------------------------------------------

/** Proiezione 2D sul piano di massima estensione della molecola. */
export function projectTo2D(atoms) {
  const n = atoms.length;
  const c = [0, 0, 0];
  atoms.forEach(a => { for (let k = 0; k < 3; k++) c[k] += a.xyz[k] / n; });
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  atoms.forEach(a => {
    const r = a.xyz.map((v, k) => v - c[k]);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) M[i][j] += r[i] * r[j];
  });
  // due autovettori principali per iterazione di potenza con deflazione
  const power = (A, avoid) => {
    let v = [0.577, 0.577, 0.577];
    if (avoid) v = [avoid[1] - avoid[2], avoid[2] - avoid[0], avoid[0] - avoid[1]];
    for (let it = 0; it < 200; it++) {
      let w = [0, 1, 2].map(i => A[i][0] * v[0] + A[i][1] * v[1] + A[i][2] * v[2]);
      if (avoid) { const d = w[0] * avoid[0] + w[1] * avoid[1] + w[2] * avoid[2]; w = w.map((x, k) => x - d * avoid[k]); }
      const l = Math.hypot(...w) || 1;
      v = w.map(x => x / l);
    }
    return v;
  };
  const e1 = power(M, null);
  let e2 = power(M, e1);
  if (Math.hypot(...e2) < 1e-6 || !Number.isFinite(e2[0])) e2 = Math.abs(e1[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  return atoms.map(a => {
    const r = a.xyz.map((v, k) => v - c[k]);
    return [r[0] * e1[0] + r[1] * e1[1] + r[2] * e1[2], r[0] * e2[0] + r[1] * e2[1] + r[2] * e2[2]];
  });
}

export function drawLewis(canvas, { atoms, bonds, orders, lone, formalCharges, resonanceCount }) {
  const { ctx, w, h } = setup(canvas);
  const ink = cssVar('--text');
  const muted = cssVar('--muted');
  const accent = cssVar('--accent');
  const pts = projectTo2D(atoms);
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (const [x, y] of pts) { xmin = Math.min(xmin, x); xmax = Math.max(xmax, x); ymin = Math.min(ymin, y); ymax = Math.max(ymax, y); }
  const pad = 34;
  const sx = (w - 2 * pad) / Math.max(xmax - xmin, 1e-6);
  const sy = (h - 2 * pad - 16) / Math.max(ymax - ymin, 1e-6);
  const s = Math.min(sx, sy, 60);
  const cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2;
  const P = pts.map(([x, y]) => [w / 2 + (x - cx) * s, (h - 16) / 2 - (y - cy) * s]);
  const labelR = 10;

  // legami
  ctx.strokeStyle = ink;
  ctx.lineWidth = 1.6;
  bonds.forEach((b, k) => {
    const [x1, y1] = P[b.a], [x2, y2] = P[b.b];
    const dx = x2 - x1, dy = y2 - y1;
    const L = Math.hypot(dx, dy) || 1;
    const ux = dx / L, uy = dy / L;
    const nx = -uy, ny = ux;
    const o = orders[k];
    const offs = o === 3 ? [-4, 0, 4] : o === 2 ? [-2.5, 2.5] : [0];
    for (const off of offs) {
      ctx.beginPath();
      ctx.moveTo(x1 + ux * labelR + nx * off, y1 + uy * labelR + ny * off);
      ctx.lineTo(x2 - ux * labelR + nx * off, y2 - uy * labelR + ny * off);
      ctx.stroke();
    }
  });
  // simboli, coppie solitarie, elettroni spaiati, cariche formali
  atoms.forEach((a, i) => {
    const [x, y] = P[i];
    const sym = ELEMENTS[a.Z - 1].symbol;
    ctx.font = `600 15px ${FONT}`;
    ctx.fillStyle = ink;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(sym, x, y);
    // direzioni dei legami
    const dirs = [];
    bonds.forEach(b => {
      if (b.a === i || b.b === i) {
        const j = b.a === i ? b.b : b.a;
        dirs.push(Math.atan2(P[j][1] - y, P[j][0] - x));
      }
    });
    const nPairs = Math.floor(lone[i] / 2);
    const nSingle = lone[i] % 2;
    const items = [...Array(nPairs).fill(2), ...Array(nSingle).fill(1)];
    const taken = [...dirs];
    for (const e of items) {
      // angolo più lontano da legami e coppie già piazzate
      let bestA = 0, bestD = -1;
      for (let k = 0; k < 72; k++) {
        const ang = k * Math.PI / 36;
        let dmin = Math.PI;
        for (const t of taken) {
          let d = Math.abs(ang - t) % (2 * Math.PI);
          if (d > Math.PI) d = 2 * Math.PI - d;
          dmin = Math.min(dmin, d);
        }
        if (dmin > bestD + 1e-9) { bestD = dmin; bestA = ang; }
      }
      taken.push(bestA);
      const r = 16;
      const px = x + Math.cos(bestA) * r, py = y + Math.sin(bestA) * r;
      ctx.fillStyle = e === 1 ? accent : ink;
      if (e === 2) {
        const qx = -Math.sin(bestA) * 3, qy = Math.cos(bestA) * 3;
        ctx.beginPath(); ctx.arc(px + qx, py + qy, 1.7, 0, 2 * Math.PI); ctx.fill();
        ctx.beginPath(); ctx.arc(px - qx, py - qy, 1.7, 0, 2 * Math.PI); ctx.fill();
      } else {
        ctx.beginPath(); ctx.arc(px, py, 2.2, 0, 2 * Math.PI); ctx.fill();
      }
    }
    const fc = formalCharges[i];
    if (fc) {
      const txt = fc > 0 ? (fc > 1 ? `${fc}+` : '+') : (fc < -1 ? `${-fc}−` : '−');
      ctx.font = `600 11px ${FONT}`;
      const tx = x + 11, ty = y - 11;
      ctx.strokeStyle = fc > 0 ? cssVar('--phase-neg') : cssVar('--phase-pos');
      ctx.fillStyle = ctx.strokeStyle;
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(tx, ty, 6.5, 0, 2 * Math.PI); ctx.stroke();
      ctx.fillText(txt, tx, ty + 0.5);
    }
  });
  ctx.font = `11px ${FONT}`;
  ctx.fillStyle = muted;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'bottom';
  ctx.fillText(resonanceCount > 1 ? `una di ${resonanceCount} strutture di risonanza equivalenti` : 'struttura di Lewis', 4, h - 2);
}

// ---------------------------------------------------------------------------
// Diagramma degli orbitali molecolari calcolati
// ---------------------------------------------------------------------------

/**
 * levels: [{ index, e (eV), occ, spin }]; restituisce le regioni cliccabili [{x0,x1,y0,y1,index,spin}].
 */
export function drawMOLevels(canvas, { levels, selected, ncoreShown, unrestricted, title }) {
  const { ctx, w, h } = setup(canvas);
  const ink = cssVar('--text');
  const muted = cssVar('--muted');
  const accent = cssVar('--accent');
  const grid = cssVar('--chart-grid');
  const pad = { l: 50, r: 10, t: 12, b: 26 };
  const es = levels.map(l => l.e);
  let emin = Math.min(...es), emax = Math.max(...es);
  const span = Math.max(emax - emin, 1);
  emin -= 0.06 * span; emax += 0.08 * span;
  const Y = (e) => pad.t + (1 - (e - emin) / (emax - emin)) * (h - pad.t - pad.b);
  // asse
  ctx.font = `11px ${MONO}`;
  ctx.fillStyle = muted;
  ctx.strokeStyle = grid;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  const step = [1, 2, 5, 10, 20, 50].find(s => span / s <= 7) ?? 100;
  for (let e = Math.ceil(emin / step) * step; e <= emax; e += step) {
    const y = Y(e);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
    ctx.fillText(fmt(e, 0), pad.l - 6, y);
  }
  if (emin < 0 && emax > 0) {
    ctx.strokeStyle = muted;
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(pad.l, Y(0)); ctx.lineTo(w - pad.r, Y(0)); ctx.stroke();
    ctx.setLineDash([]);
  }
  ctx.font = `11px ${FONT}`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'bottom';
  ctx.fillText(title ?? 'ε (eV)', pad.l, h - 4);
  if (ncoreShown) {
    ctx.textAlign = 'right';
    ctx.fillText(ncoreShown, w - pad.r, h - 4);
  }

  const regions = [];
  const columns = unrestricted ? ['α', 'β'] : [null];
  const colW = (w - pad.l - pad.r) / columns.length;
  columns.forEach((spin, ci) => {
    const list = levels.filter(l => l.spin === spin).sort((a, b) => a.e - b.e);
    // raggruppa i livelli degeneri
    const groups = [];
    for (const l of list) {
      const g = groups[groups.length - 1];
      if (g && Math.abs(g[0].e - l.e) < 0.01) g.push(l); else groups.push([l]);
    }
    const x0 = pad.l + ci * colW;
    if (spin) {
      ctx.fillStyle = muted;
      ctx.font = `600 11px ${FONT}`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      ctx.fillText(`spin ${spin}`, x0 + colW / 2, 0);
    }
    const lineW = Math.min(34, colW / 5);
    for (const g of groups) {
      const y = Y(g[0].e);
      const total = g.length * lineW + (g.length - 1) * 6;
      let x = x0 + colW * 0.42 - total / 2;
      for (const l of g) {
        const isSel = selected && selected.index === l.index && selected.spin === l.spin;
        ctx.strokeStyle = isSel ? accent : l.occ > 0 ? ink : muted;
        ctx.lineWidth = isSel ? 3 : 2;
        ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + lineW, y); ctx.stroke();
        // elettroni
        ctx.lineWidth = 1.3;
        ctx.strokeStyle = ink;
        const arrow = (ax, up) => {
          const y0 = up ? y + 7 : y - 7, y1 = up ? y - 7 : y + 7, d = up ? 1 : -1;
          ctx.beginPath(); ctx.moveTo(ax, y0); ctx.lineTo(ax, y1);
          ctx.moveTo(ax - 2.5, y1 + 4 * d); ctx.lineTo(ax, y1); ctx.lineTo(ax + 2.5, y1 + 4 * d); ctx.stroke();
        };
        if (!spin) {
          if (l.occ >= 1) arrow(x + lineW / 2 - 4, true);
          if (l.occ >= 2) arrow(x + lineW / 2 + 4, false);
        } else if (l.occ >= 1) arrow(x + lineW / 2, spin === 'α');
        regions.push({ x0: x - 3, x1: x + lineW + 3, y0: y - 8, y1: y + 8, index: l.index, spin: l.spin });
        x += lineW + 6;
      }
      const l = g[0];
      ctx.fillStyle = muted;
      ctx.font = `10.5px ${MONO}`;
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      const lbl = l.label ? `${l.label} ` : '';
      ctx.fillText(`${lbl}${fmt(l.e, 2)}`, x0 + colW * 0.42 + total / 2 + 8, y);
    }
  });
  return regions;
}

// ---------------------------------------------------------------------------
// Spettro infrarosso
// ---------------------------------------------------------------------------

export function drawIR(canvas, { modes, scale = 1, exp = [], selected = null }) {
  const { ctx, w, h } = setup(canvas);
  const ink = cssVar('--text');
  const muted = cssVar('--muted');
  const accent = cssVar('--accent');
  const grid = cssVar('--chart-grid');
  const pad = { l: 40, r: 12, t: 16, b: 32 };
  const x0 = 4000, x1 = 400;
  const X = (nu) => pad.l + (x0 - nu) / (x0 - x1) * (w - pad.l - pad.r);
  const real = modes.filter(m => m.freq > 0);
  const Imax = Math.max(1, ...real.map(m => m.ir ?? 0));
  const Y = (v) => pad.t + (1 - v / (Imax * 1.12)) * (h - pad.t - pad.b);
  ctx.font = `11px ${MONO}`;
  ctx.fillStyle = muted;
  ctx.strokeStyle = grid;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (let nu = 4000; nu >= 500; nu -= 500) {
    const x = X(nu);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, h - pad.b); ctx.stroke();
    ctx.fillText(String(nu), x, h - pad.b + 4);
  }
  ctx.font = `11px ${FONT}`;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'bottom';
  ctx.fillText('numero d\'onda (cm⁻¹)', w - pad.r, h - 2);
  ctx.save();
  ctx.translate(12, (pad.t + h - pad.b) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText('intensità (km/mol)', 0, 0);
  ctx.restore();
  // spettro con righe lorentziane (semi-larghezza 12 cm⁻¹)
  const gamma = 12;
  ctx.beginPath();
  for (let px = pad.l; px <= w - pad.r; px++) {
    const nu = x0 - (px - pad.l) / (w - pad.l - pad.r) * (x0 - x1);
    let a = 0;
    for (const m of real) {
      const d = nu - m.freq * scale;
      a += (m.ir ?? 0) * gamma * gamma / (d * d + gamma * gamma);
    }
    const y = Y(a);
    if (px === pad.l) ctx.moveTo(px, y); else ctx.lineTo(px, y);
  }
  ctx.strokeStyle = accent;
  ctx.lineWidth = 1.6;
  ctx.stroke();
  // righe dei modi (anche quelli inattivi in IR)
  real.forEach((m) => {
    const x = X(m.freq * scale);
    const isSel = selected !== null && modes[selected] === m;
    ctx.strokeStyle = isSel ? accent : muted;
    ctx.lineWidth = isSel ? 2 : 1;
    ctx.beginPath(); ctx.moveTo(x, h - pad.b); ctx.lineTo(x, h - pad.b - 6); ctx.stroke();
  });
  // fondamentali sperimentali
  ctx.fillStyle = ink;
  for (const nu of exp) {
    const x = X(nu);
    ctx.beginPath(); ctx.moveTo(x, pad.t - 2); ctx.lineTo(x - 4, pad.t - 9); ctx.lineTo(x + 4, pad.t - 9); ctx.closePath(); ctx.fill();
  }
}

// ---------------------------------------------------------------------------
// Diagramma entalpico di una reazione
// ---------------------------------------------------------------------------

export function drawEnthalpy(canvas, { left, right, values }) {
  const { ctx, w, h } = setup(canvas);
  const ink = cssVar('--text');
  const muted = cssVar('--muted');
  const accent = cssVar('--accent');
  const grid = cssVar('--chart-grid');
  const pad = { l: 56, r: 16, t: 18, b: 30 };
  const all = [0, ...values.filter(v => v.value !== null && v.value !== undefined).map(v => v.value)];
  let vmin = Math.min(...all), vmax = Math.max(...all);
  const span = Math.max(vmax - vmin, 10);
  vmin -= 0.12 * span; vmax += 0.12 * span;
  const Y = (v) => pad.t + (1 - (v - vmin) / (vmax - vmin)) * (h - pad.t - pad.b);
  ctx.font = `11px ${MONO}`;
  ctx.fillStyle = muted;
  ctx.strokeStyle = grid;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  const step = [10, 20, 50, 100, 200, 500, 1000].find(s => span / s <= 6) ?? 2000;
  for (let v = Math.ceil(vmin / step) * step; v <= vmax; v += step) {
    const y = Y(v);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
    ctx.fillText(fmt(v, 0), pad.l - 6, y);
  }
  ctx.save();
  ctx.translate(12, (pad.t + h - pad.b) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.font = `11px ${FONT}`;
  ctx.textAlign = 'center';
  ctx.fillText('H (kJ/mol)', 0, 0);
  ctx.restore();
  const xa = pad.l + 16, xb = w - pad.r - 16;
  const lw = (xb - xa) * 0.32;
  // reagenti
  ctx.strokeStyle = ink;
  ctx.lineWidth = 2.5;
  ctx.beginPath(); ctx.moveTo(xa, Y(0)); ctx.lineTo(xa + lw, Y(0)); ctx.stroke();
  ctx.font = `12px ${FONT}`;
  ctx.fillStyle = ink;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'bottom';
  wrapText(ctx, left, xa, Y(0) - 4, lw);
  // prodotti per ciascun metodo
  values.forEach((v, k) => {
    if (v.value === null || v.value === undefined) return;
    const y = Y(v.value);
    ctx.strokeStyle = v.color ?? accent;
    ctx.lineWidth = v.dashed ? 1.5 : 2.5;
    ctx.setLineDash(v.dashed ? [5, 4] : []);
    ctx.beginPath(); ctx.moveTo(xb - lw, y); ctx.lineTo(xb, y); ctx.stroke();
    ctx.setLineDash([2, 3]);
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(xa + lw, Y(0)); ctx.lineTo(xb - lw, y); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = v.color ?? accent;
    ctx.font = `11px ${MONO}`;
    ctx.textAlign = 'right';
    ctx.textBaseline = k % 2 ? 'top' : 'bottom';
    ctx.fillText(`${v.label}: ${fmt(v.value, 1)}`, xb, y + (k % 2 ? 3 : -3));
  });
  ctx.fillStyle = ink;
  ctx.font = `12px ${FONT}`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  const yp = Y(values.find(v => v.value !== null && v.value !== undefined)?.value ?? 0);
  wrapText(ctx, right, xb - lw, Math.min(h - pad.b - 14, yp + 18), lw);
}

function wrapText(ctx, text, x, y, maxW) {
  const words = text.split(' ');
  let line = '';
  const lines = [];
  for (const wd of words) {
    const t = line ? `${line} ${wd}` : wd;
    if (ctx.measureText(t).width > maxW && line) { lines.push(line); line = wd; } else line = t;
  }
  if (line) lines.push(line);
  const lh = 14;
  const base = ctx.textBaseline === 'bottom' ? y - (lines.length - 1) * lh : y;
  lines.forEach((l, i) => ctx.fillText(l, x, base + i * lh));
}

/** Barre orizzontali con segno (contributi termodinamici). */
export function drawBars(canvas, { items, unit }) {
  const { ctx, w, h } = setup(canvas);
  const ink = cssVar('--text');
  const muted = cssVar('--muted');
  const grid = cssVar('--chart-grid');
  const pad = { l: 118, r: 60, t: 10, b: 22 };
  const vmax = Math.max(1, ...items.map(i => Math.abs(i.value)));
  const X = (v) => pad.l + (w - pad.l - pad.r) / 2 + v / vmax * (w - pad.l - pad.r) / 2;
  const rowH = (h - pad.t - pad.b) / items.length;
  ctx.strokeStyle = grid;
  ctx.beginPath(); ctx.moveTo(X(0), pad.t); ctx.lineTo(X(0), h - pad.b); ctx.stroke();
  items.forEach((it, k) => {
    const y = pad.t + k * rowH + rowH * 0.2;
    const bh = rowH * 0.6;
    ctx.fillStyle = it.color ?? (it.value < 0 ? cssVar('--phase-neg') : cssVar('--phase-pos'));
    const xa = X(Math.min(0, it.value)), xb = X(Math.max(0, it.value));
    ctx.fillRect(xa, y, Math.max(1, xb - xa), bh);
    ctx.fillStyle = ink;
    ctx.font = `12px ${FONT}`;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(it.label, pad.l - 8, y + bh / 2);
    ctx.font = `11px ${MONO}`;
    ctx.fillStyle = muted;
    ctx.textAlign = it.value < 0 ? 'right' : 'left';
    ctx.fillText(fmt(it.value, 1), it.value < 0 ? xa - 4 : xb + 4, y + bh / 2);
  });
  ctx.font = `11px ${FONT}`;
  ctx.fillStyle = muted;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'bottom';
  ctx.fillText(unit, w - 4, h - 2);
}
