// Grafici 2D su canvas: distribuzione radiale, diagramma dei livelli energetici, sezione piana di ψ.

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

function setup(canvas) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

function niceStep(range, target) {
  const raw = range / target;
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const f = raw / p;
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
}

function fmt(v) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-2) {
    const e = Math.floor(Math.log10(a));
    const m = v / Math.pow(10, e);
    return `${+m.toFixed(1)}·10${superscript(e)}`;
  }
  return String(+v.toPrecision(4)).replace('.', ',');
}

export function superscript(n) {
  const map = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  return String(n).split('').map(c => map[c] ?? c).join('');
}

const FONT = '"IBM Plex Sans", system-ui, sans-serif';
const MONO = '"IBM Plex Mono", ui-monospace, monospace';

/**
 * Grafico a linee con asse x lineare o logaritmico.
 * series: [{ xs, ys, color, fill?, width?, dash? , label }]
 */
export function drawLineChart(canvas, { series, log = false, xlabel, ylabel, markers = [], signed = false, legend = true }) {
  const { ctx, w, h } = setup(canvas);
  const text = cssVar('--muted');
  const grid = cssVar('--chart-grid');
  const ink = cssVar('--text');
  const pad = { l: 46, r: 12, t: legend ? 26 : 12, b: 34 };
  const pw = w - pad.l - pad.r;
  const ph = h - pad.t - pad.b;
  if (pw <= 10 || ph <= 10 || !series.length) return;

  let xmin = Infinity, xmax = -Infinity, ymax = 0, ymin = 0;
  for (const s of series) {
    for (let i = 0; i < s.xs.length; i++) {
      if (s.ys[i] === 0 && !s.keepZero) continue;
      xmin = Math.min(xmin, s.xs[i]);
      xmax = Math.max(xmax, s.xs[i]);
      ymax = Math.max(ymax, s.ys[i]);
      ymin = Math.min(ymin, s.ys[i]);
    }
  }
  if (!signed) ymin = 0;
  if (log) xmin = Math.max(xmin, xmax * 1e-4);
  else xmin = 0;
  ymax *= 1.08;
  if (signed) { const a = Math.max(ymax, -ymin * 1.08); ymax = a; ymin = -a; }
  const X = log
    ? (x) => pad.l + (Math.log10(Math.max(x, xmin)) - Math.log10(xmin)) / (Math.log10(xmax) - Math.log10(xmin)) * pw
    : (x) => pad.l + (x - xmin) / (xmax - xmin) * pw;
  const Y = (y) => pad.t + (1 - (y - ymin) / (ymax - ymin)) * ph;

  // griglia e tacche
  ctx.font = `11px ${MONO}`;
  ctx.fillStyle = text;
  ctx.strokeStyle = grid;
  ctx.lineWidth = 1;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  if (log) {
    for (let e = Math.ceil(Math.log10(xmin)); e <= Math.log10(xmax); e++) {
      const x = X(Math.pow(10, e));
      ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
      ctx.fillText(fmt(Math.pow(10, e)), x, pad.t + ph + 5);
      for (let k = 2; k < 10; k++) {
        const xv = k * Math.pow(10, e);
        if (xv > xmax || xv < xmin) continue;
        const xx = X(xv);
        ctx.globalAlpha = 0.4;
        ctx.beginPath(); ctx.moveTo(xx, pad.t + ph); ctx.lineTo(xx, pad.t + ph - 4); ctx.stroke();
        ctx.globalAlpha = 1;
      }
    }
  } else {
    const step = niceStep(xmax - xmin, Math.max(3, Math.floor(pw / 70)));
    for (let v = 0; v <= xmax + 1e-9; v += step) {
      const x = X(v);
      ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
      ctx.fillText(fmt(v), x, pad.t + ph + 5);
    }
  }
  if (signed) {
    ctx.strokeStyle = text;
    ctx.globalAlpha = 0.5;
    ctx.beginPath(); ctx.moveTo(pad.l, Y(0)); ctx.lineTo(pad.l + pw, Y(0)); ctx.stroke();
    ctx.globalAlpha = 1;
  }
  ctx.strokeStyle = text;
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + ph); ctx.lineTo(pad.l + pw, pad.t + ph);
  ctx.stroke();

  // etichette degli assi
  ctx.fillStyle = text;
  ctx.font = `11px ${FONT}`;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'bottom';
  ctx.fillText(xlabel, pad.l + pw, h - 2);
  ctx.save();
  ctx.translate(12, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText(ylabel, 0, 0);
  ctx.restore();

  ctx.save();
  ctx.beginPath();
  ctx.rect(pad.l, pad.t - 2, pw, ph + 2);
  ctx.clip();
  for (const s of series) {
    ctx.beginPath();
    let started = false;
    for (let i = 0; i < s.xs.length; i++) {
      const x = s.xs[i];
      if (x < xmin || x > xmax) continue;
      const px = X(x);
      const py = Y(s.ys[i]);
      if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
    }
    if (s.fill) {
      ctx.save();
      ctx.lineTo(X(xmax), Y(0));
      ctx.lineTo(X(xmin), Y(0));
      ctx.closePath();
      ctx.fillStyle = s.fill;
      ctx.fill();
      ctx.restore();
      // ridisegna il contorno
      ctx.beginPath();
      started = false;
      for (let i = 0; i < s.xs.length; i++) {
        const x = s.xs[i];
        if (x < xmin || x > xmax) continue;
        if (!started) { ctx.moveTo(X(x), Y(s.ys[i])); started = true; } else ctx.lineTo(X(x), Y(s.ys[i]));
      }
    }
    ctx.strokeStyle = s.color;
    ctx.lineWidth = s.width ?? 1.6;
    ctx.setLineDash(s.dash ?? []);
    ctx.stroke();
    ctx.setLineDash([]);
  }
  // marcatori verticali (nodi, r di massima probabilità)
  ctx.font = `10.5px ${MONO}`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  for (const m of markers) {
    if (m.x < xmin || m.x > xmax) continue;
    const x = X(m.x);
    ctx.strokeStyle = m.color ?? text;
    ctx.setLineDash(m.dash ?? [3, 3]);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
    ctx.setLineDash([]);
    if (m.label) {
      ctx.fillStyle = m.color ?? text;
      ctx.fillText(m.label, x + 3, pad.t + 2 + (m.row ?? 0) * 13);
    }
  }
  ctx.restore();

  if (legend) {
    ctx.font = `11px ${FONT}`;
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';
    let lx = pad.l;
    for (const s of series) {
      if (!s.label) continue;
      const tw = ctx.measureText(s.label).width;
      if (lx + tw + 22 > w) break;
      ctx.strokeStyle = s.color;
      ctx.lineWidth = 2;
      ctx.setLineDash(s.dash ?? []);
      ctx.beginPath(); ctx.moveTo(lx, 12); ctx.lineTo(lx + 12, 12); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = ink;
      ctx.fillText(s.label, lx + 16, 12);
      lx += tw + 28;
    }
  }
}

/**
 * Diagramma dei livelli: energie di legame |ε| in scala logaritmica (eV).
 * levels: [{ label, e (eV, negativo), occ, cap, color, highlight }]
 */
export function drawLevels(canvas, levels, { title } = {}) {
  const { ctx, w, h } = setup(canvas);
  const text = cssVar('--muted');
  const ink = cssVar('--text');
  const grid = cssVar('--chart-grid');
  const accent = cssVar('--accent');
  const pad = { l: 52, r: 12, t: 14, b: 26 };
  const ph = h - pad.t - pad.b;
  const vals = levels.map(l => -l.e);
  const vmin = Math.max(Math.min(...vals) * 0.6, 0.05);
  const vmax = Math.max(...vals) * 1.6;
  const lmin = Math.log10(vmin);
  const lmax = Math.log10(vmax);
  // energia più negativa in basso
  const Y = (v) => pad.t + (Math.log10(v) - lmin) / (lmax - lmin) * ph;

  ctx.font = `11px ${MONO}`;
  ctx.fillStyle = text;
  ctx.strokeStyle = grid;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (let e = Math.ceil(lmin); e <= lmax; e++) {
    const y = Y(Math.pow(10, e));
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
    ctx.fillText(`−${fmt(Math.pow(10, e))}`, pad.l - 6, y);
  }
  ctx.font = `11px ${FONT}`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'bottom';
  ctx.fillText(title ?? 'ε (eV), scala logaritmica', pad.l, h - 4);

  // livelli con etichette separate per evitare sovrapposizioni
  const x0 = pad.l + 14;
  const barW = Math.min(90, (w - pad.l - pad.r) * 0.32);
  const labelX = x0 + barW + 16;
  const items = levels.map(l => ({ ...l, y: Y(-l.e) })).sort((a, b) => a.y - b.y);
  const minGap = 14;
  const ly = items.map(it => it.y);
  for (let i = 1; i < ly.length; i++) if (ly[i] - ly[i - 1] < minGap) ly[i] = ly[i - 1] + minGap;
  const overflow = ly.length ? ly[ly.length - 1] - (h - pad.b) : 0;
  if (overflow > 0) {
    for (let i = ly.length - 1; i >= 0; i--) {
      const limit = i === ly.length - 1 ? h - pad.b : ly[i + 1] - minGap;
      ly[i] = Math.min(ly[i], limit);
    }
  }
  items.forEach((it, i) => {
    ctx.strokeStyle = it.highlight ? accent : (it.color ?? ink);
    ctx.lineWidth = it.highlight ? 3 : 2;
    ctx.globalAlpha = it.occ > 0 || it.highlight ? 1 : 0.45;
    ctx.beginPath(); ctx.moveTo(x0, it.y); ctx.lineTo(x0 + barW, it.y); ctx.stroke();
    // riempimento proporzionale all'occupazione
    if (it.cap) {
      ctx.globalAlpha = 0.35;
      ctx.fillStyle = it.color ?? ink;
      ctx.fillRect(x0, it.y - 3, barW * (it.occ / it.cap), 6);
    }
    ctx.globalAlpha = 1;
    ctx.strokeStyle = grid;
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0 + barW + 2, it.y); ctx.lineTo(labelX - 4, ly[i]); ctx.stroke();
    ctx.fillStyle = it.highlight ? accent : ink;
    ctx.font = `${it.highlight ? 600 : 500} 11.5px ${FONT}`;
    ctx.textBaseline = 'middle';
    const occTxt = it.cap ? `${it.label}${superscript(it.occ)}` : it.label;
    ctx.fillText(occTxt, labelX, ly[i]);
    ctx.fillStyle = text;
    ctx.font = `11px ${MONO}`;
    ctx.fillText(fmt(it.e), labelX + 46, ly[i]);
  });
}

/**
 * Diagramma degli orbitali molecolari di valenza (livelli atomici calcolati ai lati, OM al centro).
 */
export function drawMODiagram(canvas, { eS, eP, levels, selected, symbol }) {
  const { ctx, w, h } = setup(canvas);
  const text = cssVar('--muted');
  const ink = cssVar('--text');
  const accent = cssVar('--accent');
  const grid = cssVar('--chart-grid');
  const pad = { t: 16, b: 30 };
  const gap = Math.max(eP - eS, 2);
  // spostamenti qualitativi (in unità del divario s–p)
  const shift = { 'sigma-s': -0.35, 'sigma-s*': 0.25, 'pi-p': -0.34, 'sigma-p': -0.5, 'pi-p*': 0.3, 'sigma-p*': 0.62 };
  const base = { 'sigma-s': eS, 'sigma-s*': eS, 'pi-p': eP, 'sigma-p': eP, 'pi-p*': eP, 'sigma-p*': eP };
  const sMix = levels.findIndex(l => l.id === 'pi-p') < levels.findIndex(l => l.id === 'sigma-p');
  if (sMix) shift['sigma-p'] = -0.06;
  const energies = levels.map(l => base[l.id] + shift[l.id] * gap);
  const emin = Math.min(eS, ...energies) - 0.15 * gap;
  const emax = Math.max(eP, ...energies) + 0.15 * gap;
  const Y = (e) => pad.t + (1 - (e - emin) / (emax - emin)) * (h - pad.t - pad.b);
  const cx = w / 2;
  const aoW = Math.min(46, w * 0.12);
  const moW = Math.min(64, w * 0.18);
  const leftX = 18;
  const rightX = w - 18 - aoW;

  const drawElectrons = (x, y, width, n, cap) => {
    const orbitals = cap / 2;
    const slot = width / orbitals;
    for (let o = 0; o < orbitals; o++) {
      const up = n > o ? 1 : 0;
      const down = n > orbitals + o ? 1 : 0;
      const ox = x + slot * (o + 0.5);
      ctx.strokeStyle = ink;
      ctx.lineWidth = 1.4;
      if (up) arrow(ctx, ox - 3, y + 7, y - 8);
      if (down) arrow(ctx, ox + 3, y - 8, y + 7);
    }
  };
  const level = (x, y, width, color, lw = 2) => {
    ctx.strokeStyle = color;
    ctx.lineWidth = lw;
    const orbitals = Math.round(width / 30) || 1;
    ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + width, y); ctx.stroke();
    return orbitals;
  };

  ctx.font = `11px ${FONT}`;
  ctx.fillStyle = text;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText(`${symbol}`, leftX + aoW / 2, h - 4);
  ctx.fillText(`${symbol}₂`, cx, h - 4);
  ctx.fillText(`${symbol}`, rightX + aoW / 2, h - 4);

  // livelli atomici (calcolati)
  for (const [e, lab] of [[eS, 's'], [eP, 'p']]) {
    const y = Y(e);
    level(leftX, y, aoW, text);
    level(rightX, y, aoW, text);
    ctx.fillStyle = text;
    ctx.font = `10.5px ${MONO}`;
    ctx.textBaseline = 'bottom';
    ctx.textAlign = 'left';
    ctx.fillText(`${lab} ${fmt(e)} eV`, leftX, y - 3);
  }
  // linee di correlazione
  ctx.setLineDash([2, 3]);
  ctx.strokeStyle = grid;
  ctx.lineWidth = 1;
  levels.forEach((lv, i) => {
    const y = Y(energies[i]);
    const ae = lv.id.includes('-s') ? eS : eP;
    ctx.beginPath();
    ctx.moveTo(leftX + aoW, Y(ae)); ctx.lineTo(cx - moW / 2, y);
    ctx.moveTo(rightX, Y(ae)); ctx.lineTo(cx + moW / 2, y);
    ctx.stroke();
  });
  ctx.setLineDash([]);
  levels.forEach((lv, i) => {
    const y = Y(energies[i]);
    const isSel = lv.id === selected;
    level(cx - moW / 2, y, moW, isSel ? accent : ink, isSel ? 3 : 2);
    drawElectrons(cx - moW / 2, y, moW, lv.electrons, lv.cap);
    ctx.fillStyle = isSel ? accent : ink;
    ctx.font = `${isSel ? 600 : 500} 12px ${FONT}`;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(lv.name, cx + moW / 2 + 6, y);
  });
}

function arrow(ctx, x, y0, y1) {
  const dir = Math.sign(y1 - y0);
  ctx.beginPath();
  ctx.moveTo(x, y0); ctx.lineTo(x, y1);
  ctx.moveTo(x - 2.5, y1 - dir * 4); ctx.lineTo(x, y1); ctx.lineTo(x + 2.5, y1 - dir * 4);
  ctx.stroke();
}

function hexToRgb(hex) {
  const m = hex.replace('#', '');
  const v = parseInt(m.length === 3 ? m.split('').map(c => c + c).join('') : m, 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

/**
 * Sezione piana: mappa di ψ (con segno) o della densità (scala logaritmica).
 * field: Float32Array n×n, riga per riga dall'alto; extent in pm (semi-lato).
 */
export function drawSlice(canvas, { field, n, signed, extentPm, axes, nuclei = [] }) {
  const { ctx, w, h } = setup(canvas);
  const size = Math.min(w, h - 22);
  const ox = (w - size) / 2;
  const oy = 0;
  const bg = hexToRgb(cssVar('--slice-bg'));
  const pos = hexToRgb(cssVar('--phase-pos'));
  const neg = hexToRgb(cssVar('--phase-neg'));
  const dens = hexToRgb(cssVar('--accent'));
  const nodeCol = hexToRgb(cssVar('--slice-node'));
  const img = ctx.createImageData(n, n);
  let max = 0;
  for (let i = 0; i < field.length; i++) max = Math.max(max, Math.abs(field[i]));
  const logMin = Math.log(max * 1e-5);
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      const v = field[j * n + i];
      let t;
      let col;
      if (signed) {
        t = Math.pow(Math.abs(v) / max, 0.5);
        col = v >= 0 ? pos : neg;
      } else {
        t = v > 0 ? Math.max(0, (Math.log(v) - logMin) / (Math.log(max) - logMin)) : 0;
        col = dens;
      }
      let rgb = [0, 1, 2].map(k => bg[k] + (col[k] - bg[k]) * t);
      if (signed) {
        // linea nodale: cambio di segno con i vicini
        const r = i < n - 1 ? field[j * n + i + 1] : v;
        const d = j < n - 1 ? field[(j + 1) * n + i] : v;
        if ((v > 0 && (r < 0 || d < 0)) || (v < 0 && (r > 0 || d > 0))) rgb = nodeCol;
      }
      const o = 4 * (j * n + i);
      img.data[o] = rgb[0];
      img.data[o + 1] = rgb[1];
      img.data[o + 2] = rgb[2];
      img.data[o + 3] = 255;
    }
  }
  const tmp = document.createElement('canvas');
  tmp.width = tmp.height = n;
  tmp.getContext('2d').putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(tmp, ox, oy, size, size);

  const text = cssVar('--muted');
  ctx.strokeStyle = cssVar('--chart-grid');
  ctx.strokeRect(ox + 0.5, oy + 0.5, size - 1, size - 1);
  ctx.fillStyle = cssVar('--nucleus');
  for (const [u, v] of nuclei) {
    const px = ox + (u / extentPm * 0.5 + 0.5) * size;
    const py = oy + (0.5 - v / extentPm * 0.5) * size;
    ctx.beginPath(); ctx.arc(px, py, 2.5, 0, 2 * Math.PI); ctx.fill();
  }
  ctx.font = `11px ${MONO}`;
  ctx.fillStyle = text;
  ctx.textBaseline = 'top';
  ctx.textAlign = 'left';
  ctx.fillText(`piano ${axes[0]}${axes[1]}`, ox, size + 5);
  ctx.textAlign = 'right';
  ctx.fillText(`lato ${fmt(2 * extentPm)} pm`, ox + size, size + 5);
  ctx.textAlign = 'center';
  ctx.fillText(`${axes[0]} →`, ox + size / 2, size + 5);
}
