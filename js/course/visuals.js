// Animated teaching diagrams. Their positions are markers, never electron paths.
export const COLORS = [
  '#5ab8e7',
  '#f0b429',
  '#b4a0ed',
  '#72c9a5',
  '#f58b76',
  '#8ecad1',
];
const sub = (s) =>
  String(s)
    .replace(/\^(\d*)([+-])/g, (_, n, sign) =>
      [...((n || '') + sign)]
        .map((c) =>
          c === '+' ? '⁺' : c === '-' ? '⁻' : '⁰¹²³⁴⁵⁶⁷⁸⁹'[Number(c)],
        )
        .join(''),
    )
    .replace(/\d/g, (c) => '₀₁₂₃₄₅₆₇₈₉'[Number(c)]);
export function canvasAnimation(
  canvas,
  render,
  { getPaused = () => false } = {},
) {
  const ctx = canvas.getContext('2d');
  if (!ctx) return () => {};
  let frame = 0,
    previous = 0,
    time = 0,
    disposed = false;
  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)');
  const resize = () => {
    const w = Math.max(240, canvas.clientWidth || 960),
      h = w * 0.38,
      dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  };
  const ro = new ResizeObserver(resize);
  ro.observe(canvas);
  resize();
  const loop = (now) => {
    if (disposed) return;
    frame = requestAnimationFrame(loop);
    const dt = previous ? Math.min((now - previous) / 1000, 0.05) : 0;
    previous = now;
    if (document.hidden || !canvas.isConnected) return;
    if (!getPaused() && !reduced.matches) time += dt;
    ctx.setTransform(canvas.width / 960, 0, 0, canvas.height / 365, 0, 0);
    const style = getComputedStyle(document.documentElement);
    const palette = {
      bg: style.getPropertyValue('--panel').trim() || '#11161e',
      text: style.getPropertyValue('--text').trim() || '#edf1f5',
      muted: style.getPropertyValue('--muted').trim() || '#93a0b4',
      line: style.getPropertyValue('--line').trim() || '#303a49',
    };
    ctx.fillStyle = palette.bg;
    ctx.fillRect(0, 0, 960, 365);
    render(ctx, time, palette);
  };
  frame = requestAnimationFrame(loop);
  return () => {
    disposed = true;
    cancelAnimationFrame(frame);
    ro.disconnect();
  };
}
function text(ctx, value, x, y, color, size = 16, align = 'left') {
  ctx.fillStyle = color;
  ctx.font = `${size}px "IBM Plex Sans", sans-serif`;
  ctx.textAlign = align;
  ctx.fillText(String(value), x, y);
}
function dot(ctx, x, y, r, color) {
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2 * Math.PI);
  ctx.fillStyle = color;
  ctx.fill();
}
export function reactionDiagram(ctx, time, p, model) {
  if (!model) {
    text(
      ctx,
      'Scegli una reazione e quantità valide',
      480,
      180,
      p.muted,
      20,
      'center',
    );
    return;
  }
  const { rx, left, right, result, fraction } = model;
  const x0 = 390,
    y0 = 85,
    bw = 180,
    bh = 210;
  ctx.strokeStyle = p.line;
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.moveTo(x0, y0);
  ctx.lineTo(x0, y0 + bh);
  ctx.quadraticCurveTo(x0, y0 + bh + 15, x0 + 15, y0 + bh + 15);
  ctx.lineTo(x0 + bw - 15, y0 + bh + 15);
  ctx.quadraticCurveTo(x0 + bw, y0 + bh + 15, x0 + bw, y0 + bh);
  ctx.lineTo(x0 + bw, y0);
  ctx.stroke();
  ctx.fillStyle = 'rgba(90,184,231,.10)';
  ctx.fillRect(x0 + 3, y0 + 60, bw - 6, bh - 63);
  const active = result.extent > 0;
  if (active && ['precip', 'deposit', 'solid'].includes(rx.visual)) {
    const h = 5 + fraction * 38;
    ctx.fillStyle =
      rx.visual === 'deposit'
        ? '#cf9165'
        : rx.id === 'idrossido'
          ? '#5ab8e7'
          : '#b9bec7';
    ctx.fillRect(x0 + 5, y0 + bh - h, bw - 10, h);
    for (let i = 0; i < 25 * fraction; i++)
      dot(
        ctx,
        x0 + 15 + ((i * 43) % 150),
        y0 + 90 + ((time * 34 + i * 17) % Math.max(20, bh - h - 90)),
        3,
        rx.id === 'idrossido' ? '#5ab8e7' : '#b9bec7',
      );
  }
  if (active && rx.visual === 'gas')
    for (let i = 0; i < Math.round(15 * fraction); i++) {
      const u = (time * 0.3 + i * 0.137) % 1,
        x = x0 + 20 + ((i * 47) % 140),
        y = y0 + bh - u * (bh + 50);
      ctx.beginPath();
      ctx.arc(x, y, 2 + u * 5, 0, 2 * Math.PI);
      ctx.strokeStyle = '#8ecad1';
      ctx.lineWidth = 1.4;
      ctx.stroke();
    }
  text(
    ctx,
    rx.visual === 'precip'
      ? 'Precipitato'
      : rx.visual === 'gas'
        ? 'Sviluppo di gas'
        : rx.visual === 'deposit'
          ? 'Deposito metallico'
          : 'Avanzamento',
    480,
    335,
    p.muted,
    14,
    'center',
  );
  const drawSide = (species, amounts, x, align) =>
    species.slice(0, 4).forEach((s, i) => {
      const y = 80 + i * 65,
        color = COLORS[i % COLORS.length];
      text(ctx, sub(s.formula), x, y, color, 22, align);
      text(
        ctx,
        `${amounts[i].toLocaleString('it-IT', { maximumSignificantDigits: 4 })} mol`,
        x,
        y + 24,
        p.text,
        14,
        align,
      );
      const dots = Math.min(
        14,
        Math.round((14 * amounts[i]) / Math.max(...amounts, 1e-12)),
      );
      for (let j = 0; j < dots; j++)
        dot(
          ctx,
          x + (align === 'right' ? -1 : 1) * (j % 7) * 14,
          y + 38 + Math.floor(j / 7) * 10,
          3,
          color,
        );
    });
  text(ctx, 'REAGENTI RESIDUI', 42, 32, p.muted, 12);
  text(ctx, 'PRODOTTI FORMATI', 650, 32, p.muted, 12);
  drawSide(left, result.remaining, 42, 'left');
  drawSide(
    right,
    result.products.map((v) => v.mol),
    650,
    'left',
  );
  text(
    ctx,
    `ξ = ${result.extent.toLocaleString('it-IT', { maximumSignificantDigits: 4 })} mol`,
    480,
    46,
    p.text,
    18,
    'center',
  );
}
export function dilutionDiagram(ctx, time, p, model) {
  if (!model) return;
  const drawBeaker = (x, volume, C, label) => {
    const fill = Math.max(
      20,
      Math.min(200, (200 * volume) / Math.max(model.Vi, model.Vf)),
    );
    ctx.strokeStyle = p.line;
    ctx.lineWidth = 3;
    ctx.strokeRect(x, 65, 240, 235);
    ctx.fillStyle = 'rgba(90,184,231,.12)';
    ctx.fillRect(x + 2, 300 - fill, 236, fill - 2);
    // Equal marker count in both vessels: dilution conserves solute, not its density.
    for (let i = 0; i < 36; i++) {
      const u = (i * 0.6180339 + time * 0.012) % 1,
        v = (i * 0.4142135 + time * 0.008) % 1;
      dot(
        ctx,
        x + 14 + u * 210,
        300 - fill + 12 + v * Math.max(1, fill - 24),
        4,
        COLORS[0],
      );
    }
    text(ctx, label, x + 120, 40, p.text, 18, 'center');
    text(
      ctx,
      `${volume} mL · ${C.toLocaleString('it-IT', { maximumSignificantDigits: 4 })} mol/L`,
      x + 120,
      330,
      p.muted,
      15,
      'center',
    );
  };
  drawBeaker(100, model.Vi, model.C, 'Prima');
  drawBeaker(610, model.Vf, model.out.finalConcentration, 'Dopo');
  text(ctx, '+ solvente', 480, 152, COLORS[1], 18, 'center');
  text(ctx, '→', 480, 198, p.text, 42, 'center');
  text(ctx, 'stesse moli', 480, 237, p.muted, 15, 'center');
}
export function scatterDiagram(ctx, time, p, model) {
  if (!model) return;
  const pts = model.path,
    xs = pts.map((v) => v.x),
    ys = pts.map((v) => v.y);
  const xmin = Math.min(0, ...xs),
    xmax = Math.max(0, ...xs),
    ymin = Math.min(0, ...ys),
    ymax = Math.max(0, ...ys);
  const scale = Math.min(
    800 / Math.max(xmax - xmin, 1),
    260 / Math.max(ymax - ymin, 1),
  );
  const ox = 480 - ((xmin + xmax) / 2) * scale,
    oy = 175 + ((ymin + ymax) / 2) * scale;
  const pos = (v) => [ox + v.x * scale, oy - v.y * scale];
  ctx.strokeStyle = COLORS[0];
  ctx.lineWidth = 2;
  ctx.beginPath();
  pts.forEach((v, i) => {
    const [x, y] = pos(v);
    i ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
  });
  ctx.stroke();
  dot(ctx, ox, oy, 7, COLORS[4]);
  text(ctx, 'nucleo', ox + 13, oy + 5, p.muted, 13);
  const u = (time / 6) % 1;
  let i = pts.findIndex((v) => v.t >= u);
  if (i < 1) i = 1;
  const a = pts[i - 1],
    b = pts[i],
    f = (u - a.t) / (b.t - a.t),
    [x, y] = pos({ x: a.x + (b.x - a.x) * f, y: a.y + (b.y - a.y) * f });
  dot(ctx, x, y, 5, COLORS[1]);
  text(
    ctx,
    `θ = ${((model.theta * 180) / Math.PI).toLocaleString('it-IT', { maximumFractionDigits: 2 })}°`,
    35,
    32,
    p.text,
    18,
  );
  text(
    ctx,
    'Particella α · potenziale coulombiano repulsivo · scala uguale sui due assi',
    35,
    345,
    p.muted,
    14,
  );
}
export function isotopeDiagram(ctx, time, p, percent) {
  const n = Math.round(percent);
  text(ctx, 'Miscela isotopica · 100 marcatori', 35, 32, p.muted, 15);
  for (let i = 0; i < 100; i++) {
    const x = 45 + (i % 20) * 29,
      y = 85 + Math.floor(i / 20) * 41;
    dot(ctx, x, y, 8, i < n ? COLORS[0] : COLORS[2]);
  }
  text(
    ctx,
    `³⁵Cl  ${percent.toLocaleString('it-IT')}%`,
    665,
    140,
    COLORS[0],
    23,
  );
  text(
    ctx,
    `³⁷Cl  ${(100 - percent).toLocaleString('it-IT', { maximumFractionDigits: 2 })}%`,
    665,
    182,
    COLORS[2],
    23,
  );
  text(
    ctx,
    'La media riguarda la popolazione, non un singolo atomo “frazionario”.',
    35,
    325,
    p.muted,
    16,
  );
}
