// Tavola periodica interattiva, colorata per blocco (s, p, d, f): il tipo di orbitale in riempimento.

import { ELEMENTS } from '../physics/elements.js';

export function blockOf(el) {
  if ((el.Z >= 57 && el.Z <= 71) || (el.Z >= 89 && el.Z <= 103)) return 'f';
  if (el.Z === 2) return 's';
  if (el.group === 1 || el.group === 2) return 's';
  if (el.group >= 13) return 'p';
  return 'd';
}

function position(el) {
  if (el.Z >= 57 && el.Z <= 71) return { row: 9, col: el.Z - 57 + 3 };
  if (el.Z >= 89 && el.Z <= 103) return { row: 10, col: el.Z - 89 + 3 };
  return { row: el.period, col: el.group };
}

export function buildPeriodicTable(container, onSelect) {
  container.innerHTML = '';
  const cells = new Map();
  for (const el of ELEMENTS) {
    const { row, col } = position(el);
    const b = document.createElement('button');
    b.type = 'button';
    b.className = `pt-cell block-${blockOf(el)}`;
    b.style.gridRow = String(row);
    b.style.gridColumn = String(col);
    b.dataset.z = String(el.Z);
    b.title = `${el.name} (Z = ${el.Z})`;
    b.setAttribute('aria-label', `${el.name}, numero atomico ${el.Z}`);
    b.innerHTML = `<span class="pt-z">${el.Z}</span><span class="pt-sym">${el.symbol}</span>`;
    b.addEventListener('click', () => onSelect(el.Z));
    container.appendChild(b);
    cells.set(el.Z, b);
  }
  // segnaposto per lantanidi e attinidi nel gruppo 3
  for (const [row, label] of [[6, '57–71'], [7, '89–103']]) {
    const s = document.createElement('div');
    s.className = 'pt-cell pt-placeholder block-f';
    s.style.gridRow = String(row);
    s.style.gridColumn = '3';
    s.innerHTML = `<span class="pt-sym">${label}</span>`;
    s.setAttribute('aria-hidden', 'true');
    container.appendChild(s);
  }
  const gap = document.createElement('div');
  gap.className = 'pt-gap';
  gap.style.gridRow = '8';
  gap.style.gridColumn = '1 / -1';
  container.appendChild(gap);

  // navigazione con le frecce della tastiera
  container.addEventListener('keydown', (ev) => {
    const cur = document.activeElement?.dataset?.z;
    if (!cur) return;
    const delta = { ArrowRight: 1, ArrowLeft: -1 }[ev.key];
    if (delta) {
      const next = cells.get(+cur + delta);
      if (next) { next.focus(); ev.preventDefault(); }
    }
  });

  return {
    select(Z) {
      for (const [z, c] of cells) c.classList.toggle('selected', z === Z);
    },
    /** Colora le caselle secondo una proprietà (mappa Z → valore); null = colori per blocco. */
    heatmap(values) {
      container.classList.toggle('heatmap', !!values);
      if (!values) {
        for (const c of cells.values()) { c.style.removeProperty('--heat'); c.style.removeProperty('--heat-ink'); c.title = c.title.split(' · ')[0]; }
        return null;
      }
      const nums = [...values.values()].filter(v => v !== null && Number.isFinite(v));
      const min = Math.min(...nums), max = Math.max(...nums);
      for (const [z, c] of cells) {
        const v = values.get(z);
        c.title = c.title.split(' · ')[0] + (v !== null && v !== undefined ? ` · ${v.toLocaleString('it-IT', { maximumFractionDigits: 2 })}` : '');
        if (v === null || v === undefined || !Number.isFinite(v)) { c.style.removeProperty('--heat'); c.style.removeProperty('--heat-ink'); continue; }
        const t = (v - min) / (max - min || 1);
        const [r, g, b] = ramp(t);
        c.style.setProperty('--heat', `rgb(${r}, ${g}, ${b})`);
        c.style.setProperty('--heat-ink', 0.299 * r + 0.587 * g + 0.114 * b > 140 ? '#10141b' : '#f4f6fa');
      }
      return { min, max };
    },
  };
}

// Scala di colori sequenziale percettivamente uniforme (approssimazione di "viridis").
const STOPS = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
export function ramp(t) {
  const x = Math.max(0, Math.min(1, t)) * (STOPS.length - 1);
  const i = Math.min(Math.floor(x), STOPS.length - 2);
  const f = x - i;
  return STOPS[i].map((c, k) => Math.round(c + (STOPS[i + 1][k] - c) * f));
}
