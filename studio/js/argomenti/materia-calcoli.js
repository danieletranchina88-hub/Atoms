// Argomento 1, parte di calcolo: unità, conversioni con l'analisi dimensionale, temperature e generatori di esercizi.
// Nessun accesso al DOM: il modulo è verificato da tests/studio.mjs.

import { parseNumber, formatSig, formatPlace, formatPlain, roundSig, placeOfSig, sigOfPlace, exponentOf, pow10 } from '../numeri.js';
import { rand, randInt, pick } from '../ui.js';

// Fattori verso un'unità di riferimento del gruppo (esatti per definizione).
export const UNITA = {
  lunghezza: { base: 'm', nome: 'Lunghezza', u: { km: 1e3, m: 1, dm: 1e-1, cm: 1e-2, mm: 1e-3, 'µm': 1e-6, nm: 1e-9, 'Å': 1e-10, pm: 1e-12 } },
  massa: { base: 'g', nome: 'Massa', u: { kg: 1e3, g: 1, mg: 1e-3, 'µg': 1e-6, ng: 1e-9 } },
  volume: { base: 'L', nome: 'Volume', u: { 'm³': 1e3, 'dm³': 1, L: 1, dL: 1e-1, mL: 1e-3, 'cm³': 1e-3, 'µL': 1e-6, 'mm³': 1e-6 } },
  densita: { base: 'g/cm³', nome: 'Densità', u: { 'g/cm³': 1, 'g/mL': 1, 'kg/L': 1, 'kg/m³': 1e-3, 'g/L': 1e-3 } },
};

export const NOTE_UNITA = {
  'Å': '1 Å (ångström) = 10⁻¹⁰ m = 100 pm: non è SI ma si usa per le distanze fra atomi.',
  'dm³': '1 dm³ = 1 L esattamente.',
  'cm³': '1 cm³ = (10⁻² m)³ = 10⁻⁶ m³ = 1 mL.',
  'm³': '1 m³ = (10 dm)³ = 10³ dm³ = 10³ L.',
  'mm³': '1 mm³ = (10⁻³ m)³ = 10⁻⁹ m³ = 1 µL.',
  'kg/m³': '1 g/cm³ = 10⁻³ kg / 10⁻⁶ m³ = 10³ kg/m³.',
  'g/L': '1 g/cm³ = 1 g/mL = 10³ g/L.',
};

// Scrive un fattore esatto: 1000 → 10³, 0,001 → 10⁻³, 3600 → 3600.
export function factorText(f) {
  const e = Math.round(Math.log10(f));
  if (Math.abs(f / 10 ** e - 1) < 1e-12) return e === 0 ? '1' : pow10(e);
  return formatPlain(f);
}

const frac = (num, den) => `<span class="frac"><span>${num}</span><span>${den}</span></span>`;
const cancel = (u) => `<span class="cancel">${u}</span>`;

// Catena di fattori di conversione: valore × (… / …) × (… / …) = risultato, con le unità che si semplificano.
export function conversionChain(valueText, from, to, group) {
  const g = UNITA[group];
  const ff = g.u[from], ft = g.u[to];
  const parts = [`${valueText} ${from === to ? from : cancel(from)}`];
  if (from !== to && from !== g.base && to !== g.base) {
    parts.push(frac(`${factorText(ff)} ${cancel(g.base)}`, `1 ${cancel(from)}`));
    parts.push(frac(`1 ${to}`, `${factorText(ft)} ${cancel(g.base)}`));
  } else if (from !== to && to === g.base) {
    parts.push(frac(`${factorText(ff)} ${to}`, `1 ${cancel(from)}`));
  } else if (from !== to) {
    parts.push(frac(`1 ${to}`, `${factorText(ft)} ${cancel(from)}`));
  }
  return parts.join(' <span class="op">×</span> ');
}

export function convert(valueText, from, to, group) {
  const p = parseNumber(valueText);
  if (!p) return null;
  const g = UNITA[group];
  const exact = p.value * g.u[from] / g.u[to];
  const sig = p.sigMax ?? 1; // zeri finali ambigui considerati misurati (il widget lo segnala)
  return { exact, sig, ambiguous: p.ambiguous, sigMin: p.sig, text: formatSig(exact, sig), chain: conversionChain(valueText.trim(), from, to, group) };
}

// Temperature: T(K) = t(°C) + 273,15 ; t(°F) = 9/5 t(°C) + 32. 273,15, 9/5 e 32 sono esatti.
export const SCALE_T = ['°C', 'K', '°F'];
export function convertTemperature(valueText, from, to) {
  const p = parseNumber(valueText);
  if (!p) return null;
  const steps = [];
  // si porta avanti la posizione dell'ultima cifra (somme) o le c.s. (prodotti);
  // gli zeri finali ambigui di un intero (es. 300 K) si considerano misurati
  let v = p.value, place = p.ambiguous ? p.exp : p.lastPlace, sig = null;
  const v0 = valueText.trim();
  if (from === to) return { exact: v, place, text: formatPlace(v, place), steps: [`${v0} ${from}`] };
  // verso °C
  if (from === 'K') {
    v = p.value - 273.15;
    steps.push(`t = T − 273,15 = ${v0} − 273,15 = ${formatPlain(v)} °C <span class="note">(somma: si tiene la posizione dell'ultima cifra di ${v0})</span>`);
  } else if (from === '°F') {
    const d = p.value - 32;
    const sd = sigOfPlace(d, place);
    v = d * 5 / 9; sig = sd;
    steps.push(`t(°F) − 32 = ${v0} − 32 = ${formatPlace(d, place)} <span class="note">(${sd} c.s.)</span>`);
    steps.push(`t = 5/9 × ${formatPlace(d, place)} = ${formatPlain(v)} °C <span class="note">(prodotto per un numero esatto: ${sd} c.s.)</span>`);
    place = placeOfSig(v, sd);
  }
  if (to === '°C') return finish(v, place, steps, '°C');
  if (to === 'K') {
    const t = v; v = t + 273.15;
    steps.push(`T = t + 273,15 = ${formatPlain(t)} + 273,15 = ${formatPlain(v)} K`);
    return finish(v, place, steps, 'K');
  }
  // verso °F
  const tC = v;
  const s = sig ?? sigOfPlace(tC, place);
  const prod = 1.8 * tC;
  const pp = placeOfSig(prod, s);
  v = prod + 32;
  steps.push(`9/5 × t = 1,8 × ${formatPlain(tC)} = ${formatPlain(prod)} <span class="note">(${s} c.s., ultima cifra: ${formatPlace(prod, pp)})</span>`);
  steps.push(`t(°F) = ${formatPlace(prod, pp)} + 32 = ${formatPlain(v)} °F`);
  return finish(v, pp, steps, '°F');
}
function finish(v, place, steps, unit) {
  const text = formatPlace(v, place);
  steps.push(`Risultato arrotondato: <b>${text} ${unit}</b>`);
  return { exact: v, place, text, steps };
}

// ——— Generatori di esercizi ———

const randSig = (a, b, n) => roundSig(rand(a, b), n);

function contaGen() {
  const n = randInt(1, 5);
  const D = [randInt(1, 9)];
  for (let i = 1; i < n; i++) D.push(Math.random() < 0.4 ? 0 : randInt(0, 9));
  const form = pick(n === 1 ? ['piccolo', 'intero', 'sci'] : ['piccolo', 'decimale', 'intero', 'sci', 'finale']);
  let s;
  if (form === 'piccolo') s = `0,${'0'.repeat(randInt(0, 3))}${D.join('')}`;
  else if (form === 'decimale') { const k = randInt(1, n - 1); s = `${D.slice(0, k).join('')},${D.slice(k).join('')}`; }
  else if (form === 'finale') { D[n - 1] = 0; const k = randInt(1, n - 1); s = `${D.slice(0, k).join('')},${D.slice(k).join('')}`; }
  else if (form === 'intero') { if (n > 1) D[n - 1] = randInt(1, 9); s = D.join(''); }
  else { const e = pick([-9, -7, -5, -4, -3, 3, 4, 5, 6, 8, 23]); s = `${D[0]}${n > 1 ? ',' + D.slice(1).join('') : ''} × ${pow10(e)}`; }
  const p = parseNumber(s);
  return {
    prompt: `<p>Quante cifre significative ha il numero <span class="num big">${s}</span>?</p>`,
    expected: { exact: p.sig, integer: true }, unit: 'c.s.',
    steps: `${digitsExplanation(p, s)}<p>Totale: <b>${p.sig}</b> cifre significative.</p>`,
    value: s,
  };
}

export function digitsExplanation(p) {
  const why = {
    'zero-iniziale': 'zeri iniziali: indicano solo dove sta la virgola, non sono significativi',
    'non-nulla': 'cifre diverse da zero: sempre significative',
    'zero-intermedio': 'zeri compresi fra cifre non nulle: significativi',
    'zero-finale': 'zeri finali con la virgola scritta: significativi, sono stati misurati',
    'zero-ambiguo': 'zeri finali di un intero senza virgola: ambigui',
  };
  const used = [...new Set(p.marks)];
  return `<ul class="tight">${used.map(r => `<li><span class="mk mk-${r}">${p.marks.filter(m => m === r).length}</span> ${why[r]}</li>`).join('')}${p.sci ? '<li>la potenza di 10 non conta: dà solo l\'ordine di grandezza</li>' : ''}</ul>`;
}

function operazioniGen() {
  if (Math.random() < 0.5) {
    const n1 = randInt(2, 5), n2 = randInt(2, 4);
    const a = randSig(1.2, 99, n1), b = randSig(1.1, 60, n2);
    const div = Math.random() < 0.5;
    const exact = div ? a / b : a * b;
    const sig = Math.min(n1, n2);
    const A = formatSig(a, n1), B = formatSig(b, n2);
    return {
      prompt: `<p>Calcola e scrivi il risultato con il numero corretto di cifre significative:</p><p class="num big">${A} ${div ? '÷' : '×'} ${B}</p>`,
      expected: { exact, sig }, unit: '',
      steps: `<p>${A} ha ${n1} c.s., ${B} ne ha ${n2}. In una ${div ? 'divisione' : 'moltiplicazione'} il risultato ha le cifre significative del dato che ne ha <b>meno</b>: ${sig}.</p><p>${A} ${div ? '÷' : '×'} ${B} = ${formatPlain(exact, 8)} → <b>${formatSig(exact, sig)}</b></p>`,
    };
  }
  const k = randInt(2, 3);
  const terms = Array.from({ length: k }, () => { const d = randInt(0, 3); return { d, v: Number(rand(1, 400).toFixed(d)) }; });
  if (terms.some(t => t.v === 0)) return operazioniGen();
  const sub = k === 2 && Math.random() < 0.4;
  if (sub && terms[0].v < terms[1].v) terms.reverse();
  const exact = sub ? terms[0].v - terms[1].v : terms.reduce((s, t) => s + t.v, 0);
  const place = -Math.min(...terms.map(t => t.d));
  const txt = terms.map(t => formatPlace(t.v, -t.d));
  const expr = sub ? `${txt[0]} − ${txt[1]}` : txt.join(' + ');
  return {
    prompt: `<p>Calcola e scrivi il risultato con il numero corretto di cifre:</p><p class="num big">${expr}</p>`,
    expected: { exact, place }, unit: '',
    steps: `<p>In una ${sub ? 'sottrazione' : 'somma'} conta la <b>posizione</b> dell'ultima cifra, non il numero di c.s.: il risultato si ferma alla cifra ${place === 0 ? 'delle unità' : place === -1 ? 'dei decimi' : place === -2 ? 'dei centesimi' : 'dei millesimi'}, quella del termine meno preciso.</p><p>${expr} = ${formatPlain(exact, 9)} → <b>${formatPlace(exact, place)}</b></p>`,
  };
}

function conversioneGen() {
  const group = pick(['lunghezza', 'massa', 'volume', 'volume']);
  const units = Object.keys(UNITA[group].u);
  let from, to;
  do { from = pick(units); to = pick(units); } while (from === to || UNITA[group].u[from] === UNITA[group].u[to]);
  const n = randInt(2, 4);
  const v = randSig(1, 999, n);
  const vt = formatSig(v, n);
  const c = convert(vt, from, to, group);
  const notes = [from, to].map(u => NOTE_UNITA[u]).filter(Boolean);
  return {
    prompt: `<p>Converti <span class="num">${vt} ${from}</span> in <b>${to}</b>.</p>`,
    expected: { exact: c.exact, sig: c.sig }, unit: to,
    steps: `${notes.length ? `<p class="note">${notes.join(' ')}</p>` : ''}<p class="chain">${c.chain} = <b>${c.text} ${to}</b></p><p class="note">I fattori di conversione sono esatti: il risultato conserva le ${c.sig} c.s. del dato.</p>`,
  };
}

const COMPOSTE = [
  { da: 'km/h', a: 'm/s', f: 1 / 3.6, catena: (x) => `${x} ${frac(cancel('km'), cancel('h'))} × ${frac(`10³ m`, `1 ${cancel('km')}`)} × ${frac(`1 ${cancel('h')}`, '3600 s')}`, range: [20, 130] },
  { da: 'm/s', a: 'km/h', f: 3.6, catena: (x) => `${x} ${frac(cancel('m'), cancel('s'))} × ${frac('1 km', `10³ ${cancel('m')}`)} × ${frac(`3600 ${cancel('s')}`, '1 h')}`, range: [3, 340] },
  { da: 'g/cm³', a: 'kg/m³', f: 1e3, catena: (x) => `${x} ${frac(cancel('g'), cancel('cm³'))} × ${frac('1 kg', `10³ ${cancel('g')}`)} × ${frac(`10⁶ ${cancel('cm³')}`, '1 m³')}`, range: [0.6, 21] },
  { da: 'kg/m³', a: 'g/cm³', f: 1e-3, catena: (x) => `${x} ${frac(cancel('kg'), cancel('m³'))} × ${frac(`10³ g`, `1 ${cancel('kg')}`)} × ${frac(`1 ${cancel('m³')}`, '10⁶ cm³')}`, range: [600, 21000] },
  { da: 'mL/min', a: 'L/h', f: 0.06, catena: (x) => `${x} ${frac(cancel('mL'), cancel('min'))} × ${frac('1 L', `10³ ${cancel('mL')}`)} × ${frac(`60 ${cancel('min')}`, '1 h')}`, range: [1, 500] },
];
function composteGen() {
  const t = pick(COMPOSTE);
  const n = randInt(2, 4);
  const v = randSig(t.range[0], t.range[1], n);
  const vt = formatSig(v, n);
  const exact = v * t.f;
  return {
    prompt: `<p>Converti <span class="num">${vt} ${t.da}</span> in <b>${t.a}</b>.</p>`,
    expected: { exact, sig: n }, unit: t.a,
    steps: `<p class="chain">${t.catena(vt)} = <b>${formatSig(exact, n)} ${t.a}</b></p><p class="note">Ogni unità si converte con il suo fattore; per cm³ e m³ il fattore è il cubo di quello delle lunghezze: (10² cm)³ = 10⁶ cm³. 3600 s/h e 60 min/h sono esatti.</p>`,
  };
}

function temperaturaGen() {
  const pairs = [['°C', 'K'], ['K', '°C'], ['°C', '°F'], ['°F', '°C'], ['K', '°F']];
  const [from, to] = pick(pairs);
  const d = randInt(0, 2);
  const base = from === '°C' ? rand(-190, 400) : from === 'K' ? rand(60, 700) : rand(-300, 700);
  let v = Number(base.toFixed(d));
  if (Math.abs(v) < 10) v += 37;
  if (d === 0 && v % 10 === 0) v += 3; // niente interi con zeri finali ambigui
  const vt = formatPlace(v, -d);
  const r = convertTemperature(vt, from, to);
  return {
    prompt: `<p>Converti <span class="num">${vt} ${from}</span> in <b>${to}</b>.</p>`,
    expected: { exact: r.exact, place: r.place }, unit: to,
    steps: `<ol class="tight">${r.steps.map(s => `<li>${s}</li>`).join('')}</ol>`,
  };
}

export const METALLI = [
  { nome: 'magnesio', d: 1.74 }, { nome: 'alluminio', d: 2.70 }, { nome: 'titanio', d: 4.51 }, { nome: 'zinco', d: 7.14 },
  { nome: 'ferro', d: 7.87 }, { nome: 'rame', d: 8.96 }, { nome: 'argento', d: 10.5 }, { nome: 'piombo', d: 11.3 }, { nome: 'oro', d: 19.3 },
];
export const LIQUIDI = [
  { nome: 'etanolo', d: '0,789' }, { nome: 'benzene', d: '0,877' }, { nome: 'acqua (20 °C)', d: '0,998' },
  { nome: 'glicerolo', d: '1,26' }, { nome: 'cloroformio', d: '1,48' }, { nome: 'mercurio', d: '13,53' },
];

function densitaGen() {
  const kind = randInt(0, 2);
  if (kind === 0) {
    const met = pick(METALLI);
    const m = Number(rand(met.d * 2.5, met.d * 9).toFixed(2));
    const V1 = Number(rand(15, 30).toFixed(1));
    const V2 = Number((V1 + m / met.d).toFixed(1));
    const dV = V2 - V1;
    const sdV = sigOfPlace(dV, -1);
    const sm = parseNumber(formatPlace(m, -2)).sig;
    const sig = Math.min(sdV, sm);
    const exact = m / roundPlaceSafe(dV);
    return {
      prompt: `<p>Un oggetto metallico di massa <span class="num">${formatPlace(m, -2)} g</span> viene immerso in un cilindro graduato: il livello dell'acqua sale da <span class="num">${formatPlace(V1, -1)} mL</span> a <span class="num">${formatPlace(V2, -1)} mL</span>. Calcola la densità del metallo.</p>`,
      expected: { exact, sig }, unit: 'g/cm³',
      steps: `<p>Volume dell'oggetto (sottrazione, si tengono i decimi): V = ${formatPlace(V2, -1)} − ${formatPlace(V1, -1)} = <b>${formatPlace(dV, -1)} mL</b> = ${formatPlace(dV, -1)} cm³ (${sdV} c.s.).</p><p>d = m / V = ${formatPlace(m, -2)} g / ${formatPlace(dV, -1)} cm³ = ${formatPlain(exact)} → <b>${formatSig(exact, sig)} g/cm³</b> (${sig} c.s., limitate dal volume).</p><p class="note">Confronta con la tabella delle densità: il valore è compatibile con il ${met.nome} (${formatSig(met.d, 3)} g/cm³). La sottrazione ha ridotto le cifre significative del volume: è il passaggio che limita la precisione.</p>`,
    };
  }
  if (kind === 1) {
    const liq = pick(LIQUIDI);
    const nV = randInt(3, 4);
    const V = randSig(12, 480, nV);
    const Vt = formatSig(V, nV);
    const pd = parseNumber(liq.d);
    const exact = pd.value * V;
    const sig = Math.min(pd.sig, nV);
    return {
      prompt: `<p>Qual è la massa di <span class="num">${Vt} mL</span> di ${liq.nome} (d = <span class="num">${liq.d} g/mL</span>)?</p>`,
      expected: { exact, sig }, unit: 'g',
      steps: `<p>m = d × V = ${liq.d} g/mL × ${Vt} mL = ${formatPlain(exact)} g → <b>${formatSig(exact, sig)} g</b></p><p class="note">La densità fa da fattore di conversione fra volume e massa. Cifre significative: il minimo fra ${pd.sig} e ${nV}.</p>`,
    };
  }
  const met = pick(METALLI);
  const m = randSig(0.2, 9.5, 3);
  const dn = parseNumber(formatSig(met.d, 3)).sig;
  const exact = m * 1000 / met.d;
  return {
    prompt: `<p>Che volume, in cm³, occupa un blocco di ${met.nome} di massa <span class="num">${formatSig(m, 3)} kg</span>? (d = <span class="num">${formatSig(met.d, 3)} g/cm³</span>)</p>`,
    expected: { exact, sig: Math.min(3, dn) }, unit: 'cm³',
    steps: `<p class="chain">${formatSig(m, 3)} ${cancel('kg')} × ${frac(`10³ ${cancel('g')}`, `1 ${cancel('kg')}`)} × ${frac('1 cm³', `${formatSig(met.d, 3)} ${cancel('g')}`)} = <b>${formatSig(exact, 3)} cm³</b></p><p class="note">Si divide per la densità: il fattore va scritto con i grammi al denominatore, così si semplificano.</p>`,
  };
}
const roundPlaceSafe = (x) => Number(x.toFixed(1));

function percentualeGen() {
  const kind = randInt(0, 2);
  if (kind === 0) {
    const pc = randSig(55, 90, 3), m = randSig(5, 95, 3);
    const exact = m * pc / 100;
    return {
      prompt: `<p>Un ottone contiene il <span class="num">${formatSig(pc, 3)}%</span> in massa di rame. Quanti grammi di rame ci sono in <span class="num">${formatSig(m, 3)} g</span> di ottone?</p>`,
      expected: { exact, sig: 3 }, unit: 'g',
      steps: `<p>La percentuale è un fattore di conversione: ${formatSig(pc, 3)} g di rame ogni 100 g di ottone (100 è esatto).</p><p class="chain">${formatSig(m, 3)} ${cancel('g ottone')} × ${frac(`${formatSig(pc, 3)} g rame`, `100 ${cancel('g ottone')}`)} = <b>${formatSig(exact, 3)} g</b></p>`,
    };
  }
  if (kind === 1) {
    const pc = randSig(0.5, 9.5, 3), m = randSig(1, 20, 3);
    const exact = m * 100 / pc;
    return {
      prompt: `<p>Una soluzione acquosa di cloruro di sodio è al <span class="num">${formatSig(pc, 3)}%</span> in massa. Quanti grammi di soluzione contengono <span class="num">${formatSig(m, 3)} g</span> di NaCl?</p>`,
      expected: { exact, sig: 3 }, unit: 'g',
      steps: `<p class="chain">${formatSig(m, 3)} ${cancel('g NaCl')} × ${frac('100 g soluzione', `${formatSig(pc, 3)} ${cancel('g NaCl')}`)} = <b>${formatSig(exact, 3)} g</b></p><p class="note">Il fattore si capovolge perché ora si parte dal soluto.</p>`,
    };
  }
  const tot = randSig(5, 40, 3), part = randSig(tot * 0.1, tot * 0.8, 3);
  const exact = part / tot * 100;
  return {
    prompt: `<p>L'analisi di un campione di <span class="num">${formatSig(tot, 3)} g</span> trova <span class="num">${formatSig(part, 3)} g</span> di carbonio. Qual è la percentuale in massa di carbonio?</p>`,
    expected: { exact, sig: 3 }, unit: '%',
    steps: `<p>% C = (massa di C / massa del campione) × 100 = (${formatSig(part, 3)} / ${formatSig(tot, 3)}) × 100 = ${formatPlain(exact)} → <b>${formatSig(exact, 3)} %</b></p>`,
  };
}

export const GENERATORI = [
  { id: 'conta', label: 'Contare le c.s.', gen: contaGen },
  { id: 'operazioni', label: 'Calcoli e c.s.', gen: operazioniGen },
  { id: 'conversione', label: 'Conversioni', gen: conversioneGen },
  { id: 'composte', label: 'Unità composte', gen: composteGen },
  { id: 'temperatura', label: 'Temperature', gen: temperaturaGen },
  { id: 'densita', label: 'Densità', gen: densitaGen },
  { id: 'percentuale', label: 'Percentuali', gen: percentualeGen },
];

export { exponentOf };
