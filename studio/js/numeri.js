// Numeri misurati: lettura come li scrive uno studente, cifre significative, arrotondamento,
// scrittura all'italiana (virgola decimale, notazione scientifica con ×10ⁿ) e correzione delle risposte.

const SUP_IN = { '⁻': '-', '⁺': '+', '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9' };
const SUP_OUT = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };

export const sup = (n) => String(n).replace('−', '-').split('').map(c => SUP_OUT[c] ?? c).join('');
export const pow10 = (e) => `10${sup(e)}`;

// Legge "0,004050", "1.5e3", "1,50 × 10^3", "1,50·10⁻³", "-2,5". Restituisce null se non è un numero.
export function parseNumber(input) {
  let s = String(input ?? '').trim().replace(/\s+/g, '').replace(/[−–]/g, '-');
  if (!s) return null;
  s = s.replace(/[⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹]+/g, m => '^' + [...m].map(c => SUP_IN[c]).join(''));
  let mant = s, exp = 0, sci = false;
  const m = s.match(/^(.*?)(?:[eE]|[x×·*]10\^?)([+-]?\d+)$/);
  if (m) { mant = m[1]; exp = parseInt(m[2], 10); sci = true; }
  mant = mant.replace(/,/g, '.');
  const mm = mant.match(/^([+-]?)(\d*)(?:(\.)(\d*))?$/);
  if (!mm || (!mm[2] && !mm[4])) return null;
  const [, sign, intPart, point, fracPart = ''] = mm;
  const value = Number(`${sign}${intPart || '0'}.${fracPart || '0'}`) * 10 ** exp;
  if (!Number.isFinite(value)) return null;
  return { value: Number(value.toPrecision(15)), sign, intPart, fracPart, hasPoint: !!point, exp, sci, ...analyse(intPart, fracPart, !!point, exp) };
}

// Classifica ogni cifra della mantissa. reason: 'zero-iniziale', 'non-nulla', 'zero-intermedio',
// 'zero-finale' (dopo la virgola o con la virgola scritta), 'zero-ambiguo' (intero senza virgola).
function analyse(intPart, fracPart, hasPoint, exp) {
  const digits = [...intPart, ...fracPart];
  const first = digits.findIndex(c => c !== '0');
  const last = digits.length - 1 - [...digits].reverse().findIndex(c => c !== '0');
  const marks = digits.map((c, i) => {
    if (first < 0 || i < first) return 'zero-iniziale';
    if (c !== '0') return 'non-nulla';
    if (i < last) return 'zero-intermedio';
    return hasPoint ? 'zero-finale' : 'zero-ambiguo';
  });
  const sig = marks.filter(r => r === 'non-nulla' || r === 'zero-intermedio' || r === 'zero-finale').length;
  const amb = marks.filter(r => r === 'zero-ambiguo').length;
  // posizione (potenza di 10) dell'ultima cifra certa: per gli zeri ambigui si prende l'ultima non nulla
  const lastDigitIndex = amb ? last : digits.length - 1;
  const lastPlace = first < 0 ? -fracPart.length + exp : (intPart.length - 1 - lastDigitIndex) + exp;
  return { marks, sig: first < 0 ? null : sig, sigMax: first < 0 ? null : sig + amb, ambiguous: amb > 0, lastPlace };
}

// Arrotondamento a n cifre significative (metà per eccesso, in valore assoluto).
export function roundSig(x, n) {
  if (x === 0 || !Number.isFinite(x)) return x;
  const e = exponentOf(x);
  return roundPlace(x, e - n + 1);
}

// Arrotondamento alla potenza di 10 indicata (place = −2 → centesimi).
export function roundPlace(x, place) {
  const k = Math.abs(x) / 10 ** place;
  const r = Math.floor(Number(k.toPrecision(12)) + 0.5);
  const v = Math.sign(x) * r * 10 ** place;
  return Number(v.toPrecision(15));
}

export function exponentOf(x) {
  if (x === 0) return 0;
  return parseInt(Math.abs(x).toExponential(12).split('e')[1], 10);
}

// Posizione dell'ultima cifra significativa di un valore scritto con n cifre significative.
export const placeOfSig = (x, n) => exponentOf(roundSig(x, n)) - n + 1;
// Numero di cifre significative di un valore la cui ultima cifra è alla posizione place.
export const sigOfPlace = (x, place) => Math.max(1, exponentOf(roundPlace(x, place)) - place + 1);

const comma = (s) => s.replace('.', ',').replace('-', '−');

// Scrive x con n cifre significative senza ambiguità: forma decimale quando possibile,
// altrimenti notazione scientifica (1500 con 3 c.s. → "1,50 × 10³").
export function formatSig(x, n, { forceSci = false } = {}) {
  if (x === 0) return n > 1 ? `0,${'0'.repeat(n - 1)}` : '0';
  const r = roundSig(x, n);
  const e = exponentOf(r);
  const dec = n - 1 - e;
  const plain = r.toFixed(Math.max(0, dec));
  const ambiguousInt = dec === 0 && /0$/.test(plain);
  if (!forceSci && e >= -3 && e <= 5 && dec >= 0 && !ambiguousInt) return comma(plain);
  const mant = (r / 10 ** e).toFixed(n - 1);
  return `${comma(mant)} × ${pow10(e)}`;
}

// Scrive x arrotondato alla posizione place (regola di addizione e sottrazione).
export function formatPlace(x, place) {
  const r = roundPlace(x, place);
  if (place <= 0) return comma(r.toFixed(-place));
  return formatSig(r, sigOfPlace(r, place));
}

// Scrittura senza regole di arrotondamento, per i passaggi intermedi (cifre di guardia).
export function formatPlain(x, maxSig = 6) {
  if (x === 0) return '0';
  const e = exponentOf(x);
  if (e < -3 || e > 6) {
    const m = Number((x / 10 ** e).toPrecision(maxSig));
    return `${comma(String(m))} × ${pow10(e)}`;
  }
  return comma(String(Number(x.toPrecision(maxSig))));
}

// Correzione di una risposta numerica.
// expected: { exact, sig } (regola di moltiplicazione/divisione) oppure { exact, place } (addizione/sottrazione)
// oppure { exact, integer: true } (risposte intere, es. un numero di cifre significative).
export function grade(input, expected) {
  const p = parseNumber(input);
  if (!p) return { status: 'invalid', msg: 'Scrivi un numero, per esempio 2,50 oppure 2,50e3 oppure 2,50 × 10^3.' };
  if (expected.integer) {
    return p.value === expected.exact ? { status: 'ok', msg: 'Corretto.' } : { status: 'wrong', msg: 'Non è corretto.' };
  }
  const place = expected.place ?? placeOfSig(expected.exact, expected.sig);
  const target = roundPlace(expected.exact, place);
  const shown = expected.place != null ? formatPlace(expected.exact, place) : formatSig(expected.exact, expected.sig);
  const unit = 10 ** place;
  const sameValue = Math.abs(p.value - target) <= unit * 1e-6 + Math.abs(target) * 1e-12;
  const close = Math.abs(p.value - expected.exact) <= unit * 0.5 * 1.0001 + Math.abs(target) * 1e-12;
  const placeOk = p.lastPlace === place || (p.ambiguous && p.lastPlace >= place && place >= p.lastPlace - (p.sigMax - p.sig));
  if (sameValue && placeOk) return { status: 'ok', msg: `Corretto: ${shown}.` };
  if (sameValue || close) {
    const want = expected.place != null ? `l'ultima cifra deve essere quella ${placeName(place)}` : `servono ${expected.sig} cifre significative`;
    const have = p.ambiguous ? 'il numero che hai scritto ha zeri finali ambigui'
      : expected.place != null ? `la tua ultima cifra è quella ${placeName(p.lastPlace)}` : `hai scritto ${p.sig} cifre significative`;
    return { status: 'almost', msg: `Il valore è giusto, ma ${want} (${have}). Risposta: ${shown}.` };
  }
  const rel = Math.abs(p.value / target);
  const powerOff = target !== 0 && Math.abs(Math.log10(rel) - Math.round(Math.log10(rel))) < 0.002 && Math.round(Math.log10(rel)) !== 0;
  return { status: 'wrong', msg: powerOff ? `Le cifre sono giuste ma l'ordine di grandezza no: sei fuori di un fattore ${pow10(Math.round(Math.log10(rel)))}. Ricontrolla i fattori di conversione.` : 'Non è corretto.' };
}

export function placeName(place) {
  const names = { 3: 'delle migliaia', 2: 'delle centinaia', 1: 'delle decine', 0: 'delle unità', '-1': 'dei decimi', '-2': 'dei centesimi', '-3': 'dei millesimi', '-4': 'dei decimillesimi' };
  return names[place] ?? `della posizione ${pow10(place)}`;
}
