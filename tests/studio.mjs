// Verifiche del sito di studio: lettura dei numeri, cifre significative, arrotondamenti, correzione delle risposte,
// conversioni e coerenza degli esercizi generati a caso.

import { parseNumber, formatSig, formatPlace, grade } from '../studio/js/numeri.js';
import { convert, convertTemperature, GENERATORI } from '../studio/js/argomenti/materia-calcoli.js';
import { DOMANDE, SEZIONI } from '../studio/js/argomenti/materia.js';
import { PROGRAMMA } from '../studio/js/programma.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};

// cifre significative secondo le regole del Petrucci
{
  const cases = [['0,004050', 4, 4], ['4006', 4, 4], ['1500', 2, 4], ['1500,', 4, 4], ['1,50 × 10^3', 3, 3], ['100,0', 4, 4],
    ['0,1020', 4, 4], ['6,022 × 10^23', 4, 4], ['1,0300×10⁴', 5, 5], ['-0.0450', 3, 3], ['2,500e3', 4, 4], ['7', 1, 1]];
  const bad = cases.filter(([s, a, b]) => { const p = parseNumber(s); return !p || p.sig !== a || p.sigMax !== b; });
  check('conteggio delle cifre significative', !bad.length, bad.map(c => c[0]).join(', ') || `${cases.length} casi`);
  check('numeri non validi rifiutati', ['abc', '', '1,2,3', '1.2.3', 'e5'].every(s => parseNumber(s) === null));
}
{
  const cases = [[formatSig(1500, 3), '1,50 × 10³'], [formatSig(16.5363, 3), '16,5'], [formatSig(0.00045, 3), '4,50 × 10⁻⁴'], [formatSig(0.0045, 3), '0,00450'],
    [formatSig(1230, 4), '1,230 × 10³'], [formatSig(298.15, 4), '298,2'], [formatSig(904.4898, 3), '904'], [formatPlace(31.123, -1), '31,1'],
    [formatPlace(10004.538, -1), '10004,5'], [formatPlace(-320.422, -2), '−320,42'], [formatSig(2.346, 3), '2,35'], [formatSig(2.344, 3), '2,34']];
  const bad = cases.filter(([a, b]) => a !== b);
  check('arrotondamento e scrittura', !bad.length, bad.map(([a, b]) => `${a} ≠ ${b}`).join('; ') || `${cases.length} casi`);
}
// correzione: valore giusto con c.s. sbagliate → "quasi"; ordine di grandezza sbagliato → errore
{
  const ok = grade('16,5', { exact: 16.5363, sig: 3 }).status === 'ok' && grade('1,65e1', { exact: 16.5363, sig: 3 }).status === 'ok';
  const almost = grade('16,54', { exact: 16.5363, sig: 3 }).status === 'almost' && grade('31,12', { exact: 31.123, place: -1 }).status === 'almost';
  const wrong = grade('1,65', { exact: 16.5363, sig: 3 }).status === 'wrong' && grade('16,6', { exact: 16.5363, sig: 3 }).status === 'wrong';
  check('correzione delle risposte', ok && almost && wrong);
}
// conversioni e temperature (esercizi svolti della pagina)
{
  const c1 = convert('2,50', 'km', 'cm', 'lunghezza').text, c2 = convert('2,50', 'cm³', 'm³', 'volume').text, c3 = convert('1,00', 'm³', 'L', 'volume').text;
  check('conversioni con unità al cubo', c1 === '2,50 × 10⁵' && c2 === '2,50 × 10⁻⁶' && c3 === '1,00 × 10³', `${c1}; ${c2}; ${c3}`);
  const t = [convertTemperature('-195,79', '°C', 'K').text, convertTemperature('-195,79', '°C', '°F').text, convertTemperature('37,0', '°C', 'K').text,
    convertTemperature('98,6', '°F', '°C').text, convertTemperature('25,0', '°C', 'K').text, convertTemperature('0,00', 'K', '°C').text];
  check('conversioni di temperatura', t.join('|') === '77,36|−320,42|310,2|37,0|298,2|−273,15', t.join(' | '));
}
// ogni esercizio generato si corregge come giusto con la propria soluzione, e con 10× come sbagliato
{
  let bad = 0, n = 0;
  for (const g of GENERATORI) for (let i = 0; i < 2000; i++) {
    const e = g.gen(), ex = e.expected; n++;
    const ans = ex.integer ? String(ex.exact) : ex.place != null ? formatPlace(ex.exact, ex.place) : formatSig(ex.exact, ex.sig);
    if (grade(ans, ex).status !== 'ok' || !e.prompt || !e.steps || /NaN|undefined|Infinity/.test(e.prompt + e.steps)) bad++;
    if (!ex.integer && grade(String(ex.exact * 10).replace('.', ','), ex).status === 'ok') bad++;
  }
  check('esercizi generati coerenti', bad === 0, `${n} esercizi, ${bad} incoerenti`);
}
// quiz e programma
{
  const bad = DOMANDE.filter(q => !(q.correct >= 0 && q.correct < q.options.length) || new Set(q.options).size !== q.options.length || !q.why);
  check('domande del quiz ben formate', !bad.length && DOMANDE.length >= 15, `${DOMANDE.length} domande`);
  check('programma completo', PROGRAMMA.length === 23 && PROGRAMMA.reduce((s, a) => s + a.ore, 0) === 80 && new Set(PROGRAMMA.map(a => a.simbolo)).size === 23 && SEZIONI.length > 10);
}

if (failures) { console.log(`\n${failures} verifiche fallite`); process.exit(1); }
console.log('\nSito di studio: tutte le verifiche superate');
