// Isolated teaching calculations. No changes to the quantum or sandbox engines.
import { ELEMENTS } from '../physics/elements.js';

export const NA = 6.02214076e23;
export const SYMBOLS = new Map(ELEMENTS.map((e) => [e.symbol, e]));
const fail = (message) => {
  throw new Error(message);
};
export function number(
  raw,
  { min = -Infinity, max = Infinity, integer = false } = {},
) {
  const s = String(raw).trim().replace(',', '.');
  const n = s === '' ? NaN : Number(s);
  if (
    !Number.isFinite(n) ||
    n < min ||
    n > max ||
    (integer && !Number.isSafeInteger(n))
  )
    fail(
      `Inserisci un numero ${integer ? 'intero ' : ''}valido${Number.isFinite(min) ? `, almeno ${min}` : ''}${Number.isFinite(max) ? ` e non oltre ${max}` : ''}.`,
    );
  return n;
}

export function parseFormula(raw) {
  let s = String(raw)
    .trim()
    .replace(/\s/g, '')
    .replace(/[₀-₉]/g, (c) => String('₀₁₂₃₄₅₆₇₈₉'.indexOf(c)))
    .replace(/[−–]/g, '-');
  if (!s || s.length > 200) fail('La formula è vuota o troppo lunga.');
  const phase = s.match(/\((aq|s|l|g)\)$/)?.[1] ?? '';
  if (phase) s = s.replace(/\((aq|s|l|g)\)$/, '');
  // Superscripts denote a charge, subscripts an atom count.
  s = s.replace(
    /([⁰¹²³⁴⁵⁶⁷⁸⁹]*[⁺⁻])$/,
    (v) =>
      '^' +
      [...v]
        .map((c) =>
          c === '⁺' ? '+' : c === '⁻' ? '-' : String('⁰¹²³⁴⁵⁶⁷⁸⁹'.indexOf(c)),
        )
        .join(''),
  );
  let charge = 0,
    body = s;
  let m = s.match(/\^(\d*)([+-])$/);
  if (m) {
    charge =
      (m[2] === '+' ? 1 : -1) *
      (m[1] ? number(m[1], { min: 1, max: 100, integer: true }) : 1);
    body = s.slice(0, m.index);
  } else if ((m = s.match(/([+-])(\d+)$/))) {
    charge =
      (m[1] === '+' ? 1 : -1) *
      number(m[2], { min: 1, max: 100, integer: true });
    body = s.slice(0, m.index);
  } else if ((m = s.match(/^([A-Z][a-z]?)(\d+)([+-])$/))) {
    body = m[1];
    charge =
      number(m[2], { min: 1, max: 100, integer: true }) *
      (m[3] === '+' ? 1 : -1);
  } else if (/[+-]$/.test(s)) {
    charge = s.endsWith('+') ? 1 : -1;
    body = s.slice(0, -1);
  }
  if (body === 'e' && charge === -1 && !phase)
    return {
      formula: 'e^-',
      body: 'e',
      atoms: {},
      charge: -1,
      phase: '',
      electron: true,
    };
  if (/[+^\-]/.test(body))
    fail('Carica non valida: usa Fe^3+, SO4^2- oppure SO4-2.');
  const atoms = {};
  const merge = (target, source, n) => {
    for (const [el, count] of Object.entries(source)) {
      target[el] = (target[el] ?? 0) + count * n;
      if (target[el] > 100000)
        fail('Formula troppo grande per questo banco didattico.');
    }
  };
  const parts = body.split(/[·:]/);
  for (let pi = 0; pi < parts.length; pi++) {
    let part = parts[pi],
      multiplier = 1;
    if (pi > 0) {
      const lead = part.match(/^\d+/);
      if (lead) {
        multiplier = number(lead[0], { min: 1, max: 10000, integer: true });
        part = part.slice(lead[0].length);
      }
    }
    let i = 0;
    const count = () => {
      const start = i;
      while (/\d/.test(part[i] ?? '')) i++;
      return i === start
        ? 1
        : number(part.slice(start, i), { min: 1, max: 10000, integer: true });
    };
    const group = (closing = '', depth = 0) => {
      if (depth > 8) fail('Troppe parentesi annidate.');
      const result = {};
      let entries = 0;
      while (i < part.length && part[i] !== closing) {
        const c = part[i];
        if (c === '(' || c === '[') {
          i++;
          const end = c === '(' ? ')' : ']';
          const inner = group(end, depth + 1);
          if (part[i] !== end)
            fail('Parentesi non chiuse o non corrispondenti.');
          i++;
          merge(result, inner, count());
        } else {
          const token = part.slice(i).match(/^[A-Z][a-z]?/);
          if (!token || !SYMBOLS.has(token[0]))
            fail(
              `Simbolo non riconosciuto vicino a «${part.slice(i)}». Maiuscole e minuscole contano.`,
            );
          i += token[0].length;
          merge(result, { [token[0]]: 1 }, count());
        }
        entries++;
      }
      if (!entries) fail('Gruppo o formula vuoti.');
      return result;
    };
    const result = group();
    if (i !== part.length) fail('Parentesi non corrispondenti.');
    merge(atoms, result, multiplier);
  }
  const formula =
    body +
    (charge
      ? `^${Math.abs(charge) === 1 ? '' : Math.abs(charge)}${charge > 0 ? '+' : '-'}`
      : '');
  return { formula, body, atoms, charge, phase, electron: false };
}

export function molarMass(raw) {
  const sp = typeof raw === 'string' ? parseFormula(raw) : raw;
  if (sp.electron)
    fail(
      'Per l’elettrone usa il bilancio di carica, non una massa molare chimica.',
    );
  let M = 0;
  for (const [el, count] of Object.entries(sp.atoms)) {
    const mass = Number(SYMBOLS.get(el).mass);
    if (!Number.isFinite(mass))
      fail(
        `${el} non ha una massa atomica standard: occorre specificare l’isotopo. Il numero tra parentesi quadre non è una massa media.`,
      );
    M += mass * count;
  }
  return M;
}

export function parseEquation(raw) {
  const parts = String(raw)
    .replace(/⟶|⟹|=>|=/g, '->')
    .split(/->|→/);
  if (parts.length !== 2)
    fail(
      'Usa una sola freccia: H2 + O2 -> H2O. Separa le specie con spazio + spazio.',
    );
  const side = (text) =>
    text
      .trim()
      .split(/\s+\+\s+/)
      .map((v) => {
        const m = v.trim().match(/^(\d+)\s*(?=[A-Z(e])/);
        const coeff = m
          ? number(m[1], { min: 1, max: 1000000, integer: true })
          : 1;
        return {
          ...parseFormula(m ? v.trim().slice(m[0].length) : v.trim()),
          coeff,
        };
      });
  const left = side(parts[0]),
    right = side(parts[1]);
  if (left.length + right.length > 16)
    fail('Usa al massimo 16 specie per questo esercizio.');
  return { left, right };
}

export function inventory(side) {
  const atoms = {};
  let charge = 0;
  for (const sp of side) {
    const c = number(sp.coeff, { min: 1, max: 1000000, integer: true });
    for (const [el, n] of Object.entries(sp.atoms))
      atoms[el] = (atoms[el] ?? 0) + n * c;
    charge += sp.charge * c;
  }
  return { atoms, charge };
}
export function checkBalance(left, right) {
  if (!left.length || !right.length)
    fail('Servono almeno un reagente e un prodotto.');
  const L = inventory(left),
    R = inventory(right);
  const rows = [
    ...new Set([...Object.keys(L.atoms), ...Object.keys(R.atoms)]),
  ].map((el) => ({ el, left: L.atoms[el] ?? 0, right: R.atoms[el] ?? 0 }));
  rows.push({ el: 'Carica', left: L.charge, right: R.charge });
  return {
    balanced: rows.every((r) => r.left === r.right),
    rows,
    left: L,
    right: R,
  };
}

// Exact rational RREF. Never round a floating nullspace into a guessed reaction.
const abs = (x) => (x < 0n ? -x : x);
const gcd = (a, b) => (b ? gcd(b, a % b) : abs(a));
const frac = (n, d = 1n) => {
  if (d < 0n) {
    n = -n;
    d = -d;
  }
  const g = gcd(n, d);
  return [n / g, d / g];
};
const add = (a, b) => frac(a[0] * b[1] + b[0] * a[1], a[1] * b[1]);
const mul = (a, b) => frac(a[0] * b[0], a[1] * b[1]);
const neg = (a) => [-a[0], a[1]];
const div = (a, b) => mul(a, frac(b[1], b[0]));
export function solveBalance(left, right) {
  const species = [...left, ...right];
  if (!left.length || !right.length || species.length > 16)
    fail('Indica da 2 a 16 specie, su entrambi i lati.');
  const els = [...new Set(species.flatMap((sp) => Object.keys(sp.atoms)))];
  const A = [
    ...els.map((el) =>
      species.map((s, i) =>
        frac(BigInt((s.atoms[el] ?? 0) * (i < left.length ? 1 : -1))),
      ),
    ),
    species.map((s, i) => frac(BigInt(s.charge * (i < left.length ? 1 : -1)))),
  ];
  const pivots = [];
  let row = 0;
  for (let col = 0; col < species.length && row < A.length; col++) {
    const p = A.findIndex((r, i) => i >= row && r[col][0] !== 0n);
    if (p < 0) continue;
    [A[row], A[p]] = [A[p], A[row]];
    const divisor = A[row][col];
    A[row] = A[row].map((v) => div(v, divisor));
    for (let i = 0; i < A.length; i++)
      if (i !== row && A[i][col][0]) {
        const factor = A[i][col];
        A[i] = A[i].map((v, j) => add(v, neg(mul(factor, A[row][j]))));
      }
    pivots.push(col);
    row++;
  }
  const free = species.map((_, i) => i).filter((i) => !pivots.includes(i));
  if (!free.length)
    fail(
      'Queste specie non ammettono un bilanciamento: controlla le formule e i prodotti.',
    );
  if (free.length !== 1)
    fail(
      'Il bilancio non è univoco. Separa le reazioni indipendenti o aggiungi vincoli chimici.',
    );
  const x = species.map(() => frac(0n));
  x[free[0]] = frac(1n);
  pivots.forEach((p, i) => {
    x[p] = neg(A[i][free[0]]);
  });
  let lcm = 1n;
  for (const [, d] of x) lcm = (lcm / gcd(lcm, d)) * d;
  let integers = x.map(([n, d]) => n * (lcm / d));
  if (integers.every((n) => n < 0n)) integers = integers.map((n) => -n);
  if (integers.some((n) => n <= 0n))
    fail(
      'Non esiste una soluzione con tutte queste specie sui lati indicati e coefficienti positivi.',
    );
  const common = integers.reduce((g, n) => gcd(g, n));
  const values = integers.map((n) =>
    number(n / common, { min: 1, max: 1000000, integer: true }),
  );
  return {
    left: left.map((s, i) => ({ ...s, coeff: values[i] })),
    right: right.map((s, i) => ({ ...s, coeff: values[left.length + i] })),
  };
}
export function equationText(left, right) {
  const side = (a) =>
    a
      .map(
        (s) =>
          `${s.coeff === 1 ? '' : s.coeff + ' '}${s.formula}${s.phase ? '(' + s.phase + ')' : ''}`,
      )
      .join(' + ');
  return side(left) + ' -> ' + side(right);
}

export function stoichiometry(
  left,
  right,
  amounts,
  fraction = 1,
  yieldPercent = 100,
) {
  if (!checkBalance(left, right).balanced)
    fail(
      'Prima bilancia atomi e carica: altrimenti la stechiometria non è definita.',
    );
  if (amounts.length !== left.length)
    fail('Indica la quantità di ogni reagente.');
  fraction = number(fraction, { min: 0, max: 1 });
  yieldPercent = number(yieldPercent, { min: 0, max: 100 });
  if (amounts.some((a) => !['g', 'mol'].includes(a.unit)))
    fail('Usa quantità in g oppure mol.');
  const initial = amounts.map((a, i) =>
    a.unit === 'g'
      ? number(a.value, { min: 0, max: 1e9 }) / molarMass(left[i])
      : number(a.value, { min: 0, max: 1e9 }),
  );
  const capacities = initial.map((n, i) => n / left[i].coeff);
  const maximum = Math.min(...capacities),
    extent = maximum * fraction;
  const limiting = capacities
    .map((v, i) =>
      Math.abs(v - maximum) <= 1e-9 * Math.max(1, maximum) ? i : -1,
    )
    .filter((i) => i >= 0);
  const remaining = initial.map((n, i) =>
    Math.max(0, n - left[i].coeff * extent),
  );
  const products = right.map((s) => ({
    formula: s.formula,
    mol: s.coeff * extent,
    theoreticalMol: s.coeff * maximum,
    grams: s.electron ? null : s.coeff * extent * molarMass(s),
    theoreticalGrams: s.electron ? null : s.coeff * maximum * molarMass(s),
    isolatedGrams: s.electron
      ? null
      : (s.coeff * maximum * molarMass(s) * yieldPercent) / 100,
  }));
  return { initial, remaining, maximum, extent, limiting, products };
}

export function oxidation(raw) {
  const sp = typeof raw === 'string' ? parseFormula(raw) : raw;
  if (sp.electron)
    fail(
      'L’elettrone ha carica -1; non è un elemento con un numero di ossidazione.',
    );
  const els = Object.keys(sp.atoms),
    states = {};
  if (els.length === 1)
    return {
      sp,
      states: { [els[0]]: sp.charge / sp.atoms[els[0]] },
      note: sp.charge
        ? 'La somma dei numeri di ossidazione è la carica dello ione.'
        : 'Una specie elementare ha numero di ossidazione 0.',
    };
  const fixed = {
    Li: 1,
    Na: 1,
    K: 1,
    Rb: 1,
    Cs: 1,
    Be: 2,
    Mg: 2,
    Ca: 2,
    Sr: 2,
    Ba: 2,
    Al: 3,
    Zn: 2,
    Ag: 1,
    F: -1,
  };
  els.forEach((el) => {
    if (el in fixed) states[el] = fixed[el];
  });
  let note =
    'La somma nᵢ × n.o.ᵢ deve uguagliare la carica. H normalmente +1, O normalmente -2.';
  if (sp.atoms.H) {
    const hydride =
      els.length === 2 &&
      els.some(
        (el) =>
          el !== 'H' &&
          ['Li', 'Na', 'K', 'Rb', 'Cs', 'Ca', 'Ba', 'Mg'].includes(el),
      );
    states.H = hydride ? -1 : 1;
    if (hydride) note = 'Idruro metallico: H ha numero di ossidazione -1.';
  }
  if (sp.atoms.O) {
    states.O = ['H2O2', 'Na2O2', 'K2O2', 'BaO2'].includes(sp.body)
      ? -1
      : sp.body === 'KO2'
        ? -0.5
        : -2;
    if (states.O !== -2)
      note =
        states.O === -1
          ? 'Perossido: O ha numero di ossidazione -1.'
          : 'Superossido: numero di ossidazione medio di O = -1/2.';
    if (els.length === 2 && sp.atoms.F) {
      delete states.O;
      note =
        'Il fluoro resta -1 anche nei composti con ossigeno; ricava O dal bilancio.';
    }
  }
  for (const hal of ['Cl', 'Br', 'I'])
    if (sp.atoms[hal] && !sp.atoms.O && !sp.atoms.F && els.length === 2)
      states[hal] = -1;
  if (/(?:SO4|\(SO4\)\d*)$/.test(sp.body)) {
    states.S = 6;
    note += ' Gruppo solfato noto: S è +6.';
  }
  if (/(?:NO3|\(NO3\)\d*)$/.test(sp.body)) {
    states.N = 5;
    note += ' Gruppo nitrato noto: N è +5.';
  }
  if (sp.body === '(NH4)2Cr2O7') {
    states.N = -3;
    states.Cr = 6;
    note = 'Composto noto: ammonio NH₄⁺ (N −3) e dicromato Cr₂O₇²⁻ (Cr +6).';
  }
  const unknown = els.filter((el) => !(el in states));
  const assigned = () =>
    els
      .filter((el) => el in states)
      .reduce((a, el) => a + sp.atoms[el] * states[el], 0);
  if (unknown.length === 1) {
    states[unknown[0]] = (sp.charge - assigned()) / sp.atoms[unknown[0]];
    if (sp.atoms[unknown[0]] > 1)
      note += ` Per ${unknown[0]} il bilancio determina un valore medio: gli stati dei singoli siti possono richiedere la struttura.`;
  } else if (unknown.length)
    return {
      sp,
      states,
      unresolved: unknown,
      note: 'Più incognite: la formula da sola non basta. Servono struttura o ulteriori stati noti.',
    };
  if (Math.abs(assigned() - sp.charge) > 1e-8)
    fail(
      'Le regole standard non sono compatibili con questa formula. Serve conoscere la struttura o un’eccezione.',
    );
  if (els.some((el) => !Number.isInteger(states[el])))
    note +=
      ' Valore medio: non implica che ogni atomo abbia uno stato frazionario (possibile valenza mista).';
  return { sp, states, note };
}

export function redoxLedger(left, right) {
  if (!checkBalance(left, right).balanced)
    fail('Bilancia prima la reazione completa, compresa la carica.');
  const bins = (side) => {
    const result = new Map();
    side
      .filter((s) => !s.electron)
      .forEach((sp) => {
        const ox = oxidation(sp);
        if (ox.unresolved)
          fail(
            `Stati non determinati in ${sp.formula}: specifica una reazione con stati noti.`,
          );
        for (const [el, n] of Object.entries(sp.atoms)) {
          if (!result.has(el)) result.set(el, new Map());
          const map = result.get(el);
          map.set(ox.states[el], (map.get(ox.states[el]) ?? 0) + n * sp.coeff);
        }
      });
    return result;
  };
  const L = bins(left),
    R = bins(right),
    changes = [];
  for (const [el, froms] of L) {
    const tos = R.get(el);
    for (const [state, n] of froms) {
      const same = Math.min(n, tos?.get(state) ?? 0);
      froms.set(state, n - same);
      if (tos?.has(state)) tos.set(state, tos.get(state) - same);
    }
    const from = [...froms].filter(([, n]) => n > 0),
      to = [...(tos ?? [])].filter(([, n]) => n > 0);
    if (from.length > 1 && to.length > 1)
      fail(
        `Più percorsi redox per ${el}: il bilancio globale non identifica un trasferimento unico.`,
      );
    for (const f of from)
      for (const t of to) {
        const n = Math.min(f[1], t[1]);
        if (!n) continue;
        changes.push({
          el,
          from: f[0],
          to: t[0],
          atoms: n,
          electrons: (t[0] - f[0]) * n,
        });
        f[1] -= n;
        t[1] -= n;
      }
  }
  const lost = changes
    .filter((c) => c.electrons > 0)
    .reduce((s, c) => s + c.electrons, 0);
  const gained = -changes
    .filter((c) => c.electrons < 0)
    .reduce((s, c) => s + c.electrons, 0);
  const external = [...left, ...right].some((s) => s.electron);
  return {
    changes,
    lost,
    gained,
    external,
    balanced: external || Math.abs(lost - gained) < 1e-9,
  };
}

export function ionicFormula(cation, anion) {
  const c = parseFormula(cation),
    a = parseFormula(anion);
  if (c.charge <= 0 || a.charge >= 0)
    fail('Scegli un catione positivo e un anione negativo.');
  const g = Number(gcd(BigInt(c.charge), BigInt(-a.charge)));
  const nc = -a.charge / g,
    na = c.charge / g;
  const group = (sp, n) =>
    n === 1
      ? sp.body
      : `${Object.keys(sp.atoms).length > 1 || Object.values(sp.atoms).some((v) => v > 1) ? '(' + sp.body + ')' : sp.body}${n}`;
  return {
    formula: group(c, nc) + group(a, na),
    cationCount: nc,
    anionCount: na,
  };
}

export function dilution(concentration, initialMl, finalMl) {
  const C = number(concentration, { min: 0, max: 1000 }),
    Vi = number(initialMl, { min: 0.001, max: 1e9 }),
    Vf = number(finalMl, { min: 0.001, max: 1e9 });
  if (Vf < Vi)
    fail(
      'Una diluizione richiede V finale ≥ V iniziale. Per concentrare occorre un altro processo.',
    );
  return {
    mol: (C * Vi) / 1000,
    finalConcentration: (C * Vi) / Vf,
    addedMl: Vf - Vi,
  };
}
export function isotopeMean(m1, m2, percent) {
  const a = number(m1, { min: 0.001 }),
    b = number(m2, { min: 0.001 }),
    x = number(percent, { min: 0, max: 100 }) / 100;
  return a * x + b * (1 - x);
}
export function significantFigures(raw) {
  const s = String(raw).trim().replace(',', '.');
  if (!/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(s))
    fail('Scrivi un numero, anche in notazione scientifica: 1.20e3.');
  let mantissa = s.split(/e/i)[0].replace(/^[+-]/, '');
  const decimal = mantissa.includes('.');
  mantissa = mantissa.replace('.', '').replace(/^0+/, '');
  if (!decimal && !/e/i.test(s)) mantissa = mantissa.replace(/0+$/, '');
  return mantissa.length || 1;
}

export function empiricalFormula(entries) {
  if (entries.length < 2 || entries.length > 8)
    fail('Indica da 2 a 8 elementi distinti.');
  if (new Set(entries.map((e) => e.el)).size !== entries.length)
    fail('Ogni elemento va inserito una sola volta.');
  const mol = entries.map(
    ({ el, mass }) => number(mass, { min: 1e-12, max: 1e9 }) / molarMass(el),
  );
  const smallest = Math.min(...mol),
    ratios = mol.map((n) => n / smallest);
  for (let factor = 1; factor <= 12; factor++) {
    const rounded = ratios.map((r) => Math.round(r * factor));
    const error = Math.max(
      ...rounded.map(
        (n, i) => Math.abs(n - ratios[i] * factor) / (ratios[i] * factor),
      ),
    );
    if (error < 0.02 && rounded.every((n) => n > 0 && n <= 1000)) {
      const divisor = Number(rounded.reduce((g, n) => gcd(g, BigInt(n)), 0n));
      const counts = rounded.map((n) => n / divisor);
      return {
        mol,
        ratios,
        error,
        formula: entries
          .map((e, i) => e.el + (counts[i] === 1 ? '' : counts[i]))
          .join(''),
      };
    }
  }
  fail(
    'Nessun piccolo rapporto intero compatibile entro il 2%. Verifica le masse e la precisione dei dati.',
  );
}

export function chromatography(K, phaseRatio = 1, front = 1) {
  K = number(K, { min: 0, max: 10000 });
  phaseRatio = number(phaseRatio, { min: 0.001, max: 10000 });
  front = number(front, { min: 0, max: 1 });
  const Rf = 1 / (1 + K * phaseRatio);
  return { Rf, fractionAlongPaper: front * Rf };
}
export function measurementStatistics(values) {
  if (values.length < 2 || values.length > 1000)
    fail('Servono da 2 a 1000 misure.');
  const data = values.map((v) => number(v)),
    mean = data.reduce((s, v) => s + v, 0) / data.length;
  const sd = Math.sqrt(
    data.reduce((s, v) => s + (v - mean) ** 2, 0) / (data.length - 1),
  );
  return { mean, sd, count: data.length };
}
export function waterPhase(celsius) {
  const T = number(celsius, { min: -273.15, max: 1000 });
  if (Math.abs(T) < 1e-8) return 'solido + liquido';
  if (Math.abs(T - 100) < 1e-8) return 'liquido + gas';
  return T < 0 ? 'solido' : T < 100 ? 'liquido' : 'gas';
}
export function measuredCalculation(aRaw, bRaw, operation) {
  const a = number(aRaw),
    b = number(bRaw);
  const sa = significantFigures(aRaw),
    sb = significantFigures(bRaw);
  if (!['multiply', 'divide', 'add', 'subtract'].includes(operation))
    fail('Operazione non riconosciuta.');
  if (operation === 'divide' && b === 0) fail('Non puoi dividere per zero.');
  const value =
    operation === 'multiply'
      ? a * b
      : operation === 'divide'
        ? a / b
        : operation === 'add'
          ? a + b
          : a - b;
  if (!Number.isFinite(value))
    fail('Il risultato supera il campo numerico di questo calcolo.');
  if (operation === 'multiply' || operation === 'divide') {
    const precision = Math.min(15, sa, sb);
    return {
      value,
      reported: value.toPrecision(precision),
      rule: `${precision} cifre significative, come il dato meno preciso.`,
    };
  }
  const places = (raw) => {
    const [mantissa, exp = '0'] = String(raw)
      .trim()
      .replace(',', '.')
      .split(/e/i);
    return (mantissa.split('.')[1]?.length ?? 0) - Number(exp);
  };
  const decimals = Math.min(places(aRaw), places(bRaw));
  if (decimals > 50 || decimals < -50)
    fail('Usa esponenti e precisioni fra -50 e 50 per questo esempio.');
  const rounded =
    decimals >= 0
      ? value.toFixed(decimals)
      : String(Math.round(value / 10 ** -decimals) * 10 ** -decimals);
  return {
    value,
    reported: rounded,
    rule: `Arrotondamento alla posizione decimale ${decimals}, come il dato con minore precisione decimale.`,
  };
}

// Repulsive Coulomb hyperbola; a = k/(2E), e = sqrt(1+(b/a)^2).
// t is proportional to e*sinh(u)+u, so playback follows physical time.
export function rutherford({ Z = 79, E = 5, b = 20, points = 321 } = {}) {
  Z = number(Z, { min: 1, max: 118, integer: true });
  E = number(E, { min: 0.01, max: 1000 });
  b = number(b, { min: 0, max: 10000 });
  const k = 2 * Z * 1.439964548,
    a = k / (2 * E),
    eccentricity = Math.hypot(1, b / a);
  const theta = b === 0 ? Math.PI : 2 * Math.atan2(a, b),
    closest = a * (1 + eccentricity);
  if (b === 0) {
    const uMax = 2.8,
      maxTime = Math.sinh(uMax) + uMax;
    const path = Array.from({ length: points }, (_, i) => {
      const u = uMax * (1 - (2 * i) / (points - 1));
      return {
        x: -a * (1 + Math.cosh(u)),
        y: 0,
        t: (maxTime - Math.sinh(u) - u) / (2 * maxTime),
      };
    });
    return { theta, closest, a, eccentricity, path };
  }
  const rot = Math.PI - Math.atan2(b, a),
    co = Math.cos(rot),
    si = Math.sin(rot),
    uMax = 3.2;
  const maxTime = eccentricity * Math.sinh(uMax) + uMax;
  const path = Array.from({ length: points }, (_, i) => {
    const u = uMax * (1 - (2 * i) / (points - 1)),
      x = a * (eccentricity + Math.cosh(u)),
      y = b * Math.sinh(u);
    return {
      x: x * co - y * si,
      y: x * si + y * co,
      t: (maxTime - eccentricity * Math.sinh(u) - u) / (2 * maxTime),
    };
  });
  return { theta, closest, a, eccentricity, path };
}
