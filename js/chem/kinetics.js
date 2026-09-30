// Cinetica chimica: meccanismi scritti liberamente, legge di azione di massa, Arrhenius, integratore rigido.
//
// Sintassi di una riga del meccanismo (una reazione elementare o un passo con legge cinetica nota):
//   2 N2O5 -> 4 NO2 + O2 ; A=4.94e13 Ea=103
//   A + B <=> C          ; k=1e3 kr=10
//   O + O2 + N2 -> O3 + N2 ; k0=5.6e-34 n=-2.6        (k = k0 (T/300)^n, forma usata da IUPAC e JPL)
//   O + O3 -> 2 O2       ; A=8.0e-12 EaR=2060          (E/R in kelvin)
//   Bi210 -> Po210       ; t12=5.012 d                 (tempo di dimezzamento: s, min, h, d, a)
//   2 NO + O2 -> 2 NO2   ; k=... ord.NO=2 ord.O2=1     (ordini diversi dalla stechiometria)
// Parametri: k (costante), A + Ea (kJ/mol) o EaR (K), n per il fattore (T/300)^n; per la reazione inversa
// gli stessi con il suffisso r (kr, Ar, Ear, EaRr, nr). Le righe "costanti: X, Y" tengono fisse le specie X e Y
// (reagenti in grande eccesso, terzo corpo M, fotoni). Il resto della riga dopo # è un commento.
//
// Velocità: r = k(T) Π c_i^o_i con o_i = coefficiente stechiometrico salvo indicazione; dc_j/dt = Σ ν_j r.
// Integratore: metodo di Rosenbrock ROS2 (Verwer, Spee, Blom, Hundsdorfer, SIAM J. Sci. Comput. 20, 1456,
// 1999), L-stabile, del secondo ordine, con stima dell'errore dalla soluzione incorporata del primo ordine e
// jacobiana analitica: adatto a sistemi rigidi come il ciclo di Chapman o l'Oregonator.

export const R_GAS = 8.314462618; // J/(mol K), CODATA 2018 (esatta)
export const KB = 1.380649e-23; // J/K, esatta
export const H_PLANCK = 6.62607015e-34; // J s, esatta
export const N_A = 6.02214076e23; // 1/mol, esatta

const TIME_UNITS = { s: 1, min: 60, h: 3600, d: 86400, a: 365.25 * 86400, y: 365.25 * 86400, anni: 365.25 * 86400, giorni: 86400, ore: 3600 };

const parseNum = (s) => {
  const v = Number(String(s).replace(',', '.').replace(/×10\^?/, 'e').replace(/·10\^?/, 'e'));
  return Number.isFinite(v) ? v : NaN;
};

function parseSide(side, errs, line) {
  const out = [];
  for (let term of side.split(/\s\+\s|^\+\s|\s\+$/)) {
    term = term.trim();
    if (!term || /^(hν|hv|hnu|luce)$/i.test(term)) continue;
    const m = term.match(/^(\d+\/\d+|\d+(?:[.,]\d+)?(?=\s|[A-Z(\[]))?\s*(.+)$/);
    let coef = 1;
    if (m[1]) coef = m[1].includes('/') ? m[1].split('/').reduce((a, b) => +a / +b) : parseNum(m[1]);
    const sp = m[2].trim();
    if (!/^[A-Za-z(\[][\w()[\]+\-−*'·.]*$/.test(sp)) { errs.push(`riga ${line}: specie non valida "${sp}"`); continue; }
    if (!(coef > 0)) { errs.push(`riga ${line}: coefficiente non valido per ${sp}`); continue; }
    const prev = out.find(x => x[0] === sp);
    if (prev) prev[1] += coef; else out.push([sp, coef]);
  }
  return out;
}

/**
 * Legge un meccanismo. Restituisce { reactions, species, fixed, errors, units }.
 * reactions: [{ text, reactants:[[sp,ν]], products:[[sp,ν]], orders:{sp:o}, fwd:{...}, rev:{...}|null }]
 */
export function parseMechanism(text) {
  const reactions = [], errors = [], fixed = new Set();
  let units = 'M';
  const lines = text.split('\n');
  lines.forEach((raw, idx) => {
    const line = idx + 1;
    const s = raw.replace(/#.*$/, '').trim();
    if (!s) return;
    const cm = s.match(/^(costanti|fisse|fissi|constant)\s*:\s*(.*)$/i);
    if (cm) { for (const x of cm[2].split(/[,\s]+/).filter(Boolean)) fixed.add(x); return; }
    const um = s.match(/^unit[àa]\s*:\s*(.*)$/i);
    if (um) { units = /molec/i.test(um[1]) ? 'molec' : 'M'; return; }
    const [eqPart, paramPart = ''] = s.split(';').length > 1 ? [s.slice(0, s.indexOf(';')), s.slice(s.indexOf(';') + 1)] : [s, ''];
    const arrow = eqPart.match(/<=>|⇌|<->|->|→|=>/);
    if (!arrow) { errors.push(`riga ${line}: manca la freccia (->, <=>)`); return; }
    const reversible = ['<=>', '⇌', '<->'].includes(arrow[0]);
    const [lhs, rhs] = [eqPart.slice(0, arrow.index), eqPart.slice(arrow.index + arrow[0].length)];
    const reactants = parseSide(lhs, errors, line), products = parseSide(rhs, errors, line);
    if (!reactants.length) { errors.push(`riga ${line}: nessun reagente`); return; }
    const p = {};
    const orders = {};
    for (const tok of paramPart.replace(/;/g, ' ').split(/\s+/).filter(Boolean)) {
      const kv = tok.match(/^([A-Za-z_.0-9()\[\]+\-]+?)\s*=\s*(.+)$/);
      if (!kv) {
        // unità di tempo dopo t12 (t12=5.012 d)
        if (p.__t12 !== undefined && TIME_UNITS[tok] !== undefined) { p.t12 = p.__t12 * TIME_UNITS[tok]; delete p.__t12; continue; }
        errors.push(`riga ${line}: parametro non riconosciuto "${tok}"`); continue;
      }
      const key = kv[1], val = parseNum(kv[2]);
      if (!Number.isFinite(val)) { errors.push(`riga ${line}: valore non valido per ${key}`); continue; }
      if (key.startsWith('ord.')) { orders[key.slice(4)] = val; continue; }
      if (key === 't12') { p.__t12 = val; continue; }
      p[key] = val;
    }
    if (p.__t12 !== undefined) { p.t12 = p.__t12; delete p.__t12; }
    const fwd = rateSpec(p, '', errors, line);
    const rev = reversible ? rateSpec(p, 'r', errors, line) : null;
    if (!fwd) return;
    if (reversible && !rev) { errors.push(`riga ${line}: per una reazione reversibile servono anche i parametri della reazione inversa (kr, oppure Ar ed Ear)`); return; }
    reactions.push({ text: eqPart.trim(), reactants, products, orders, fwd, rev, line });
  });
  const species = [];
  for (const r of reactions) for (const [sp] of [...r.reactants, ...r.products]) if (!species.includes(sp)) species.push(sp);
  for (const f of fixed) if (!species.includes(f)) errors.push(`la specie costante ${f} non compare nel meccanismo`);
  return { reactions, species, fixed, errors, units };
}

function rateSpec(p, suf, errors, line) {
  const g = (k) => p[k + suf];
  if (suf === '' && p.t12 !== undefined) return { k: Math.LN2 / p.t12, kind: 't12', t12: p.t12 };
  if (g('k') !== undefined) return { k: g('k'), kind: 'k' };
  if (g('k0') !== undefined) return { A: g('k0'), n: g('n') ?? 0, EaR: 0, kind: 'Tn' };
  if (g('A') !== undefined) {
    if (g('Ea') !== undefined) return { A: g('A'), EaR: g('Ea') * 1000 / R_GAS, n: g('n') ?? 0, kind: 'arr' };
    if (g('EaR') !== undefined) return { A: g('A'), EaR: g('EaR'), n: g('n') ?? 0, kind: 'arr' };
    return { A: g('A'), EaR: 0, n: g('n') ?? 0, kind: 'arr' };
  }
  if (suf === '') errors.push(`riga ${line}: manca la costante (k, A ed Ea, k0 e n, oppure t12)`);
  return null;
}

/** Costante cinetica alla temperatura T: k = A (T/300)^n exp(−Ea/RT). */
export function rateConstant(spec, T) {
  if (!spec) return 0;
  if (spec.k !== undefined) return spec.k;
  return spec.A * (spec.n ? Math.pow(T / 300, spec.n) : 1) * Math.exp(-spec.EaR / T);
}

/** Modello numerico: specie variabili, stechiometria, ordini. */
export function buildModel(mech, T) {
  const vars = mech.species.filter(s => !mech.fixed.has(s));
  const vi = new Map(vars.map((s, i) => [s, i]));
  const steps = [];
  for (const r of mech.reactions) {
    const net = new Map();
    for (const [s, v] of r.reactants) net.set(s, (net.get(s) ?? 0) - v);
    for (const [s, v] of r.products) net.set(s, (net.get(s) ?? 0) + v);
    const nu = [...net].filter(([s, v]) => vi.has(s) && Math.abs(v) > 1e-14).map(([s, v]) => [vi.get(s), v]);
    const ordF = r.reactants.map(([s, v]) => [s, r.orders[s] ?? v]);
    for (const [s, o] of Object.entries(r.orders)) if (!ordF.some(x => x[0] === s)) ordF.push([s, o]);
    steps.push({ nu, ord: ordF.filter(([, o]) => o !== 0), k: rateConstant(r.fwd, T) });
    if (r.rev) steps.push({ nu: nu.map(([i, v]) => [i, -v]), ord: r.products.map(([s, v]) => [s, v]), k: rateConstant(r.rev, T) });
  }
  return { vars, vi, steps };
}

/**
 * Integra dc/dt = f(c). c0: {specie: concentrazione}; fisse: valori costanti presi da c0.
 * Restituisce { t:[...], c:[[...]], vars, stats }.
 */
export function simulate(mech, c0, T, tEnd, opts = {}) {
  const { vars, vi, steps } = buildModel(mech, T);
  const n = vars.length;
  const fixedVal = new Map([...mech.fixed].map(s => [s, c0[s] ?? 0]));
  // ogni fattore di concentrazione: indice della variabile o valore fisso
  const fac = steps.map(st => st.ord.map(([s, o]) => (vi.has(s) ? { i: vi.get(s), o } : { v: Math.pow(Math.max(fixedVal.get(s) ?? 0, 0), o), o: 0 })));
  const rates = new Float64Array(steps.length);
  const f = (y, out) => {
    out.fill(0);
    for (let j = 0; j < steps.length; j++) {
      let r = steps[j].k;
      for (const q of fac[j]) r *= q.i === undefined ? q.v : Math.pow(Math.max(y[q.i], 0), q.o);
      rates[j] = r;
      for (const [i, v] of steps[j].nu) out[i] += v * r;
    }
  };
  const jac = (y, J) => {
    for (const row of J) row.fill(0);
    for (let j = 0; j < steps.length; j++) {
      for (let a = 0; a < fac[j].length; a++) {
        const qa = fac[j][a];
        if (qa.i === undefined) continue;
        const ya = Math.max(y[qa.i], 0);
        // ∂r/∂y_a = k o_a y_a^(o_a−1) Π_{b≠a} y_b^o_b
        let d = steps[j].k * qa.o * (qa.o === 1 ? 1 : ya > 0 ? Math.pow(ya, qa.o - 1) : 0);
        for (let b = 0; b < fac[j].length; b++) {
          if (b === a) continue;
          const qb = fac[j][b];
          d *= qb.i === undefined ? qb.v : Math.pow(Math.max(y[qb.i], 0), qb.o);
        }
        for (const [i, v] of steps[j].nu) J[i][qa.i] += v * d;
      }
    }
  };
  const y = Float64Array.from(vars, s => c0[s] ?? 0);
  const scale = Math.max(1e-300, ...vars.map(s => Math.abs(c0[s] ?? 0)), ...[...fixedVal.values()].map(Math.abs));
  const rtol = opts.rtol ?? 1e-6, atol = opts.atol ?? 1e-12 * scale;
  const gamma = 1 + 1 / Math.SQRT2;
  const J = Array.from({ length: n }, () => new Float64Array(n));
  const A = Array.from({ length: n }, () => new Float64Array(n));
  const f0 = new Float64Array(n), f1 = new Float64Array(n), k1 = new Float64Array(n), k2 = new Float64Array(n), y1 = new Float64Array(n), yn = new Float64Array(n);
  const tOut = [0], cOut = [Array.from(y)];
  const maxPoints = opts.maxPoints ?? 20000;
  let t = 0;
  f(y, f0);
  let h = opts.h0 ?? Math.min(tEnd / 1000, 1e-3 * (scale / Math.max(1e-300, Math.max(...f0.map(Math.abs)))) || tEnd / 1000);
  h = Math.max(h, tEnd * 1e-14);
  let accepted = 0, rejected = 0, evals = 0;
  while (t < tEnd && accepted + rejected < (opts.maxSteps ?? 400000)) {
    if (t + h > tEnd) h = tEnd - t;
    f(y, f0); jac(y, J); evals++;
    for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) A[i][k] = (i === k ? 1 : 0) - gamma * h * J[i][k];
    const LU = luDecompose(A);
    if (!LU) { h *= 0.25; rejected++; continue; }
    k1.set(f0); luSolve(LU, k1);
    for (let i = 0; i < n; i++) y1[i] = y[i] + h * k1[i];
    f(y1, f1);
    for (let i = 0; i < n; i++) k2[i] = f1[i] - 2 * k1[i];
    luSolve(LU, k2);
    let err = 0;
    for (let i = 0; i < n; i++) {
      yn[i] = y[i] + 1.5 * h * k1[i] + 0.5 * h * k2[i];
      const e = 0.5 * h * (k1[i] + k2[i]); // differenza con la soluzione del primo ordine y + h k1
      const sc = atol + rtol * Math.max(Math.abs(y[i]), Math.abs(yn[i]));
      err += (e / sc) ** 2;
    }
    err = Math.sqrt(err / Math.max(n, 1));
    if (!Number.isFinite(err)) { h *= 0.25; rejected++; continue; }
    if (err <= 1) {
      t += h;
      y.set(yn);
      accepted++;
      tOut.push(t); cOut.push(Array.from(y));
      // si tengono tutti i passi accettati (densi dove la soluzione cambia in fretta); oltre il limite si dimezzano
      if (tOut.length > maxPoints) {
        for (let i = 1, j = 1; i < tOut.length; i += 2, j++) { tOut[j] = tOut[i]; cOut[j] = cOut[i]; }
        tOut.length = cOut.length = Math.ceil(tOut.length / 2);
        if (tOut[tOut.length - 1] !== t) { tOut.push(t); cOut.push(Array.from(y)); }
      }
    } else rejected++;
    h *= Math.min(4, Math.max(0.2, 0.9 / Math.sqrt(Math.max(err, 1e-10))));
  }
  return { t: tOut, c: cOut, vars, finished: t >= tEnd * (1 - 1e-12), stats: { accepted, rejected, evals } };
}

function luDecompose(A0) {
  const n = A0.length;
  const A = A0.map(r => Float64Array.from(r));
  const p = Array.from({ length: n }, (_, i) => i);
  for (let k = 0; k < n; k++) {
    let m = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[m][k])) m = i;
    if (Math.abs(A[m][k]) < 1e-300) return null;
    if (m !== k) { [A[k], A[m]] = [A[m], A[k]]; [p[k], p[m]] = [p[m], p[k]]; }
    for (let i = k + 1; i < n; i++) {
      const f = A[i][k] / A[k][k];
      A[i][k] = f;
      if (f !== 0) for (let j = k + 1; j < n; j++) A[i][j] -= f * A[k][j];
    }
  }
  return { A, p };
}
function luSolve({ A, p }, b) {
  const n = A.length;
  const x = Array.from(p, i => b[i]);
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) x[i] -= A[i][j] * x[j];
  for (let i = n - 1; i >= 0; i--) { for (let j = i + 1; j < n; j++) x[i] -= A[i][j] * x[j]; x[i] /= A[i][i]; }
  for (let i = 0; i < n; i++) b[i] = x[i];
}

// ---------------------------------------------------------------------------
// Analisi
// ---------------------------------------------------------------------------

/** Regressione lineare y = a + b x con coefficiente di determinazione R². */
export function linfit(xs, ys) {
  const n = xs.length;
  if (n < 3) return null;
  let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; syy += ys[i] * ys[i]; }
  const d = n * sxx - sx * sx;
  if (d === 0) return null;
  const b = (n * sxy - sx * sy) / d, a = (sy - b * sx) / n;
  const ssTot = syy - sy * sy / n;
  let ssRes = 0;
  for (let i = 0; i < n; i++) ssRes += (ys[i] - a - b * xs[i]) ** 2;
  return { a, b, r2: ssTot > 0 ? 1 - ssRes / ssTot : 1 };
}

/**
 * Leggi integrate per una specie che si consuma: [A] (ordine 0), ln[A] (ordine 1), 1/[A] (ordine 2) contro t,
 * sui dati fino al 90 % di conversione. La legge con R² più vicino a 1 indica l'ordine apparente.
 */
export function orderAnalysis(t, c) {
  const c0 = c[0];
  if (!(c0 > 0)) return null;
  const xs = [], ys = [];
  for (let i = 0; i < t.length; i++) { if (c[i] < 0.1 * c0) break; if (c[i] > 0) { xs.push(t[i]); ys.push(c[i]); } }
  if (xs.length < 5 || ys[ys.length - 1] > 0.98 * c0) return null;
  const fits = [
    { order: 0, label: '[A] = [A]₀ − k t', fit: linfit(xs, ys), k: (f) => -f.b },
    { order: 1, label: 'ln[A] = ln[A]₀ − k t', fit: linfit(xs, ys.map(Math.log)), k: (f) => -f.b },
    { order: 2, label: '1/[A] = 1/[A]₀ + k t', fit: linfit(xs, ys.map(v => 1 / v)), k: (f) => f.b },
  ].filter(x => x.fit);
  for (const x of fits) x.kobs = x.k(x.fit);
  const best = fits.reduce((a, b) => (1 - b.fit.r2 < 1 - a.fit.r2 ? b : a));
  // tempo di dimezzamento misurato (interpolazione lineare)
  let t12 = null;
  for (let i = 1; i < t.length; i++) if (c[i] <= c0 / 2 && c[i - 1] > c0 / 2) { t12 = t[i - 1] + (t[i] - t[i - 1]) * (c[i - 1] - c0 / 2) / (c[i - 1] - c[i]); break; }
  return { fits, best, t12 };
}

/**
 * Parametri di attivazione dalla teoria dello stato di transizione (equazione di Eyring):
 * k = (k_B T / h) exp(ΔS‡/R) exp(−ΔH‡/RT), con ΔH‡ = Ea − RT per una reazione unimolecolare (Atkins).
 * Vale per costanti del primo ordine in s⁻¹.
 */
export function eyring(spec, T) {
  if (!spec || spec.kind !== 'arr' || spec.n) return null;
  const Ea = spec.EaR * R_GAS;
  const dH = Ea - R_GAS * T;
  const dS = R_GAS * (Math.log(spec.A * H_PLANCK / (KB * T)) - 1);
  return { dH, dS, dG: dH - T * dS };
}

/** Bilancio di elementi e cariche per le specie che sono formule chimiche (null se non verificabile). */
export function balance(r) {
  const count = (side) => {
    const tot = {};
    for (const [sp, v] of side) {
      const f = parseFormula(sp);
      if (!f) return null;
      for (const [e, k] of Object.entries(f)) tot[e] = (tot[e] ?? 0) + v * k;
    }
    return tot;
  };
  const a = count(r.reactants), b = count(r.products);
  if (!a || !b) return null;
  const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
  const bad = [...keys].filter(k => Math.abs((a[k] ?? 0) - (b[k] ?? 0)) > 1e-9);
  return { ok: bad.length === 0, bad };
}

const ELEMENTS = new Set('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu'.split(' '));
/** Formula → conteggio atomi (carica come elemento fittizio "q"); null se non è una formula. Gli isotopi (Po210) non sono formule. */
export function parseFormula(sp) {
  let s = sp.replace(/−/g, '-');
  let q = 0;
  const cm = s.match(/(\d*)([+-])$/);
  if (cm) { q = (cm[2] === '+' ? 1 : -1) * (cm[1] ? +cm[1] : 1); s = s.slice(0, -cm[0].length); }
  const stack = [{}];
  let i = 0;
  while (i < s.length) {
    const ch = s[i];
    if (ch === '(' || ch === '[') { stack.push({}); i++; continue; }
    if (ch === ')' || ch === ']') {
      i++;
      let d = ''; while (i < s.length && /\d/.test(s[i])) d += s[i++];
      const g = stack.pop(); if (!stack.length) return null;
      const top = stack[stack.length - 1];
      for (const [e, k] of Object.entries(g)) top[e] = (top[e] ?? 0) + k * (d ? +d : 1);
      continue;
    }
    const m = s.slice(i).match(/^([A-Z][a-z]?)(\d*)/);
    if (!m || !ELEMENTS.has(m[1])) return null;
    const top = stack[stack.length - 1];
    top[m[1]] = (top[m[1]] ?? 0) + (m[2] ? +m[2] : 1);
    i += m[0].length;
  }
  if (stack.length !== 1) return null;
  const out = stack[0];
  if (!Object.keys(out).length) return null;
  if (q) out.q = q;
  return out;
}

// ---------------------------------------------------------------------------
// Esperienze con costanti di letteratura
// ---------------------------------------------------------------------------

// Atmosfera standard USA 1976 a 30 km: T = 226,51 K, p = 1197 Pa → n = p/(k_B T)
const N30 = 1197.0 / (KB * 226.51) / 1e6; // molecole/cm³

export const KINETICS_PRESETS = [
  {
    id: 'n2o5', name: 'Decomposizione di N₂O₅ (primo ordine)',
    mech: `# Atkins, Physical Chemistry, tabella 22.4 (dati di Nicholas; Frost e Pearson)
# 2 N2O5 → 4 NO2 + O2 con −d[N2O5]/dt = k[N2O5]
N2O5 -> 2 NO2 + 1/2 O2 ; A=4.94e13 Ea=103`,
    c0: { N2O5: 0.01 }, T: 298.15, tEnd: 3 * 86400, focus: 'N2O5',
    text: 'Reazione del primo ordine: il tempo di dimezzamento t½ = ln 2 / k non dipende dalla concentrazione iniziale. Atkins riporta k(25 °C) = 3,38·10⁻⁵ s⁻¹ (t½ = 5,70 h); i parametri di Arrhenius della stessa tabella danno un valore vicino.',
  },
  {
    id: 'cyclopropane', name: 'Isomerizzazione del ciclopropano',
    mech: `# Atkins, tabella 22.4: A = 1,58·10¹⁵ s⁻¹, Ea = 272 kJ/mol (limite ad alta pressione)
C3H6cyclo -> C3H6 ; A=1.58e15 Ea=272`,
    c0: { C3H6cyclo: 0.01 }, T: 773.15, tEnd: 5400, focus: 'C3H6cyclo',
    text: 'Una reazione unimolecolare in fase gas. A 500 °C Atkins (tabella 22.1) misura k = 6,71·10⁻⁴ s⁻¹ (t½ = 17,2 min): confronta con quello che dà Arrhenius. Nel grafico di Arrhenius si vede quanto conta l\'energia di attivazione alta: 50 K in più moltiplicano k per quasi 10.',
  },
  {
    id: 'sn2', name: 'Sostituzione SN2: C₂H₅Br + OH⁻ in acqua',
    mech: `# Atkins, tabella 22.4 (reazioni in soluzione): A = 4,30·10¹¹ dm³ mol⁻¹ s⁻¹, Ea = 89,5 kJ/mol
C2H5Br + OH- -> C2H5OH + Br- ; A=4.30e11 Ea=89.5`,
    c0: { 'C2H5Br': 0.1, 'OH-': 0.1 }, T: 333.15, tEnd: 4 * 86400, focus: 'C2H5Br',
    text: 'Reazione del secondo ordine (primo in ciascun reagente). Con concentrazioni iniziali uguali vale 1/[A] = 1/[A]₀ + k t e il tempo di dimezzamento 1/(k[A]₀) dipende dalla concentrazione. Prova a cambiare [OH⁻]₀: con un grande eccesso di base la cinetica diventa di pseudo-primo ordine.',
  },
  {
    id: 'sn1', name: 'Solvolisi SN1 del cloruro di terz-butile',
    mech: `# Atkins, tabella 22.4: (CH3)3CCl in acqua, A = 7,1·10¹⁶ s⁻¹, Ea = 100 kJ/mol
(CH3)3CCl + H2O -> (CH3)3COH + H+ + Cl- ; A=7.1e16 Ea=100 ord.H2O=0
costanti: H2O`,
    c0: { '(CH3)3CCl': 0.01, H2O: 55.3 }, T: 298.15, tEnd: 600, focus: '(CH3)3CCl',
    text: 'Lo stadio lento è la ionizzazione (CH₃)₃CCl → (CH₃)₃C⁺ + Cl⁻: la velocità non dipende dall\'acqua (ordine zero, ord.H2O=0) anche se l\'acqua è un reagente. Confronta con gli altri solventi della tabella di Atkins (metanolo: A = 2,3·10¹³ s⁻¹, Ea = 107 kJ/mol).',
  },
  {
    id: 'decay', name: 'Decadimento radioattivo in serie: ²¹⁰Bi → ²¹⁰Po → ²⁰⁶Pb',
    mech: `# Tempi di dimezzamento: 210Bi 5,012 d (β⁻), 210Po 138,376 d (α) (NNDC/ENSDF)
Bi-210 -> Po-210 ; t12=5.012 d
Po-210 -> Pb-206 ; t12=138.376 d`,
    c0: { 'Bi-210': 1 }, T: 298.15, tEnd: 400 * 86400, focus: 'Bi-210', relative: true,
    text: 'Due primi ordini in serie. La soluzione esatta (Bateman) è [Po](t) = [Bi]₀ k₁/(k₂−k₁)(e^(−k₁t) − e^(−k₂t)): la simulazione deve coincidere. La costante di decadimento non dipende dalla temperatura.',
  },
  {
    id: 'chapman', name: 'Ozono stratosferico: ciclo di Chapman (30 km)',
    mech: `# Costanti termiche: IUPAC (Atkinson et al., Atmos. Chem. Phys. 4, 1461, 2004)
# Fotolisi: ordini di grandezza a 30 km (Jacob, Introduction to Atmospheric Chemistry, cap. 10)
unità: molecole/cm3
costanti: O2 N2
O2 -> 2 O ; k=1e-11              # O2 + hν (λ < 242 nm)
O + O2 + N2 -> O3 + N2 ; k0=5.6e-34 n=-2.6
O + O2 + O2 -> O3 + O2 ; k0=6.0e-34 n=-2.6
O3 -> O + O2 ; k=1e-3            # O3 + hν
O + O3 -> 2 O2 ; A=8.0e-12 EaR=2060`,
    c0: { O2: 0.2095 * N30, N2: 0.7808 * N30, O: 0, O3: 0 }, T: 226.51, tEnd: 3e7, focus: 'O3', logT: true,
    text: `A 30 km (atmosfera standard: 226,5 K, n = ${(N30 / 1e17).toFixed(2).replace('.', ',')}·10¹⁷ molecole/cm³) l\'ossigeno atomico vive centesimi di secondo, l\'ozono mesi: un sistema rigido che richiede un integratore implicito. Lo stato stazionario vale [O₃] = [O₂] √(J₁ k₂ [M] / (J₃ k₄)). Il meccanismo di Chapman sovrastima l\'ozono osservato (dell\'ordine di 3·10¹² molecole/cm³ a 30 km): mancano i cicli catalitici di NOx, HOx e ClOx.`,
  },
  {
    id: 'oregonator', name: 'Reazione oscillante di Belousov–Zhabotinsky (Oregonator)',
    mech: `# Modello di Field e Noyes (J. Chem. Phys. 60, 1877, 1974), costanti di Field e Försterling (1986) con [H⁺] ≈ 0,8 M
# X = HBrO2, Y = Br-, Z = Ce(IV), A = BrO3- (costante), B = acido malonico (costante)
# Passi globali: H⁺ e H₂O sono inclusi nelle costanti, quindi i bilanci di H e O non tornano
costanti: BrO3- MA
BrO3- + Br- -> HBrO2 + HOBr ; k=1.28
HBrO2 + Br- -> 2 HOBr ; k=2.4e6
BrO3- + HBrO2 -> 2 HBrO2 + 2 Ce4+ ; k=33.6
2 HBrO2 -> BrO3- + HOBr ; k=3e3
MA + Ce4+ -> 1/2 Br- ; k=1          # passo globale: f = 1 (parametro fenomenologico)`,
    c0: { 'BrO3-': 0.06, MA: 0.02, 'HBrO2': 0, 'Br-': 1e-6, 'Ce4+': 0, HOBr: 0 }, T: 298.15, tEnd: 2000, focus: 'Ce4+', logY: true,
    text: 'Un ciclo di retroazione: HBrO₂ si autocatalizza (terzo passo) finché lo ione bromuro, rigenerato dall\'acido malonico, lo spegne. Ne nasce un ciclo limite: Ce(IV) e Br⁻ oscillano per ordini di grandezza. Il quinto passo riassume la chimica organica ed è fenomenologico (f e k_c sono parametri del modello, non costanti misurate).',
  },
  {
    id: 'mm', name: 'Enzima: meccanismo di Michaelis–Menten (modello)',
    mech: `# Modello con parametri da scegliere: nessun enzima reale
E + S <=> ES ; k=1e6 kr=1e2
ES -> E + P ; k=10`,
    c0: { E: 1e-6, S: 1e-3, ES: 0, P: 0 }, T: 298.15, tEnd: 300, focus: 'S',
    text: 'Per [E]₀ ≪ [S]₀ il complesso ES raggiunge presto uno stato quasi stazionario e v = k_cat[E]₀[S]/(K_M + [S]) con K_M = (k₋₁ + k_cat)/k₁. I numeri sono un modello da modificare, non i dati di un enzima reale.',
  },
  {
    id: 'empty', name: 'Meccanismo libero',
    mech: `# Scrivi qui le tue reazioni, una per riga. Esempi di sintassi:
# A + B -> C ; k=0.5                  costante fissa
# A -> B ; A=1e13 Ea=100              Arrhenius (Ea in kJ/mol)
# A <=> B ; k=2 kr=1                  reversibile
# X -> Y ; t12=3.5 h                  tempo di dimezzamento
# 2 A -> B ; k=1 ord.A=1              ordine diverso dalla stechiometria
# costanti: M                         specie a concentrazione fissa
A + B -> C ; k=0.5`,
    c0: { A: 1, B: 1, C: 0 }, T: 298.15, tEnd: 10, focus: 'A',
    text: 'Qualunque meccanismo: la velocità di ogni passo segue la legge di azione di massa e la sua costante la legge di Arrhenius, se dai A ed Ea.',
  },
];
