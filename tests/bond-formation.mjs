// "Come nasce un legame": scansioni HF/UHF/FCI delle molecole biatomiche contro valori pubblicati e contro le
// costanti spettroscopiche sperimentali (Huber–Herzberg, NIST).

import { DIATOMICS, experimentalWell, morseCurve, reducedMass, CM1_EV } from '../js/chem/diatomicData.js';
import { makeContext, scanHF, scanGrid, curveSummary, curveMinimum, derivative, virialSplit, geometry, HARTREE_EV, ANGSTROM_TO_BOHR } from '../js/chem/bondFormation.js';
import { runHF } from '../js/chem/hf.js';
import { fci2 } from '../js/chem/bondFormation.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};
const mol = (id) => DIATOMICS.find(d => d.id === id);
const f = (v, d = 4) => v.toFixed(d);

// 1. H₂ in STO-3G a R = 1,4 bohr: Szabo e Ostlund, Modern Quantum Chemistry, §3.5.2 e §4.1 (E_RHF = −1,1167 Eh,
//    E_FCI = E_RHF − 0,0206 = −1,1373 Eh)
{
  const atoms = geometry([1, 1], 1.4);
  const r = runHF(atoms, { basis: 'STO-3G', conv: 1e-12 });
  const fc = fci2({ n: r.n, C: r.Ca, H: r.H, eri: r.eri, Enuc: r.Enuc });
  check('H₂ STO-3G (R = 1,4 bohr): RHF e FCI come Szabo–Ostlund', Math.abs(r.energy + 1.1167) < 2e-4 && Math.abs(fc.energy + 1.1373) < 2e-4,
    `RHF ${f(r.energy, 5)} Eh (libro −1,1167), FCI ${f(fc.energy, 5)} Eh (libro −1,1373), occupazioni naturali ${fc.natOcc.slice(0, 2).map(v => f(v, 4)).join(', ')}`);
}

// 2. Dati sperimentali coerenti con la teoria esatta: De = D0 + G(0) contro Kołos–Wolniewicz/Pachucki (H₂) e la
//    soluzione esatta di H₂⁺ (Bates, Ledsham e Stewart 1953)
{
  const h2 = experimentalWell(mol('H2')).De, h2p = experimentalWell(mol('H2+')).De;
  check('De sperimentale di H₂ e H₂⁺ = valori teorici esatti', Math.abs(h2 - 38292.9 * CM1_EV) < 0.002 && Math.abs(h2p - 0.1026342 * HARTREE_EV) < 0.002,
    `H₂ ${f(h2, 4)} eV (teoria 4,7477), H₂⁺ ${f(h2p, 4)} eV (teoria ${f(0.1026342 * HARTREE_EV, 4)})`);
  // la curva di Morse ha il minimo −De in re e curvatura che restituisce ωe
  const m = mol('N2'), M = morseCurve(m), h = 1e-4;
  const k = (M.V(M.re + h) - 2 * M.V(M.re) + M.V(M.re - h)) / (h * h) / HARTREE_EV / ANGSTROM_TO_BOHR ** 2;
  const we = Math.sqrt(k / (reducedMass(m.Z) * 1822.888486209)) * 219474.6313632;
  check('Morse di N₂: minimo −De in re e ωe dalla curvatura', Math.abs(M.V(M.re) + M.De) < 1e-12 && Math.abs(we - m.exp.we) < 0.5, `ωe = ${f(we, 2)} cm⁻¹ (dato ${m.exp.we})`);
}

// 3. H₂: FCI dissocia negli atomi, RHF no, UHF sì (con ⟨S²⟩ → 1); FCI vicino all'esperimento
{
  const m = mol('H2'), ctx = makeContext(m, '6-31G**');
  const grid = scanGrid(m);
  const sc = (k) => scanHF(ctx, k, { grid });
  const rhf = sc('rhf'), uhf = sc('uhf'), fci = sc('fci');
  const far = (pts) => (pts[0].E - ctx.Eatoms) * HARTREE_EV;
  check('H₂ a 4 Å: FCI e UHF danno due atomi H, RHF resta molto più in alto', Math.abs(far(fci)) < 0.01 && Math.abs(far(uhf)) < 0.01 && far(rhf) > 5 && uhf[0].S2 > 0.95,
    `E − 2E(H): FCI ${f(far(fci), 4)}, UHF ${f(far(uhf), 4)} (⟨S²⟩ = ${f(uhf[0].S2, 3)}), RHF ${f(far(rhf), 2)} eV`);
  const s = curveSummary(m, fci.map(p => p.R), fci.map(p => p.E), ctx.Eatoms), ex = experimentalWell(m);
  check('H₂ FCI/6-31G**: re, De e ωe vicini all\'esperimento', Math.abs(s.re - ex.re) < 0.01 && Math.abs(s.De - ex.De) < 0.25 && Math.abs(s.we - ex.we) < 250,
    `re ${f(s.re, 4)} Å (esp. ${ex.re}), De ${f(s.De, 3)} eV (esp. ${f(ex.De, 3)}), ωe ${f(s.we, 0)} cm⁻¹ (esp. ${ex.we})`);
  // FCI è sempre sotto RHF (principio variazionale) e la correlazione cresce allungando il legame
  const corr = fci.map((p, i) => rhf[i].E - p.E);
  check('FCI sotto RHF a ogni distanza, correlazione crescente con R', corr.every(c => c > 0) && corr[0] > corr[corr.length - 1] * 5,
    `E_RHF − E_FCI da ${f(corr[corr.length - 1] * HARTREE_EV, 3)} eV (R corta) a ${f(corr[0] * HARTREE_EV, 3)} eV (4 Å)`);
  // teorema del viriale per la funzione d'onda RHF nel suo minimo: −⟨V⟩/⟨T⟩ = 2 a meno dell'incompletezza della base
  // (le gaussiane non sono riscalate: qualche millesimo)
  const sr = curveSummary(m, rhf.map(p => p.R), rhf.map(p => p.E), ctx.Eatoms);
  const [atMin] = scanHF(ctx, 'rhf', { grid: [sr.re] });
  check('H₂ RHF: rapporto viriale −V/T ≈ 2 nel minimo', Math.abs(-atMin.V / atMin.T - 2) < 0.006, `−V/T = ${f(-atMin.V / atMin.T, 5)} a ${f(sr.re, 4)} Å`);
  const iMin = rhf.reduce((b, p, i) => (p.E < rhf[b].E ? i : b), 0);
  // composizione: 1σg metà 1s di A e metà 1s di B, legante; 1σu antilegante
  const p = rhf[iMin], s1 = p.mos[0], s2 = p.mos[1];
  const w = (o, atom) => o.comp.filter(c => c.atom === atom && c.label === '1s').reduce((a, c) => a + c.w, 0);
  check('H₂: 1σg = ½(1s_A + 1s_B) legante, 1σu antilegante', s1.label === '1σg' && s1.character === 'legante' && Math.abs(w(s1, 0) - w(s1, 1)) < 1e-6 && w(s1, 0) > 0.45 && s2.label === '1σu' && s2.character === 'antilegante',
    `${s1.label}: ${f(w(s1, 0), 3)} + ${f(w(s1, 1), 3)} da 1s; ${s2.label} ${s2.character}`);
}

// 4. Teorema del viriale applicato a una curva: T = −E − R dE/dR, V = 2E + R dE/dR sommano a E
{
  const M = morseCurve(mol('HCl'));
  const R = Array.from({ length: 300 }, (_, k) => M.re + (k - 10) * 0.005), E = R.map(r => M.V(r) / HARTREE_EV);
  const d = derivative(R, E), tv = virialSplit(R, E, d);
  const iMin = 10;
  const ok = tv.every((x, i) => Math.abs(x.T + x.V - E[i]) < 1e-12) && Math.abs(tv[iMin].T + E[iMin]) < 2e-3 * Math.abs(E[iMin]);
  check('viriale su una curva: T + V = E e, nel minimo, T = −E', ok, `nel minimo T = ${f(tv[iMin].T * HARTREE_EV, 3)} eV, −E = ${f(-E[iMin] * HARTREE_EV, 3)} eV`);
}

// 5. N₂ e F₂ con Hartree–Fock: risultati noti della letteratura
{
  const n2 = mol('N2'), ctx = makeContext(n2, '6-31G**');
  const grid = scanGrid(n2, { near: 12, far: 6 });
  const rhf = scanHF(ctx, 'rhf', { grid });
  const s = curveSummary(n2, rhf.map(p => p.R), rhf.map(p => p.E), ctx.Eatoms), ex = experimentalWell(n2);
  // HF accorcia il legame e lo indebolisce: re ≈ 1,07–1,08 Å, De circa metà del valore vero, ωe troppo alta
  check('N₂ RHF: legame più corto, più debole e più rigido di quello vero', s.re < ex.re && s.re > 1.06 && s.De > 3.5 && s.De < 0.6 * ex.De && s.we > ex.we,
    `re ${f(s.re, 4)} Å (esp. ${ex.re}), De ${f(s.De, 2)} eV (esp. ${f(ex.De, 2)}), ωe ${f(s.we, 0)} (esp. ${ex.we})`);
  const p = rhf.reduce((b, q) => (Math.abs(q.R - ex.re) < Math.abs(b.R - ex.re) ? q : b));
  const occ = p.mos.filter(o => o.occ > 1.5).map(o => o.label);
  check('N₂: configurazione 1σg² 1σu² 2σg² 2σu² 3σg² 1πu⁴, ordine di legame ≈ 3', ['1σg', '1σu', '2σg', '2σu', '3σg', '1πu'].every(l => occ.includes(l)) && occ.filter(l => l === '1πu').length === 2 && p.bondOrder > 2.6,
    `occupati: ${occ.join(' ')}; Mayer ${f(p.bondOrder, 2)}`);
  const f2 = mol('F2'), cf = makeContext(f2, '6-31G**');
  const rf = scanHF(cf, 'rhf', { grid: scanGrid(f2, { near: 10, far: 4 }) });
  const sf = curveSummary(f2, rf.map(q => q.R), rf.map(q => q.E), cf.Eatoms);
  check('F₂ RHF non legato rispetto agli atomi (De < 0)', sf.De < 0, `minimo RHF a ${f(sf.re, 3)} Å, ${f(sf.De, 2)} eV rispetto a 2 F (esperimento +${f(experimentalWell(f2).De, 2)} eV)`);
}

// 6. LiH: legame polare, la coppia sta più sull'idrogeno
{
  const m = mol('LiH'), ctx = makeContext(m, '6-31G**');
  const pts = scanHF(ctx, 'rhf', { grid: [m.exp.re] });
  const p = pts[0], bond = p.mos.find(o => o.label === '2σ');
  const wH = bond.comp.filter(c => c.atom === 1).reduce((a, c) => a + c.w, 0), wLi = bond.comp.filter(c => c.atom === 0).reduce((a, c) => a + c.w, 0);
  check('LiH: orbitale di legame più sull\'idrogeno, carica Li positiva', wH > wLi && p.charges[0] > 0.1, `2σ: ${f(100 * wH, 0)} % su H, ${f(100 * wLi, 0)} % su Li; q(Li) = ${f(p.charges[0], 2)}`);
}

// 7. UHF di singoletto segue gli atomi separati anche per HCl (spin opposti su H e Cl)
{
  const m = mol('HCl'), ctx = makeContext(m, '6-31G**');
  const grid = scanGrid(m, { near: 14, far: 6 });
  const uhf = scanHF(ctx, 'uhf', { grid }), rhf = scanHF(ctx, 'rhf', { grid });
  const farU = (uhf[0].E - ctx.Eatoms) * HARTREE_EV;
  const below = uhf.every((p, i) => p.E <= rhf[i].E + 1e-7);
  const s = curveSummary(m, uhf.map(p => p.R), uhf.map(p => p.E), ctx.Eatoms);
  check('HCl UHF: dissocia in H + Cl, mai sopra RHF, minimo come RHF', Math.abs(farU) < 0.01 && below && Math.abs(s.re - curveMinimum(rhf.map(p => p.R), rhf.map(p => p.E)).re) < 0.01,
    `E(${f(uhf[0].R, 2)} Å) − E(atomi) = ${f(farU, 4)} eV, spin su Cl ${f(uhf[0].spinA, 2)}, re ${f(s.re, 3)} Å`);
}

if (failures) { console.log(`\n${failures} verifiche sulla formazione del legame fallite`); process.exit(1); }
console.log('\nFormazione del legame: tutte le verifiche superate.');
