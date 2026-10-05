// Orbitali molecolari localizzati di Pipek–Mezey (J. Pipek, P. G. Mezey, J. Chem. Phys. 90, 4916, 1989) e analisi
// dell'ibridazione. Gli orbitali canonici occupati sono delocalizzati su tutta la molecola; una rotazione unitaria fra
// di loro non cambia né la densità né l'energia, e quella che massimizza Σ_i Σ_A (Q_A^ii)² (popolazioni di Mulliken)
// produce gli orbitali "chimici": legami σ e π a due centri, doppietti solitari, legami a più centri. Pipek–Mezey,
// a differenza di Foster–Boys, non mescola σ e π (niente "legami a banana").
// Per ogni orbitale e ogni centro si calcola la composizione s/p/d della sua parte su quell'atomo: il rapporto
// p/s è l'esponente n dell'ibrido spⁿ (sp³ ⇒ n = 3).

/**
 * Ct: coefficienti degli orbitali da localizzare, per righe (Ct[k·n + μ]); S: sovrapposizione n×n;
 * aoAtom, aoL: atomo e momento angolare di ogni funzione; pos: coordinate (per classificare σ/π), qualsiasi unità.
 * Restituisce { C (righe localizzate), orbitals: [{ centers, kind, label, hybrids }], sweeps, converged }.
 */
export function pipekMezey(Ct, m, n, S, aoAtom, aoL, nat, pos, { maxSweeps = 60, tol = 1e-8 } = {}) {
  const C = Float64Array.from(Ct.subarray(0, m * n));
  // SC_k = S c_k, ruotato insieme a c_k
  const SC = new Float64Array(m * n);
  for (let k = 0; k < m; k++) for (let i = 0; i < n; i++) {
    let v = 0; const ri = i * n, rk = k * n;
    for (let j = 0; j < n; j++) v += S[ri + j] * C[rk + j];
    SC[rk + i] = v;
  }
  const first = new Int32Array(nat + 1);
  for (let mu = 0; mu < n; mu++) first[aoAtom[mu] + 1] = mu + 1;
  for (let A = 1; A <= nat; A++) first[A] = Math.max(first[A], first[A - 1]);
  // popolazioni Q_A^kk per saltare le coppie lontane (nessun atomo in comune con peso apprezzabile)
  const pops = (k) => {
    const q = new Float64Array(nat), rk = k * n;
    for (let mu = 0; mu < n; mu++) q[aoAtom[mu]] += C[rk + mu] * SC[rk + mu];
    return q;
  };
  let converged = false, sweeps = 0;
  for (sweeps = 1; sweeps <= maxSweeps; sweeps++) {
    let change = 0;
    const P = Array.from({ length: m }, (_, k) => pops(k));
    for (let i = 0; i < m; i++) for (let j = 0; j < i; j++) {
      let share = 0;
      for (let A = 0; A < nat; A++) share += Math.min(Math.abs(P[i][A]), Math.abs(P[j][A]));
      if (share < 1e-3) continue;
      // A_ij e B_ij della rotazione di Jacobi ottimale
      let a = 0, b = 0;
      const ri = i * n, rj = j * n;
      for (let A = 0; A < nat; A++) {
        let qij = 0, qii = 0, qjj = 0;
        for (let mu = first[A]; mu < first[A + 1]; mu++) {
          qij += 0.5 * (C[ri + mu] * SC[rj + mu] + C[rj + mu] * SC[ri + mu]);
          qii += C[ri + mu] * SC[ri + mu];
          qjj += C[rj + mu] * SC[rj + mu];
        }
        const d = qii - qjj;
        a += qij * qij - 0.25 * d * d;
        b += qij * d;
      }
      const h = Math.hypot(a, b);
      if (h + a < tol) continue;
      const g = 0.25 * Math.atan2(b, -a);
      const c = Math.cos(g), s = Math.sin(g);
      for (let mu = 0; mu < n; mu++) {
        const x = C[ri + mu], y = C[rj + mu];
        C[ri + mu] = c * x + s * y; C[rj + mu] = -s * x + c * y;
        const u = SC[ri + mu], v = SC[rj + mu];
        SC[ri + mu] = c * u + s * v; SC[rj + mu] = -s * u + c * v;
      }
      P[i] = pops(i); P[j] = pops(j);
      change = Math.max(change, h + a);
    }
    if (change < tol) { converged = true; break; }
  }
  // analisi: popolazioni per atomo e per momento angolare
  const orbitals = [];
  for (let k = 0; k < m; k++) {
    const rk = k * n;
    const byAtom = new Float64Array(nat), byL = new Float64Array(nat * 3);
    for (let mu = 0; mu < n; mu++) {
      const q = C[rk + mu] * SC[rk + mu];
      byAtom[aoAtom[mu]] += q;
      byL[aoAtom[mu] * 3 + Math.min(aoL[mu], 2)] += q;
    }
    const centers = [];
    for (let A = 0; A < nat; A++) if (byAtom[A] > 0.08) centers.push({ atom: A, pop: byAtom[A] });
    centers.sort((x, y) => y.pop - x.pop);
    // ibrido su ciascun centro: frazioni s, p, d e direzione del vettore p (per σ/π)
    const hybrid = (A) => {
      const s = Math.max(byL[3 * A], 0), p = Math.max(byL[3 * A + 1], 0), d = Math.max(byL[3 * A + 2], 0);
      const tot = s + p + d || 1;
      // componenti p in ordine (y, z, x): direzione dell'ibrido (una sola shell p per atomo nelle basi minime)
      let px = 0, py = 0, pz = 0, cnt = 0;
      for (let mu = first[A]; mu < first[A + 1]; mu++) if (aoL[mu] === 1) {
        const off = cnt++ % 3;
        if (off === 0) py += C[rk + mu]; else if (off === 1) pz += C[rk + mu]; else px += C[rk + mu];
      }
      return { atom: A, s: s / tot, p: p / tot, d: d / tot, n: s > 1e-3 ? p / s : Infinity, dir: [px, py, pz] };
    };
    const top = centers[0]?.pop ?? 0, two = top + (centers[1]?.pop ?? 0);
    let kind, label;
    const hybrids = centers.slice(0, 3).map(c => ({ ...hybrid(c.atom), pop: c.pop }));
    if (top >= 0.85 || centers.length === 1) {
      kind = 'lone'; label = 'doppietto solitario';
    } else if (two >= 0.8) {
      const [A, B] = [centers[0].atom, centers[1].atom];
      const ax = [pos[3 * B] - pos[3 * A], pos[3 * B + 1] - pos[3 * A + 1], pos[3 * B + 2] - pos[3 * A + 2]];
      const al = Math.hypot(...ax) || 1;
      // carattere σ (s più p lungo l'asse) contro carattere π (p perpendicolare all'asse), pesato sui due centri
      let sig = 0, pi = 0;
      for (const h of hybrids.slice(0, 2)) {
        const pl = h.dir[0] * h.dir[0] + h.dir[1] * h.dir[1] + h.dir[2] * h.dir[2];
        const along = pl > 1e-12 ? ((h.dir[0] * ax[0] + h.dir[1] * ax[1] + h.dir[2] * ax[2]) / al) ** 2 / pl : 1;
        sig += h.pop * (h.s + h.d + h.p * along); pi += h.pop * h.p * (1 - along);
      }
      kind = pi > sig ? 'pi' : 'sigma';
      label = kind === 'pi' ? (two < 0.9 ? 'legame π (in parte delocalizzato)' : 'legame π') : 'legame σ';
    } else {
      // orbitale a più centri: π delocalizzato se su ogni centro è quasi tutto p, altrimenti legame a più centri
      const allP = hybrids.every(h => h.s < 0.1 && h.d < 0.2);
      kind = allP ? 'pi' : 'multi';
      label = allP ? `π delocalizzato su ${centers.length} centri` : `legame a ${centers.length} centri`;
    }
    orbitals.push({ index: k, centers, kind, label, hybrids, pops: byAtom });
  }
  return { C, orbitals, sweeps, converged };
}

/** Etichetta dell'ibrido: "s", "p", "sp²·¹"… dal rapporto p/s (e d/s se c'è). */
export function hybridLabel(h) {
  if (h.s > 0.95) return 's';
  if (h.s < 0.04 && h.d < 0.1) return 'p';
  if (h.d > 0.25 && h.s < 0.05) return 'd';
  const n = h.p / Math.max(h.s, 1e-9), dn = h.d / Math.max(h.s, 1e-9);
  const sup = (x) => x.toFixed(1).replace('.', ',');
  return `sp${n < 0.15 ? '' : `^${sup(n)}`}${dn > 0.15 ? `d^${sup(dn)}` : ''}`;
}

/**
 * Ibridazione di un atomo: media pesata degli ibridi dei suoi legami σ e doppietti solitari (non dei π, che usano
 * orbitali p puri). Restituisce { n, label, sigma, lone, pi } con n = p/s medio.
 */
export function atomHybridization(orbitals, A) {
  let s = 0, p = 0, sigma = 0, lone = 0, pi = 0;
  for (const o of orbitals) {
    const h = o.hybrids.find(x => x.atom === A);
    if (!h) continue;
    if (o.kind === 'pi') { pi++; continue; }
    if (o.kind === 'lone' && o.centers[0].atom !== A) continue;
    if (o.kind === 'lone') lone++; else sigma++;
    s += h.s; p += h.p;
  }
  if (!s && !p) return null;
  const n = p / Math.max(s, 1e-9);
  // nomi canonici quando il rapporto è vicino a 1, 2, 3
  const named = n < 0.4 ? 's' : Math.abs(n - 1) < 0.35 ? 'sp' : Math.abs(n - 2) < 0.5 ? 'sp²' : Math.abs(n - 3) < 0.8 ? 'sp³' : `sp^${n.toFixed(1).replace('.', ',')}`;
  return { n, label: named, sigma, lone, pi };
}
