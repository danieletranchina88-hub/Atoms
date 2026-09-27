// Modello completo di un atomo: calcolo autoconsistente + grandezze derivate.

import { runSCF, radialStats, enclosingRadius } from './scf.js';
import {
  groundStateConfiguration, slaterZeff, slaterNStar, hundTerm, subshellLabel,
} from './configuration.js';

export const HARTREE_EV = 27.211386245988;
export const BOHR_PM = 52.917721090;

/** Sottolivello da cui si rimuove l'elettrone nella prima ionizzazione: n massimo, poi l massimo. */
export function ionizedSubshell(config) {
  return config.reduce((a, b) => (b.n > a.n || (b.n === a.n && b.l > a.l) ? b : a));
}

/** Energia di prima ionizzazione come differenza di energie totali (ΔSCF): E(catione) − E(atomo). */
export function ionizationEnergy(scf) {
  if (scf.Z === 1) {
    // Il catione H⁺ non ha elettroni: E = 0.
    return { hartree: -scf.energy.total, eV: -scf.energy.total * HARTREE_EV, converged: true, from: { n: 1, l: 0 } };
  }
  const from = ionizedSubshell(scf.config);
  const cationCfg = scf.config
    .map(s => (s === from ? { ...s, occ: s.occ - 1 } : { ...s }))
    .filter(s => s.occ > 0);
  const cation = runSCF(scf.Z, cationCfg);
  const de = cation.energy.total - scf.energy.total;
  return { hartree: de, eV: de * HARTREE_EV, converged: cation.converged, from };
}

/** Calcola l'atomo con numero atomico Z e restituisce un oggetto serializzabile. */
export function computeAtom(Z) {
  const config = groundStateConfiguration(Z);
  const t0 = Date.now();
  const scf = runSCF(Z, config);
  const ion = ionizationEnergy(scf);
  const elapsed = Date.now() - t0;

  const orbitals = scf.orbitals.map(o => {
    const stats = radialStats(scf.grid, o.u);
    return {
      n: o.n,
      l: o.l,
      label: subshellLabel(o.n, o.l),
      occ: o.occ,
      e: o.e,
      u: o.u,
      rAvg: stats.rAvg,
      rMaxProb: stats.rMaxProb,
      nodes: stats.nodes,
      r90: enclosingRadius(scf.grid, o.u, 0.9),
      zeffSlater: slaterZeff(Z, config, o.n, o.l),
      nStar: slaterNStar(o.n),
    };
  });

  // Raggio che contiene il 90% della carica elettronica totale.
  const { grid, rho } = scf;
  const u2 = new Float64Array(grid.N);
  for (let i = 0; i < grid.N; i++) u2[i] = Math.sqrt(4 * Math.PI * rho[i]) * grid.r[i];
  const r90 = enclosingRadius(grid, u2, 0.9);

  return {
    Z,
    config,
    term: hundTerm(config),
    grid: { Z, N: grid.N, h: grid.h, xmin: grid.xmin },
    orbitals,
    V: scf.V,
    rho,
    energy: scf.energy,
    ionization: ion,
    r90,
    iterations: scf.iterations,
    converged: scf.converged,
    elapsed,
  };
}
