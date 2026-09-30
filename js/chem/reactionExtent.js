// Ideal gas mixture, fixed T and V, one balanced reaction. Standard pressure 1 atm
// matches thermochemistry() in vibrations.js. No rate law is inferred from ΔG°.
const R = 8.314462618; // J mol−1 K−1
export function extentState(species, xi, { T = 298.15, volumeL = 10, lnK = 0 } = {}) {
  if (!Number.isFinite(xi) || !Number.isFinite(T) || T <= 0 || !Number.isFinite(volumeL) || volumeL <= 0 || !Number.isFinite(lnK))
    throw new Error('Stato termodinamico non valido.');
  if (!species.length || species.some(s => !Number.isFinite(s.nu) || !Number.isFinite(s.n0) || s.n0 < 0))
    throw new Error('Quantità o stechiometria non valida.');
  const amounts = species.map(s => ({ ...s, n: s.n0 + s.nu * xi }));
  if (amounts.some(s => s.n < -1e-12)) throw new Error('Avanzamento oltre il reagente limitante.');
  let lnQ = 0;
  for (const s of amounts) {
    s.n = Math.max(0, s.n);
    s.pressurePa = s.n * R * T / (volumeL * 1e-3);
    if (s.nu) lnQ += s.nu * Math.log(s.pressurePa / 101325);
  }
  return { amounts, lnQ, dG: R * T * (lnQ - lnK) / 1000,
    pressureBar: amounts.reduce((sum, s) => sum + s.pressurePa, 0) / 1e5 };
}
export function extentBounds(species) {
  const lo = Math.max(...species.filter(s => s.nu > 0).map(s => -s.n0 / s.nu));
  const hi = Math.min(...species.filter(s => s.nu < 0).map(s => s.n0 / -s.nu));
  if (!Number.isFinite(lo) || !Number.isFinite(hi) || lo > hi) throw new Error('Intervallo di avanzamento non valido.');
  return { lo, hi };
}
export function equilibriumExtent(species, options) {
  let { lo, hi } = extentBounds(species);
  if (hi === lo) return { xi: lo, resolved: false };
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2;
    if (mid === lo || mid === hi) break;
    const state = extentState(species, mid, options);
    if (state.dG > 0) hi = mid; else lo = mid;
  }
  const xi = (lo + hi) / 2, state = extentState(species, xi, options);
  return { xi, resolved: Number.isFinite(state.dG) && Math.abs(state.dG) < 1e-5 };
}
