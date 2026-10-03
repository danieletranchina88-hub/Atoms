import { pairLJ, shiftedEnergy } from '../js/chem/ljNoble.js';

let failures = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
  if (!ok) failures++;
};

const ar = pairLJ(18, 18);
const rMin = ar.sigma * 2 ** (1 / 6);
check('coppia Ar', ar.sigma === 3.405 && ar.epsilon > 0, `σ=${ar.sigma} ε=${ar.epsilon.toExponential(3)} eV`);
check('minimo vicino a σ·2^(1/6)', Math.abs(shiftedEnergy(rMin, ar.sigma, ar.epsilon) + ar.epsilon) / ar.epsilon < 0.08, `U(rMin)=${shiftedEnergy(rMin, ar.sigma, ar.epsilon).toExponential(3)}`);
check('taglio nullo', shiftedEnergy(2.5 * ar.sigma, ar.sigma, ar.epsilon) === 0);
check('mixing Ne–Ar', Math.abs(pairLJ(10, 18).sigma - 0.5 * (2.789 + 3.405)) < 1e-12);

if (failures) process.exit(1);
