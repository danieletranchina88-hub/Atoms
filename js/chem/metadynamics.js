// Metadinamica Well-Tempered (Barducci, Bussi, Parrinello - PRL 100, 020603, 2008)
// Implementazione rigorosa per l'esplorazione di reazioni chimiche nella sandbox.
//
// Variabili Collettive (CV) supportate:
//   - Distanza tra due atomi: s = |R_i - R_j|
//   - Angolo tra tre atomi: s = angle(R_i, R_j, R_k)
//   - Numero di coordinazione: s = sum_j (1 - (r_ij/r0)^n) / (1 - (r_ij/r0)^m)
//
// Il bias si accumula come somma di gaussiane nello spazio delle CV:
//   V(s,t) = sum_{t'<t} W(t') * exp(-sum_k (s_k - s_k(t'))^2 / (2*sigma_k^2))
//
// con decadimento well-tempered dell'altezza:
//   W(t') = W0 * exp(-V(s(t'),t') / (kB * DeltaT))
//   DeltaT = (biasFactor - 1) * T
//
// La forza di bias è calcolata analiticamente come gradiente:
//   F_bias_i = -grad_{R_i} V(s,t) = -sum_k (dV/ds_k) * (ds_k/dR_i)
//
// La superficie di energia libera si ricostruisce asintoticamente:
//   F(s) = -(biasFactor / (biasFactor - 1)) * V(s) + cost
//
// Riferimenti:
//   [1] A. Barducci, G. Bussi, M. Parrinello, Phys. Rev. Lett. 100, 020603 (2008)
//   [2] A. Laio, M. Parrinello, Proc. Natl. Acad. Sci. 99, 12562 (2002)
//   [3] G. Bussi, A. Laio, Rev. Mod. Phys. 92, 045006 (2020)

export class CollectiveVariable {
  /**
   * @param {'distance'|'angle'|'coordination'} type - Tipo di variabile collettiva
   * @param {number[]} indices - Indici atomici coinvolti
   * @param {object} params - Parametri aggiuntivi (sigma, r0, n, m)
   */
  constructor(type, indices, params = {}) {
    this.type = type;
    this.indices = indices;
    this.params = params;
    this.sigma = params.sigma || 0.1; // Larghezza gaussiana in unità della CV
  }

  /**
   * Valuta il valore della CV data la configurazione atomica.
   * @param {number[]} Z - Numeri atomici
   * @param {Float64Array} pos - Posizioni atomiche (x0,y0,z0,x1,y1,z1,...)
   * @returns {number} Valore della CV
   */
  evaluate(Z, pos) {
    if (this.type === 'distance') {
      const [i, j] = this.indices;
      const dx = pos[3*i] - pos[3*j];
      const dy = pos[3*i+1] - pos[3*j+1];
      const dz = pos[3*i+2] - pos[3*j+2];
      return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }
    if (this.type === 'angle') {
      const [i, j, k] = this.indices;
      const v1x = pos[3*i]-pos[3*j], v1y = pos[3*i+1]-pos[3*j+1], v1z = pos[3*i+2]-pos[3*j+2];
      const v2x = pos[3*k]-pos[3*j], v2y = pos[3*k+1]-pos[3*j+1], v2z = pos[3*k+2]-pos[3*j+2];
      const dot = v1x*v2x + v1y*v2y + v1z*v2z;
      const n1 = Math.sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
      const n2 = Math.sqrt(v2x*v2x + v2y*v2y + v2z*v2z);
      return Math.acos(Math.max(-1, Math.min(1, dot/(n1*n2))));
    }
    if (this.type === 'coordination') {
      const [i] = this.indices;
      const {r0, n, m} = {r0: 1.8, n: 6, m: 12, ...this.params};
      let cn = 0;
      for (let j = 0; j < Z.length; j++) {
        if (j === i) continue;
        const dx = pos[3*i] - pos[3*j];
        const dy = pos[3*i+1] - pos[3*j+1];
        const dz = pos[3*i+2] - pos[3*j+2];
        const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-10) continue;
        const x = r/r0;
        // Evita divisione per zero quando x^m → 1
        const xn = Math.pow(x, n);
        const xm = Math.pow(x, m);
        const denom = 1 - xm;
        if (Math.abs(denom) < 1e-10) {
          // Limite per x → 1: (1-x^n)/(1-x^m) → n/m
          cn += n / m;
        } else {
          cn += (1 - xn) / denom;
        }
      }
      return cn;
    }
    throw new Error(`CV type non supportato: ${this.type}`);
  }

  /**
   * Calcola il gradiente analitico ds/dR per ogni atomo.
   * @param {number[]} Z - Numeri atomici
   * @param {Float64Array} pos - Posizioni atomiche
   * @returns {Float64Array} Gradiente (dimensione 3*N)
   */
  gradient(Z, pos) {
    const N = Z.length;
    const grad = new Float64Array(3 * N);

    if (this.type === 'distance') {
      const [i, j] = this.indices;
      const dx = pos[3*i] - pos[3*j];
      const dy = pos[3*i+1] - pos[3*j+1];
      const dz = pos[3*i+2] - pos[3*j+2];
      const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
      if (r < 1e-10) return grad;
      const inv_r = 1.0 / r;
      grad[3*i]   =  dx * inv_r;
      grad[3*i+1] =  dy * inv_r;
      grad[3*i+2] =  dz * inv_r;
      grad[3*j]   = -dx * inv_r;
      grad[3*j+1] = -dy * inv_r;
      grad[3*j+2] = -dz * inv_r;
    }

    if (this.type === 'angle') {
      const [i, j, k] = this.indices;
      const v1 = [pos[3*i]-pos[3*j], pos[3*i+1]-pos[3*j+1], pos[3*i+2]-pos[3*j+2]];
      const v2 = [pos[3*k]-pos[3*j], pos[3*k+1]-pos[3*j+1], pos[3*k+2]-pos[3*j+2]];
      const n1 = Math.sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
      const n2 = Math.sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
      
      if (n1 < 1e-10 || n2 < 1e-10) return grad;
      
      const dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
      const cos_t = Math.max(-1, Math.min(1, dot / (n1 * n2)));
      const sin_t = Math.sqrt(Math.max(1e-10, 1 - cos_t*cos_t));

      // Formula corretta per il gradiente dell'angolo
      // θ = arccos(v1·v2 / (|v1||v2|))
      // dθ/dR = -1/sin(θ) * d(cos θ)/dR
      
      const n1_sq = n1 * n1;
      const n2_sq = n2 * n2;
      const n1n2 = n1 * n2;
      const inv_sin = 1.0 / sin_t;

      // d(cos θ)/dR_i = (v2/(n1*n2) - cos(θ)*v1/n1²)
      // d(cos θ)/dR_k = (v1/(n1*n2) - cos(θ)*v2/n2²)
      // d(cos θ)/dR_j = -d(cos θ)/dR_i - d(cos θ)/dR_k
      
      for (let a = 0; a < 3; a++) {
        const dcos_dRi = v2[a] / n1n2 - cos_t * v1[a] / n1_sq;
        const dcos_dRk = v1[a] / n1n2 - cos_t * v2[a] / n2_sq;
        const dcos_dRj = -dcos_dRi - dcos_dRk;
        
        // dθ/dR = -1/sin(θ) * d(cos θ)/dR
        grad[3*i+a] = -inv_sin * dcos_dRi;
        grad[3*k+a] = -inv_sin * dcos_dRk;
        grad[3*j+a] = -inv_sin * dcos_dRj;
      }
    }

    if (this.type === 'coordination') {
      const [i] = this.indices;
      const {r0, n, m} = {r0: 1.8, n: 6, m: 12, ...this.params};
      for (let j = 0; j < N; j++) {
        if (j === i) continue;
        const dx = pos[3*i] - pos[3*j];
        const dy = pos[3*i+1] - pos[3*j+1];
        const dz = pos[3*i+2] - pos[3*j+2];
        const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-10) continue;
        
        const x = r / r0;
        const xn = Math.pow(x, n);
        const xm = Math.pow(x, m);
        const denom = 1 - xm;
        
        // Evita instabilità numerica quando x ≈ 1
        if (Math.abs(denom) < 1e-8) {
          // Per x ≈ 1, la derivata è circa (n-m)/(m*r0)
          const dcn_dr = (n - m) / (m * r0);
          const inv_r = 1.0 / r;
          for (let a = 0; a < 3; a++) {
            const d = [dx, dy, dz][a];
            grad[3*i+a] += dcn_dr * d * inv_r;
            grad[3*j+a] -= dcn_dr * d * inv_r;
          }
        } else {
          // Derivata completa: d/dr [(1-x^n)/(1-x^m)]
          // = (1/r0) * [n*x^(n-1)*(1-x^m) - (1-x^n)*m*x^(m-1)] / (1-x^m)^2
          const numerator = n * Math.pow(x, n-1) * (1 - xm) - (1 - xn) * m * Math.pow(x, m-1);
          const dcn_dr = numerator / (r0 * denom * denom);
          const inv_r = 1.0 / r;
          for (let a = 0; a < 3; a++) {
            const d = [dx, dy, dz][a];
            grad[3*i+a] += dcn_dr * d * inv_r;
            grad[3*j+a] -= dcn_dr * d * inv_r;
          }
        }
      }
    }

    return grad;
  }
}

export class WellTemperedMetaDynamics {
  /**
   * @param {CollectiveVariable[]} cvs - Variabili collettive
   * @param {object} options - Parametri del bias
   * @param {number} options.W0 - Altezza iniziale gaussiana (Hartree)
   * @param {number} options.biasFactor - Fattore gamma = T_bias/T (tipico: 5-20)
   * @param {number} options.T - Temperatura (K)
   * @param {number} options.depositionStride - Deposita ogni N step MD
   */
  constructor(cvs, options = {}) {
    this.cvs = cvs;
    this.W0 = options.W0 || 0.001;
    this.sigma = cvs.map(cv => cv.sigma);
    this.biasFactor = options.biasFactor || 10;
    this.T = options.T || 300;
    this.kB = 3.166808578e-6; // Boltzmann in Hartree/K
    this.DeltaT = (this.biasFactor - 1) * this.T;
    this.depositionStride = options.depositionStride || 100;

    this.history = []; // [{s: Float64Array, weight: number}]
    this.timestep = 0;
    this.totalBias = 0; // Bias accumulato per diagnostica
  }

  /**
   * Valuta il potenziale di bias V(s) dato il vettore CV corrente.
   * @param {Float64Array|number[]} s - Valori correnti delle CV
   * @returns {number} V(s) in Hartree
   */
  evaluateBias(s) {
    let V = 0;
    for (let h = 0; h < this.history.length; h++) {
      const hist = this.history[h];
      let exp_arg = 0;
      for (let k = 0; k < this.cvs.length; k++) {
        const ds = s[k] - hist.s[k];
        exp_arg -= (ds * ds) / (2 * this.sigma[k] * this.sigma[k]);
      }
      V += hist.weight * Math.exp(exp_arg);
    }
    return V;
  }

  /**
   * Calcola il gradiente del bias rispetto alle CV: dV/ds_k.
   * @param {Float64Array|number[]} s - Valori correnti delle CV
   * @returns {Float64Array} dV/ds
   */
  biasGradientCV(s) {
    const dVds = new Float64Array(this.cvs.length);
    for (let h = 0; h < this.history.length; h++) {
      const hist = this.history[h];
      let gauss = hist.weight;
      for (let k = 0; k < this.cvs.length; k++) {
        const ds = s[k] - hist.s[k];
        gauss *= Math.exp(-(ds * ds) / (2 * this.sigma[k] * this.sigma[k]));
      }
      for (let k = 0; k < this.cvs.length; k++) {
        const ds = s[k] - hist.s[k];
        dVds[k] -= gauss * ds / (this.sigma[k] * this.sigma[k]);
      }
    }
    return dVds;
  }

  /**
   * Calcola e aggiunge le forze di bias alle forze atomiche.
   * Questa è la funzione principale da chiamare ad ogni step MD.
   *
   * @param {number[]} Z - Numeri atomici
   * @param {Float64Array} pos - Posizioni atomiche
   * @param {Float64Array} forces - Forze atomiche (modificate in-place)
   * @returns {number} Energia di bias corrente V(s)
   */
  computeBiasForces(Z, pos, forces) {
    // Valuta CV correnti
    const s = new Float64Array(this.cvs.length);
    for (let k = 0; k < this.cvs.length; k++) {
      s[k] = this.cvs[k].evaluate(Z, pos);
    }

    // Deposita nuova gaussiana se è il momento
    if (this.timestep > 0 && this.timestep % this.depositionStride === 0) {
      const V_curr = this.evaluateBias(s);
      const weight = this.W0 * Math.exp(-V_curr / (this.kB * this.DeltaT));
      this.history.push({ s: Float64Array.from(s), weight });
      this.totalBias += weight;
    }

    // Calcola gradiente del bias rispetto alle CV
    const dVds = this.biasGradientCV(s);

    // Propaga alle forze atomiche: F_bias_i = -sum_k (dV/ds_k) * (ds_k/dR_i)
    for (let k = 0; k < this.cvs.length; k++) {
      if (Math.abs(dVds[k]) < 1e-15) continue;
      const grad = this.cvs[k].gradient(Z, pos);
      for (let i = 0; i < forces.length; i++) {
        forces[i] -= dVds[k] * grad[i];
      }
    }

    this.timestep++;
    return this.evaluateBias(s);
  }

  /**
   * Ricostruisci la superficie di energia libera (FES) dalla storia del bias.
   * Formula: F(s) = -(gamma / (gamma - 1)) * V(s) + cost
   * Supportato solo per 1 CV (per semplicità di visualizzazione).
   *
   * @param {number} gridSize - Numero di punti della griglia
   * @param {number[]} range - [s_min, s_max] opzionale
   * @returns {{s: number[], F: number[]}} Griglia FES
   */
  reconstructFES(gridSize = 100, range = null) {
    if (this.cvs.length !== 1) {
      throw new Error('FES 1D supportata solo per 1 CV. Per 2 CV usare reconstructFES2D().');
    }

    if (this.history.length === 0) {
      // Nessuna gaussiana depositata, ritorna FES piatta
      const cv = this.cvs[0];
      const s_min = range ? range[0] : -2;
      const s_max = range ? range[1] : 2;
      const ds = (s_max - s_min) / (gridSize - 1);
      const s_grid = [];
      const F_grid = [];
      for (let i = 0; i < gridSize; i++) {
        s_grid.push(s_min + i * ds);
        F_grid.push(0);
      }
      return { s: s_grid, F: F_grid };
    }

    const cv = this.cvs[0];
    let s_min, s_max;

    if (range) {
      [s_min, s_max] = range;
    } else {
      const values = this.history.map(h => h.s[0]);
      s_min = Math.min(...values) - 3 * cv.sigma;
      s_max = Math.max(...values) + 3 * cv.sigma;
    }

    const ds = (s_max - s_min) / (gridSize - 1);
    const s_grid = [];
    const F_grid = [];
    const prefactor = this.biasFactor / (this.biasFactor - 1);

    for (let i = 0; i < gridSize; i++) {
      const s_val = s_min + i * ds;
      const V = this.evaluateBias([s_val]);
      s_grid.push(s_val);
      F_grid.push(-prefactor * V);
    }

    // Normalizza: F_min = 0
    const F_min = Math.min(...F_grid);
    for (let i = 0; i < gridSize; i++) F_grid[i] -= F_min;

    return { s: s_grid, F: F_grid };
  }

  /**
   * Ricostruisci FES 2D per due variabili collettive.
   * @param {number} gridSize - Punti per asse
   * @returns {{s1: number[], s2: number[], F: number[][]}}
   */
  reconstructFES2D(gridSize = 50) {
    if (this.cvs.length !== 2) {
      throw new Error('FES 2D richiede esattamente 2 CV.');
    }

    if (this.history.length === 0) {
      throw new Error('Nessuna gaussiana depositata per FES 2D.');
    }

    const cv1 = this.cvs[0], cv2 = this.cvs[1];
    const values1 = this.history.map(h => h.s[0]);
    const values2 = this.history.map(h => h.s[1]);

    const s1_min = Math.min(...values1) - 3 * cv1.sigma;
    const s1_max = Math.max(...values1) + 3 * cv1.sigma;
    const s2_min = Math.min(...values2) - 3 * cv2.sigma;
    const s2_max = Math.max(...values2) + 3 * cv2.sigma;

    const ds1 = (s1_max - s1_min) / (gridSize - 1);
    const ds2 = (s2_max - s2_min) / (gridSize - 1);
    const prefactor = this.biasFactor / (this.biasFactor - 1);

    const s1_grid = [], s2_grid = [];
    const F_grid = [];

    for (let i = 0; i < gridSize; i++) s1_grid.push(s1_min + i * ds1);
    for (let j = 0; j < gridSize; j++) s2_grid.push(s2_min + j * ds2);

    let F_min = Infinity;
    for (let i = 0; i < gridSize; i++) {
      F_grid[i] = new Float64Array(gridSize);
      for (let j = 0; j < gridSize; j++) {
        const V = this.evaluateBias([s1_grid[i], s2_grid[j]]);
        F_grid[i][j] = -prefactor * V;
        F_min = Math.min(F_min, F_grid[i][j]);
      }
    }

    // Normalizza
    for (let i = 0; i < gridSize; i++)
      for (let j = 0; j < gridSize; j++)
        F_grid[i][j] -= F_min;

    return { s1: s1_grid, s2: s2_grid, F: F_grid };
  }

  /**
   * Diagnostiche di convergenza.
   * @returns {object} Statistiche del bias
   */
  diagnostics() {
    const nGaussians = this.history.length;
    const recentWeights = this.history.slice(-Math.min(20, nGaussians)).map(h => h.weight);
    const earlyWeights = this.history.slice(0, Math.min(20, nGaussians)).map(h => h.weight);

    const avg = arr => arr.length > 0 ? arr.reduce((a, b) => a + b, 0) / arr.length : 0;

    return {
      nGaussians,
      timestep: this.timestep,
      avgWeightEarly: avg(earlyWeights),
      avgWeightRecent: avg(recentWeights),
      totalBias: this.totalBias,
      convergenceRatio: avg(earlyWeights) > 0 ? avg(recentWeights) / avg(earlyWeights) : 1
    };
  }

  /**
   * Reset completo del bias (per nuova simulazione).
   */
  reset() {
    this.history = [];
    this.timestep = 0;
    this.totalBias = 0;
  }
}

// Utility: stima automatica di sigma per una CV di distanza
export function estimateSigma(cvs, Z, pos, nSamples = 200) {
  // Stima la fluttuazione termica della CV per impostare sigma ~ delta_s
  const samples = cvs.map(() => []);

  for (let i = 0; i < nSamples; i++) {
    const noisy_pos = Float64Array.from(pos);
    for (let k = 0; k < noisy_pos.length; k++) {
      noisy_pos[k] += (Math.random() - 0.5) * 0.02;
    }
    for (let k = 0; k < cvs.length; k++) {
      samples[k].push(cvs[k].evaluate(Z, noisy_pos));
    }
  }

  return samples.map(s => {
    const mean = s.reduce((a, b) => a + b, 0) / s.length;
    const variance = s.reduce((a, b) => a + (b - mean) ** 2, 0) / s.length;
    return Math.max(0.05, Math.sqrt(variance));
  });
}
