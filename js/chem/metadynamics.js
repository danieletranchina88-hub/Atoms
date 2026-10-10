// Metadinamica Well-Tempered (Barducci, Bussi, Parrinello - PRL 2008)
// Implementazione rigorosa per l'esplorazione di reazioni chimiche nella sandbox
//
// Variabili Collettive (CV) supportate:
//   - Distanza tra due atomi: s = |R_i - R_j|
//   - Angolo tra tre atomi: s = angle(R_i, R_j, R_k)
//   - Numero di coordinazione: s = sum_j (1 - (r_ij/r0)^n) / (1 - (r_ij/r0)^m)
//
// Il bias si accumula come somma di gaussiane nello spazio delle CV:
//   V(s,t) = sum_{t'<t} W0 * exp(-V(s(t'),t')/(kB*DeltaT)) * exp(-sum_k (s_k - s_k(t'))^2 / (2*sigma_k^2))
//
// La forza di bias è calcolata analiticamente come gradiente:
//   F_bias = -grad_R V(s,t)

export class CollectiveVariable {
  constructor(type, indices, params = {}) {
    this.type = type; // 'distance', 'angle', 'coordination'
    this.indices = indices;
    this.params = params;
    this.sigma = params.sigma || 0.1; // Larghezza gaussiana in unità della CV
  }

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
      const v1 = [pos[3*i]-pos[3*j], pos[3*i+1]-pos[3*j+1], pos[3*i+2]-pos[3*j+2]];
      const v2 = [pos[3*k]-pos[3*j], pos[3*k+1]-pos[3*j+1], pos[3*k+2]-pos[3*j+2]];
      const dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
      const n1 = Math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2);
      const n2 = Math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2);
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
        const x = r/r0;
        cn += (1 - x**n) / (1 - x**m);
      }
      return cn;
    }
    throw new Error(`CV type non supportato: ${this.type}`);
  }

  gradient(Z, pos) {
    // Restituisce gradiente dCV/dR come array di dimensione 3*N
    const N = Z.length;
    const grad = new Float64Array(3 * N);
    
    if (this.type === 'distance') {
      const [i, j] = this.indices;
      const dx = pos[3*i] - pos[3*j];
      const dy = pos[3*i+1] - pos[3*j+1];
      const dz = pos[3*i+2] - pos[3*j+2];
      const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
      if (r < 1e-10) return grad;
      grad[3*i] = dx/r; grad[3*i+1] = dy/r; grad[3*i+2] = dz/r;
      grad[3*j] = -dx/r; grad[3*j+1] = -dy/r; grad[3*j+2] = -dz/r;
    }
    if (this.type === 'angle') {
      const [i, j, k] = this.indices;
      const v1 = [pos[3*i]-pos[3*j], pos[3*i+1]-pos[3*j+1], pos[3*i+2]-pos[3*j+2]];
      const v2 = [pos[3*k]-pos[3*j], pos[3*k+1]-pos[3*j+1], pos[3*k+2]-pos[3*j+2]];
      const n1 = Math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2);
      const n2 = Math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2);
      const dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
      const cos_t = dot/(n1*n2);
      const sin_t = Math.sqrt(Math.max(0, 1 - cos_t*cos_t));
      if (sin_t < 1e-10) return grad;
      
      // d(angle)/dR_i = (v2/(n1*n2) - cos_t*v1/n1^2) / sin_t
      for (let a = 0; a < 3; a++) {
        grad[3*i+a] = (v2[a]/(n1*n2) - cos_t*v1[a]/(n1**2)) / sin_t;
        grad[3*k+a] = (v1[a]/(n1*n2) - cos_t*v2[a]/(n2**2)) / sin_t;
        grad[3*j+a] = -grad[3*i+a] - grad[3*k+a];
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
        const x = r/r0;
        const xmn = x**(m-n);
        const deriv = (n*xmn - m) / ((1 - x**m)**2) * (n * x**(n-1)) / r0;
        for (let a = 0; a < 3; a++) {
          const d = [dx, dy, dz][a];
          grad[3*i+a] += deriv * d/r;
          grad[3*j+a] -= deriv * d/r;
        }
      }
    }
    return grad;
  }
}

export class WellTemperedMetaDynamics {
  constructor(cvs, options = {}) {
    this.cvs = cvs;
    this.W0 = options.W0 || 0.001; // Altezza iniziale gaussiana (Hartree)
    this.sigma = cvs.map(cv => cv.sigma); // Larghezze gaussiane
    this.biasFactor = options.biasFactor || 10; // gamma = T_bias/T
    this.T = options.T || 300; // Temperatura (K)
    this.kB = 3.166808e-6; // Costante di Boltzmann in Hartree/K
    this.DeltaT = (this.biasFactor - 1) * this.T;
    
    this.history = []; // Array di {s: Float64Array, weight: number}
    this.timestep = 0;
    this.depositionStride = options.depositionStride || 100; // Deposita ogni N step
  }

  // Calcola il potenziale di bias dato lo stato corrente
  evaluateBias(s) {
    let V = 0;
    for (const h of this.history) {
      let exp_arg = 0;
      for (let k = 0; k < this.cvs.length; k++) {
        exp_arg -= (s[k] - h.s[k])**2 / (2 * this.sigma[k]**2);
      }
      V += h.weight * Math.exp(exp_arg);
    }
    return V;
  }

  // Calcola gradiente del bias rispetto alle CV: dV/ds_k
  biasGradient(s) {
    const dVds = new Float64Array(this.cvs.length);
    for (const h of this.history) {
      let gauss = 1;
      const diff = new Float64Array(this.cvs.length);
      for (let k = 0; k < this.cvs.length; k++) {
        diff[k] = s[k] - h.s[k];
        gauss *= Math.exp(-diff[k]**2 / (2 * this.sigma[k]**2));
      }
      for (let k = 0; k < this.cvs.length; k++) {
        dVds[k] += h.weight * gauss * (-diff[k] / this.sigma[k]**2);
      }
    }
    return dVds;
  }

  // Calcola forze di bias da aggiungere alle forze atomiche
  computeBiasForces(Z, pos, forces) {
    // Valuta CV correnti
    const s = new Float64Array(this.cvs.length);
    for (let k = 0; k < this.cvs.length; k++) {
      s[k] = this.cvs[k].evaluate(Z, pos);
    }
    
    // Deposita nuova gaussiana se è il momento
    if (this.timestep % this.depositionStride === 0) {
      const V_curr = this.evaluateBias(s);
      const weight = this.W0 * Math.exp(-V_curr / (this.kB * this.DeltaT));
      this.history.push({s: s.slice(), weight});
    }
    
    // Calcola gradiente del bias rispetto alle CV
    const dVds = this.biasGradient(s);
    
    // Propaga alle forze atomiche: F_bias_i = -sum_k (dV/ds_k) * (ds_k/dR_i)
    for (let k = 0; k < this.cvs.length; k++) {
      const grad = this.cvs[k].gradient(Z, pos);
      for (let i = 0; i < 3 * Z.length; i++) {
        forces[i] -= dVds[k] * grad[i];
      }
    }
    
    this.timestep++;
    return this.evaluateBias(s);
  }

  // Ricostruisci superficie di energia libera dalla storia del bias
  // F(s) ≈ -V(s) * (T + DeltaT) / DeltaT + cost
  reconstructFES(gridSize = 50) {
    if (this.cvs.length !== 1) {
      throw new Error('FES reconstruction supportata solo per 1 CV');
    }
    
    const cv = this.cvs[0];
    const s_min = Math.min(...this.history.map(h => h.s[0])) - 3*cv.sigma;
    const s_max = Math.max(...this.history.map(h => h.s[0])) + 3*cv.sigma;
    const ds = (s_max - s_min) / (gridSize - 1);
    
    const grid = new Float64Array(gridSize);
    for (let i = 0; i < gridSize; i++) {
      const s = s_min + i * ds;
      const V = this.evaluateBias([s]);
      grid[i] = -V * (this.T + this.DeltaT) / this.DeltaT;
    }
    
    // Normalizza: F(s_min) = 0
    const F_min = Math.min(...grid);
    for (let i = 0; i < gridSize; i++) grid[i] -= F_min;
    
    return {s: Array.from({length: gridSize}, (_, i) => s_min + i*ds), F: grid};
  }
}
