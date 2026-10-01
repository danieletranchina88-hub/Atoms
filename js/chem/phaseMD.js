// Reduced-unit Lennard–Jones molecular dynamics, NOT a reactive force field.
// U_FS(r)=U_LJ(r)-U_LJ(rc)+(r-rc)F_LJ(rc); F_FS=F_LJ(r)-F_LJ(rc).
// NIST SRS: linear force shift, rc=2.5 sigma. No tail corrections.
export const ARGON = { sigma: 3.405, epsilonK: 119.8, mass: 39.948, tauPs: 2.156 }; // approximate LJ mapping
export const RC = 2.5;
const fc = 24 * (2 / RC ** 13 - 1 / RC ** 7);
const uc = 4 * (RC ** -12 - RC ** -6);
export function ljPair(r) {
  if (r >= RC) return { u: 0, f: 0 };
  if (!(r > 0)) throw new Error('Atomi sovrapposti.');
  const s6 = r ** -6;
  return { u: 4 * (s6 * s6 - s6) - uc + (r - RC) * fc, f: 24 * (2 * s6 * s6 - s6) / r - fc };
}
export class PhaseMD {
  constructor({ density = 0.95, temperature = 0.35, seed = 2026, cells = 4 } = {}) {
    if (!Number.isFinite(density) || density < .02 || density > 1.1 || !Number.isFinite(temperature) || temperature <= 0 || temperature > 5) throw new Error('Densità o temperatura iniziali non valide.');
    if (![3,4,5].includes(cells)) throw new Error('Numero di celle non supportato.');
    this.cells = cells; this.n = 4 * cells ** 3; this.dof = 3 * this.n - 3;
    this.seed = seed >>> 0; this.dt = .002; this.time = 0; this.heat = 0; this.work = 0;
    this.target = temperature; this.thermostat = true;
    this.L = Math.cbrt(this.n / density);
    if (!(this.L > 2 * RC)) throw new Error('La scatola deve superare 2 rc.');
    this.x = new Float64Array(3 * this.n); this.v = new Float64Array(this.x.length); this.f = new Float64Array(this.x.length);
    const basis = [[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]];
    let k = 0;
    for (let i=0;i<cells;i++) for(let j=0;j<cells;j++) for(let l=0;l<cells;l++) for(const b of basis) {
      for(const [axis,t] of [i,j,l].entries()) this.x[k++] = (t + b[axis] + .25) * this.L / cells;
    }
    for(let i=0;i<this.v.length;i++) this.v[i] = this.normal();
    this.removeCOM();
    const scale = Math.sqrt(temperature / this.temperature());
    for(let i=0;i<this.v.length;i++) this.v[i] *= scale;
    this.reference = this.x.slice(); this.forces(); this.initialEnergy = this.U + this.kinetic();
  }
  random() { let t = this.seed = (this.seed + 0x6D2B79F5) >>> 0; t = Math.imul(t ^ t >>> 15, t | 1); t ^= t + Math.imul(t ^ t >>> 7, t | 61); return ((t ^ t >>> 14) >>> 0) / 4294967296; }
  normal() { return Math.sqrt(-2 * Math.log(Math.max(this.random(), 1e-15))) * Math.cos(2 * Math.PI * this.random()); }
  removeCOM() { for(let a=0;a<3;a++){ let m=0;for(let i=a;i<this.v.length;i+=3)m+=this.v[i];m/=this.n;for(let i=a;i<this.v.length;i+=3)this.v[i]-=m; } }
  kinetic() { return this.v.reduce((s,v)=>s+.5*v*v,0); }
  temperature() { return 2*this.kinetic()/this.dof; }
  distance(i,j) { const d=[0,0,0]; for(let a=0;a<3;a++){let v=this.x[3*i+a]-this.x[3*j+a];d[a]=v-this.L*Math.round(v/this.L);}return d; }
  forces() {
    this.f.fill(0); this.U=0; this.virial=0;
    for(let i=0;i<this.n;i++)for(let j=i+1;j<this.n;j++){
      const d=this.distance(i,j),r2=d[0]**2+d[1]**2+d[2]**2;
      if(r2>=RC*RC)continue;
      const r=Math.sqrt(r2),p=ljPair(r);this.U+=p.u;this.virial+=r*p.f;
      for(let a=0;a<3;a++){const f=p.f*d[a]/r;this.f[3*i+a]+=f;this.f[3*j+a]-=f;}
    }
  }
  setDensity(rho) {
    if (!Number.isFinite(rho) || rho < .02 || rho > 1.1) throw new Error('Densità fuori intervallo.');
    const newL=Math.cbrt(this.n/rho);if(newL<=2*RC)throw new Error('L deve superare 2 rc.');
    const oldU=this.U,scale=newL/this.L;this.x=this.x.map(v=>v*scale);this.L=newL;
    this.forces();this.work+=this.U-oldU;this.reference=this.x.slice();
  }
  step(count=1) {
    if (!Number.isFinite(this.target) || this.target < 0 || this.target > 5 || !(this.dt > 0 && this.dt <= .004)) throw new Error('Parametri MD non validi.');
    for(let k=0;k<count;k++){
      // Transactional step: never silently clip forces, speed or displacement.
      const saved={x:this.x.slice(),v:this.v.slice(),heat:this.heat,seed:this.seed};
      try {
        const dt=this.dt;
        for(let i=0;i<this.x.length;i++){this.v[i]+=.5*dt*this.f[i];this.x[i]+=.5*dt*this.v[i];}
        if(this.thermostat){
          const K=this.kinetic(),c=Math.exp(-dt),noise=Math.sqrt(this.target*(1-c*c));
          for(let i=0;i<this.v.length;i++)this.v[i]=c*this.v[i]+noise*this.normal();
          this.removeCOM();this.heat+=this.kinetic()-K;
        }
        for(let i=0;i<this.x.length;i++)this.x[i]+=.5*dt*this.v[i];
        this.forces();
        for(let i=0;i<this.v.length;i++)this.v[i]+=.5*dt*this.f[i];
        if(!Number.isFinite(this.U) || this.x.some((v,i)=>Math.abs(v-saved.x[i])>.15))throw new Error('Passo instabile: riduci Δt o la compressione. Stato ripristinato.');
        this.time+=dt;
      } catch(e){this.x=saved.x;this.v=saved.v;this.heat=saved.heat;this.seed=saved.seed;this.forces();throw e;}
    }
    return this.stats();
  }
  stats() {
    const K=this.kinetic(),T=2*K/this.dof,V=this.L**3;
    const msd=this.x.reduce((s,v,i)=>s+(v-this.reference[i])**2,0)/this.n;
    // Bragg coherence at the initial FCC (200), averaged along three axes.
    let order=0;
    for(let a=0;a<3;a++){let re=0,im=0;for(let i=a;i<this.x.length;i+=3){const q=4*Math.PI*this.cells*this.x[i]/this.L;re+=Math.cos(q);im+=Math.sin(q);}order+=(re*re+im*im)/(3*this.n*this.n);}
    return {T,K,U:this.U,E:K+this.U,pressure:(2*K+this.virial)/(3*V),density:this.n/V,order,msd,time:this.time,heat:this.heat,work:this.work,drift:K+this.U-this.initialEnergy-this.heat-this.work};
  }
  rdf(bins=60) {
    const max=this.L/2,dr=max/bins,counts=new Float64Array(bins);
    for(let i=0;i<this.n;i++)for(let j=i+1;j<this.n;j++){const d=this.distance(i,j),r=Math.hypot(...d);if(r<max)counts[Math.floor(r/dr)]+=2;}
    return {r:Array.from(counts,(_,i)=>(i+.5)*dr),g:Array.from(counts,(c,i)=>c/((this.n*(this.n-1)/this.L**3)*4*Math.PI/3*((i+1)**3-i**3)*dr**3))};
  }
}
