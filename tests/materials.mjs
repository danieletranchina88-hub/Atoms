import assert from 'node:assert/strict';
import { PhaseMD, ljPair, RC } from '../js/chem/phaseMD.js';
import { Beaker, equilibrate, phaseByName } from '../js/chem/aqueous.js';
import { CATIONS, ANIONS, saltRecipe, solutionOf, quickShelf, beakerPresets, WATER, METAL_SOLIDS, metalRecipe } from '../js/chem/beakerReagents.js';
import { malachiteCrystal, MALACHITE } from '../js/chem/crystals.js';
import { sampleSpecies } from '../js/chem/particleSample.js';
const near=(a,b,tol,msg)=>assert.ok(Math.abs(a-b)<tol,`${msg}: ${a} vs ${b}`);
for(const r of [.97,1.15,1.6,2.49]) {
  const h=1e-6;near(ljPair(r).f,-(ljPair(r+h).u-ljPair(r-h).u)/(2*h),1e-6,'LJ force is -dU/dr');
}
assert.deepEqual(ljPair(RC),{u:0,f:0});near(ljPair(RC-1e-7).f,0,1e-7,'force continuity');
const md=new PhaseMD();md.thermostat=false;const initial=md.stats();md.step(1500);
assert.ok(Math.abs(md.stats().drift/md.n)<1e-4,'NVE energy drift per atom');
for(let a=0;a<3;a++)near(md.v.filter((_,i)=>i%3===a).reduce((s,x)=>s+x,0),0,1e-10,'momentum');
const energy=md.U;md.x=md.x.map(x=>x+md.L);md.forces();near(md.U,energy,1e-9,'periodic translational invariance');
const bath=new PhaseMD();bath.step(1500);assert.ok(Math.abs(bath.stats().drift/bath.n)<1e-4,'NVT energy minus thermostat heat');
const prior=bath.stats();bath.setDensity(.9);near(bath.stats().E-prior.E,bath.work-prior.work,1e-9,'volume work');
const melt=new PhaseMD({density:.85,temperature:1.5});melt.step(3000);assert.ok(melt.stats().order<.25 && melt.stats().msd>.1,'FCC lattice melts dynamically');
const one=malachiteCrystal([1,1,1]);assert.equal(one.atoms.length,40);
const counts={};for(const a of one.atoms)counts[a.element]=(counts[a.element]??0)+1;
assert.deepEqual(counts,{Cu:8,O:20,C:4,H:8});
const [a,b,c,,beta]=MALACHITE.cell;near(a*b*c*Math.sin(beta*Math.PI/180),364.347,.002,'experimental unit cell volume');
for(let i=0;i<one.atoms.length;i++)for(let j=i+1;j<one.atoms.length;j++)assert.ok(Math.hypot(...one.atoms[i].position.map((x,k)=>x-one.atoms[j].position[k]))>.65,'no duplicated crystallographic sites');
assert.equal(sampleSpecies([{c:1},{c:2}],120).reduce((n,s)=>n+s.count,0),120);
const sol=(cat,an,c)=>solutionOf(saltRecipe(CATIONS.findIndex(x=>x[0]===cat),ANIONS.findIndex(x=>x[0]===an)),c);
const CuCl2=sol('Cu+2','Cl-',.1),NaHCO3=sol('Na+','HCO3-',.5),HCl=sol('H+','Cl-',1),NaOH=sol('Na+','OH-',1),HNO3=sol('H+','NO3-',1),KI=sol('K+','I-',.1);
for(const p of beakerPresets()){const b=new Beaker();if(p.thermostat)b.setThermostat(true,298.15);for(const [r,n]of p.steps)b.add(r,n);if(b.st.eq)assert.ok(b.st.eq.converged,`preset ${p.id}`);}
const bkr=new Beaker();bkr.setThermostat(true,298.15);bkr.add(WATER,80);bkr.add(CuCl2,10);bkr.add(NaHCO3,5);
const eq=bkr.st.eq,n=eq.solids.Malachite;assert.ok(n>.00049 && n<=.0005,'Cu stoichiometry limits malachite');
near(eq.SI.Malachite,0,1e-6,'malachite saturation');
const cu=eq.activity('Cu+2'),carbonate=eq.activity('CO3-2'),h=eq.activity('H+');
near(Math.log10(cu*cu*carbonate/h**2),phaseByName('Malachite').lk,1e-6,'malachite mass action includes H+');
bkr.add(HCl,10);assert.ok(!bkr.st.eq.solids.Malachite,'acid dissolves basic carbonate');
assert.ok(bkr.undo());near(bkr.st.eq.solids.Malachite,n,1e-12,'undo restores exact previous phase quantity');
const copper=new Beaker();copper.setThermostat(true,298.15);copper.add(WATER,90);copper.add(HCl,10);copper.add(metalRecipe(METAL_SOLIDS.find(m=>m.id==='Cu')),.1);
assert.ok(copper.st.eq.solids.Cumetal>.999*.1/63.546,'copper stays metallic in dilute nonoxidizing HCl');assert.ok(copper.st.eq.species['Cu+2']<1e-10,'no macroscopic Cu(II) production');
const before=JSON.stringify(copper.st.totals),logLength=copper.log.length;assert.throws(()=>copper.add(HNO3,1),/NOₓ/);assert.equal(JSON.stringify(copper.st.totals),before);assert.equal(copper.log.length,logLength);
assert.equal(copper.add(NaOH,NaN),null);assert.equal(copper.add(NaOH,Infinity),null);
assert.throws(()=>copper.setThermostat(true,Infinity));assert.equal(JSON.stringify(copper.st.totals),before);
const iodide=new Beaker();iodide.add(CuCl2,1);assert.throws(()=>iodide.add(KI,1),/I₂/);
let accepted=0,declined=0;
for(const r of quickShelf().flatMap(([,r])=>r).filter(r=>r.kind==='solution'))for(const second of [NaHCO3,HCl,NaOH]){
  const b=new Beaker();b.setThermostat(true,298.15);b.add(WATER,100);b.add(r,.5);
  try{b.add(second,.5);assert.ok(b.st.eq.converged);near(b.st.eq.chargeResidual,0,1e-10,'charge conservation');accepted++;}
  catch(e){assert.match(e.message,/non convergente|non modellata|database/);declined++;}
}
assert.throws(()=>equilibrate({},0));
console.log(`✓ Materials: force gradients, energy/heat/work, melting, crystallography, equilibria, rollback; ${accepted} mixtures converged, ${declined} explicitly declined`);
