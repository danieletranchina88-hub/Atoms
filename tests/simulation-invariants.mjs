import assert from 'node:assert/strict';
import { Simulation, MV2 } from '../js/chem/md.js';
import { HALF_REACTIONS } from '../js/chem/labData.js';
import { redoxCell, electrodePotential, halfStoichiometry } from '../js/chem/electrochemistry.js';
import { extentState, equilibriumExtent } from '../js/chem/reactionExtent.js';
import { makeGrid, solveRadial } from '../js/physics/scf.js';
import { RadialFunction } from '../js/physics/wavefunction.js';
const close = (a,b,t=1e-9) => assert.ok(Math.abs(a-b)<t, `${a} != ${b}`);
const h = id => HALF_REACTIONS.find(h => h.id === id);
const momentum = sim => [0,1,2].map(c=>sim.Z.reduce((sum,_,i)=>sum+sim.mass[i]*sim.vel[3*i+c],0));

// Energy impulse must add hc/lambda at any initial radial velocity, preserving momentum.
for (const v of [-0.04, 0, 0.04]) {
  const sim = new Simulation({dt:0.1});
  sim.addAtoms([{Z:1,pos:[-0.37,0,0],vel:[v,0.002,0]},{Z:1,pos:[0.37,0,0],vel:[-v,0,0.001]}]);
  sim.forces(); const k=sim.kinetic(), p=momentum(sim);
  const hit=sim.photon(400); assert.ok(hit);
  close(sim.kinetic()-k,1239.841984/400); momentum(sim).forEach((x,i)=>close(x,p[i]));
}

// NVE integration convergence is measured against the same physical trajectory duration.
function drift(dt) {
  const sim=new Simulation({dt,box:100}); sim.thermostat=false;
  sim.addAtoms([{Z:1,pos:[-0.41,0,0]},{Z:1,pos:[0.41,0,0]}]); sim.forces(); sim.resetMeasurements();
  let max=0; for(let t=0;t<20/dt;t++){sim.step();max=Math.max(max,Math.abs(sim.diagnostics().drift));}
  return max;
}
const d1=drift(0.1),d2=drift(0.05);
assert.ok(d2<d1*0.35,`Verlet convergence: ${d1}, ${d2}`);
assert.ok(d2<0.0002);
console.log('NVE: max drift eV',d1,d2);

const sim=new Simulation({dt:0.1,box:10}); sim.thermostat=false;
sim.addAtoms([{Z:18,pos:[4.8,0,0],vel:[0.001,0,0]}]); sim.forces(); sim.resetMeasurements();
sim.changeBox(8); close(sim.diagnostics().drift,0);
sim.setGrab({i:0,target:[4.7,0,0]}); close(sim.diagnostics().drift,0);
sim.setGrab({i:0,target:[4.5,0,0]}); close(sim.diagnostics().drift,0);
sim.setGrab(null); close(sim.diagnostics().drift,0);
sim.editInventory(()=>sim.addAtoms([{Z:18,pos:[-3,0,0]}]));close(sim.diagnostics().drift,0);
sim.vel.fill(0);sim.resetMeasurements();sim.thermostat=true;sim.bussi();close(sim.diagnostics().drift,0);
sim.clear();assert.equal(sim.N,0);assert.equal(sim.totalEnergy(),0);assert.equal(sim.measurePressure(),0);assert.equal(sim.frags.length,0);

// Unsafe steps abort without a hidden thermostat or velocity clipping.
const fast=new Simulation({dt:2});fast.addAtoms([{Z:1,pos:[0,0,0],vel:[1,0,0]}]);fast.forces();
assert.throws(()=>fast.step(),/instabile/);close(fast.vel[0],1);close(fast.pos[0],0);close(fast.time,0);
assert.throws(()=>new Simulation({dt:-1}));
assert.throws(()=>fast.photon(0));

// Seed determinism for reproducible initial conditions and CSVR samples.
const seeded=()=>{const s=new Simulation({seed:42,dt:0.1});s.addAtoms([{Z:18,pos:[0,0,0]}]);s.rethermalize();for(let i=0;i<50;i++)s.step();return [...s.pos,...s.vel];};
assert.deepEqual(seeded(),seeded());

// Nernst: Daniell, pH slope, stoichiometric powers, soluble reduced species, reversal.
close(redoxCell(h('Zn'),h('Cu')).E,1.10);
close(redoxCell(h('Cu'),h('Zn')).E,-1.10);
const slope=8.314462618*298.15/96485.33212*Math.LN10;
close(electrodePotential(h('H2'),{'H⁺':1e-7,'H₂':1}).E,-7*slope);
close(electrodePotential(h('O2'),{'H⁺':1e-7,'O₂':1}).E,1.23-7*slope);
close(electrodePotential(h('MnO4'),{'H⁺':0.1}).E,1.51-8/5*slope);
close(electrodePotential(h('Fe3'),{'Fe³⁺':0.1,'Fe²⁺':1}).E,0.77-slope);
close(electrodePotential(h('Cl2'),{'Cl⁻':0.1,'Cl₂':1}).E,1.36+slope);
const conc=redoxCell(h('Cu'),h('Cu'),{leftActivities:{'Cu²⁺':1},rightActivities:{'Cu²⁺':1e-3}});
assert.ok(conc.reverse);close(conc.E,-1.5*slope);
const cell=redoxCell(h('Zn'),h('MnO4'),{rightActivities:{'H⁺':1e-3}});
close(cell.n,10);close(cell.E,cell.E0-8.314462618*298.15/(cell.n*96485.33212)*cell.lnQ);
assert.throws(()=>electrodePotential(h('Zn'),{'Zn²⁺':0}));
for (const half of HALF_REACTIONS) assert.ok(halfStoichiometry(half).length >= 2);

// Extent: H2 + I2 <=> 2HI with supplied K=4 has xi=0.5 at equimolar feed.
const species=[{id:'H2',nu:-1,n0:1},{id:'I2',nu:-1,n0:1},{id:'HI',nu:2,n0:0}];
const opts={lnK:Math.log(4),T:300,volumeL:10};
const eq=equilibriumExtent(species,opts);close(eq.xi,0.5);assert.ok(eq.resolved);
const at=extentState(species,eq.xi,opts);close(at.dG,0);close(at.amounts[0].n*2+at.amounts[2].n,2);
assert.throws(()=>extentState(species,1.1,opts));
// Compression favors fewer gas particles: N2O4 <=> 2NO2, same supplied K.
const diss=[{nu:-1,n0:1},{nu:2,n0:0}];
assert.ok(equilibriumExtent(diss,{...opts,volumeL:1}).xi < equilibriumExtent(diss,opts).xi);

// Hydrogen exact Coulomb reference for cumulative radial probability.
const grid=makeGrid(1), V=Float64Array.from(grid.r,r=>-1/r);
const radial=solveRadial(grid,V,1,0,undefined,1);
const r=new RadialFunction(grid,radial.u,0);
close(r.probabilityWithin(1),1-5*Math.exp(-2),3e-5);
close(r.probabilityWithin(0),0);close(r.probabilityWithin(1e6),1);
console.log('Simulation invariants: passed (energy, momentum, convergence, state reset, Nernst, extent, radial probability).');
