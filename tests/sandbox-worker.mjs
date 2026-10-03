import assert from 'node:assert/strict';
import { Worker } from 'node:worker_threads';
import { PRESETS } from '../js/chem/sandboxData.js';
const url=new URL('../js/chem/sandboxWorker.js',import.meta.url).href;
const worker=new Worker(`const {parentPort}=require('node:worker_threads');
globalThis.onmessage=null;globalThis.postMessage=(m)=>parentPort.postMessage(m);
import(${JSON.stringify(url)}).then(()=>parentPort.on('message',data=>globalThis.onmessage({data})));`,{eval:true});
function wait(predicate){return new Promise((resolve,reject)=>{const timer=setTimeout(()=>{cleanup();reject(new Error('worker timeout'));},30000);const msg=m=>{if(m.type==='error'){cleanup();reject(new Error(m.text));}else if(predicate(m)){cleanup();resolve(m);}};const cleanup=()=>{clearTimeout(timer);worker.off('message',msg);};worker.on('message',msg);});}
try {
  await wait(m=>m.type==='frame');
  let next=wait(m=>m.type==='frame'&&m.N===3&&m.stats.forceField==='hf');
  worker.postMessage({type:'preset',preset:PRESETS.find(p=>p.id==='aimd-h3')});
  let frame=await next;assert.equal(frame.stats.hf.converged,true);assert.equal(frame.stats.paused,true);
  next=wait(m=>m.type==='frame'&&m.stats.t>0);worker.postMessage({type:'step',n:10});frame=await next;
  assert.ok(Number.isFinite(frame.stats.diagnostics.drift));assert.ok(Math.abs(frame.stats.diagnostics.drift)<0.02);
  next=wait(m=>m.type==='frame'&&m.N===0);worker.postMessage({type:'clear'});frame=await next;
  assert.equal(frame.stats.Etot,0);assert.equal(frame.stats.nMol,0);assert.equal(frame.stats.P,0);assert.equal(frame.stats.forceField,'reactive');
  next=wait(m=>m.type==='frame'&&m.N===2);worker.postMessage({type:'add',smiles:'[H][H]',count:1});frame=await next;
  assert.ok(frame.bonds.length>0);assert.equal(frame.stats.nMol,1);
  // il passo scelto non viene aumentato da solo: con idrogeno e Δt ≥ 0,5 fs l'energia non si conserva
  next=wait(m=>m.type==='frame'&&m.N===120&&m.stats.t>300);
  worker.postMessage({type:'preset',preset:PRESETS.find(p=>p.id==='water')});worker.postMessage({type:'set',paused:false,thermostat:false});
  frame=await next;assert.ok(frame.stats.dt<=0.2+1e-12,`dt ${frame.stats.dt}`);assert.ok(Math.abs(frame.stats.diagnostics.drift)<0.1,`drift ${frame.stats.diagnostics.drift}`);
  console.log('Sandbox worker: UHF step, clear, immediate molecular census and fixed time step passed.');
} finally {await worker.terminate();}
