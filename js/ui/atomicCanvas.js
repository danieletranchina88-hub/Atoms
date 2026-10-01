// Orthographic projection of actual Cartesian coordinates; no generated trajectories.
const COLORS={H:'#eeeeee',C:'#718096',O:'#eb635a',Cu:'#cf8b59',Ar:'#80bfc1',Na:'#b5a0f4',Cl:'#83c779'};
export function drawAtoms(canvas, atoms, {bonds=[],yaw=.5,pitch=.35,labels=false,slice=false,span=null,unit='Å'}={}) {
  const w=canvas.clientWidth,h=canvas.clientHeight,dpr=Math.min(devicePixelRatio||1,2);
  if(!w||!h)return [];
  if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);}
  const g=canvas.getContext('2d');g.setTransform(dpr,0,0,dpr,0,0);g.clearRect(0,0,w,h);
  if(!atoms.length)return [];
  const low=[0,1,2].map(a=>Math.min(...atoms.map(t=>t.position[a]))),high=[0,1,2].map(a=>Math.max(...atoms.map(t=>t.position[a])));
  const center=low.map((v,a)=>(v+high[a])/2),extent=span??Math.max(...high.map((v,a)=>v-low[a]))+2;
  const scale=Math.min(w,h)*.64/extent;
  const pts=atoms.map((atom,i)=>{
    const [x,y,z]=atom.position.map((v,a)=>v-center[a]);
    const X=x*Math.cos(yaw)+z*Math.sin(yaw),Z=-x*Math.sin(yaw)+z*Math.cos(yaw);
    return {...atom,i,x:w/2+X*scale,y:h/2+(y*Math.cos(pitch)-Z*Math.sin(pitch))*scale,z:y*Math.sin(pitch)+Z*Math.cos(pitch),hidden:slice&&z>0};
  });
  g.strokeStyle='#8596a777';g.lineWidth=2;
  for(const [i,j] of bonds){const a=pts[i],b=pts[j];if(a.hidden||b.hidden)continue;g.beginPath();g.moveTo(a.x,a.y);g.lineTo(b.x,b.y);g.stroke();}
  for(const p of [...pts].sort((a,b)=>a.z-b.z)){
    if(p.hidden)continue;
    const r=Math.max(2,(p.radius??({H:.22,C:.34,O:.38,Cu:.53}[p.element]??.38))*scale);
    p.screenRadius=r;
    const color=p.color??COLORS[p.element]??'#b39add';
    const grad=g.createRadialGradient(p.x-r*.3,p.y-r*.3,r*.05,p.x,p.y,r);grad.addColorStop(0,'#ffffff');grad.addColorStop(.35,color);grad.addColorStop(1,color);
    g.fillStyle=grad;g.beginPath();g.arc(p.x,p.y,r,0,Math.PI*2);g.fill();g.strokeStyle='#15212c66';g.lineWidth=.6;g.stroke();
    if(labels){g.fillStyle=getComputedStyle(document.documentElement).getPropertyValue('--text');g.font='11px sans-serif';g.textAlign='center';g.fillText(p.label??p.element,p.x,p.y-r-3);}
  }
  const bar=extent>15?5:1;g.strokeStyle=getComputedStyle(document.documentElement).getPropertyValue('--text');g.fillStyle=g.strokeStyle;g.lineWidth=2;
  g.beginPath();g.moveTo(20,h-24);g.lineTo(20+bar*scale,h-24);g.stroke();g.font='12px monospace';g.textAlign='left';g.fillText(`${bar} ${unit}`,20,h-32);
  return pts.filter(p=>!p.hidden);
}
export function bindRotation(canvas,state,redraw,onPick=()=>{}) {
  let start=null,previous=null;
  const down=e=>{start=[e.clientX,e.clientY];previous=start;canvas.setPointerCapture(e.pointerId);};
  const move=e=>{if(!previous)return;state.yaw+=(e.clientX-previous[0])*.008;state.pitch+=(e.clientY-previous[1])*.008;previous=[e.clientX,e.clientY];redraw();};
  const up=e=>{if(start&&Math.hypot(e.clientX-start[0],e.clientY-start[1])<5){const r=canvas.getBoundingClientRect();onPick(e.clientX-r.left,e.clientY-r.top);}start=previous=null;};
  canvas.addEventListener('pointerdown',down);canvas.addEventListener('pointermove',move);canvas.addEventListener('pointerup',up);canvas.addEventListener('pointercancel',up);
  return ()=>{canvas.removeEventListener('pointerdown',down);canvas.removeEventListener('pointermove',move);canvas.removeEventListener('pointerup',up);canvas.removeEventListener('pointercancel',up);};
}
