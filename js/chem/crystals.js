// Experimental neutron diffraction coordinates, NOT generated lattice positions.
// Zigan, Joswig, Schuster & Mason (1977), Z. Kristallogr. 145, 412–426.
// AMCSD 0010795; primary paper https://rruff.info/uploads/ZK145_421.pdf
export const MALACHITE = {
  name: 'Malachite · Cu₂CO₃(OH)₂', spaceGroup: 'P 1 21/a 1',
  cell: [9.502,11.974,3.240,90,98.75,90], Z: 4,
  source: 'https://rruff.info/uploads/ZK145_421.pdf',
  sites: [
    ['Cu',.49814,.28793,.89250],['Cu',.23242,.39331,.38800],
    ['O',.13150,.13646,.34170],['O',.33325,.23591,.45000],['O',.33412,.05622,.63080],
    ['O',.09403,.35155,.91910],['O',.37725,.41615,.85980],['C',.26622,.14075,.47270],
    ['H',.01670,.40480,.84440],['H',.41050,.49180,.82960],
  ],
};
export function malachiteCrystal(repeats=[1,1,2]) {
  if (repeats.length!==3 || repeats.some(n=>!Number.isInteger(n)||n<1||n>4))throw new Error('Ripetizioni non valide.');
  const [a,b,c,,beta]=MALACHITE.cell,angle=beta*Math.PI/180;
  const atoms=[];
  for(let i=0;i<repeats[0];i++)for(let j=0;j<repeats[1];j++)for(let k=0;k<repeats[2];k++){
    for(const [element,x,y,z] of MALACHITE.sites){
      for(const f of [[x,y,z],[.5+x,.5-y,z],[.5-x,.5+y,-z],[-x,-y,-z]]){
        const [u,v,w]=f.map(t=>((t%1)+1)%1);
        atoms.push({element,position:[a*(u+i)+c*Math.cos(angle)*(w+k),b*(v+j),c*Math.sin(angle)*(w+k)],fractional:[u,v,w]});
      }
    }
  }
  const bonds=[];
  for(let i=0;i<atoms.length;i++)for(let j=i+1;j<atoms.length;j++){
    const pair=[atoms[i].element,atoms[j].element].sort().join('-');
    const cutoff={'H-O':1.15,'C-O':1.55,'Cu-O':2.55}[pair];
    if(cutoff && Math.hypot(...atoms[i].position.map((v,a)=>v-atoms[j].position[a]))<cutoff)bonds.push([i,j]);
  }
  return {atoms,bonds,cell:MALACHITE.cell,repeats};
}
