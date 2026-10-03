// Largest-remainder sampling: glyphs are statistical species representatives, not trajectories.
export function sampleSpecies(species, budget=120) {
  const rows=species.filter(s=>s.c>0 && Number.isFinite(s.c));
  const total=rows.reduce((n,s)=>n+s.c,0);if(!total)return [];
  const out=rows.map(s=>({...s,exact:s.c/total*budget,count:Math.floor(s.c/total*budget)}));
  let left=budget-out.reduce((n,s)=>n+s.count,0);
  for(const row of [...out].sort((a,b)=>(b.exact-b.count)-(a.exact-a.count)))if(left-->0)row.count++;
  return out.filter(s=>s.count>0);
}
