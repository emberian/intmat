/* Geometry of recorded observations only; no simulation or temporal interpolation. */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.LivingConstructionGeometry=api;})(typeof globalThis!=='undefined'?globalThis:this,()=>{
  function periodicSegments(a,b,extent,origin=[0,0,0]){
    const delta=b.map((v,d)=>{const x=v-a[d];return x-Math.round(x/extent[d])*extent[d];}),segments=[];
    for(let ox=-1;ox<=1;ox++)for(let oy=-1;oy<=1;oy++)for(let oz=-1;oz<=1;oz++){
      const start=a.map((v,d)=>v+[ox,oy,oz][d]*extent[d]);let lower=0,upper=1;
      for(let d=0;d<3;d++){
        if(Math.abs(delta[d])<1e-14){if(start[d]<origin[d]||start[d]>origin[d]+extent[d])upper=-1;}
        else{const x=(origin[d]-start[d])/delta[d],y=(origin[d]+extent[d]-start[d])/delta[d];lower=Math.max(lower,Math.min(x,y));upper=Math.min(upper,Math.max(x,y));}
      }
      if(lower<upper)segments.push([start.map((v,d)=>v+lower*delta[d]),start.map((v,d)=>v+upper*delta[d])]);
    }
    return segments;
  }
  function positions(frames,id,throughTime=Infinity){return frames.flatMap((f,index)=>{const c=f.cells.find(c=>c.id===id);return c&&f.time<=throughTime?[{time:f.time,index,position:c.position,radius:c.radius}]:[];});}
  function network(frame){
    const adj=new Map(frame.cells.map(c=>[c.id,new Set()]));for(const b of frame.bonds||[]){if(adj.has(b.a_id)&&adj.has(b.b_id)){adj.get(b.a_id).add(b.b_id);adj.get(b.b_id).add(b.a_id);}}
    const remaining=new Set(adj.keys()),sizes=[];while(remaining.size){const seed=remaining.values().next().value,stack=[seed];remaining.delete(seed);let size=0;while(stack.length){const id=stack.pop();size++;for(const next of adj.get(id))if(remaining.delete(next))stack.push(next);}sizes.push(size);}
    return {bonds:(frame.bonds||[]).length,components:sizes.length,largest:Math.max(0,...sizes),degree:id=>adj.get(id)?.size??0};
  }
  return {periodicSegments,positions,network};
});
