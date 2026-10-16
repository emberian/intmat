// Actual resident Julia integration: owns an isolated Playwright browser only.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/ember/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const assert=require('node:assert/strict'),fs=require('node:fs/promises'),path=require('node:path');
(async()=>{const url=process.argv[2]||'http://127.0.0.1:8807/',out=path.resolve(process.argv[3]||'outputs/native-workbench-ui-review-20261016');await fs.mkdir(out,{recursive:true});
 const deadline=Date.now()+400000;while(true){try{const r=await fetch(new URL('/api/status',url));if(r.ok)break;}catch{}if(Date.now()>deadline)throw Error('Gateway readiness deadline');await new Promise(r=>setTimeout(r,250));}
 const browser=await chromium.launch({headless:true});try{
  const page=await browser.newPage({viewport:{width:1640,height:1420}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.goto(url);await page.locator('#advance').waitFor({state:'visible'});
  await page.waitForFunction(()=>!document.querySelector('#advance').disabled,{},{timeout:600000});
  const frames=[page.frame({url:/branches\/off\/index.html/}),page.frame({url:/branches\/sensing\/index.html/})];
  assert(frames.every(Boolean),'Both embedded native viewers exist');
  for(const frame of frames)await frame.waitForFunction(()=>window.nativeObservatory?.metadata?.population>0,{},{timeout:180000});
  const read=async()=>await(await page.request.get(new URL('/api/status',url).href)).json();
  const command=async(action,branch='both',steps=16)=>{
   await page.selectOption('#branch',branch);if(action==='step')await page.fill('#steps',String(steps));
   await page.click(`[data-action="${action}"]`);
   if(action==='stop')await page.waitForFunction(()=>document.querySelector('#status').textContent.startsWith('stopped'),{},{timeout:180000});
   else await page.waitForFunction(()=>!document.querySelector('#advance').disabled&&document.querySelector('#message').textContent==='completed',{},{timeout:180000});
   return read();
  };
  const initial=await read();assert.equal(initial.status,'ready');assert.equal(initial.branches.off.step,0);assert.equal(initial.branches.sensing.step,0);
  const paired=await command('step');assert.equal(paired.branches.off.step,16);assert.equal(paired.branches.sensing.step,16);
  let receptorIntervention=null;
  if(process.env.WORKBENCH_INTERVENTIONS==='1'){
   for(const [name,value] of [['raw','.5'],['bound','.03'],['catalyst','.006']])await page.fill('#scale-'+name,value);
   const applied=await command('receptor_scales','sensing');
   assert.equal(applied.branches.sensing.step,16);assert.equal(applied.branches.off.step,16);
   assert.equal(applied.branches.sensing.interventions.length,1);
   assert.equal(applied.branches.sensing.material_parameters.material_scales[0],.5);
   const index=await(await page.request.get(new URL('/branches/sensing/observatory.json',url).href)).json();
   const eventFrames=index.frames.slice(-2);assert(eventFrames[0].file.includes('-before/'));assert(eventFrames[1].file.includes('-after/'));
   const metas=await Promise.all(eventFrames.map(async f=>await(await page.request.get(new URL('/branches/sensing/'+f.file,url).href)).json()));
   assert.equal(metas[0].binary_sha256,metas[1].binary_sha256);assert.notEqual(eventFrames[0].source_sha256,eventFrames[1].source_sha256);
   const sense=frames[1];await sense.waitForFunction(()=>nativeObservatory.metadata?.intervention?.phase==='after',{},{timeout:20000});
   await sense.click('#previous');await sense.waitForFunction(()=>nativeObservatory.metadata?.intervention?.phase==='before');
   await sense.click('#next');await sense.waitForFunction(()=>nativeObservatory.metadata?.intervention?.phase==='after');
   await sense.check('#follow');
   for(const frame of frames)await frame.evaluate(()=>window.scrollTo(0,0));
   await page.screenshot({path:path.join(out,'receptor-intervention.png'),fullPage:true});
   receptorIntervention={applied,eventFrames,metas};
  }
  const separate=await command('step','sensing',8);assert.equal(separate.branches.off.step,16);assert.equal(separate.branches.sensing.step,24);
  const saved=await command('checkpoint','sensing');assert.equal(saved.branches.sensing.latest_checkpoint.step,24);
  for(const [i,frame] of frames.entries())await frame.waitForFunction(step=>window.nativeObservatory?.metadata?.step===step,i===0?16:24,{timeout:30000});
  const views=await Promise.all(frames.map(frame=>frame.evaluate(()=>({step:nativeObservatory.metadata.step,population:nativeObservatory.metadata.population,shape:nativeObservatory.metadata.shape,fieldValues:nativeObservatory.arrays.fields.length,hiddenValues:nativeObservatory.arrays.hidden.length,error:document.querySelector('#error').textContent,stamp:document.querySelector('#stamp').textContent}))));
  assert(views.every(v=>v.population===4&&v.fieldValues===8**3*9&&v.hiddenValues===32&&v.error===''));
  const denied=await page.request.get(new URL('/'+saved.branches.sensing.latest_checkpoint.file,url).href);assert.equal(denied.status(),404);
  const oversized=await page.request.post(new URL('/api/control',url).href,{data:{action:'step',branch:'both',steps:257}});assert.equal(oversized.status(),400);
  await page.screenshot({path:path.join(out,'workbench.png'),fullPage:true});
  const stopped=await command('stop');assert.equal(stopped.status,'stopped');assert.equal(stopped.branches.off.step,16);assert.equal(stopped.branches.sensing.step,24);
  assert.deepEqual(errors,[]);
  const result={initial,paired,separate,saved,stopped,views,receptorIntervention,browserErrors:errors,checkpointHttpStatus:denied.status(),oversizedHttpStatus:oversized.status()};
  await fs.writeFile(path.join(out,'browser-results.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify({status:stopped.status,views,errors,output:out}));
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
