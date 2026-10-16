// Fresh actual Julia session only: current-state fork, compare, intervene, stop.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/ember/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const assert=require('node:assert/strict'),fs=require('node:fs/promises'),path=require('node:path');
(async()=>{const url=process.argv[2]||'http://127.0.0.1:8810/',out=path.resolve(process.argv[3]||'outputs/native-workbench-fork-review-20261016');await fs.mkdir(out,{recursive:true});const browser=await chromium.launch({headless:true});try{
 const page=await browser.newPage({viewport:{width:1640,height:1420}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.goto(url);await page.waitForFunction(()=>!document.querySelector('#advance').disabled,{},{timeout:600000});
 const read=async()=>await(await page.request.get(new URL('/api/status',url).href)).json();
 const command=async(action,branch,steps=4)=>{await page.selectOption('#branch',branch);if(action==='step')await page.fill('#steps',String(steps));await page.click(`[data-action="${action}"]`);await page.waitForFunction(action=>action==='stop'?document.querySelector('#status').textContent.startsWith('stopped'):!document.querySelector('#advance').disabled&&document.querySelector('#message').textContent==='completed',action,{timeout:120000});return read()};
 const initial=await read();assert.equal(initial.commands_completed,0);
 await command('step','sensing',4);const forked=await command('fork','sensing');
 assert.equal(forked.branches.fork1.step,4);assert.equal(forked.branches.fork1.parent.branch,'sensing');assert.equal(await page.locator('#view-left').inputValue(),'sensing');assert.equal(await page.locator('#view-right').inputValue(),'fork1');
 const parentAdvanced=await command('step','sensing',4);const matched=await command('step','fork1',4);
 const getMeta=async b=>await(await page.request.get(new URL('/'+b.latest_observation.file,url).href)).json();
 const [parentMeta,childMeta]=await Promise.all([getMeta(matched.branches.sensing),getMeta(matched.branches.fork1)]);assert.equal(parentMeta.binary_sha256,childMeta.binary_sha256);
 for(const [n,value]of[['raw','.5'],['bound','.03'],['catalyst','.006']])await page.fill('#scale-'+n,value);
 const changed=await command('receptor_scales','fork1');assert.equal(changed.branches.fork1.material_parameters.material_scales[0],.5);assert.equal(changed.branches.sensing.material_parameters.material_scales[0],4);assert.equal(changed.branches.fork1.step,8);
 for(const name of ['sensing','fork1']){const frame=page.frame({url:new RegExp('/branches/'+name+'/index.html')});await frame.waitForFunction(()=>nativeObservatory?.metadata?.step===8,{},{timeout:20000});await frame.evaluate(()=>window.scrollTo(0,0));}
 await page.screenshot({path:path.join(out,'fork-comparison.png'),fullPage:true});
 const fourth=await command('fork','off');assert.equal(Object.keys(fourth.branches).length,4);assert(await page.locator('[data-action="fork"]').isDisabled());
 const denied=await page.request.post(new URL('/api/control',url).href,{data:{action:'fork',branch:'sensing'}});assert.equal(denied.status(),202);const receipt=await denied.json();await page.waitForFunction(async id=>(await(await fetch('/api/ack')).json()).id===id,receipt.id);const capAck=await(await page.request.get(new URL('/api/ack',url).href)).json();assert.equal(capAck.accepted,false);
 const stopped=await command('stop','sensing');assert.equal(stopped.status,'stopped');assert.deepEqual(errors,[]);
 const result={initial,forked,parentAdvanced,matched,parentMeta,childMeta,changed,fourth,capAck,stopped,browserErrors:errors};await fs.writeFile(path.join(out,'browser-results.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify({status:stopped.status,branches:Object.keys(stopped.branches),parityBinary:parentMeta.binary_sha256,errors}));
 }finally{await browser.close()}})().catch(e=>{console.error(e);process.exitCode=1});
