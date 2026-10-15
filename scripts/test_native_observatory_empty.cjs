// In-memory empty-population fixture only: no recorded observations are changed.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/ember/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const assert=require('node:assert/strict'),crypto=require('node:crypto');
(async()=>{const browser=await chromium.launch({headless:true});try{
 const page=await browser.newPage();const errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.goto(process.argv[2]||'http://127.0.0.1:8769/outputs/native-garden-observatory-4b747919ec56/');
 await page.waitForFunction(()=>window.nativeObservatory?.metadata?.population>0);
 const original=await page.evaluate(()=>({metadata:nativeObservatory.metadata,entry:nativeObservatory.entry,material:document.querySelector('#material').textContent}));
 const metadata=structuredClone(original.metadata);const response=await page.request.get(new URL(metadata.binary,new URL(original.entry.file,page.url())).href);const binary=await response.body();
 metadata.population=0;metadata.binary='empty-state.bin';metadata.binary_sha256=crypto.createHash('sha256').update(binary).digest('hex');
 for(const [name,spec] of Object.entries(metadata.arrays))if(name!=='fields'){spec.bytes=0;spec.shape[0]=0;spec.offset=0;}
 await page.route('**/empty-frame.json',route=>route.fulfill({json:metadata}));
 await page.route('**/empty-state.bin',route=>route.fulfill({body:binary,contentType:'application/octet-stream'}));
 await page.evaluate(async()=>{document.querySelector('#material').textContent='STALE';await nativeObservatory.load({file:'empty-frame.json',source_sha256:'empty-population-test-fixture',time:0,side:32});});
 const state=await page.evaluate(()=>({population:nativeObservatory.metadata.population,graph:document.querySelector('#graph').textContent,material:document.querySelector('#material').textContent,cell:document.querySelector('#cell').childElementCount,neural:document.querySelector('#neural').childElementCount,bars:document.querySelector('#bars').childElementCount,error:document.querySelector('#error').textContent}));
 assert.equal(state.population,0);assert.equal(state.cell,0);assert.equal(state.neural,0);assert.equal(state.bars,0);assert.equal(state.material,original.material);assert.match(state.graph,/Components0Largest component0Reciprocal bonds0/);assert.equal(state.error,'');assert.deepEqual(errors,[]);console.log('Empty population clears cell/neural/bars and refreshes zero graph + actual field totals.');
}finally{await browser.close();}})().catch(error=>{console.error(error);process.exitCode=1;});
