// Read-only browser gate: actual frozen histories plus explicit UI edge fixtures.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/ember/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const assert=require('node:assert/strict'),fs=require('node:fs/promises'),path=require('node:path');
(async()=>{const url=process.argv[2]||'http://127.0.0.1:8810/',out=path.resolve(process.argv[3]||'outputs/native-cell-inspector-review-20261016');await fs.mkdir(out,{recursive:true});const browser=await chromium.launch({headless:true});try{
 const page=await browser.newPage({viewport:{width:1640,height:1420}}),errors=[],writes=[];page.on('pageerror',e=>errors.push(e.message));page.on('request',r=>{if(r.method()!=='GET')writes.push(r.url())});
 await page.goto(url);await page.waitForFunction(()=>document.querySelector('#status').textContent.startsWith('stopped'));
 await page.selectOption('#view-left','sensing');await page.selectOption('#view-right','fork1');await page.click('#cell-comparison > summary');
 await page.waitForFunction(()=>document.querySelector('#cell-comparison').dataset.comparisonStatus==='matched');
 const panel=page.locator('#cell-comparison');const text=await panel.innerText();assert(text.includes('Cross-branch identity is unverified'));assert(text.includes('after receptor change #5'));assert(text.includes('controller commands · Δ norm 0.0000'));
 assert.equal(await panel.locator('[data-role="right-phase"] option').count(),3);
 await panel.locator('[data-role="right-phase"]').selectOption('1');await page.waitForFunction(()=>document.querySelector('#cell-comparison [data-role="context"]').textContent.includes('before receptor change #5'));
 await panel.locator('[data-role="right-phase"]').selectOption('2');await panel.locator('[data-role="history"]').check();await panel.locator('svg circle').first().waitFor();assert.equal(await panel.locator('svg circle').count(),4);
 await panel.locator('[data-role="values"] summary').filter({hasText:'hidden'}).click();
 await panel.screenshot({path:path.join(out,'cell-comparison.png')});
 await page.selectOption('#view-right','off');await page.waitForFunction(()=>document.querySelector('#cell-comparison [data-role="status"]').textContent.includes('Latest captures differ'));
 const historical=await panel.locator('[data-role="status"]').textContent();assert(historical.includes('earlier shared capture'));
 await panel.locator('[data-role="id"]').evaluate(el=>{el.add(new Option('ID99 · explicit absent-ID test fixture','99'));el.value='99';el.dispatchEvent(new Event('change'))});await page.waitForFunction(()=>document.querySelector('#cell-comparison').dataset.comparisonStatus==='absent-id');
 // Explicit empty-history fixture checks that no nearest-time pairing occurs.
 await page.route('**/branches/off/observatory.json',r=>r.fulfill({json:{frames:[]}}));await panel.locator('[data-role="refresh"]').click();await page.waitForFunction(()=>document.querySelector('#cell-comparison').dataset.comparisonStatus==='unmatched-time');
 assert.equal(await panel.locator('[data-role="values"]').textContent(),'');assert.deepEqual(errors,[]);assert.deepEqual(writes,[]);
 const result={actual:'sensing/fork1 share exact times .04/.08; same-time before/after phases reachable; current state Δ0 with changed scales',actualMatchedDotCount:4,historical,fixtures:['absent numeric ID99','empty history/no exact match'],browserErrors:errors,writeRequests:writes};await fs.writeFile(path.join(out,'browser-results.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result));
 }finally{await browser.close()}})().catch(e=>{console.error(e);process.exitCode=1});
