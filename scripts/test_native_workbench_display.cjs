// Read-only UI regression. Startup is an explicit UI fixture; displayed frames
// come from the actual session behind the gateway, and no command is submitted.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/ember/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const assert=require('node:assert/strict'),fs=require('node:fs/promises'),path=require('node:path');
(async()=>{const url=process.argv[2]||'http://127.0.0.1:8807/',out=path.resolve(process.argv[3]||'outputs/native-workbench-ui-review-20261016');await fs.mkdir(out,{recursive:true});
 const browser=await chromium.launch({headless:true});try{
  const page=await browser.newPage({viewport:{width:1640,height:1420}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.route('**/api/status',route=>route.fulfill({json:{status:'starting',phase:'Loading Julia and preparing paired native states'}}));
  await page.goto(url);await page.waitForFunction(()=>document.querySelector('#status').textContent.startsWith('Initializing Julia'));
  assert(await page.locator('#advance').isDisabled());assert.equal(await page.locator('#clock-left').textContent(),'initializing');
  await page.screenshot({path:path.join(out,'starting-fixture.png'),fullPage:true});
  await page.unroute('**/api/status');await page.waitForFunction(()=>document.querySelector('#status').textContent.startsWith('stopped'));
  const frames=[page.frame({url:/branches\/off\/index.html/}),page.frame({url:/branches\/sensing\/index.html/})];
  for(const frame of frames){await frame.waitForFunction(()=>window.nativeObservatory?.metadata?.population>0);assert(await frame.locator('#status').isVisible());assert.equal(await frame.locator('h1').isVisible(),false);assert.equal(await frame.locator('#error').textContent(),'');}
  await page.screenshot({path:path.join(out,'workbench-compact.png'),fullPage:true});
  assert.deepEqual(errors,[]);console.log(JSON.stringify({startupFixture:'disabled controls, no JS errors',actualSession:'stopped, both native viewers populated',compactHeader:'status visible, hero hidden',errors}));
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
