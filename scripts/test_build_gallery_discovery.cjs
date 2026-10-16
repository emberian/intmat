#!/usr/bin/env node
// Real generated-page checks, with dated/undated and missing-viewer fixtures.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {execFileSync} = require('node:child_process');
const {pathToFileURL} = require('node:url');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
(async()=>{
  const directory=fs.mkdtempSync(path.join(os.tmpdir(),'interior-gallery-discovery-'));
  const fixtures=[
    ['collection-20261016/new-light','New light','waves',true],
    ['collection-20261015/old-phase','Old phase','phase',true],
    ['refinement-20261015/refined','Refined','phase',false],
    ['collection-20260230/invalid-date','Invalid date','writing',false],
    ['collection-20261016-study/dated-suffix','Dated suffix','mechanics',false],
    ['legacy','Legacy','ecology',false],
  ];
  for(const [id,title,family,viewer] of fixtures){const p=path.join(directory,id);fs.mkdirSync(p,{recursive:true});fs.writeFileSync(path.join(p,'run.json'),JSON.stringify({title,family,viewer:'viewer.html'}));fs.writeFileSync(path.join(p,'preview.png'),Buffer.from('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+jR1kAAAAASUVORK5CYII=','base64'));if(viewer)fs.writeFileSync(path.join(p,'viewer.html'),'<title>Actual recording</title>');}
  execFileSync('python3',[path.join(__dirname,'build_gallery.py'),directory]);
  const html=fs.readFileSync(path.join(directory,'index.html'),'utf8'),data=JSON.parse(html.match(/<script id="specimens" type="application\/json">(.*?)<\/script>/s)[1]);
  const byId=Object.fromEntries(data.map(s=>[s.id,s]));
  assert.equal(data[0].id,'collection-20261016/new-light');
  assert.equal(byId['collection-20261016-study/dated-suffix'].batch_date,'2026-10-16');
  assert.equal(byId['collection-20260230/invalid-date'].batch_date,null);
  assert.equal(byId.legacy.batch_date,null);
  assert.equal(byId.legacy.viewer,null);
  assert.equal(data.length,6);
  const browser=await chromium.launch({headless:true});
  try{
    const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
    const url=pathToFileURL(path.join(directory,'index.html')).href;
    await page.goto(url+'#legacy');assert.equal(await page.locator('#name').textContent(),'Legacy');
    assert.equal(await page.locator('.tile.active').getAttribute('aria-label'),'Select Legacy');
    assert.equal(await page.locator('.recording').count(),2);
    await page.selectOption('#sort','family');
    const first=await page.locator('.tile').first().getAttribute('aria-label');assert.match(first,/Old phase|Refined/);
    assert.equal(await page.locator('.tile.active').getAttribute('aria-label'),'Select Legacy');
    await page.selectOption('#format','interactive');assert.equal(await page.locator('.tile:visible').count(),2);
    assert.equal(await page.locator('#selection-note').isVisible(),true);
    await page.selectOption('#family','phase');assert.equal(await page.locator('.tile:visible').count(),1);
    await page.locator('.tile:visible').focus();await page.keyboard.press('Enter');
    assert.equal(await page.locator('#name').textContent(),'Old phase');assert.equal(await page.locator('#selection-note').isVisible(),false);
    await page.selectOption('#sort','latest');assert.equal(await page.locator('.tile:visible').count(),1);
    await page.fill('#search','no such art');assert.equal(await page.locator('#no-results').isVisible(),true);assert.equal(await page.locator('.tile:visible').count(),0);
    await page.fill('#search','');await page.selectOption('#family','');
    await page.locator('.tile:visible').first().click();assert.equal(await page.locator('#name').textContent(),'New light');
    await page.goBack();await page.waitForFunction(()=>document.getElementById('name').textContent==='Old phase');
    await page.evaluate(()=>location.hash='collection-20261016%2Fnew-light');await page.waitForFunction(()=>document.getElementById('name').textContent==='New light');
    await page.setViewportSize({width:390,height:844});
    assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),true);
    await page.selectOption('#format','all');assert.equal(await page.locator('.tile:visible').count(),6);
    assert.equal(await page.locator('#batch').textContent(),'Batch 2026-10-16');
    assert.deepEqual(errors,[]);
    let navigationChecks=0;
    if(process.argv[2]){
      await page.goto(pathToFileURL(path.resolve(process.argv[2],'index.html')).href);
      const expected=['animations/','laboratory/','arena/arena.html','metabolic/arena.html','recurrent/','cellular/','constructed/','https://github.com/emberian/intmat','https://github.com/emberian/intmat/blob/dev/docs/alife-adaptation-results.md'];
      assert.deepEqual((await page.locator('.header-actions a').evaluateAll(nodes=>nodes.map(n=>n.getAttribute('href')))).sort(),expected.sort());
      assert.equal(await page.locator('.header-actions>.primary').count(),2);
      assert.equal(await page.locator('.project-menu').getAttribute('open'),null);
      assert.equal(await page.locator('#sort').isVisible(),true);
      assert.equal(await page.locator('#hero').evaluate(n=>n.getBoundingClientRect().top<550),true);
      await page.locator('.project-menu summary').focus();await page.keyboard.press('Enter');
      assert.equal(await page.locator('.project-links a:visible').count(),7);
      assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),true);
      await page.keyboard.press('Enter');assert.equal(await page.locator('.project-menu').getAttribute('open'),null);
      assert.deepEqual(errors,[]);navigationChecks=9;
    }
    console.log(JSON.stringify({checks:24+navigationChecks,fixtures:6,browser_errors:errors,fixture_directory:directory}));
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
