// Compare old/new rendered pixels and exercise playback, optical controls and mobile layout.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/ember/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const crypto=require('crypto');
const fs=require('fs'),pathUtil=require('path');
const backupDirectory=process.argv[2];
if(!backupDirectory)throw Error('Usage: node scripts/test_viewer_storage_browser.cjs PRIVATE_BACKUP_DIRECTORY');
const repository=pathUtil.resolve(__dirname,'..');
const names=['electrical','collection-20261015/resonant-reliquary','choir','neural'];
(async()=>{const browser=await chromium.launch({headless:true});try{
for(const [i,name] of names.entries()){
const backup=fs.readdirSync(backupDirectory).find(n=>n.startsWith(i+'-'));
let hashes=[];
for(const path of [backupDirectory+'/'+backup,repository+'/site/collection/'+name+'/viewer.html']){
const page=await browser.newPage({viewport:{width:800,height:600}});let errors=[];page.on('pageerror',e=>errors.push(e.message));const begin=performance.now();await page.goto('file://'+path);await page.waitForFunction(()=>window.choirViewer);await page.evaluate(()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))));
let info=await page.evaluate(()=>{const c=document.querySelector('canvas'),gl=c.getContext('webgl2'),pixels=new Uint8Array(c.width*c.height*4);gl.readPixels(0,0,c.width,c.height,gl.RGBA,gl.UNSIGNED_BYTE,pixels);return {pixels:Array.from(pixels),error:document.querySelector('#error').textContent,gl:gl.getError(),frame:choirViewer.frame};});
let hash=crypto.createHash('sha256').update(Buffer.from(info.pixels)).digest('hex');hashes.push(hash);console.log(name,path.startsWith(backupDirectory)?'old':'gzip',JSON.stringify({bytes:fs.statSync(path).size,ready_ms:Math.round(performance.now()-begin),hash,error:info.error,gl:info.gl,frame:info.frame}));
if(!path.startsWith(backupDirectory)){
await page.locator('#time').evaluate(e=>{e.value=0;e.dispatchEvent(new Event('input'));});if(await page.evaluate(()=>choirViewer.frame)!==0)throw Error('scrub');
await page.locator('#mode').selectOption('1');await page.locator('#slice').evaluate(e=>{e.value='.5';e.dispatchEvent(new Event('input'));});
const b=await page.locator('canvas').boundingBox();await page.mouse.move(b.x+100,b.y+100);await page.mouse.down();await page.mouse.move(b.x+180,b.y+150);await page.mouse.up();await page.mouse.wheel(0,100);
await page.locator('#play').click();await page.waitForFunction(()=>choirViewer.frame>0);await page.locator('#play').click();await page.setViewportSize({width:390,height:844});await page.locator('#mode').selectOption('0');await page.screenshot({path:'/private/tmp/intmat-'+name.replaceAll('/','-')+'-mobile.png',fullPage:true,timeout:120000});
if(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth))throw Error('mobile overflow');if(errors.length)throw Error(errors.join(';'));console.log(name,'desktop scrub/slice/drag/zoom/play + mobile pass');
}await page.close();
}if(hashes[0]!==hashes[1])throw Error(name+' pixel mismatch');
}}finally{await browser.close();}})().catch(e=>{console.error(e);process.exitCode=1;});
