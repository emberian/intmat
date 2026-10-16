// Exact static recording transport. The publisher supplies decoded-byte digests.
// A host may serve a .gz file directly or decode it via HTTP Content-Encoding.
(()=>{'use strict';
const base=new URL('.',document.currentScript.src),manifest=__RECORDINGS__;
const records=new Map(manifest.recordings.map(r=>[new URL(r.compressed_path,base).href,r]));
window.recordingFetch=async function recordingFetch(input,options){
 const url=new URL(input,document.baseURI),record=records.get(url.href);
 if(!record)throw Error('This recording has no published integrity record.');
 const response=await fetch(url,options);if(!response.ok)return response;
 let bytes=new Uint8Array(await response.arrayBuffer());
 if(bytes[0]===0x1f&&bytes[1]===0x8b){
  if(typeof DecompressionStream==='undefined')throw Error('This recording needs gzip decompression. Please use a current Chrome, Firefox, Safari or Edge browser.');
  bytes=new Uint8Array(await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
 }
 if(bytes.byteLength!==record.original_bytes)throw Error('Recording byte count mismatch. Please reload the page.');
 if(!globalThis.crypto?.subtle)throw Error('Recording integrity checks require a secure browser context (HTTPS).');
 const digest=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),v=>v.toString(16).padStart(2,'0')).join('');
 if(digest!==record.sha256)throw Error('Recording checksum mismatch. Please reload the page.');
 return new Response(bytes,{status:response.status,statusText:response.statusText,headers:{'Content-Type':'application/json'}});
};
})();
