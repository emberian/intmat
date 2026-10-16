const assert=require('node:assert/strict');
const C=require('../web/native_cell_compare.js');let checks=0;
const check=f=>{f();checks++};
function frame({time=1,ids=[7,9],parents=[1,2],energy=[1.25,2.5],hidden=[1,10,2,20],omit=[],generation=[1,1]}={}){
 const values={ids:new Int32Array(ids),parents:new Int32Array(parents),lineages:new Int32Array(ids.map(()=>1)),generations:new Int32Array(generation),energy:new Float64Array(energy),hidden:new Float32Array(hidden)};
 const metadata={time,step:10,endian:'little',population:ids.length,arrays:{}},chunks=[];let offset=0;
 for(const[name,data]of Object.entries(values)){if(omit.includes(name))continue;const bytes=Buffer.from(data.buffer);metadata.arrays[name]={offset,bytes:bytes.length,type:data instanceof Float64Array?'f64':data instanceof Int32Array?'i32':'f32',shape:name==='hidden'?[ids.length,ids.length?hidden.length/ids.length:0]:[ids.length]};chunks.push(bytes);offset+=bytes.length}
 const buffer=Buffer.concat(chunks);return C.decode(metadata,buffer.buffer.slice(buffer.byteOffset,buffer.byteOffset+buffer.byteLength));
}
const a=frame(),b=frame({energy:[1.5,2.75],hidden:[4,10,6,20]});
check(()=>assert.deepEqual(C.row(a,'hidden',7),[1,2]));
check(()=>assert.deepEqual(C.row(a,'hidden',9),[10,20]));
check(()=>assert.deepEqual(C.compare(a,b,7).rows.find(r=>r.name==='hidden').delta,[3,4]));
check(()=>assert.equal(C.compare(a,b,7).rows.find(r=>r.name==='energy').delta[0],.25));
check(()=>assert.deepEqual(C.sharedTimes([{time:.1},{time:.1},{time:.2}],[{time:.1},{time:.2+Number.EPSILON}]),[.1]));
check(()=>assert.equal(C.compare(a,frame({time:1+Number.EPSILON}),7).status,'unmatched-time'));
check(()=>assert.equal(C.compare(a,frame({ids:[9,10]}),7).status,'missing-id'));
check(()=>assert.equal(C.compare(a,b,100).status,'absent-id'));
check(()=>assert.equal(C.compare(a,frame({parents:[8,2]}),7).status,'ancestry-mismatch'));
check(()=>assert(C.compare(a,frame({parents:[8,2]}),7).rows.every(r=>r.delta===null)));
check(()=>assert.equal(C.compare(a,frame({omit:['parents']}),7).ancestryKnown,false));
check(()=>assert.equal(C.compare(a,frame({omit:['hidden']}),7).rows.find(r=>r.name==='hidden').delta,null));
check(()=>assert.equal(C.compare(a,frame({hidden:[NaN,1,2,3]}),7).rows.find(r=>r.name==='hidden').delta,null));
const events=[{time:1,event:'before'},{time:1,event:'after'}];
check(()=>assert.equal(C.atTime(events,1).length,2));
check(()=>assert.equal(C.latest(events,1).event,'after'));
check(()=>assert.throws(()=>frame({ids:[7,7]}),/Duplicate/));
check(()=>assert.throws(()=>C.decode({endian:'big'},new ArrayBuffer(0)),/byte order/));
check(()=>assert.throws(()=>C.decode({endian:'little',arrays:{ids:{type:'i32',shape:[2],bytes:8,offset:0}}},new ArrayBuffer(4)),/Invalid cell array/));
check(()=>assert.equal(C.compare(frame({ids:[],parents:[],energy:[],hidden:[],generation:[]}),a,7).status,'missing-id'));
console.log(`${checks} native cell comparison checks passed`);
