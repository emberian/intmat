#!/usr/bin/env python3
"""Package actual living-construction-v1 captures; no resimulation or frame invention."""
import argparse
import base64
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import struct

FIELDS = ('raw_a', 'raw_b', 'catalyst_a', 'catalyst_b', 'bound_a', 'bound_b', 'free_a', 'free_b')
LEDGERS = ('material', 'potential', 'work', 'heat', 'supplied', 'removed', 'supplied_potential', 'removed_potential', 'material_residual', 'potential_residual')

def require(test, message):
    if not test: raise ValueError(message)

def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)

def public(value):
    if isinstance(value, dict): return {k:public(v) for k,v in value.items()}
    if isinstance(value, list): return [public(v) for v in value]
    if isinstance(value, str) and value.startswith(('/Users/', '/home/', '/tmp/', '/private/', 'file://')):
        return Path(value).name
    return value

def pack(record, source_sha256):
    require(record.get('schema') == 'living-construction-v1', 'Expected living-construction-v1')
    m = record['metadata']; shape=m['shape']; extent=m['extent']; origin=m['origin']
    require(len(shape)==3 and all(type(v) is int and 0<v<=128 for v in shape), 'Invalid grid shape')
    require(len(extent)==len(origin)==3 and all(finite(v) and v>0 for v in extent) and all(finite(v) for v in origin), 'Invalid physical domain')
    require(m['field_order']=='x-fast' and m['material_units']=='amount per voxel', 'Unknown field order or units')
    require(finite(m['dt']) and m['dt']>0, 'Invalid physical step')
    require(all(finite(m[k]) and m[k]>=0 for k in ('horizon','actor_rate','sigma','inheritance_fraction')) and m['inheritance_fraction']<=1, 'Invalid learning parameters')
    require(type(m['birth_transfer_adjusted']) is bool, 'Unknown birth transfer adjustment')
    require(all(len(m[k])==2 and all(finite(v) and v>=0 for v in m[k]) for k in ('source_rates','source_budgets')), 'Invalid finite source parameters')
    frames=record['frames']; require(0<len(frames)<=1024, 'Expected 1–1024 completed captures')
    count=math.prod(shape); maximum=[0.0]*len(FIELDS); binary=bytearray(); output=[]; identities={}; previous=-math.inf
    for index,f in enumerate(frames):
        require(finite(f['time']) and f['time']>=previous, 'Captures must be ordered by actual time'); previous=f['time']
        ids=set()
        for c in f['cells']:
            require(type(c['id']) is int and 0<c['id']<=2**53-1 and c['id'] not in ids, 'Duplicate or invalid cell ID')
            ids.add(c['id']); identity=tuple(c[k] for k in ('parent_id','lineage','generation'))
            require(all(type(v) is int and 0<=v<=2**53-1 for v in identity), 'Invalid ancestry')
            require(c['id'] not in identities or identities[c['id']]==identity, 'A cell ID changed its recorded ancestry')
            identities[c['id']]=identity
            require(len(c['position'])==3 and len(c['stores'])==2, 'Wrong cell vector size')
            require(all(finite(v) for v in c['position']+c['stores']+[c[k] for k in ('radius','mass','energy','learned_l2','learning_work')]), 'Nonfinite cell observation')
            require(c['radius']>=0, 'Negative physical radius')
        cells_by_id={c['id']:c for c in f['cells']};seen_bonds=set()
        for bond in f.get('bonds',[]):
            a,b=bond['a_id'],bond['b_id'];pair=tuple(sorted((a,b)))
            require(type(a) is int and type(b) is int and a in cells_by_id and b in cells_by_id and a!=b and pair not in seen_bonds, 'Invalid or duplicate live bond IDs')
            seen_bonds.add(pair)
            require(finite(bond['rest_length']) and bond['rest_length']>0, 'Invalid physical bond rest length')
            require(bond['a_position']==cells_by_id[a]['position'] and bond['b_position']==cells_by_id[b]['position'], 'Bond endpoint differs from its recorded cell')
        require(all(finite(f['ledger'][key]) for key in LEDGERS), 'Missing/nonfinite ledger observation')
        require(type(f['births']) is int and type(f['deaths']) is int and min(f['births'],f['deaths'])>=0, 'Invalid cumulative births/deaths')
        row={k:public(v) for k,v in f.items() if k!='material'}; row['field_offset']=len(binary)//4
        for channel,key in enumerate(FIELDS):
            values=f['material'][key]
            require(len(values)==count and all(finite(v) for v in values), f'Invalid {key} field in capture {index}')
            maximum[channel]=max(maximum[channel],max(values))
            try: binary.extend(struct.pack('<'+str(count)+'f', *values))
            except (OverflowError, struct.error) as exc: raise ValueError('Field cannot be represented as Float32') from exc
        output.append(row)
    events=record.get('birth_events',[])
    for event in events:
        require(all(finite(event[k]) for k in ('time','parent_id','child_id','child_usable_before_copy','copied_l1','copy_work','adjusted_pending_count')), 'Incomplete birth event')
        require(type(event['copy_funded']) is bool, 'Birth copy_funded must be Boolean')
        require(all(type(event[k]) is int and 0<event[k]<=2**53-1 for k in ('parent_id','child_id')) and event['parent_id']!=event['child_id'], 'Invalid birth identity')
    payload=dict(schema='living-construction-viewer-v1',source_schema=record['schema'],source_sha256=source_sha256,
        metadata=public(m),frames=output,birth_events=public(events),
        fields=dict(names=FIELDS,encoding='float32-le/gzip/base64',maximum=maximum,bytes=len(binary),sha256=hashlib.sha256(binary).hexdigest(),data=base64.b64encode(gzip.compress(binary,mtime=0)).decode()))
    for key in ('title','description','synthetic','source'):
        if key in record:payload[key]=public(record[key])
    return payload

def validate_packed(payload):
    """Rebuild from the compact repository artifact without trusting its offsets/hash."""
    require(payload.get('schema')=='living-construction-viewer-v1' and payload.get('source_schema')=='living-construction-v1', 'Unknown packed schema')
    require(re.fullmatch('[a-f0-9]{64}',payload['source_sha256']) is not None, 'Invalid source digest')
    field=payload['fields'];require(tuple(field['names'])==FIELDS and field['encoding']=='float32-le/gzip/base64', 'Unsupported packed field encoding')
    shape=payload['metadata']['shape'];require(len(shape)==3 and all(type(v) is int and 0<v<=128 for v in shape), 'Invalid packed shape')
    count=math.prod(shape);size=count*4*len(FIELDS)*len(payload['frames'])
    require(0<size<=128*1024*1024 and field['bytes']==size, 'Invalid packed byte count')
    compressed=base64.b64decode(field['data'],validate=True)
    # Bound decompression to the declared bytes; reject trailing or oversized payloads.
    import zlib
    decoder=zlib.decompressobj(16+zlib.MAX_WBITS);binary=decoder.decompress(compressed,size+1)
    require(len(binary)==size and decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail, 'Invalid compressed field size')
    require(hashlib.sha256(binary).hexdigest()==field['sha256'], 'Packed field digest mismatch')
    frames=[]
    for index,f in enumerate(payload['frames']):
        offset=index*count*len(FIELDS);require(f['field_offset']==offset, 'Invalid frame field offset')
        row={k:v for k,v in f.items() if k!='field_offset'}
        row['material']={key:list(struct.unpack_from('<'+str(count)+'f',binary,4*(offset+channel*count))) for channel,key in enumerate(FIELDS)}
        frames.append(row)
    raw=dict(schema='living-construction-v1',metadata=payload['metadata'],frames=frames,birth_events=payload['birth_events'])
    checked=pack(raw,payload['source_sha256'])
    require(len(field['maximum'])==len(FIELDS) and all(finite(v) and math.isclose(v,w,rel_tol=1e-6,abs_tol=1e-30) for v,w in zip(field['maximum'],checked['fields']['maximum'])), 'Invalid fixed field scales')
    return public(payload)

def build(source, destination, max_bytes=1_000_000, packed_output=None, title=None, description=None):
    raw=source.read_bytes();record=json.loads(raw)
    payload=validate_packed(record) if record.get('schema')=='living-construction-viewer-v1' else pack(record,hashlib.sha256(raw).hexdigest())
    if title is not None:payload['title']=title
    if description is not None:payload['description']=description
    encoded=json.dumps(payload,separators=(',',':'),allow_nan=False).replace('</','<\\/')
    require(not re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:',encoded), 'Private/local reference in public payload')
    template=Path(__file__).resolve().parent.parent/'web/living_construction_viewer.html'
    geometry=template.with_name('living_construction_geometry.js').read_text()
    html=template.read_text().replace('__DATA__',encoded).replace('__GEOMETRY__',geometry)
    size=len(html.encode());require(size<=max_bytes,f'Viewer is {size} bytes, above explicit {max_bytes} byte budget; select fewer actual captures upstream or explicitly change budget')
    destination.mkdir(parents=True,exist_ok=True);(destination/'index.html').write_text(html)
    if packed_output:
        packed_output.parent.mkdir(parents=True,exist_ok=True);packed_output.write_text(json.dumps(payload,separators=(',',':'),allow_nan=False)+'\n')
    print(f'{len(payload["frames"])} actual captures; {size} bytes; {destination}/index.html')
    return payload,size

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('destination',type=Path);p.add_argument('--max-bytes',type=int,default=1_000_000);p.add_argument('--packed-output',type=Path);p.add_argument('--title');p.add_argument('--description');a=p.parse_args();build(a.source,a.destination,a.max_bytes,a.packed_output,a.title,a.description)
if __name__=='__main__':main()
