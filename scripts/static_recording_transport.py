"""Losslessly package explicit legacy JSON recordings in a Pages staging tree.

Repository recordings remain unchanged. Only their staged copies become gzip;
all browser references are rewritten before the uncompressed staging copies go.
"""
import gzip,hashlib,json,pathlib,re
RECORDINGS=('arena/replay.json','metabolic/replay.json',
    'cellular/broad-early.json','cellular/broad-mature.json','cellular/broad-history.json',
    'cellular/mosaic-early.json','cellular/mosaic-mature.json','cellular/mosaic-history.json')

def stage_recordings(destination,template=None):
    destination=pathlib.Path(destination)
    template=pathlib.Path(template) if template else pathlib.Path(__file__).resolve().parents[1]/'web/recording_transport.js'
    records=[];prepared=[]
    for relative in RECORDINGS:
        source=destination/relative;raw=source.read_bytes()
        if re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:',raw.decode('utf-8')):raise ValueError(f'Private/local reference in recording: {relative}')
        packed=gzip.compress(raw,compresslevel=9,mtime=0)
        restored=gzip.decompress(packed)
        if restored!=raw:raise ValueError(f'Recording transport changed bytes: {relative}')
        original=json.loads(raw)
        if json.loads(restored)!=original:raise ValueError(f'Recording transport changed state: {relative}')
        frames=original if isinstance(original,list) else original.get('frames',[])
        record=dict(original_path=relative,compressed_path=relative+'.gz',original_bytes=len(raw),compressed_bytes=len(packed),sha256=hashlib.sha256(raw).hexdigest(),gzip_sha256=hashlib.sha256(packed).hexdigest(),frames=len(frames))
        records.append(record);prepared.append((source,packed))
    # Refuse unfamiliar loader source instead of leaving a broken published URL.
    pages={}
    for section in ('arena','metabolic'):
        path=destination/section/'arena.html';text=path.read_text();old="await fetch('replay.json')"
        if text.count(old)!=1:raise ValueError(f'Unexpected recording loader: {path}')
        text=text.replace(old,"await recordingFetch('replay.json.gz')",1)
        pages[path]=text
    path=destination/'cellular/index.html';text=path.read_text();old="await fetch($('dataset').value)"
    if text.count(old)!=1 or text.count("'mosaic-early.json'")!=1:raise ValueError('Unexpected cellular recording loader')
    pages[path]=text.replace(old,"await recordingFetch($('dataset').value)",1).replace("'mosaic-early.json'","'mosaic-early.json.gz'",1)
    catalog_path=destination/'cellular/catalog.json';catalog=json.loads(catalog_path.read_text())
    expected={pathlib.PurePosixPath(r).name for r in RECORDINGS if r.startswith('cellular/')}
    if {r['file'] for r in catalog}!=expected:raise ValueError('Cellular catalog does not match explicit recording transport set')
    for record in catalog:record['file']+='.gz'
    manifest=dict(schema='static-recording-transport-v1',encoding='gzip of exact original UTF-8 JSON bytes',recordings=records)
    script=template.read_text()
    if script.count('__RECORDINGS__')!=1:raise ValueError('Unexpected recording transport template')
    script=script.replace('__RECORDINGS__',json.dumps(manifest,separators=(',',':')))
    for source,packed in prepared:
        target=source.with_suffix(source.suffix+'.gz');target.write_bytes(packed)
    for path,text in pages.items():
        first=text.find('<script>')
        if first<0:raise ValueError(f'No inline bootstrap in {path}')
        text=text[:first]+'<script src="../recording-transport.js"></script>'+text[first:]
        path.write_text(text)
    catalog_path.write_text(json.dumps(catalog,indent=2)+'\n')
    (destination/'recording-transport.js').write_text(script)
    (destination/'recording-transport.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # These are our freshly copied staging files, not source recordings or art.
    for source,_ in prepared:source.unlink()
    saved=sum(r['original_bytes']-r['compressed_bytes'] for r in records)
    print(f'Exact recording transport: {len(records)} files, {saved:,} bytes freed before shared loader/manifest overhead')
    return manifest
