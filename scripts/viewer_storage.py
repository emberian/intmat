"""Lossless, deterministic storage for the viewer's already-quantized fields.

Only field transport changes. Scales, frames, shaders and decoded bytes stay
identical. A single gzip stream shares the dictionary across channels/frames.
"""
import base64
import gzip
import json
import re

RECORDING = re.compile(r'(<script id="recording" type="application/json">)(.*?)(</script>)', re.S)
BOOTSTRAP = """if(data.fields_gzip){
if(typeof DecompressionStream==='undefined')throw Error('This recording needs gzip decompression. Please use a current browser (Chrome, Firefox, Safari or Edge).');
const packed=Uint8Array.from(atob(data.fields_gzip),c=>c.charCodeAt(0));
const stream=new Blob([packed]).stream().pipeThrough(new DecompressionStream('gzip'));
data.fields=JSON.parse(await new Response(stream).text());delete data.fields_gzip;
}
"""

def unpack(data):
    """Return the original base64 field arrays, for byte-level verification."""
    if 'fields_gzip' in data:
        return json.loads(gzip.decompress(base64.b64decode(data['fields_gzip'])))
    return data['fields']

def compress_html(text):
    """Upgrade a complete existing viewer, retaining all custom source edits."""
    match = RECORDING.search(text)
    if not match:
        raise ValueError('No embedded recording found')
    data = json.loads(match[2])
    if 'fields_gzip' in data:
        return text
    original = data.pop('fields')
    packed = gzip.compress(json.dumps(original,separators=(',',':')).encode(),compresslevel=9,mtime=0)
    data['fields_gzip'] = base64.b64encode(packed).decode()
    assert unpack(data) == original
    payload = json.dumps(data,separators=(',',':')).replace('</','<\\/')
    text = text[:match.start(2)] + payload + text[match.end(2):]
    start = "const data=JSON.parse(document.getElementById('recording').textContent)"
    if text.count(start) != 1 or text.count('if(!gl)') != 1:
        raise ValueError('Unrecognized viewer bootstrap; refusing to alter it')
    text = text.replace(start,'(async()=>{\n'+start,1)
    text = text.replace('if(!gl)',BOOTSTRAP+'if(!gl)',1)
    end = text.rfind('</script>')
    text = text[:end]+"\n})().catch(e=>{const box=document.getElementById('error');box.style.display='block';box.textContent=e.message;console.error(e);});\n"+text[end:]
    return text

if __name__ == '__main__':
    import argparse
    import pathlib
    parser = argparse.ArgumentParser(description='Losslessly compress explicit existing viewer HTML files; preserve originals in a private backup directory.')
    parser.add_argument('--backup', required=True, type=pathlib.Path)
    parser.add_argument('viewers', nargs='+', type=pathlib.Path)
    args = parser.parse_args()
    args.backup.mkdir(parents=True,exist_ok=True)
    for index,path in enumerate(args.viewers):
        old = path.read_text()
        new = compress_html(old)
        if new == old:
            print(f'{path}: already compressed')
            continue
        backup = args.backup / f'{index}-{path.parent.name}-{path.name}'
        with backup.open('x') as handle:
            handle.write(old)
        temporary = path.with_suffix('.html.tmp')
        temporary.write_text(new)
        temporary.replace(path)
        print(f'{path}: {len(old.encode())} -> {len(new.encode())} bytes; original {backup}')
