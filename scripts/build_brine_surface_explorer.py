"""Run in Blender Python: package a display surface from immutable native brine fields."""
import argparse,gzip,hashlib,json,math,sys,time
from pathlib import Path
import numpy as np
import openvdb as vdb

def require(ok,message):
    if not ok:raise ValueError(message)

def sample(field,points):
    shape=np.array(field.shape);q=np.clip((points+1)*shape/2-.5,0,shape-1);lo=np.floor(q).astype(int);hi=np.minimum(lo+1,shape-1);weight=q-lo;result=np.zeros(len(points))
    for i in (0,1):
        for j in (0,1):
            for k in (0,1):
                indices=[(hi if bit else lo)[:,d] for d,bit in enumerate((i,j,k))]
                w=np.prod(np.column_stack([weight[:,d] if bit else 1-weight[:,d] for d,bit in enumerate((i,j,k))]),axis=1)
                result+=field[tuple(indices)]*w
    return result

def quantize(values):
    lo=np.min(values,axis=0);hi=np.max(values,axis=0);span=hi-lo
    encoded=np.rint((values-lo)/np.where(span>0,span,1)*65535).clip(0,65535).astype('<u2')
    restored=lo+encoded.astype(float)/65535*span
    return encoded,lo.tolist(),span.tolist(),np.max(np.abs(values-restored),axis=0).tolist()

def build(epoch,destination,iso=.30,adaptivity=.02,max_bytes=8_000_000):
    started=time.time();report_path=epoch/'run.json';report=json.loads(report_path.read_text());inventory=json.loads((epoch/'verified-inventory.json').read_text())
    require(report.get('success') is True and not report.get('export_failure'),'Source epoch did not complete cleanly')
    frame=report['frames'][-1];shape=tuple(frame['shape']);require(frame['kind']=='final' and len(shape)==3,'Use completed native final fields, not replay')
    require(iso>0 and 0<=adaptivity<=1,'Invalid display extraction settings')
    expected={r['path']:r for r in inventory['files']};fields={};hashes={};dtypes={};stats={}
    for name in ('mineral','phi','a','b'):
        dtype={'Float32':'<f4','Float64':'<f8'}[frame['dtypes'][name]];relative=frame['directory']+'/'+name+('.f64' if dtype=='<f8' else '.f32');file=epoch/relative
        require(file.stat().st_size==math.prod(shape)*np.dtype(dtype).itemsize==expected[relative]['bytes'],'Source field size mismatch')
        hashes[name]=hashlib.file_digest(file.open('rb'),'sha256').hexdigest();require(hashes[name]==expected[relative]['sha256'],'Source field digest mismatch: '+name)
        fields[name]=np.memmap(file,mode='r',dtype=dtype,shape=shape,order='F');dtypes[name]=frame['dtypes'][name];stats[name]=(file.stat().st_size,file.stat().st_mtime_ns)
    grid=vdb.DoubleGrid() if dtypes['mineral']=='Float64' else vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(fields['mineral']))
    points,triangles,quads=grid.convertToPolygons(isovalue=iso,adaptivity=adaptivity);del grid
    vertices=(np.asarray(points,dtype=np.float64)+.5)*2/np.array(shape)-1
    triangles=np.asarray(triangles,dtype=np.uint32).reshape(-1,3);quads=np.asarray(quads,dtype=np.uint32).reshape(-1,4)
    faces=np.concatenate((triangles,quads[:,[0,1,2]],quads[:,[0,2,3]]),axis=0)
    require(len(vertices)>0 and len(faces)>0 and int(faces.max())<len(vertices),'Invalid extracted topology')
    normals=np.zeros_like(vertices);cross=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]])
    for i in range(3):np.add.at(normals,faces[:,i],cross)
    length=np.linalg.norm(normals,axis=1);require(np.all(length>0),'Zero-area vertex normal');normals/=length[:,None]
    sampled=np.column_stack([sample(fields[name],vertices) for name in ('phi','a','b')]);require(np.isfinite(sampled).all(),'Nonfinite sampled fields')
    positions,pmin,pspan,perror=quantize(vertices);optics,smin,sspan,serror=quantize(sampled)
    nencoded=np.rint(normals*32767).clip(-32767,32767).astype('<i2')
    packed=bytearray();arrays={}
    for name,array in [('position',positions),('normal',nencoded),('fields',optics),('triangles',faces.astype('<u4'))]:
        while len(packed)%4:packed.append(0)
        arrays[name]=dict(offset=len(packed),bytes=array.nbytes,type=str(array.dtype),shape=list(array.shape));packed.extend(array.tobytes(order='C'))
    compressed=gzip.compress(packed,compresslevel=9,mtime=0)
    extent=report['recipe'].get('extent',24.);extent=[float(extent)]*3 if isinstance(extent,(float,int)) else extent
    metadata=dict(schema='brine-surface-explorer-v1',title='The chamber the brine kept',source=dict(native_shape=shape,extent=extent,array_order='Julia column-major, x-fast',field_dtypes=dtypes,field_sha256=hashes,run_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),recipe_sha256=report['recipe_sha256'],base_commit=report.get('base_commit'),sources=report.get('sources'),step=frame['step'],time=frame['time']),
        surface=dict(isovalue=iso,adaptivity=adaptivity,native_grid_used_directly=True,grid_type='DoubleGrid' if dtypes['mineral']=='Float64' else 'FloatGrid',vertices=len(vertices),source_triangles=len(triangles),source_quads=len(quads),display_triangles=len(faces),clipping='Client-side uncapped plane; complete display surface retained',sampling='Trilinear samples from original native phi/A/B at extracted vertices; edge-clamped to native voxel centers',field_names=['phi','a','b'],sample_ranges=dict(phi=[float(sampled[:,0].min()),float(sampled[:,0].max())],ratio=[float((sampled[:,1]/(sampled[:,1]+sampled[:,2]+1e-12)).min()),float((sampled[:,1]/(sampled[:,1]+sampled[:,2]+1e-12)).max())],total=[float((sampled[:,1]+sampled[:,2]).min()),float((sampled[:,1]+sampled[:,2]).max())]),position_min=pmin,position_span=pspan,position_max_error_model_units=(np.array(perror)*np.array(extent)/2).tolist(),field_min=smin,field_span=sspan,field_max_error=serror,normal_quantization='Signed16 components; renormalized in shader',simplification='OpenVDB surface adaptivity only; no scalar downsampling or mesh decimation'),
        binary=dict(file='mesh.bin.gz',encoding='gzip',bytes=len(packed),gzip_bytes=len(compressed),sha256=hashlib.sha256(packed).hexdigest(),gzip_sha256=hashlib.sha256(compressed).hexdigest(),arrays=arrays),
        builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),display=dict(default_clip_axis=0,default_clip_model_position=13.2,default_keep='less',optical_mapping='Chosen optical colors and lighting from saved phase and reagents; not an oxidation or optical-field simulation'))
    template=Path(__file__).resolve().parent.parent/'web/brine_surface_explorer.html';html=template.read_text();payload=json.dumps(metadata,separators=(',',':'),allow_nan=False)+'\n'
    size=len(compressed)+len(payload.encode())+len(html.encode());require(size<=max_bytes,f'Exhibit {size} bytes exceeds explicit {max_bytes} byte budget; choose and record a different display adaptivity or budget')
    for name in fields:
        suffix='.f64' if dtypes[name]=='Float64' else '.f32';file=epoch/frame['directory']/(name+suffix);require(stats[name]==(file.stat().st_size,file.stat().st_mtime_ns),'Source changed during extraction')
    destination.mkdir(parents=True,exist_ok=True);(destination/'mesh.bin.gz').write_bytes(compressed);(destination/'mesh.json').write_text(payload);(destination/'index.html').write_text(html)
    print(json.dumps(dict(vertices=len(vertices),triangles=len(faces),compressed_bytes=len(compressed),exhibit_bytes=size,adaptivity=adaptivity,position_error_model_units=metadata['surface']['position_max_error_model_units'],elapsed_seconds=time.time()-started,destination=str(destination))),flush=True)
    return metadata

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('epoch',type=Path);parser.add_argument('destination',type=Path);parser.add_argument('--iso',type=float,default=.30);parser.add_argument('--adaptivity',type=float,default=.02);parser.add_argument('--max-bytes',type=int,default=8_000_000)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else None);build(args.epoch,args.destination,args.iso,args.adaptivity,args.max_bytes)
if __name__=='__main__':main()
