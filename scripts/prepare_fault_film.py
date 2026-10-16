"""Prepare the stored-state film and contact sheets for explicit visual review."""
import argparse,hashlib,json
from pathlib import Path
from PIL import Image,ImageDraw
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);p.add_argument('frames',type=Path);a=p.parse_args()
record=json.loads((a.frames/'frames.json').read_text());verify=json.loads((a.capture/'capture-verification.json').read_text());assert len(record['frames'])==len(verify['frames'])==40
images=[];frames=[]
for r,v in zip(record['frames'],verify['frames']):
 assert (r['step'],r['time'],r['source_sha256'])==(v['step'],v['time'],v['sha256'])
 path=a.frames/r['image'];images.append(path)
 frames.append({'file':path.name,'time':r['time'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
sequence={'id':'fault-remembers-motion','title':'The fault remembers','description':'A laminated material takes tension and shear. Conductive dose marks its sheets as a fault opens through seven layers. Forty actual recorded states.','fps':8,'poster_frame_index':39,'source':{'study':'studies/fault_remembers.jl','render_study':'scripts/render_fault_motion.py','capture_verification_sha256':hashlib.sha256((a.capture/'capture-verification.json').read_bytes()).hexdigest(),'render_record_sha256':hashlib.sha256((a.frames/'frames.json').read_bytes()).hexdigest(),'same_approved_final_state':True,'temporal_interpolation':False,'actual_samples':40,'physical_time_interval':[0,1.015],'global_dose_optical_max':.04419005268498425,'camera':'fixed throughout','native_nodes':5670,'native_bonds':63209},'frames':frames}
(a.frames/'sequence.json').write_text(json.dumps(sequence,indent=2)+'\n')
for start in range(0,40,10):
 sheet=Image.new('RGB',(1280,552),(12,18,24));draw=ImageDraw.Draw(sheet)
 for local,path in enumerate(images[start:start+10]):
  i=start+local;x=(local%5)*256;y=(local//5)*276
  im=Image.open(path).convert('RGB');im.thumbnail((256,256));sheet.paste(im,(x,y));r=record['frames'][i]
  draw.text((x+8,y+258),f"{i+1}/40   step {r['step']}   t {r['time']:.4f}",fill=(235,225,204))
 sheet.save(a.frames/f'contact-{start//10}.jpg',quality=92)
print('Prepared40frame sequence and4contact sheets; visual approval still required.')
