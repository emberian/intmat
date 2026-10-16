"""Encode forty actual observations and prepare decoded review contact sheets."""
import argparse,hashlib,json,subprocess
from pathlib import Path
from PIL import Image,ImageDraw
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);p.add_argument('frames',type=Path);p.add_argument('destination',type=Path);a=p.parse_args();record=json.loads((a.frames/'frames.json').read_text());verify=json.loads((a.capture/'capture-verification.json').read_text());assert sorted(map(int,record['frames']))==list(range(40));a.destination.mkdir(parents=True,exist_ok=True)
frames=[]
for i in range(40):
 r=record['frames'][str(i)];v=verify['frames'][i+1];assert (r['step'],r['time'])==(v['step'],v['time']);assert r['image_sha256']==hashlib.sha256((a.frames/r['image']).read_bytes()).hexdigest()
 for k,value in v['fields_sha256'].items():assert r['source_fields'][k+'.f32']==value
 frames.append({'video_time':i/8,'simulation_time':r['time'],'step':r['step'],'source_image_sha256':r['image_sha256']})
film=a.destination/'film.mp4';subprocess.run(['ffmpeg','-y','-loglevel','error','-framerate','8','-i',str(a.frames/'%05d.png'),'-frames:v','40','-c:v','libx264','-preset','slow','-crf','20','-pix_fmt','yuv420p','-movflags','+faststart',str(film)],check=True)
Image.open(a.frames/'00039.png').convert('RGB').save(a.destination/'poster.webp','WEBP',quality=95,method=6)
probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-count_frames','-show_entries','stream=width,height,nb_read_frames,r_frame_rate,duration','-of','json',str(film)]));stream=probe['streams'][0];assert stream['nb_read_frames']=='40' and stream['width']==record['size']==stream['height']
meta={'schema':'interior-animation','id':'wave-memory-motion-128','title':'Where the note goes','description':'Three coherent packets write and remodel finite mineral. Their changing interfaces slow subsequent waves; warm optical light follows actual instantaneous pressure. Forty stored native Julia observations.','source':{'study':'studies/wave_memory.jl','render_study':'scripts/render_wave_memory_frame.py','native_shape':[128]*3,'actual_samples':40,'temporal_interpolation':False,'restart_checkpoints':False,'camera':'fixed throughout','capture_verification_sha256':hashlib.sha256((a.capture/'capture-verification.json').read_bytes()).hexdigest(),'render_record_sha256':hashlib.sha256((a.frames/'frames.json').read_bytes()).hexdigest(),'renderer_sha256':record['renderer_sha256'],'physical_time_interval':[frames[0]['simulation_time'],frames[-1]['simulation_time']],'mineral_isovalue':.3,'pressure_optics':'Emission strength .3*min(abs(pressure),2), same fixed transfer for every frame; assigned optical surrogate.','cutaway':'Exact plane clipping, open mineral y=-.04 and precursor shell y=.12.','model_qualification':'Finite donor-limited conversion/dissolution; no latent-energy debit, no calibrated acoustic or crystallization claim.'},'playback':{'fps':8,'duration_seconds':5.,'width':record['size'],'height':record['size'],'rendered_width':record['size'],'rendered_height':record['size'],'samples':record['samples'],'interpolation':'none; one video frame per stored simulation observation'},'frames':frames,'assets':{'video':'film.mp4','poster':'poster.webp'},'sha256':{'film.mp4':hashlib.sha256(film.read_bytes()).hexdigest(),'poster.webp':hashlib.sha256((a.destination/'poster.webp').read_bytes()).hexdigest()},'review':{'status':'pending'}}
meta['video_sha256']=meta['sha256']['film.mp4'];meta['poster_sha256']=meta['sha256']['poster.webp'];meta['poster_frame_index']=39;meta['source']['interpretation']='Prepared pressure packets and a finite ellipsoidal support are authored. Computed pressure converts precursor to slower mineral; finite dissolution returns mineral to precursor. Warm light is an assigned instantaneous-pressure surrogate. No calibrated crystallization or closed energy accounting claim. The optical cut planes expose actual material interfaces.'
(a.destination/'film.json').write_text(json.dumps(meta,indent=2))
decoded=a.frames/'decoded';decoded.mkdir(exist_ok=True);subprocess.run(['ffmpeg','-y','-loglevel','error','-i',str(film),str(decoded/'%05d.png')],check=True)
assert len(list(decoded.glob('*.png')))==40
for block in range(4):
 sheet=Image.new('RGB',(1280,552),(12,18,24));draw=ImageDraw.Draw(sheet)
 for local in range(10):
  i=block*10+local;x=(local%5)*256;y=(local//5)*276;im=Image.open(decoded/f'{i+1:05d}.png').convert('RGB');im.thumbnail((256,256));sheet.paste(im,(x,y));draw.text((x+8,y+258),f"{i+1}/40  t={frames[i]['simulation_time']:.2f}",fill=(235,225,204))
 sheet.save(a.frames/f'contact-{block}.jpg',quality=93)
print(json.dumps({'files':{q.name:q.stat().st_size for q in a.destination.iterdir() if q.is_file()},'decoded_frames':40,'review':'pending'},indent=2))
