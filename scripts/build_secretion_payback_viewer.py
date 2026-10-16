#!/usr/bin/env python3
"""Build a recorded secretion-cost instrument from completed native Julia exports."""
import argparse,hashlib,json,math,pathlib,re
FIELDS=('raw','catalyst','bound','free')
def pack(report_path,frame_path):
    source=json.loads(report_path.read_text());frames=json.loads(frame_path.read_text())
    assert source['schema']=='native-secretion-payback-export-v1'
    reports=[]
    for report in source['reports']:
        assert report['schema']=='native-secretion-payback-v2'
        results=[]
        for result in report['results']:
            trace=[]
            assert len(result['control'])==len(result['pulse'])==len(result['differences'])
            for control,pulse,difference in zip(result['control'],result['pulse'],result['differences']):
                assert control['time']==pulse['time']==difference['time']
                delta=[a-b for a,b in zip(pulse['usable'],control['usable'])]
                assert len(delta)==4
                assert abs(delta[0]-difference['private'])<1e-12
                assert abs(sum(delta[1:])-difference['neighbors'])<1e-12
                assert abs(sum(delta)-difference['total'])<1e-12
                for branch in (control,pulse):
                    assert abs(branch['material'])<1e-8 and abs(branch['potential'])<1e-8
                trace.append(dict(time=control['time'],producer=delta[0],neighbors=sum(delta[1:]),total=sum(delta),
                    secreted=sum(pulse[k][0] for k in ('secreted_a','secreted_b')),catalyst=sum(pulse['enzyme']),
                    extra_work=sum(pulse['work'])-sum(control['work']),
                    producer_usable=[control['usable'][0],pulse['usable'][0]],neighbor_usable=[sum(control['usable'][1:]),sum(pulse['usable'][1:])]))
            assert all(a['time']<b['time'] for a,b in zip(trace,trace[1:]))
            # The export retains original samples; report the last sampled increase,
            # not an invented exact event time inside the sampling interval.
            increases=[v['time'] for old,v in zip(trace,trace[1:]) if v['secreted']>old['secreted']+1e-12]
            results.append(dict(raw=result['raw'],trace=trace,last_sampled_secretion_increase=max(increases,default=0),private_payback=result['private_payback'],neighbor_benefit=result['neighbor_benefit'],group_payback=result['group_payback'],source=result['source'],habitat_source=result['habitat_source']))
        assert {r['raw'] for r in results}=={0,2}
        reports.append(dict(seed=report['seed'],horizon=report['horizon'],duration=report['duration'],investment=report['investment'],allocation=report['allocation'],dt=report['dt'],driver_sha256=report['driver_sha256'],full_report_sha256=report['full_report_sha256'],results=results))
    assert len({r['seed'] for r in reports})==len(reports)
    for branch in ('control','pulse'):
        for frame in frames[branch]:
            count=math.prod(frame['shape']);assert frame['shape']==[8,8,8]
            assert len(frame['positions'])==12 and len(frame['radii'])==4
            for name in FIELDS:assert len(frame[name])==count and all(math.isfinite(v) and v>=0 for v in frame[name])
    assert [f['time'] for f in frames['control']]==[f['time'] for f in frames['pulse']]
    captured=next(v for r in source['reports'] if r['seed']==211 for v in r['results'] if v['raw']==2)
    for branch in ('control','pulse'):
        for frame in frames[branch]:
            observation=next(v for v in captured[branch] if abs(v['time']-frame['time'])<1e-10)
            volume=math.prod(frame['extent'])/math.prod(frame['shape'])
            for field,key in (('raw','raw'),('catalyst','enzyme'),('bound','bound'),('free','free')):
                assert abs(sum(frame[field])*volume-sum(observation[key]))<1e-8, (branch,frame['time'],field)
    data=dict(schema='secretion-payback-explorer-v1',sampling=source['sampling'],reports=reports,
        fields=dict(seed=211,raw=2,storage='x fastest; positions are four x, four y, four z coordinates',frames=frames),
        sources=dict(report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),frames_sha256=hashlib.sha256(frame_path.read_bytes()).hexdigest()))
    text=json.dumps(data,separators=(',',':'),allow_nan=False)
    assert not re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:',text)
    return data,text

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=pathlib.Path,default=pathlib.Path('docs/figures/native-secretion-payback.json'));p.add_argument('--frames',type=pathlib.Path,default=pathlib.Path('docs/figures/native-secretion-material-frames.json'));p.add_argument('--destination',type=pathlib.Path,default=pathlib.Path('site/laboratory/secretion-payback'));a=p.parse_args()
    data,payload=pack(a.report,a.frames)
    html=pathlib.Path(__file__).with_name('secretion_payback_viewer.html').read_text().replace('__DATA__',payload.replace('</','<\\/'))
    a.destination.mkdir(parents=True,exist_ok=True)
    (a.destination/'index.html').write_text(html)
    (a.destination/'data.json').write_text(payload+'\n')
    print(f'{len(data["reports"])} banks; {sum(len(x["trace"]) for r in data["reports"] for x in r["results"])} exact trajectory observations; {len(html.encode())} HTML bytes')
if __name__=='__main__':main()
