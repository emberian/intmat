#!/usr/bin/env python3
"""Package completed construction-authority evidence; no simulation or inferred frames."""
import argparse,hashlib,json,math,pathlib,re
TRACE_KEYS=('time','ids','usable','learning_work','work','heat','supplied','secretion_a','secretion_b','cumulative_requests','enzyme','bound','energy','stores','reserve_action','material','potential')
def observation(v):
    assert v['ids']==[1,2,3,4]
    assert abs(v['material'])<1e-8 and abs(v['potential'])<1e-8
    assert len(v['cumulative_requests'])==8 and len(v['secretion_a'])==len(v['secretion_b'])==4
    assert sum(v['secretion_a'])+sum(v['secretion_b'])<=sum(v['cumulative_requests'])+1e-6
    if 'all_residents_valid' in v:assert v['all_residents_valid']
    return {k:v[k] for k in TRACE_KEYS if k in v}
def pack(root):
    paths={key:root/f'native-construction-{key}.json' for key in ('zero-authority','fixed-reserve','authority','credit-gpu')}
    source={k:json.loads(p.read_text()) for k,p in paths.items()};pilots={}
    for key in ('zero-authority','fixed-reserve'):
        report=source[key];arms={}
        for arm,value in report['results'].items():
            trace=[observation(v) for v in value['trace']];assert [v['time'] for v in trace]==list(range(2,25,2))
            branches={}
            for name,b in value['branches'].items():
                assert name in ('learned_continued','erased_continued','learned_shifted','erased_shifted') and b['ids']==[1,2,3,4]
                assert b['time']==28 and len(b['return_value'])==4
                assert abs(b['material'])<1e-8 and abs(b['potential'])<1e-8
                branches[name]=dict(time=b['time'],return_value=b['return_value'])
            assert len(branches)==4
            arms[arm]=dict(trace=trace,branches=branches)
        pilots[key]=dict(seed=report['seed'],noise_seed=report['noise_seed'],training=report['training'],test_time=report['test_time'],arms=arms,driver_sha256=report['driver_sha256'],source=report['source'],cellular_source=report['cellular_source'])
    authority=[dict(reserve=r['reserve'],trace=[observation(v) for v in r['trace']]) for r in source['authority']]
    assert [r['reserve'] for r in authority]==[None,.2]
    assert all([v['time'] for v in r['trace']]==[2,4,6,8] for r in authority)
    gpu=source['credit-gpu'];gate=dict(cases=[dict(feed=r['feed'],dt=r['dt'],errors=r['errors'],residual=r['residual']) for r in gpu['results']],study_sha256=gpu['study_sha256'],source=gpu['source'],habitat_source=gpu['habitat_source'])
    data=dict(schema='construction-credit-explorer-v1',authority=authority,pilots=pilots,gpu=gate,source_hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()})
    payload=json.dumps(data,separators=(',',':'),allow_nan=False);assert not re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:',payload)
    return data,payload

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=pathlib.Path,default=pathlib.Path('docs/figures'));p.add_argument('--destination',type=pathlib.Path,default=pathlib.Path('site/laboratory/construction-credit'));a=p.parse_args();data,payload=pack(a.source)
    html=pathlib.Path(__file__).with_name('construction_credit_viewer.html').read_text().replace('__DATA__',payload.replace('</','<\\/'));a.destination.mkdir(parents=True,exist_ok=True)
    (a.destination/'index.html').write_text(html);print(f'{len(html.encode())} HTML bytes; {len(payload)} data bytes; completed paired authority and two pilots')
if __name__=='__main__':main()
