#!/usr/bin/env python3
"""Publish completed episodic action-credit reports as a compact recorded instrument.

No simulation, learned tensor or private source path is published. All three arms
and every training episode/frozen test must exist before a bank is accepted.
"""
import argparse,hashlib,json,math,pathlib
ARMS=('credit','no_credit','reversed_credit')
def compact_bank(directory):
    report_path=directory/'report.json';report=json.loads(report_path.read_text())
    assert report['schema']=='native-action-credit-v1'
    count=report['episodes'];before=report['before'];assert len(before)==4
    arms={}
    def result(value):
        for key in ('material_residual','potential_residual'):assert abs(value[key])<1e-5
        reward=value['reward'];assert len(reward)==4 and all(math.isfinite(v) for v in reward)
        return dict(reward=reward,material_residual=value['material_residual'],potential_residual=value['potential_residual'])
    before=[result(v) for v in before]
    for name in ARMS:
        arm=json.loads((directory/(name+'.json')).read_text());history=arm['history'];assert len(history)==count and len(arm['after'])==4
        traces=[]
        for trial,value in enumerate(history,1):
            assert value['trial']==trial and value['context']==1+(trial-1)%2
            traces.append(dict(trial=trial,context=value['context'],**result(value)))
        after=[result(v) for v in arm['after']]
        deltas=[[after[j]['reward'][cell]-before[j]['reward'][cell] for cell in range(4)] for j in (2,3)]
        assert all(abs(a-b)<1e-9 for row,expected in zip(deltas,arm['heldout_delta']) for a,b in zip(row,expected))
        assert abs(sum(sum(row) for row in deltas)/8-report['summary'][name]['heldout_mean_delta'])<1e-9
        arms[name]=dict(history=traces,after=after,heldout_delta=deltas,mark_max=arm['mark_max'])
    return dict(seed=report['seed'],noise_seed=report['noise_seed'],episodes=count,horizon=report['horizon'],rate=report['rate'],sigma=report['sigma'],contexts=report['contexts'],before=before,arms=arms,report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),driver_sha256=report['driver_sha256'],factory_sha256=report['factory_sha256'],source_pin=report['source_pin'],cellular_source_pin=report['cellular_source_pin'])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('reports',type=pathlib.Path,nargs='+');p.add_argument('--destination',type=pathlib.Path,default=pathlib.Path('site/laboratory/action-credit'));a=p.parse_args()
    banks=[compact_bank(d) for d in a.reports];assert len({b['seed'] for b in banks})==len(banks)
    for bank in banks[1:]:
        for key in ('episodes','horizon','rate','sigma','contexts','source_pin','cellular_source_pin'):
            assert bank[key]==banks[0][key], f'Banks differ in {key}'
    data=dict(schema='action-credit-explorer-v1',banks=banks,scope='Externally reset independent physical episodes; only action readout offsets persist. Frozen tests use unseen food mixtures and no exploration noise. This is study-only episodic training, not autonomous lifetime learning.',payoff='Change in battery plus efficiency-weighted stored food chemical potential during one physical episode. All plotted differences subtract the original untrained policy in the same bank, cell and food mixture.',physical_grid=[8,8,8],cells_per_bank=4,action_rows=[4,5,7],growth_division=False,physical_sources=False)
    payload=json.dumps(data,separators=(',',':'),allow_nan=False)
    import re
    assert not re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:',payload), 'Private reference in compact metadata'
    template=pathlib.Path(__file__).with_name('action_credit_viewer.html').read_text()
    text=template.replace('__DATA__',payload.replace('</','<\\/'))
    a.destination.mkdir(parents=True,exist_ok=True)
    (a.destination/'index.html').write_text(text)
    (a.destination/'data.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(f'{len(banks)} completed banks; {sum(b["episodes"] for b in banks)*3} actual training episodes; {len(text.encode())} HTML bytes')
if __name__=='__main__':main()
