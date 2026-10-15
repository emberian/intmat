#!/usr/bin/env python3
"""Static scientific plot of measured paid-patch assay rows; no simulation."""
import json, math, sys, hashlib
from pathlib import Path
source, target = map(Path, sys.argv[1:3])
assert not target.exists(), "fresh figure required"
data=json.loads(source.read_text()); arms=data["arms"]
styles={"immobile_full":("#256d82","immobile raw · matrix",""),"mobile_full":("#bf4936","mobile raw · matrix",""),"immobile_catalytic":("#256d82","immobile raw · catalyst only","7 4"),"mobile_catalytic":("#bf4936","mobile raw · catalyst only","7 4")}
assert set(arms)==set(styles), "four complete arms required"
W,H=1200,870
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">', '<rect width="1200" height="870" fill="#fbfaf6"/>', '<style>text{font-family:Arial,sans-serif;fill:#283238}.title{font-size:27px;font-weight:600}.label{font-size:16px}.tick{font-size:13px}.note{font-size:14px}</style>', '<text x="48" y="43" class="title">A paid patch: replenishment and retention</text>', '<text x="48" y="69" class="note">Same catalyst dose and costs · actual native Julia field laws · stationary builder after one paid pulse</text>']
for n,(key,(color,label,dash)) in enumerate(styles.items()):
    x=48+(n%2)*540;y=100+(n//2)*25
    s += [f'<path d="M{x},{y}h42" stroke="{color}" stroke-width="3" stroke-dasharray="{dash}"/>',f'<text x="{x+52}" y="{y+5}" class="label">{label}</text>']
panels=[("Raw converted → soluble food", "total mass, both species",lambda r:sum(r["audit"]["solubilized"])),("Food remaining reversibly bound", "total mass, both species",lambda r:sum(r["audit"]["bound"])),("Free-food exposure at builder", "sum concentration × time; read-only sensor",lambda r:sum(r["exposure"][0::2])),("Free-food exposure 0.4 away", "sum concentration × time; read-only sensor",lambda r:sum(r["exposure"][1::2]))]
for index,(title,unit,value) in enumerate(panels):
    x=85+(index%2)*580;y=207+(index//2)*300;pw,ph=485,205
    ymax=max(value(r) for a in arms.values() for r in a["rows"])*1.08
    xmax=max(a["rows"][-1]["time"] for a in arms.values())
    s += [f'<text x="{x}" y="{y-33}" class="label" font-weight="600">{title}</text>',f'<text x="{x}" y="{y-12}" class="tick">{unit}</text>']
    for j in range(5):
        v=ymax*j/4;yy=y+ph-ph*j/4
        label=f"{v:.3g}"
        s += [f'<path d="M{x},{yy}h{pw}" stroke="#dce1dc"/>',f'<text x="{x-9}" y="{yy+4}" text-anchor="end" class="tick">{label}</text>']
    for t in (0,20,40,60):
        xx=x+pw*t/xmax
        s += [f'<text x="{xx}" y="{y+ph+23}" text-anchor="middle" class="tick">{t}</text>']
    s += [f'<text x="{x+pw/2}" y="{y+ph+43}" text-anchor="middle" class="tick">model time</text>']
    for key,(color,label,dash) in styles.items():
        points=" ".join(f"{x+pw*r['time']/xmax:.2f},{y+ph-ph*value(r)/ymax:.2f}" for r in arms[key]["rows"])
        s += [f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.7" stroke-dasharray="{dash}"/>']
s += ['<text x="48" y="819" class="note">No uptake, motion, maintenance or further construction after the pulse. Exposure is opportunity, not captured food or fitness.</text>', '<text x="48" y="843" class="note">Matrix-off retains body porosity. Conversion curves overlap within each raw-transport pair; binding and retention differ.</text>', '</svg>']
target.write_text("\n".join(s))
target.with_suffix(".identity.json").write_text(json.dumps({"report_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),"plot_driver_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"source":data["source"],"habitat_source":data["habitat_source"],"assay_driver_sha256":data["driver_sha256"]},indent=2))
