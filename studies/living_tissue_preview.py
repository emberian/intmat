"""Render actual captured sections; this file contains no simulation algorithm."""
import json
from pathlib import Path
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFilter

def preview(root):
    root = Path(root)
    run = json.loads((root / "run.json").read_text())
    dims = tuple(run["dimensions"])
    selected = [run["frame_steps"][round((len(run["frame_steps"])-1)*f)] for f in (0,.2,.4,.6,.8,1)]
    sheet = Image.new("RGB", (3*600, 2*640), "#071015")
    for index, step in enumerate(selected):
        directory = root / "frames" / f"{step:07d}"
        def field(name):
            raw = np.fromfile(directory / f"{name}.f32", dtype="<f4").reshape(dims, order="F")
            return raw[:, :, (dims[2]-1)//2].T[::-1]
        b, m, x = (field(name) for name in ("biomass", "matrix", "excitation"))
        # Fixed exposure for every time: violet living bulk, amber matrix, and
        # cyan positive excitation. This is an optical choice, not emitted light.
        bulk = 1-np.exp(-2.7*b)
        matrix = 1-np.exp(-8*m)
        active = np.clip(x,0,1)*bulk
        rgb = np.zeros((*b.shape,3), dtype=np.float32)+np.array([.015,.027,.039])
        rgb += bulk[...,None]*np.array([.22,.16,.33])
        rgb += matrix[...,None]*np.array([.70,.35,.10])
        rgb += active[...,None]*np.array([.10,.90,1.10])
        rgb = np.clip(rgb,0,1)
        tile = Image.fromarray(np.uint8(rgb**.8*255)).resize((560,560),Image.Resampling.BICUBIC)
        left, top = (index%3)*600+20, (index//3)*640+20
        sheet.paste(tile,(left,top))
        draw=ImageDraw.Draw(sheet)
        draw.text((left,top+578),f"t = {step*run['dt']:.1f}  ·  actual middle section",fill="#b6c4c8")
    path=root/"sections.png";sheet.save(path)
    return path

if __name__ == "__main__":
    print(preview(sys.argv[1]))
