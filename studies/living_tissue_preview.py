"""Render actual captured sections; this file contains no simulation algorithm."""
import json
from pathlib import Path
import sys
import argparse
import shutil
import subprocess
import numpy as np
from PIL import Image, ImageDraw, ImageFont

def read_tissue_fields(directory, dims):
    def field(name):
        raw = np.fromfile(directory / f"{name}.f32", dtype="<f4").reshape(dims, order="F")
        return raw[:, :, (dims[2]-1)//2].T[::-1]
    return tuple(field(name) for name in ("biomass", "matrix", "excitation"))

def tissue_rgb(fields, signed=False):
    b, m, x = fields
    rgb = np.zeros((*b.shape, 3), dtype=np.float32)
    if signed:
        positive, negative = np.maximum(0, x), np.maximum(0, -x)
        rgb += np.array([.02, .035, .05])
        rgb += positive[...,None]*np.array([.24, .82, .86])
        rgb += negative[...,None]*np.array([.52, .22, .57])
    else:
        bulk = 1-np.exp(-2.7*b)
        matrix = 1-np.exp(-8*m)
        active = np.clip(x,0,1)*bulk
        rgb += np.array([.015,.027,.039])
        rgb += bulk[...,None]*np.array([.22,.16,.33])
        rgb += matrix[...,None]*np.array([.70,.35,.10])
        rgb += active[...,None]*np.array([.10,.90,1.10])
    return np.uint8(np.clip(rgb,0,1)**.8*255)

def preview(root):
    root = Path(root)
    run = json.loads((root / "run.json").read_text())
    dims = tuple(run["dimensions"])
    selected = [run["frame_steps"][round((len(run["frame_steps"])-1)*f)] for f in (0,.2,.4,.6,.8,1)]
    sheet = Image.new("RGB", (3*600, 2*640), "#071015")
    for index, step in enumerate(selected):
        directory = root / "frames" / f"{step:07d}"
        fields = read_tissue_fields(directory, dims)
        # Fixed exposure for every time: violet living bulk, amber matrix, and
        # cyan positive excitation. This is an optical choice, not emitted light.
        tile = Image.fromarray(tissue_rgb(fields)).resize((560,560),Image.Resampling.BICUBIC)
        left, top = (index%3)*600+20, (index//3)*640+20
        sheet.paste(tile,(left,top))
        draw=ImageDraw.Draw(sheet)
        draw.text((left,top+578),f"t = {step*run['dt']:.1f}  ·  actual middle section",fill="#b6c4c8")
    path=root/"sections.png";sheet.save(path)
    return path

def movie(root, fps=10, size=576):
    root = Path(root)
    run = json.loads((root / "run.json").read_text())
    dims = tuple(run["dimensions"])
    destination = root / "activity-film.mp4"
    if destination.exists():
        raise FileExistsError(destination)
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        raise RuntimeError("ffmpeg is required for the recorded film")
    width, height = 2*size+64, size+112
    title_font = ImageFont.load_default(size=23)
    label_font = ImageFont.load_default(size=14)
    command = [ffmpeg, "-n", "-hide_banner", "-loglevel", "warning", "-f", "rawvideo",
               "-pix_fmt", "rgb24", "-s", f"{width}x{height}", "-r", str(fps), "-i", "-",
               "-an", "-c:v", "libx264", "-threads", "1", "-crf", "19", "-preset", "fast",
               "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(destination)]
    with (root / "movie-ffmpeg.log").open("x") as log:
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stderr=log)
        try:
            for step in run["frame_steps"]:
                fields = read_tissue_fields(root / "frames" / f"{step:07d}", dims)
                frame = Image.new("RGB", (width, height), "#071015")
                draw = ImageDraw.Draw(frame)
                draw.text((20, 13), run["title"], font=title_font, fill="#d2e5e8")
                draw.text((width-190, 18), f"t = {step*run['dt']:.1f}", font=label_font, fill="#9fc9d0")
                for panel, signed in enumerate((False, True)):
                    tile = Image.fromarray(tissue_rgb(fields, signed=signed)).resize((size,size), Image.Resampling.BICUBIC)
                    frame.paste(tile, (20+panel*(size+24), 54))
                draw.text((20, size+66), "Tissue / built matrix / positive activity", font=label_font, fill="#94aeb8")
                draw.text((size+44, size+66), "Fast state: + cyan / - violet", font=label_font, fill="#94aeb8")
                view = "recorded sheet" if dims[2]==1 else "recorded middle section"
                draw.text((20, size+90), f"{dims[0]} x {dims[1]} samples · {view} · fixed exposure", font=label_font, fill="#5c7c88")
                process.stdin.write(frame.tobytes())
        finally:
            process.stdin.close()
            code = process.wait()
        if code:
            raise RuntimeError(f"ffmpeg failed; see {root/'movie-ffmpeg.log'}")
    return destination

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root");parser.add_argument("--film",action="store_true")
    arguments=parser.parse_args()
    print(movie(arguments.root) if arguments.film else preview(arguments.root))
