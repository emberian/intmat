"""Package recorded tissue sections into an offline, view-only HTML player."""
import argparse
import base64
import json
from pathlib import Path


def package_playback(root, destination=None):
    root = Path(root)
    run = json.loads((root / "run.json").read_text())
    section = run["section"]
    width, height, channels = section["shape"]
    if channels != 3 or section["array_order"] != "F":
        raise ValueError("expected three channel-major section planes")
    frames = []
    history = {sample["step"]: sample for sample in run["history"]}
    for step in run["frame_steps"]:
        data = (root / "frames" / f"{step:07d}" / "section.u8").read_bytes()
        if len(data) != width * height * channels:
            raise ValueError(f"malformed section at step {step}")
        metrics = {key: history[step].get("physical_" + key, history[step][key]) for key in
                   ("biomass", "matrix", "precursor", "spent", "mean_permeability")}
        if "reserve" in history[step]:
            metrics["reserve"] = history[step].get("physical_reserve", history[step]["reserve"])
        frames.append({"step": step, "time": history[step]["time"],
                       "data": base64.b64encode(data).decode("ascii"),
                       "metrics": metrics})
    payload = {"title": run["title"], "description": run["description"],
               "dimensions": run["dimensions"], "extent": run["extent"],
               "section": section, "frames": frames, "probe_site": run["probe_site"],
               "source_sha256": run["source_sha256"], "dt": run["dt"],
               "initialization": run.get("initialization", "Explicit seeded tissue geometry."),
               "inventory_units": "integrated concentration × physical area/volume" if "physical_material" in run["history"][0] else "native concentration sums",
               "budget": run.get("budget", run.get("material_convention", "Finite initial material."))}
    template = Path(__file__).with_suffix(".html").read_text()
    encoded = json.dumps(payload, separators=(",", ":")).replace("<", "\\u003c")
    page = template.replace("/*__TISSUE_DATA__*/", encoded)
    destination = Path(destination) if destination else root / "playback.html"
    with destination.open("x") as stream:
        stream.write(page)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root")
    parser.add_argument("--output")
    arguments = parser.parse_args()
    print(package_playback(arguments.root, arguments.output))
