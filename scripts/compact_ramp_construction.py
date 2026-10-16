#!/usr/bin/env python3
"""Validate and compact one bounded continuous-observation ramp experiment."""
import json
import pathlib
import statistics
import sys

root = pathlib.Path(__file__).resolve().parents[1]
source = root / sys.argv[1]
data = json.loads((source / "report.json").read_text())
arms = data["arms"]
max_norm_error = 0.0
max_closure = 0.0
funding = []
for arm in arms.values():
    for row in arm["history"] + arm["probes"]:
        f = row["features"]  # JSON3 matrix serialization is column-major.
        max_norm_error = max(max_norm_error, *(abs(sum(f[i + 4*k]**2 for k in range(3))-2) for i in range(4)))
        max_closure = max(max_closure, abs(row["material"]), abs(row["potential"]))
        funding.append(sum(row["funded"]))
assert max_norm_error < 1e-6
assert max_closure < 1e-10
assert max(funding) - min(funding) < 1e-7
for a, b in zip(arms["history"]["history"], arms["instantaneous"]["history"]):
    assert a["context"] == b["context"]
    assert a["balance"] == b["balance"]
    assert a["raw_total"] == b["raw_total"]
    assert list(map(abs,a["features"])) == list(map(abs,b["features"]))
baseline = {r["context"]: statistics.mean(r["reward"]) for r in arms["no_credit"]["probes"] if not r["swap"]}
compact = {k: data[k] for k in ("schema","seed","noise_seed","episodes","rate","contexts","scope","features","schedule","source","cellular_source","driver_sha256","elapsed_seconds")}
compact["critic_scale"] = data.get("critic_scale", 1.0)
compact["validation"] = dict(max_norm_error=max_norm_error, max_closure=max_closure, funded_spread=max(funding)-min(funding))
compact["arms"] = {}
for name, arm in arms.items():
    probes = [{**r, "gain": statistics.mean(r["reward"])-baseline[r["context"]]} for r in arm["probes"]]
    compact["arms"][name] = dict(probes=probes, learned=arm["learned"], total_learning_work=arm["total_learning_work"], total_update_l1=arm["total_update_l1"])
    print(name, {"training_gain":statistics.mean(r["gain"] for r in probes if r["context"]<=2 and not r["swap"]), "heldout_gain":statistics.mean(r["gain"] for r in probes if r["context"]>2 and not r["swap"]), "swap_gain":statistics.mean(r["gain"] for r in probes if r["context"]<=2 and r["swap"]), "paid":arm["total_learning_work"]})
destination = root / sys.argv[2]
destination.write_text(json.dumps(compact, indent=2)+"\n")
print(compact["validation"])
