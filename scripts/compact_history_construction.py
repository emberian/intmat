#!/usr/bin/env python3
"""Validate and compact the bounded history-allocation study, no simulation."""
import hashlib
import json
import statistics
import sys
from pathlib import Path


def read(path):
    data = Path(path).read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def compact(path):
    report, sha = read(path)
    arms = report["arms"]
    unsigned = arms["unsigned_history"]["history"]
    signed = arms["history"]["history"]
    assert len(signed) == len(unsigned) == report["episodes"]
    for first, second in zip(signed, unsigned):
        assert first["trial"] == second["trial"] and first["context"] == second["context"]
        assert [abs(x) for x in first["feature"]] == second["feature"]
    baseline = {(p["context"], p["decision"]): statistics.mean(p["reward"])
                for p in arms["no_credit"]["probes"] if not p["swap"]}
    result = {}
    for name, arm in arms.items():
        probes = []
        for p in arm["probes"]:
            mean = statistics.mean(p["reward"])
            probes.append({"context": p["context"], "decision": p["decision"], "swap": p["swap"],
                           "mean_return": mean, "per_resident_return": p["reward"],
                           "delta_vs_matched_no_credit": mean - baseline[p["context"], p["decision"]],
                           "feature": p["feature"], "drive": p["drive"], "funded": p["funded"]})
        result[name] = {"probes": probes, "learned_column_major_4_by_2": arm["learned"],
                        "total_learning_work": arm["total_learning_work"],
                        "total_update_l1": arm["total_update_l1"]}
    records = [p for a in arms.values() for key in ("history", "probes") for p in a[key]]
    residual = max(abs(p[k]) for p in records for k in ("material", "potential"))
    assert residual < 1e-10
    # Allocation changes composition, with matched total paid construction.
    totals = [sum(p["funded"]) for p in records]
    spread = max(totals) - min(totals)
    assert spread < 1e-7
    return {"seed": report["seed"], "noise_seed": report["noise_seed"], "episodes": report["episodes"],
            "contexts": report["contexts"], "rate": report["rate"], "source_report_sha256": sha,
            "driver_sha256": report["driver_sha256"], "source": report["source"],
            "cellular_source": report["cellular_source"], "exact_paired_feature_norm_match": True,
            "funded_construction_total_spread": spread, "max_conservation_residual": residual,
            "elapsed_seconds_excluding_setup": report["elapsed_seconds"], "arms": result}


if __name__ == "__main__":
    destination, *sources = sys.argv[1:]
    assert sources
    result = {"schema": "history-dependent-construction-compact-v1", "runs": [compact(p) for p in sources]}
    Path(destination).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"runs": len(result["runs"]), "noise_seeds": [r["noise_seed"] for r in result["runs"]]}))
