#!/usr/bin/env python3
"""Summarise sim_results.csv per scenario into sim_summary.csv (mean and SD across replicates)."""
import csv, numpy as np
rows = list(csv.DictReader(open("sim_results.csv")))
keys = ["S4","M1","M2","S3","K1","K2","lost"]
out = []
for sc in ["constant","expansion","contraction"]:
    r = [x for x in rows if x["scenario"] == sc]
    v = {k: np.array([float(x[k]) for x in r]) for k in keys + ["M2_to_K1","M1_lost","M1_to_K2"]}
    rec = {"scenario": sc, "reps": len(r),
           "D4_mean": round(float(np.nanmean([float(x["D4"]) for x in r])), 3),
           "lost_frac": round(v["lost"].sum()/v["S4"].sum(), 4),
           "M2_to_K1_frac": round(v["M2_to_K1"].sum()/v["M2"].sum(), 4),
           "M1_lost_frac": round(v["M1_lost"].sum()/v["M1"].sum(), 4),
           "M1_to_K2_frac": round(v["M1_to_K2"].sum()/v["M1"].sum(), 4),
           "K1_share_3seq": round(v["K1"].sum()/v["S3"].sum(), 4)}
    for k in keys:
        rec[k+"_mean"] = round(float(v[k].mean()), 2); rec[k+"_sd"] = round(float(v[k].std()), 2)
    out.append(rec)
with open("sim_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
