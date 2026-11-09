#!/usr/bin/env python3
"""Pull the DT tables out of a NetPerturb HTML report into TSV files.

Usage: extract_report_tables.py <netperturb_report.html> <out_dir>

Writes cell_counts.tsv (8.1 Cells per identity), identity_qc.tsv (9.1) and
network_qc.tsv (9.2), dq_flags.tsv and runtime.tsv when the report carries them, and lists every other
table it found so nothing is silently skipped.
"""
import html, json, os, re, sys

src, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
text = open(src, encoding="utf8", errors="ignore").read()

wanted = {"Cells in input": "cell_counts.tsv", "median_counts": "identity_qc.tsv", "hub_share_top1pct": "network_qc.tsv",
          "Flagged at": "dq_flags.tsv", "Cached": "runtime.tsv"}
for m in re.finditer(r'<script type="application/json" data-for="[^"]+">(.*?)</script>', text, re.S):
    try:
        x = json.loads(m.group(1))["x"]
    except (ValueError, KeyError):
        continue
    if not isinstance(x, dict) or "container" not in x or "data" not in x:
        continue
    cols = [html.unescape(re.sub("<[^>]+>", "", c)) for c in re.findall(r"<th>(.*?)</th>", x["container"])]
    data = x["data"]
    if len(cols) != len(data):
        continue
    name = next((f for k, f in wanted.items() if k in cols), None)
    if name is None:
        print("skipped table:", ", ".join(cols[:6]), "...", file=sys.stderr)
        continue
    fmt = lambda v: "" if v is None else str(v)
    with open(os.path.join(out_dir, name), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in zip(*data):
            fh.write("\t".join(fmt(v) for v in row) + "\n")
    print("wrote", name, len(data[0]), "rows")
