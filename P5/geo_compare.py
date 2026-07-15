#!/usr/bin/env python3
# P5 geometry audit comparator: constructed (container) vs GDML-loaded (native).
import sys, collections

def load(fn):
    agg = collections.defaultdict(lambda: [0, 0.0])  # (name, material) -> [count, sum cv]
    for line in open(fn):
        if not line.startswith("GEOAUD ") or line.startswith("GEOAUD-END"):
            continue
        body = line[7:].strip()
        try:
            name, mat, cv = [x.strip() for x in body.split("|")]
        except ValueError:
            continue
        k = (name, mat)
        agg[k][0] += 1
        try:
            agg[k][1] += float(cv)
        except ValueError:
            pass
    return agg

ref = load("geo_ref.txt")
sa = load("geo_sa.txt")
only_ref = sorted(set(ref) - set(sa))
only_sa = sorted(set(sa) - set(ref))
print(f"GEOCMP keys: ref {len(ref)} sa {len(sa)} | only-ref {len(only_ref)} only-sa {len(only_sa)}")
for k in only_ref[:12]:
    print(f"GEOCMP only-REF: {k[0]} [{k[1]}] x{ref[k][0]}")
for k in only_sa[:12]:
    print(f"GEOCMP only-SA:  {k[0]} [{k[1]}] x{sa[k][0]}")
ncnt = nvol = 0
worst = []
for k in sorted(set(ref) & set(sa)):
    rc, rv = ref[k]
    sc, sv = sa[k]
    if rc != sc:
        ncnt += 1
        if ncnt <= 8:
            print(f"GEOCMP count-mismatch: {k[0]} [{k[1]}] ref x{rc} sa x{sc}")
    if rv > 0 and sv > 0:
        rel = abs(rv - sv) / max(rv, sv)
        if rel > 0.02:
            nvol += 1
            worst.append((rel, k, rv, sv))
worst.sort(reverse=True)
for rel, k, rv, sv in worst[:8]:
    print(f"GEOCMP volume-diff {rel*100:.1f}%: {k[0]} [{k[1]}] ref {rv:.4g} sa {sv:.4g} cm3")
print(f"GEOCMP SUMMARY: shared {len(set(ref) & set(sa))} | count-mismatch {ncnt} | volume>2% {nvol}")
