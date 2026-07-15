#!/usr/bin/env python3
# v4.0: thin fit on the full census (fired target 8200), accept list + mbd file
import math, random
kept = {}
for l in open('v40_census.txt'):
    f, e, k = l.split()
    kept[(int(f), int(e))] = max(1.0, float(k))
fired = {}
for i in range(10):
    for l in open(f'/home/rog/sPHENIX/3D_ClusterFindingML/P5/angantyr/fired_{i}.txt'):
        p = l.split()
        if p[0] == 'FIRED': fired[(i, int(p[1]) - 1)] = int(p[2])
keys = sorted(kept)
ks = [kept[k] for k in keys]
fs = [fired.get(k, 0) for k in keys]
a = 3.0
lo, hi = 100.0, 60000.0
for _ in range(60):
    k0 = 0.5 * (lo + hi)
    w = [min(1.0, (k0 / k) ** a) for k in ks]
    mf = sum(wi * ki for wi, ki, fi in zip(w, ks, fs) if fi) / max(1e-9, sum(wi for wi, fi in zip(w, fs) if fi))
    if mf > 8200: hi = k0
    else: lo = k0
rng = random.Random(20260717)
acc = {}
for key, k in zip(keys, ks):
    if rng.random() < min(1.0, (k0 / k) ** a):
        acc.setdefault(key[0], []).append(key[1])
tot = sum(len(v) for v in acc.values())
print(f"THINFIT alpha 3.0 k0 {k0:.0f} | accepted {tot}/{len(ks)}")
with open('v40_accept_list.txt', 'w') as fo:
    for f in sorted(acc):
        for j, e in enumerate(acc[f]):
            fo.write(f"{f} {e} {j} {fired.get((f,e),0)}\n")
with open('v40_mbd.txt', 'w') as fo:
    fo.write("# file event_index_in_file north south fired\n")
    for f in sorted(acc):
        for j, e in enumerate(acc[f]):
            fl = fired.get((f, e), 0)
            fo.write(f"ang_run_{f}.dat {j} {fl} {fl} {fl}\n")
