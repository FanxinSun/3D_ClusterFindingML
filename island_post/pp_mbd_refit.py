#!/usr/bin/env python3
# pp_mbd_refit.py — MBD-proxy floor calibration (v5.0, 2026-07-22).
#
# STATUS: PARKED (2026-07-22). Anchor-fit converges at floor f=0.529 ->
# eps 0.749, UNPHYSICAL for pp MBD (collaboration ~0.55): the flat-floor
# model buys bias reduction only by overpaying efficiency. Also the fired
# event is ~1 of ~50 collisions per window -> its +15% excess is ~0.5% of
# window content; the window-critical quantity (uniform mean) matches real
# at +1.9% with NO refit. Production therefore runs the GEOMETRIC flags +
# physical RSPEC (rmbd deciles / 0.519); the fired-mean +15.2% is a
# DECLARED residual measured at the step anchor. Revisit only if the
# acceptance suite shows step/window tension that a bias model would fix.
#
# Gate finding: raw geometric proxy (>=1 charged in 3.51<|eta|<4.61, both
# arms) gives fired-mean 9445 kept px (+15.2% vs the rate-free real anchor
# 8200) while the UNIFORM mean matches real to +1.9% — the excess is trigger
# BIAS (x1.285 proxy vs x1.137 real), not generator multiplicity. The real
# MBD fires on softer events (secondary feed-in to the arms); model this as a
# per-arm firing floor f: P(arm) = 1 if n_charged>=1 else f. One parameter,
# fit so the EXPECTED fired-mean equals the anchor; flags then realized with
# a seeded RNG into pp_mbd.txt (fired column only; n/s columns untouched, so
# the refit is repeatable/reversible from the same file).
import random, sys

ANCHOR = 8200.0
SEED = 20260722

rows = []   # (file, idx, nN, nS, line-of-header?) preserved order
for L in open("pp_mbd.txt"):
    if L.startswith("#"):
        continue
    f_, e, n, s, fl = L.split()
    rows.append([f_, int(e), int(n), int(s)])
cnt = {}
for L in open("pp_census.txt"):
    c, e, n = L.split()
    cnt[(c, int(e))] = int(n)
px = [cnt.get((r[0].split("_")[-1].split(".")[0], r[1]), 0) for r in rows]

def stats(f):
    sp = spx = 0.0
    for r, p in zip(rows, px):
        P = (1.0 if r[2] > 0 else f) * (1.0 if r[3] > 0 else f)
        sp += P; spx += P * p
    return sp / len(rows), spx / sp   # eps, expected fired-mean

lo, hi = 0.0, 1.0
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if stats(mid)[1] > ANCHOR: lo = mid   # more floor -> softer sample -> lower mean
    else: hi = mid
f = 0.5 * (lo + hi)
eps, fm = stats(f)
u = sum(px) / len(px)
print(f"floor f = {f:.4f} | eps {eps:.4f} (geometric was 0.519) | "
      f"expected fired-mean {fm:.0f} vs anchor {ANCHOR:.0f} | bias x{fm/u:.3f} (real x1.137)")

rng = random.Random(SEED)
nf = 0
out = ["# file event_index_in_file north south fired  "
       f"(fired = floor-refit f={f:.4f} seed {SEED}, anchor {ANCHOR:.0f}; see pp_mbd_refit.py)\n"]
for r in rows:
    fired = int((r[2] > 0 or rng.random() < f) and (r[3] > 0 or rng.random() < f))
    nf += fired
    out.append(f"{r[0]} {r[1]} {r[2]} {r[3]} {fired}\n")
open("pp_mbd.txt", "w").writelines(out)
rf = [p for r, p, L in zip(rows, px, out[1:]) if int(L.split()[4])]
print(f"realized: fired {nf}/{len(rows)} ({nf/len(rows):.4f}) | "
      f"fired-mean {sum(rf)/max(1,len(rf)):.0f} | pp_mbd.txt rewritten")
print("RSPEC deciles (real rmbd deciles / eps): " +
      ",".join(f"{d/eps:.0f}:1" for d in (320, 341, 350, 352, 353, 354, 355, 356, 358, 362)))
