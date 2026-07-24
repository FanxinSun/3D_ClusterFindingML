#!/usr/bin/env python3
# env_solve.py — envelope node solver (v5.1 protocol, scripted 2026-07-23).
# f_k = (real_k - trig_k) / (flat_k - trig_k) per coarse bin; nodes 0..7 from
# bins 0..7, node 8 = avg(bin8, bin9); normalized to unit MEAN over the 9 node
# values (envelope = shape only; content scale stays visible in acceptance).
# usage: env_solve.py "<real 10 vals>" "<trig 10 vals>" "<flat 10 vals>"
# prints the composer ENV string. Validated against the hand-derived ENV51.
import sys

real = [float(x) for x in sys.argv[1].split()]
trig = [float(x) for x in sys.argv[2].split()]
flat = [float(x) for x in sys.argv[3].split()]
f = [(r - t) / (fl - t) for r, t, fl in zip(real, trig, flat)]
nodes = f[:8] + [0.5 * (f[8] + f[9])]
m = sum(nodes) / len(nodes)
nodes = [n / m for n in nodes]
ts = ["5.80", "11.10", "16.40", "21.70", "27.00", "32.30", "37.60", "42.90", "47.96"]
print(",".join(f"{t}:{n:.3f}" for t, n in zip(ts, nodes)))
