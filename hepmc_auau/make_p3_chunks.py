#!/usr/bin/env python3
# make_p3_chunks.py — regenerate p3_chunk_{a,b,c}.dat from the two committed
# AuAu HepMC files (auau200_signal.dat + auau200_pileup.dat, 95 events total).
#
# Provenance (2026-07-09 session, P3 wave 1 — pre-species-pivot): the three
# chunks fed the parallel single-collision G4 jobs p3s_{a,b,c} (15 events each)
# that validated the frame-composer architecture on AuAu, before run 79507 was
# re-identified as p+Au. Interleaved split: a = events 0,3,6,..., b = 1,4,7,...,
# c = 2,5,8,... of signal-then-pileup concatenation (32/32/31 events).
# Logic below is VERBATIM the original session code -> byte-identical output
# (sha256 in island_post/production_manifest.txt).

import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def read_events(fn):
    head, events, foot = [], [], []
    cur = None
    for line in open(fn):
        if line.startswith('E '):
            if cur is not None:
                events.append(cur)
            cur = [line]
        elif cur is None:
            head.append(line)
        elif line.startswith('HepMC::IO_GenEvent-END'):
            events.append(cur)
            cur = None
            foot.append(line)
        else:
            cur.append(line)
    if cur is not None:
        events.append(cur)
    if not foot:
        foot = ['HepMC::IO_GenEvent-END_EVENT_LISTING\n']
    return head, events, foot


h1, e1, f1 = read_events('auau200_signal.dat')
h2, e2, f2 = read_events('auau200_pileup.dat')
allev = e1 + e2
print(f'total events available: {len(allev)}')
chunks = [allev[0::3], allev[1::3], allev[2::3]]
for i, (name, ch) in enumerate(zip('abc', chunks)):
    with open(f'p3_chunk_{name}.dat', 'w') as o:
        o.writelines(h1)
        for ev in ch:
            o.writelines(ev)
        o.writelines(f1)
    print(f'p3_chunk_{name}.dat: {len(ch)} events')
