#!/usr/bin/env python3
# make_chunks.py — regenerate pau_chunk_{a,b,c,d}.dat from pau200_lib.dat (committed).
#
# Provenance (2026-07-10 session): the four chunk files fed the G4 library jobs
# pauL_{a,b,c,d} (75 events each -> raw_lib_pau{La,Lb,Lc,Ld}.root, 300 collisions total).
#   a = even events of pau200_lib.dat (150), b = odd events (150)
#   c = events 76-150 of a, d = events 76-150 of b
#     (Fun4All reads sequentially; wave 1 consumed the first 75 of a/b, so the
#      remaining 75 were staged as separate files for jobs c/d.)
# The split logic below is VERBATIM the original session code, so output is
# byte-identical to the originals (sha256 in island_post/production_manifest.txt).
# pau200_lib.dat itself: sHIJING, pau_lib.xml, SEED 20260713, 300 events.

import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def read_events(fn):
    head, events, cur = [], [], None
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
        else:
            cur.append(line)
    if cur is not None:
        events.append(cur)
    return head, events


foot = ['HepMC::IO_GenEvent-END_EVENT_LISTING\n']

# stage 1: even/odd interleave of the 300-event library
h, ev = read_events('pau200_lib.dat')
for name, ch in zip('ab', [ev[0::2], ev[1::2]]):
    with open(f'pau_chunk_{name}.dat', 'w') as o:
        o.writelines(h)
        [o.writelines(e) for e in ch]
        o.writelines(foot)
    print(f'pau_chunk_{name}.dat: {len(ch)} events')

# stage 2: library extension = events 76-150 of each chunk, as new files
for src, dst in [('pau_chunk_a.dat', 'pau_chunk_c.dat'), ('pau_chunk_b.dat', 'pau_chunk_d.dat')]:
    h, ev = read_events(src)
    with open(dst, 'w') as o:
        o.writelines(h)
        [o.writelines(e) for e in ev[75:]]
        o.writelines(foot)
    print(f'{dst}: {len(ev[75:])} events')
