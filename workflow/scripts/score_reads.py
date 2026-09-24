#!/usr/bin/env python3
"""Emit per-read length, mean quality and window quality, matching Filtlong."""
import sys

WINDOW = 250

def mean_quality(q):
    return 100.0 * sum(q) / len(q)

def window_quality(q, w=WINDOW):
    if len(q) <= w:
        return mean_quality(q)
    s = sum(q[:w])
    lo = s
    for j in range(w, len(q)):
        s += q[j] - q[j - w]
        lo = min(lo, s)
    lo /= w
    return 0.0 if lo < 0.5 / w else 100.0 * lo

print("read_id\tlength\tmean_q\twindow_q")
with open(sys.argv[1]) as fh:
    for i, line in enumerate(fh):
        if i % 4 == 0:
            name = line[1:].split()[0]
        elif i % 4 == 3:
            q = [1.0 - 10.0 ** (-(ord(c) - 33) / 10.0) for c in line.rstrip("\n")]
            print(f"{name}\t{len(q)}\t{mean_quality(q):.6f}\t{window_quality(q):.6f}")
