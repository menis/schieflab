#!/anaconda3/bin/python

import sys

remap = {}
with open(sys.argv[2]) as haystack:
    for line in haystack:
        parts = line.strip().split()
        if len(parts) > 2:
            remap[parts[0]] = parts[2]

with open(sys.argv[1]) as needle:
    for line in needle:
        l = line.strip()
        if l in remap:
            print(remap[l])
        elif len(l) > 0:
            print("UNKNOWN")
