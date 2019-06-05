#!/anaconda3/bin/python

import glob
from collections import OrderedDict

ab1s = glob.glob("*.ab1")



volunteers = OrderedDict()

for filename in ab1s:
    fields = filename.split("_")
    volunteer = str(fields[0])
    timepoint = str(fields[1])

    if volunteer not in volunteers:
        volunteers[volunteer] = [timepoint]
    else:
        volunteers[volunteer].append(timepoint)

print("### Summary of the ab1 files in this directory")
print (len(ab1s),"sequences were found")
print("")
for volunteer in volunteers:
    print (len(glob.glob(volunteer+"*.ab1")), "sequences for", "GW volunteer" if "G00158" in volunteer else "FH volunteer", volunteer, "at timepoints", " ".join(set(volunteers[volunteer])))

