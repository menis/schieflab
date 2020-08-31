import sys
import csv

if '-f' in sys.argv:
    TargetFile = sys.argv[sys.argv.index('-f') + 1]
    with open(TargetFile, 'r') as f:
        rows = csv.reader(f, delimiter=',')
        uniqs = []
        first = True
        for row in rows:
            if first:
                uniqs = [""] * len(row)
                first = False
            i = 0
            for pos in row:
                if (i == 0):
                    # first pos is the name of the mAb
                    pass
                elif(pos.strip() != "~"):
                    for aa in pos:
                        if(aa not in uniqs[i]):
                            uniqs[i] = uniqs[i] + aa.upper()
                i = i + 1

        print (','.join(uniqs))

