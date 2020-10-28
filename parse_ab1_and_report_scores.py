import sys, os
import csv
from Bio.Seq import Seq
from Bio import SeqIO
import statistics
import numpy as np
import glob
from io import StringIO

# taken from https://stackoverflow.com/a/43099751
def geo_mean_overflow(iterable):
    a = np.log(iterable)
    return np.exp(a.sum()/len(a))

def fuzzyfind(needle, haystack):
    n = len(needle)
    h = len(haystack)
    for position in range(h):
        if (position+n) > h:
            return -1
        else:
            trim = haystack[position:position+n]
            if(len(trim) != len(needle)):
                exit("length mismatch")

            #Compare
            equal = True
            for i in range(len(needle)):
                n_c = needle[i]
                h_c = trim[i]
                if n_c.lower() == 'n' or h_c.lower() == 'n':
                    continue
                else:
                    if n_c.lower() != h_c.lower():
                        equal = False
                        break
            if equal:
                return position

    #if you made it this far there are no matches
    return -1




stats = {}
fastas = {}
if '-ab1s' in sys.argv:
    TargetPath = sys.argv[sys.argv.index('-ab1s') + 1]
    files = glob.glob(TargetPath+"/**/*.ab1", recursive=True)

    for ab1 in files:
        handle = StringIO("")
        handle_qual = StringIO("")
        of = os.path.basename(ab1)
        oroot = os.path.dirname(ab1)

        # convert the file
        dump = SeqIO.convert(os.path.join(
            oroot, of), "abi-trim", handle, "fasta-2line")
        dump = SeqIO.convert(os.path.join(
            oroot, of), "abi-trim", handle_qual, "qual")

        dump = handle.seek(0)
        dump = handle_qual.seek(0)
        seqid_full = handle.readline().replace(">", "").replace("\n", "")
        dump = handle_qual.readline()
        sequence = handle.readline().replace("\n", "")
        qualdata = ""
        for line in handle_qual:
            qualdata += line.replace("\n", " ")

        #recast to an array
        qualdata = [int(q) for q in qualdata.split()]

        if seqid_full not in fastas:
            fastas[seqid_full] = sequence
            stats[seqid_full] = qualdata
        else:
            exit("Duplicate ab1s found. Names must be unique.")
else:
    exit("Please specify the path to the ab1 files using the -f flag.")

if '-filetoupdate' in sys.argv:
    TargetFile = sys.argv[sys.argv.index('-filetoupdate') + 1]
    basename = os.path.splitext(os.path.basename(TargetFile))[0]
    headerfound = False
    headers = ""
    with open(basename+"_w_stats.csv", 'w') as of:
        spamwriter = csv.writer(of)
        with open(TargetFile) as f:
            reader = csv.reader(f)
            for row in reader:
                if not headerfound:
                    headers = row
                    headers.append("H_Length")
                    headers.append("H_Mean")
                    headers.append("H_Geomean")
                    headers.append("H_below20")
                    headers.append("H_below30")
                    headers.append("H_below40")
                    headers.append("L_Length")
                    headers.append("L_Mean")
                    headers.append("L_Geomean")
                    headers.append("L_below20")
                    headers.append("L_below30")
                    headers.append("L_below40")

                    headerfound = True
                    spamwriter.writerow(headers)
                else:
                    H_below20 = ""
                    H_below30 = ""
                    H_below40 = ""
                    H_length = ""
                    H_mean = ""
                    H_geomean = ""
                    L_below20 = ""
                    L_below30 = ""
                    L_below40 = ""
                    L_length = ""
                    L_mean = ""
                    L_geomean = ""

                    for chain in ["heavy", "light"]:
                        seq_id_to_match = ""
                        vdj = ""
                        if chain == "heavy":
                            if row[7] != "":
                                seq_id_to_match = row[0]
                                vdj = row[32]
                            else:
                                continue
                        else:
                            if row[41] != "":
                                seq_id_to_match = row[0]
                                vdj = row[62]
                            else:
                                continue

                        # Calculate offsets
                        # needle  = pairs_and_unpairs[seq_id_to_match]
                        needle = vdj
                        haystack = []
                        for k, v in fastas.items():
                            if k.find(seq_id_to_match) > -1:
                                haystack.append((k, v))
                        if len(haystack) == 0:
                            print("Seq is not found in fasta file for " + seq_id_to_match)
                            spamwriter.writerow(row)
                            continue

                        start = -1
                        end = -1
                        keytouse = ""
                        flip = False
                        for k, v in haystack:
                            start = fuzzyfind(needle, v)
                            # Some ab1 are already reverse build in a hack to
                            # check the other orientation. This is something strange.
                            if start == -1:
                                fastaseqflipped = str(Seq(v).reverse_complement())
                                start = fuzzyfind(needle, fastaseqflipped)
                                if start > -1:
                                    flip = True

                            if start > -1:
                                end = start + len(needle)
                                keytouse = k
                                break

                        if start < 0:
                            print("Seq in VDJ nt doesn't match any in the fasta file for " + seq_id_to_match)
                            spamwriter.writerow(row)
                            continue

                        statsdatatoslice = stats[keytouse]
                        if not flip:
                            focusarea = statsdatatoslice[start:end]
                        else:
                            focusarea = statsdatatoslice[::-1][start:end]

                        if len(focusarea) < 1:
                            print("Can't find the VDJ nt in the stats file for " + seq_id_to_match)
                            spamwriter.writerow(row)
                            continue

                        # Generate stats
                        if chain == "heavy":
                            H_below20 = len([s for s in focusarea if s < 20])
                            H_below30 = len([s for s in focusarea if s < 30])
                            H_below40 = len([s for s in focusarea if s < 40])
                            H_length = len(focusarea)
                            H_mean = statistics.mean(focusarea)
                            H_geomean = geo_mean_overflow(focusarea)
                        else:
                            L_below20 = len([s for s in focusarea if s < 20])
                            L_below30 = len([s for s in focusarea if s < 30])
                            L_below40 = len([s for s in focusarea if s < 40])
                            L_length = len(focusarea)
                            L_mean = statistics.mean(focusarea)
                            L_geomean = geo_mean_overflow(focusarea)
                    # Make a new line and save
                    appended_line = row + [H_length, H_mean, H_geomean, H_below20, H_below30, H_below40] + [L_length,
                                                                                                            L_mean,
                                                                                                            L_geomean,
                                                                                                            L_below20,
                                                                                                            L_below30,
                                                                                                            L_below40]
                    spamwriter.writerow(appended_line)
else:
    exit("Please specify the file to update. Use the -filetoupdate flag.")

