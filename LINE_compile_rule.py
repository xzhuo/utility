# Find out how LINEs are merged in repeatmasker

import sys
from intervaltree import IntervalTree, Interval
import re
import json

outfile = sys.argv[1]  # Repeatmasker.org output
alignfile = sys.argv[2]  # Repeatmasker.org align file
chrtree = {}

# put RMalign entries on the intervaltree first:
with open(alignfile, "r") as RMalign:
    for line in RMalign:
        if re.match(r'\d+\s',line):
            line = line.split()
            if not len(line) == 15:
                print("RMalin file split wrong?")
            chr = line[4]
            start = int(line[5])
            end = int(line[6])
            nameclass = line[9]
            names = nameclass.split('#')
            name = names[0]
            TEclass = names[1]
            RepEnd = int(line[11])
            score = float(line[1])
            strand = line[8]
            RepStart = line[10]
            ID = line[-1]
            RepLeft = int(line[12].translate({ord('('): None, ord(')'): None}))
            if strand == 'C':
                strand = "-"
                RepStart = line[12]
                RepLeft = int(line[10].translate({ord('('): None, ord(')'): None}))
            length = end - start + 1
            RepLength = RepLeft + RepEnd
            if re.match(r'LINE',TEclass):  # only LINEs are included here
                if chr in chrtree:
                    chrtree[chr][start:end] = (name, ID)
                else:
                    chrtree[chr] = IntervalTree()
                    chrtree[chr][start:end] = (name, ID)


LINEs = {}
with open(outfile, "r") as RMout:
    for line in RMout:
        if not (line == [] or line[:5] == '   SW' or line[:5] == 'score'):
            line = line.split()
            chr = line[4]
            start = int(line[5])
            end = int(line[6])
            name = line[9]
            TEclass = line[10]
            RepEnd = int(line[12])
            score = float(line[1])
            strand = line[8]
            RepStart = line[11]
            ID = line[-1]
            RepLeft = int(line[13].translate({ord('('): None, ord(')'): None}))
            if strand == 'C':
                strand = "-"
                RepStart = line[13]
                RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))
            length = end - start + 1
            RepLength = RepLeft + RepEnd
            if re.match(r'LINE', TEclass):  # only LINEs are included here
                for frag in chrtree[chr][start:end]:
                    if frag.data[1] == ID:
                        if name not in LINEs:
                            LINEs[name] = {}
                        if frag.data[0] not in LINEs[name]:
                            LINEs[name][frag.data[0]] = 1
                        else:
                            LINEs[name][frag.data[0]] += 1


print(json.dumps(LINEs))
