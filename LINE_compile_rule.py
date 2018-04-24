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
        if re.match(r'\d+\s', line):
            line = line.split()
            if len(line) < 14:
                print("RMalin file split wrong?")
            chr = line[4]
            start = int(line[5])
            end = int(line[6])
            score = float(line[1])

            # for Repeatmasker open 4.0+:
            ID = line[-1]
            if line[8] == 'C':
                strand = line[8]
                nameclass = line[9]
                names = nameclass.split('#')
                name = names[0]
                TEclass = names[1]
                RepStart = line[12]
                RepEnd = int(line[11])
                RepLeft = int(line[10].translate({ord('('): None, ord(')'): None}))
            else:
                strand = '+'
                nameclass = line[8]
                names = nameclass.split('#')
                name = names[0]
                TEclass = names[1]
                RepStart = line[9]
                RepEnd = int(line[10])
                RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))

            # for Repeatmasker open 3.0:
            # strand = line[8]
            # name = line[9]
            # TEclass = line[10]
            # RepEnd = int(line[12])
            # if strand == '+':
            #     RepStart = line[11]
            #     RepLeft = int(line[13].translate({ord('('): None, ord(')'): None}))
            # if strand == 'C':
            #     RepStart = line[13]
            #     RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))


            length = end - start + 1
            RepLength = RepLeft + RepEnd
            if re.match(r'LINE', TEclass):  # only LINEs are included here

                # for Repeatmakser 4.0+, we can use ID, but ID is not available for 3.0.
                # for 4.0:
                # if chr in chrtree:
                #     chrtree[chr][start:end] = (name, ID)
                # else:
                #     chrtree[chr] = IntervalTree()
                #     chrtree[chr][start:end] = (name, ID)

                # for 3.0:
                if chr in chrtree:
                    chrtree[chr][start:end] = (name, TEclass, strand)
                else:
                    chrtree[chr] = IntervalTree()
                    chrtree[chr][start:end] = (name, TEclass, strand)

print("done with RMalignment!")

LINEs = {}
with open(outfile, "r") as RMout:
    for line in RMout:
        line = line.split()
        if not (line == [] or line[0] == 'SW' or line[0] == 'score'):
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
                RepStart = line[13]
                RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))
            length = end - start + 1
            RepLength = RepLeft + RepEnd
            if re.match(r'LINE', TEclass):  # only LINEs are included here
                if chr in chrtree:
                    for frag in chrtree[chr].search(start, end):  # envelope search here.
                        # if frag.data[1] == ID: #  can't use ID for RM3.0
                        if frag.data[1] == TEclass and frag.data[2] == strand:  # use strand and TEclass as condition for RM3.0:
                            if name not in LINEs:
                                LINEs[name] = {}
                            if frag.data[0][-5:] == '_5end':
                                if "5end" not in LINEs[name]:
                                    LINEs[name]["5end"] = {}
                                if frag.data[0] not in LINEs[name]["5end"]:
                                    LINEs[name]["5end"][frag.data[0]] = {}
                                    LINEs[name]["5end"][frag.data[0]]['cp'] = 1
                                    LINEs[name]["5end"][frag.data[0]]['len'] = length
                                else:
                                    LINEs[name]["5end"][frag.data[0]]['cp'] += 1
                                    LINEs[name]["5end"][frag.data[0]]['len'] += length
                            elif frag.data[0][-5:] == '_orf2':
                                if "orf2" not in LINEs[name]:
                                    LINEs[name]["orf2"] = {}
                                if frag.data[0] not in LINEs[name]["orf2"]:
                                    LINEs[name]["orf2"][frag.data[0]] = {}
                                    LINEs[name]["orf2"][frag.data[0]]['cp'] = 1
                                    LINEs[name]["orf2"][frag.data[0]]['len'] = length
                                else:
                                    LINEs[name]["orf2"][frag.data[0]]['cp'] += 1
                                    LINEs[name]["orf2"][frag.data[0]]['len'] += length
                            elif frag.data[0][-5:] == '_3end':
                                if "3end" not in LINEs[name]:
                                    LINEs[name]["3end"] = {}
                                if frag.data[0] not in LINEs[name]["3end"]:
                                    LINEs[name]["3end"][frag.data[0]] = {}
                                    LINEs[name]["3end"][frag.data[0]]['cp'] = 1
                                    LINEs[name]["3end"][frag.data[0]]['len'] = length
                                else:
                                    LINEs[name]["3end"][frag.data[0]]['cp'] += 1
                                    LINEs[name]["3end"][frag.data[0]]['len'] += length
                            else:
                                if "other" not in LINEs[name]:
                                    LINEs[name]["other"] = {}
                                if frag.data[0] not in LINEs[name]["other"]:
                                    LINEs[name]["other"][frag.data[0]] = {}
                                    LINEs[name]["other"][frag.data[0]]['cp'] = 1
                                    LINEs[name]["other"][frag.data[0]]['len'] = length
                                else:
                                    LINEs[name]["other"][frag.data[0]]['cp'] += 1
                                    LINEs[name]["other"][frag.data[0]]['len'] += length
    print("done with RM output!")

#print(json.dumps(LINEs, sort_keys=True, indent=4))

sortedLINEs = {}
for name in LINEs:
    sortedLINEs[name] = {}
    for frag in LINEs[name]:
        sortedLINEs[name][frag] = {}
        for te in LINEs[name][frag]:
            if len(sortedLINEs[name][frag]) == 0:
                sortedLINEs[name][frag]['te'] = te
                sortedLINEs[name][frag]['cp'] = LINEs[name][frag][te]['cp']
                sortedLINEs[name][frag]['len'] = LINEs[name][frag][te]['len']
            elif LINEs[name][frag][te]['len'] > sortedLINEs[name][frag]['len']:
                sortedLINEs[name][frag]['te'] = te
                sortedLINEs[name][frag]['cp'] = LINEs[name][frag][te]['cp']
                sortedLINEs[name][frag]['len'] = LINEs[name][frag][te]['len']

print(json.dumps(sortedLINEs, sort_keys=True, indent=4))
