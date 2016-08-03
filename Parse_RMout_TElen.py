import sys
import os
import statistics
file = sys.argv[1]
# file = "/Users/xiaoyu/project/mm10_TE_RM/mm10_RM.fa.out"
repdict = {}  # like empty hash table in perl
infile = open(file, "r")
line = infile.readline()
for line in infile:
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
        RepLeft = int(line[13].translate({ord('('): None, ord(')'): None}))
        if strand == 'C':
            strand = "-"
            RepStart = line[13]
            RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))
        length = end - start + 1
        RepLength = RepLeft + RepEnd
        if name in repdict:
            if TEclass != repdict[name]['class']:
                print("something wrong here with %s! %s or %s?\n\n" % (name, TEclass, repdict[name]['class']))
            repdict[name]['length'].append(RepLength)  # append value to the existing list in the dict of dict
        else:
            repdict[name] = {'length': [RepLength]}  # create the dict of dict
            repdict[name]['class'] = TEclass
# the constucted dict structure is
# {repname1:{'length':[length1],'div':[divergence1],...}, repname2:{'length':[length2],'div':[divergence2],...},...}
infile.close()

for name in repdict:
    repdict[name]['median'] = statistics.median(repdict[name]['length'])
    repdict[name]['same'] = 0
    repdict[name]['diff'] = 0
    for i in repdict[name]['length']:
        if (i == repdict[name]['median']):
            repdict[name]['same'] += 1
        else:
            repdict[name]['diff'] += 1

outfile = open(file+'.length', 'w')
for repname in repdict:
    outfile.write("%s\t%s\t%d\t%d\t%d\n" % (repname, repdict[repname]['class'], repdict[repname]['median'], repdict[repname]['same'], repdict[repname]['diff']))

outfile.close()
exit()
