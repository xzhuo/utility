# Find out how LINEs are merged in repeatmasker

import sys
import re

alignfile = sys.argv[1]  # Repeatmasker.org align file
TE = sys.argv[2]

outfile = open(alignfile+'.'+TE+'.out', 'w')
# write to bed file:
with open(alignfile, "r") as RMalign:
    for line in RMalign:
        if re.match(r'\d+\s', line):
            linelist = line.split()
            if len(line) < 14:
                print("RMalin file split wrong?")
            score = int(round(float(linelist[1])))
            chrom = linelist[4]
            start = int(linelist[5])
            end = int(linelist[6])
            if linelist[8] == 'C':
                strand = '-'
                names = linelist[9]
                name, TEclass = names.split('#')
            else:
                strand = '+'
                names = linelist[8]
                name, TEclass = names.split('#')
            if name == TE:
                outfile.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (chrom, start, end, name, score, strand))

            # or write everything:
            # outfile.write(line)
