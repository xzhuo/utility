# Find out how LINEs are merged in repeatmasker

import sys
import re

alignfile = sys.argv[1]  # Repeatmasker.org align file
outfile = open(alignfile+'.out', 'w')
# put RMalign entries on the intervaltree first:
with open(alignfile, "r") as RMalign:
    for line in RMalign:
        if re.match(r'\d+\s', line):
            linelist = line.split()
            if len(line) < 14:
                print("RMalin file split wrong?")
            outfile.write(line)
