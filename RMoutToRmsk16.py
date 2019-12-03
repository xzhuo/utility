import sys
import gzip
file = sys.argv[1]
# with open(file, "rU") as infile:
with gzip.open(file, 'rt') as infile:
    with open(file[:-6] + 'bed', 'w') as outfile:  # remove .out.gz at the end of file name.
        # line = infile.readline()
        for line in infile:
            line = line.split()
            if not (line == [] or line[0] == 'SW' or line[0] == 'score'):
                chrom = line[4]
                start = int(line[5])
                end = int(line[6])
                name = line[9]
                score = 0
                strand = line[8]
                swScore = int(line[0])
                milliDiv = int(float(line[1]) * 10)
                milliDel = int(float(line[2]) * 10)
                milliIns = int(float(line[2]) * 10)
                genoLeft = -int(line[7].translate({ord('('): None, ord(')'): None}))
                repName = line[9]
                repS = line[10].split("/")
                if len(repS) == 2:
                    repClass, repFamily = repS
                else:
                    repClass = repS[0]
                    repFamily = repS[0]

                repEnd = int(line[12])
                if strand == '+':
                    repStart = int(line[11])
                    repLeft = -int(line[13].translate({ord('('): None, ord(')'): None}))
                else:
                    strand = "-"
                    repStart = -int(line[11].translate({ord('('): None, ord(')'): None}))
                    repLeft = int(line[13])
                out = [chrom, start, end, name, score, strand, swScore, milliDiv, milliDel, milliIns, genoLeft, repClass, repFamily, repStart, repEnd, repLeft]
                outfile.write('{}\n'.format('\t'.join(map(str, out))))
exit()
