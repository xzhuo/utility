import sys
import gzip


def read_line(infile, outfile):
    # line = infile.readline()
    for line in infile:
        line = line.split()
        if not (line == [] or line[0] == 'SW' or line[0] == 'score'):
            chrom = line[5]
            start = int(line[6])
            end = int(line[7])
            name = line[10]
            score = 0
            strand = line[9]
            swScore = int(line[1])
            milliDiv = int(float(line[2]) * 10)
            milliDel = int(float(line[3]) * 10)
            milliIns = int(float(line[4]) * 10)
            genoLeft = int(line[8].translate({ord('('): None, ord(')'): None}))
            repName = line[10]
            repClass = line[11]
            repFamily = line[12]

            repStart = int(line[13])
            repEnd = abs(int(line[14]))
            repLeft = int(line[15])
            out = [chrom, start, end, name, score, strand, swScore, milliDiv, milliDel, milliIns, genoLeft, repClass, repFamily, repStart, repEnd, repLeft]
            outfile.write('{}\n'.format('\t'.join(map(str, out))))


file = sys.argv[1]
if file[-3:] == '.gz':
    with gzip.open(file, 'rt') as infile:
        with open(file[:-6] + 'rm.bed', 'w') as outfile:  # remove .out.gz at the end of file name.
            read_line(infile, outfile)
else:
    with open(file, "r") as infile:
        with open(file[:-3] + 'rm.bed', 'w') as outfile:  # remove .out at the end of file name.
            read_line(infile, outfile)
exit()
