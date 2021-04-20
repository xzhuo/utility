import sys
import gzip
import codecs
import tarfile

def read_line(infile, outfile):
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
            swScore = abs(int(line[0]))
            milliDiv = abs(int(float(line[1]) * 10))
            milliDel = abs(int(float(line[2]) * 10))
            milliIns = abs(int(float(line[3]) * 10))
            genoLeft = -int(line[7].translate({ord('('): None, ord(')'): None}))
            repName = line[9]
            repS = line[10].split("/")
            if len(repS) == 2:
                repClass, repFamily = repS
            else:
                repClass = repS[0]
                repFamily = repS[0]

            repEnd = abs(int(line[12]))
            if strand == '+':
                repStart = int(line[11])
                repLeft = -int(line[13].translate({ord('('): None, ord(')'): None}))
            else:
                strand = "-"
                repStart = -int(line[11].translate({ord('('): None, ord(')'): None}))
                repLeft = int(line[13])
            out = [chrom, start, end, name, score, strand, swScore, milliDiv, milliDel, milliIns, genoLeft, repClass, repFamily, repStart, repEnd, repLeft]
            outfile.write('{}\n'.format('\t'.join(map(str, out))))


file = sys.argv[1]
if file[-7:] == '.tar.gz':
    tar = tarfile.open(file, "r:gz")
    with open(file[:-10] + 'rm.bed', 'w') as outfile:
        for member in tar.getmembers():
            f = tar.extractfile(member)
            if f is not None:
                f = codecs.getreader("utf-8")(f) # convert input from binary to text
                # read_line(f.readlines(), outfile)
                # read_line(f.read(), outfile)
                read_line(f, outfile)

elif file[-3:] == '.gz':
    with gzip.open(file, 'rt') as infile:
        with open(file[:-6] + 'rm.bed', 'w') as outfile:  # remove .out.gz at the end of file name.
            read_line(infile, outfile)
else:
    with open(file, "r") as infile:
        with open(file[:-3] + 'rm.bed', 'w') as outfile:  # remove .out at the end of file name.
            read_line(infile, outfile)
exit()
