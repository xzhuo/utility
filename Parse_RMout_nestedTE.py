import sys
import os

file = sys.argv[1]
TE = sys.argv[2]
# file = "/Users/xiaoyu/project/mm10_TE_RM/mm10_RM.fa.out"

infile = open(file, "r")
line = infile.readline()
alllines = []  # an array or array saving all records in memory.
for line in infile:
    line = line.split()
    if not (line == [] or line[0] == 'SW' or line[0] == 'score'):
        alllines.append(line)

infile.close()
last_TE = {}
outfile = open(file + '.' + TE + '.bed', 'w')
for index, line in enumerate(alllines):
    chrom = line[4]
    start = int(line[5])
    end = int(line[6])
    TEname = line[9]
    TEclass = line[10]
    RepEnd = int(line[12])
    score = float(line[1])
    strand = line[8]
    RepStart = line[11]
    RepID = line[14]
    RepLeft = int(line[13].translate({ord('('): None, ord(')'): None}))
    if strand == 'C':
        strand = "-"
        RepStart = line[13]
        RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))
    length = end - start + 1
    RepLength = RepLeft + RepEnd
    if TEname == TE:
        curr_TE = {
            "chrStart": start,
            "chrEnd": end,
            "chr": chrom,
            "TEname": TEname,
            "strand": strand,
            "repStart": RepStart,
            "repEnd": RepEnd,
            "repID": RepID,
            "rownum": index,
            "SWscore": score,
            "nested": [],
        }
        if bool(last_TE) and last_TE["chr"] == curr_TE["chr"] and last_TE["repID"] == curr_TE["repID"]:
            # stuff here...
            last_TE["chrEnd"] = curr_TE["chrEnd"]
            if curr_TE["strand"] == "+":
                last_TE["repEnd"] = curr_TE["repEnd"]
            else:
                last_TE["repStart"] = curr_TE["repStart"]
            last_TE["SWscore"] += curr_TE["SWscore"]
            if curr_TE["rownum"] - last_TE["rownum"] > 1:
                for nested_rownum in range(last_TE["rownum"], curr_TE["rownum"]):
                    TEwithStrand = alllines[nested_rownum]['TEname'].alllines[nested_rownum]['strand']
                    last_TE['nested'].append(TEwithStrand)

        else:
            if bool(last_TE):
                # print last_TE in bed format
                allnested = ",".join(last_TE['nested'])
                outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % (last_TE["chr"], last_TE["chrStart"], last_TE["chrEnd"], last_TE["TEname"], last_TE["SWscore"], last_TE["strand"], allnested)

# after loop print last_TE in bed format if bool(last_TE)
if bool(last_TE):
    # print last_TE in bed format
    outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % (last_TE["chr"], last_TE["chrStart"], last_TE["chrEnd"], last_TE["TEname"], last_TE["SWscore"], last_TE["strand"], allnested)

outfile.close()
exit()
