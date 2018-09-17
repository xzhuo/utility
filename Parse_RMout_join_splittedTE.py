import os
import argparse
# import ipdb


def _get_args():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--file',
        '-f',
        action="store",
        dest="file",
        help='The input rmsk out file.',
    )
    parser.add_argument(
        '--TE',
        '-t',
        action="store",
        dest="te",
        help='The target TE subfamily name.',
    )
    parser.add_argument(
        '--left',
        '-l',
        action="store",
        dest="left",
        help="number of missing bps tolerated at the end of TE.",
    )
    parser.add_argument(
        '--right',
        '-r',
        action="store",
        dest="right",
        help='number of missing bps tolerated at the start of TE.',
    )
    return parser.parse_args()


def main():
    args = _get_args()
    file = args.file
    TE = args.te
    TEleft_cut = int(args.left)
    # file = "/Users/xiaoyu/project/mm10_TE_RM/mm10_RM.fa.out"

    infile = open(file, "r")
    line = infile.readline()
    all_lines = []  # an array or array saving all records in memory.
    for line in infile:
        line = line.split()
        if not (line == [] or line[0] == 'SW' or line[0] == 'score') and line[9] == TE:
            all_lines.append(line)

    infile.close()
    last_TE = {}
    outfile = open(os.path.basename(file)[:-6] + TE + '.bed', 'w')  # remove fa.out and add TE.bed to the outfile
    for index, line in enumerate(all_lines):
        chrom = line[4]
        start = int(line[5])
        end = int(line[6])
        TEname = line[9]
        TEclass = line[10]
        RepEnd = int(line[12])
        score = float(line[0])
        strand = line[8]
        RepID = line[-1]  # some rows does not have repID. It is not perfect but works for MER57E3 in hg38.fa.out (2 rows without repID, use the last column instead)
        if strand == '+':
            RepStart = int(line[11])
            RepLeft = int(line[13].translate({ord('('): None, ord(')'): None}))
        else:
            RepStart = int(line[13])
            RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))
        length = end - start + 1
        RepLength = RepLeft + RepEnd
        curr_TE = {
            "chrStart": start,
            "chrEnd": end,
            "chr": chrom,
            "TEname": TEname,
            "strand": strand,
            "repStart": RepStart,
            "repEnd": RepEnd,
            "repLeft": RepLeft,
            "repID": RepID,
            "rownum": index,
            "SWscore": score
        }
        if bool(last_TE) and last_TE["chr"] == curr_TE["chr"] and last_TE["strand"] == curr_TE["strand"] and last_TE["repID"] == curr_TE["repID"]:
            # stuff here...
            last_TE["chrEnd"] = curr_TE["chrEnd"]
            if curr_TE["strand"] == "+":
                last_TE["repEnd"] = curr_TE["repEnd"]
            else:
                last_TE["repStart"] = curr_TE["repStart"]
            last_TE["SWscore"] += curr_TE["SWscore"]

        else:
            if bool(last_TE) and last_TE["repLeft"] < TEleft_cut:
                # print last_TE in bed format
                if not (args.right and last_TE["repStart"] - 1 > int(args.right)):  # don't print only if args.right exists and repStart bigger than it.
                    outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n" % (last_TE["chr"], last_TE["chrStart"], last_TE["chrEnd"], last_TE["TEname"], last_TE["SWscore"], last_TE["strand"], last_TE["repStart"], last_TE["repEnd"], last_TE["repLeft"]))

            last_TE = curr_TE

    # after loop print last_TE in bed format if bool(last_TE)
    if bool(last_TE) and last_TE["repLeft"] < TEleft_cut:
        # print last_TE in bed format
        if not (args.right and last_TE["repStart"] - 1 > int(args.right)):  # don't print only if args.right exists and repStart bigger than it.
            outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n" % (last_TE["chr"], last_TE["chrStart"], last_TE["chrEnd"], last_TE["TEname"], last_TE["SWscore"], last_TE["strand"], last_TE["repStart"], last_TE["repEnd"], last_TE["repLeft"]))

    outfile.close()
    exit()
