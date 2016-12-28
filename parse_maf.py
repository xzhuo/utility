import sys
import argparse
import warnings
import re
from collections import OrderedDict


def main():
    args = _get_args()
    genomes = args.filter.split(",")  # The list to store all included genomes
    curr_block = OrderedDict()
    last_block = OrderedDict()
    with open(args.maf, 'r') as Fh:
        with open(args.out, 'w') as Out:
            for num, line in enumerate(Fh, 1):
                if line.startswith('#'):
                    Out.write(line)
                    continue
                linelist = line.split()
                if linelist.pop(0) == 'a':  # new alignmentblock
                    last_block = merge_blocks(last_block, curr_block, genomes, Out)
                    curr_block = OrderedDict()
                    for item in linelist:
                        try:
                            (key, value) = item.split("=")
                            curr_block[key] = value
                        except RuntimeWarning:
                            print("Found abnormal a line in %d!" % (num))
                if linelist.pop(0) == 's':
                    (assembly, chrom) = linelist.pop(0).split(".")
                    if assembly not in genomes:
                        continue
                    (start, length, strand, chrlenth, seq) = linelist
                    curr_block[assembly] = {}
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['start'] = start
                    curr_block[assembly]['length'] = length
                    curr_block[assembly]['strand'] = strand
                    curr_block[assembly]['chrlenth'] = chrlenth
                    curr_block[assembly]['seq'] = seq
                if linelist.pop(0) == 'i':
                    (assembly, chrom) = linelist.pop(0).split(".")
                    if assembly not in genomes:
                        continue
                    (leftStatus, leftCount, rightStatus, rightCount) = linelist
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['leftStatus'] = leftStatus
                    curr_block[assembly]['leftCount'] = leftCount
                    curr_block[assembly]['rightStatus'] = rightStatus
                    curr_block[assembly]['rightCount'] = rightCount
                if linelist.pop(0) == 'e':
                    (assembly, chrom) = linelist.pop(0).split(".")
                    if assembly not in genomes:
                        continue
                    (start, length, strand, chrlenth, gapStatus) = linelist
                    curr_block[assembly] = {}
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['start'] = start
                    curr_block[assembly]['length'] = length
                    curr_block[assembly]['strand'] = strand
                    curr_block[assembly]['chrlenth'] = chrlenth
                    curr_block[assembly]['gapStatus'] = gapStatus
                if linelist.pop(0) == 'q':
                    (assembly, chrom) = linelist.pop(0).split(".")
                    quality = linelist[0]
                    curr_block[assembly]['quality'] = quality

            else:
                last_block = merge_blocks(last_block, curr_block, genomes, Out)
                print_block(last_block, Out)


def _get_args():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--inmaf',
        '-i',
        action="store",
        dest="maf",
        help='The input maf file.',
    )
    parser.add_argument(
        '--outmaf',
        '-o',
        action="store",
        dest="out",
        help='The onput maf file.',
    )
    parser.add_argument(
        '--filter',
        '-f',
        action="store",
        dest="filter",
        help="the genomes included in the output, comma seperated. The reference genome must be the first in the list.",
    )
    parser.add_argument(
        '--reference',
        '-r',
        action="store",
        dest="ref",
        help='the reference genome used to anchor the alignment,',
    )


def merge_blocks(last_block, curr_block, genomes, Out):
    complete = 1  # if curr_block contains all the species in genomes
    for assembly in genomes:
        if assembly not in curr_block:
            complete = 0
            break

    if complete:
        merge = 1  # if all species in curr_block are continue.
        for assembly in genomes:
            if curr_block[assembly]['leftStatus'] != 'C' or curr_block[assembly]['leftCount'] != 0:
                merge = 0
                break
        if merge:  # all species in curr_block are continue, merge curr_block into last_block and return it.
            last_block['score'] += curr_block['score']
            for key in ('length', 'seq', 'quality'):
                if key in last_block[assembly]:
                    last_block[assembly][key] += curr_block[assembly][key]
            for key in ('rightStatus', "rightCount"):
                if key in last_block[assembly]:
                    last_block[assembly][key] = curr_block[assembly][key]
            return(last_block)
        else:  # not all speceis are continue in the curr_block, print last_block, and return curr_block as next last_block.
            print_block(last_block, Out)
            return(curr_block)

    else:  # skip curr_block and print last_block because some species in missing in the curr_block.
        print_block(last_block, Out)
        curr_block.clear()
        return (curr_block)


def print_block(block, Out):
    refList = []
    alnList = []
    gapList = []
    if bool(block):
        Out.write("a\tscore= %.6f" % (block['score']))
        for key in block:
            if 'seq' in block[key] and 'leftStatus' not in block[key]:
                refList.append(key)
            if 'seq' in block[key] and 'leftStatus' in block[key]:
                alnList.append(key)
            if 'gapStatus' in block[key]:
                gapList.append(key)
        if len(refList) > 1: print("missing i line?")
        for key in refList:
            Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['seq']))
        for key in alnList:
            Out.write("i\t%s.%s\t%s\t%d\t%s\t%d" % (key, block[key]['chrom'], block[key]['leftStatus'], block[key]['leftCount'], block[key]['rightStatus'], block[key]['rightCount']))
        for key in gapList:
            Out.write("e\t%s.%s\t%d\t%d\t%s\t%d\t%s" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['gapStatus']))
