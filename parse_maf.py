import sys
import argparse
import warnings
from collections import OrderedDict

def main():
    args = _get_args()
    genomes = args.assemblies.split(",")  # The list to store all included genomes
    curr_block = OrderedDict()
    last_block = OrderedDict()
    with open(args.maf, 'r') as Fh:
        with open(args.out, 'w') as Out:
            for num, line in enumerate(Fh, 1):
                if line.startswith('#'):
                    Out.write(line + "\n")
                    continue
                linelist = line.split()
                if len(linelist) < 2:
                    continue
                lead = linelist.pop(0)
                if lead == 'a':  # new alignmentblock
                    if _can_merge(last_block, curr_block):
                        last_block = merge_blocks(last_block, curr_block, genomes, Out)
                    else:
                        print_block(last_block, Out)
                        last_block = curr_block
                    curr_block = OrderedDict()
                    for item in linelist:
                        try:
                            (key, value) = item.split("=")
                            curr_block[key] = float(value)
                        except ValueError:
                            print("Found abnormal a line in %d!" % (num))
                if lead == 's':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    if assembly not in genomes:
                        continue
                    (start, length, strand, chrlenth, seq) = linelist
                    curr_block[assembly] = {}
                    curr_block[assembly]['aln'] = 1  # it is an alignment seq in the block
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['start'] = int(start)
                    curr_block[assembly]['length'] = int(length)
                    curr_block[assembly]['strand'] = strand
                    curr_block[assembly]['chrlenth'] = int(chrlenth)
                    curr_block[assembly]['seq'] = seq
                if lead == 'i':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    if assembly not in genomes:
                        continue
                    (leftStatus, leftCount, rightStatus, rightCount) = linelist
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['leftStatus'] = leftStatus
                    curr_block[assembly]['leftCount'] = int(leftCount)
                    curr_block[assembly]['rightStatus'] = rightStatus
                    curr_block[assembly]['rightCount'] = int(rightCount)
                if lead == 'e':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    if assembly not in genomes:
                        continue
                    (start, length, strand, chrlenth, gapStatus) = linelist
                    curr_block[assembly] = {}
                    curr_block[assembly]['aln'] = 0  # it is not an alignment seq in the block
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['start'] = int(start)
                    curr_block[assembly]['length'] = int(length)
                    curr_block[assembly]['strand'] = strand
                    curr_block[assembly]['chrlenth'] = int(chrlenth)
                    curr_block[assembly]['gapStatus'] = gapStatus
                if lead == 'q':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
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
        '--assemblies',
        '-a',
        action="store",
        dest="assemblies",
        help="the genomes included in the output, comma seperated. The reference genome must be the first in the list.",
    )
    parser.add_argument(
        '--reference',
        '-r',
        action="store",
        dest="ref",
        help='the reference genome used to anchor the alignment,',
    )
    parser.add_argument(
        '--directory',
        '-d',
        action="store",
        dest="dir",
        help='The directory with all the maf files,',
    )
    parser.add_argument(
        '--format',
        '-f',
        action="store",
        dest="format",
        help='The output format, could be maf or bed.',
    )
    return parser.parse_args()


def _is_complete(curr_block, genomes):
    complete = 1  # if curr_block contains all the species in genomes
    for assembly in genomes:
        if assembly not in curr_block:
            complete = 0
            break
    return complete


def _is_continue(last_block, curr_block, genomes):
    # leftStatus C or I <10
    # or gap in the same species with C
    # or gap in the same species with I?
    continued = 1
    if len(last_block) == 0:
        continued = 0
    else:
        for assembly in genomes[1:]:
            if curr_block[assembly]['aln'] == 1:
                if curr_block[assembly]['leftStatus'] != 'C' or curr_block[assembly]['leftCount'] != 0:  # worry about small indels later
                    continued = 0
                    break
            elif last_block[assembly]['aln'] == 0:
                if curr_block[assembly]['gapStatus'] != 'C' or last_block[assembly]['gapStatus'] != 'C':  # now only merge 'C' gap status
                    continued = 0
                    break
            else:
                continued = 0
                break
    return continued


def compare_blocks(last_block, curr_block, genomes, Out):
    status = ()
    anchor = genomes[0]
    if (curr_block[anchor]['chrom'] == last_block[anchor]['chrom'] and
        curr_block[anchor]['start'] == last_block[anchor]['start'] + last_block[anchor]['length'] and
            curr_block[anchor]['strand'] == last_block[anchor]['strand']):
                status[anchor] = 'c'  # continue
                for assembly in genomes[1:]:
                    if curr_block[assembly]['aln'] == 1 and last_block[assembly]['aln'] == 1:
                        if curr_block[assembly]['leftStatus'] == 'C' and last_block[assembly]['rightStatus'] == 'C':
                            if curr_block[assembly]['leftCount'] > 0 or last_block[assembly]['rightCount'] > 0:
                                print("error C not associated with 0 at chrom start!")
                            else:
                                status[assembly] = 'c'
                        elif curr_block[assembly]['leftStatus'] == 'i' and last_block[assembly]['rightStatus'] == 'i':

                        
    else:
        status[anchor] = 'b'  # break








def merge_blocks(last_block, curr_block, genomes, Out):
    if _is_complete(curr_block, genomes):  # if the curr_block contains all the species in alignment
        if _is_continue(last_block, curr_block, genomes):  # all species in curr_block are continue, merge curr_block into last_block and return it.
            last_block['score'] += curr_block['score']
            for assembly in genomes:
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
        Out.write("a\tscore= %.6f\n" % (block['score']))
        for key in block:
            if key == 'score':
                continue
            if 'seq' in block[key] and 'leftStatus' not in block[key]:
                refList.append(key)
            if 'seq' in block[key] and 'leftStatus' in block[key]:
                alnList.append(key)
            if 'gapStatus' in block[key]:
                gapList.append(key)
        if len(refList) > 1: print("missing i line?")
        for key in refList:
            Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['seq']))
        for key in alnList:
            Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['seq']))
            Out.write("i\t%s.%s\t%s\t%d\t%s\t%d\n" % (key, block[key]['chrom'], block[key]['leftStatus'], block[key]['leftCount'], block[key]['rightStatus'], block[key]['rightCount']))
        for key in gapList:
            Out.write("e\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['gapStatus']))
        Out.write("\n")

if __name__ == '__main__':
    main()
