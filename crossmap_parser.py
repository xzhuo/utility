'''
A script to parse crossmap.py results for SV pr0ject
'''
import sys
import argparse
import ipdb
import copy

def _get_args():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--input',
        '-i',
        action="store",
        dest="input",
        help='The input crossmap output file.',
    )
    parser.add_argument(
        '--distance',
        '-d',
        action="store",
        dest="distance",
        default=50,
        help='The distance used to combine small indels into a large region. Default is 50bp',
    )
    return parser.parse_args()


def main():
    args = _get_args()
    with open(args.input, 'r') as Fh:
        regions = []
        last_region
        for line in Fh:
            line = line.rstrip()
            split = line.split()[6]
            if split == '->' or split.startswith('split.1'):
                regions.append(copy.deepcopy(last_region))
                region.digest(line)


class region:
    def __init__:
        self.frags = []

    def digest(self, line):
        list = line.split()
        self.from_chr = list[0]
        self.from_start = list[1]
        self.from_end = list[2]
        self.from_summit = list[3]
        self.from_signal = list[4]
        self.from_strand = list[5]
        frag = frag.digest_frag(list)
        self.frags.append(frag)


class frag:
    def digest_frag(self,line):
        list = line.split()
        self.frag_to_chr = list[7]
        self.frag_to_start = list[8]
        self.frag_to_end = list[9]
        self.frag_to_strand = list[12]
        if list[6] == '->':
            self.frag_from_chr = list[0]
            self.frag_from_start = list[1]
            self.frag_from_end = list[2]
            self.frag_from_strand = list[5]
        else:
            from_frag_list = list[6].split(":")
            self.frag_from_chr = from_frag_list[1]
            self.frag_from_start = from_frag_list[2]
            self.frag_from_end = from_frag_list[3]
            self.frag_from_strand = from_frag_list[4]