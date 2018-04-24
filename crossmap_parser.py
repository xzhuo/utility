'''
A script to parse crossmap.py results for SV project
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


def can_merge(frag1, frag2, distance):
    return (frag1.from_chr == frag2.from_chr and
            frag1.from_strand == frag2.from_strand and
            abs(max(frag1.from_start, frag2.from_start) - min(frag1.from_end, frag2.from_end) < distance) and
            frag1.to_chr == frag2.to_chr and
            frag1.to_strand == frag2.to_strand and
            abs(max(frag1.to_start, frag2.to_start) - min(frag1.to_end, frag2.to_end)) < distance)


def main():
    args = _get_args()
    with open(args.input, 'r') as Fh:
        regions = []
        last_region = Region.new()
        for line in Fh:
            line = line.rstrip()
            split = line.split()[6]
            if split != 'Fail':
                if split == '->' or split.startswith('split.1:'):
                    if len(last_region.frags) > 0:
                        regions.append(copy.deepcopy(last_region))
                    last_region = Region.init(line)
                else:
                    frag = Frag.digest_line(line)
                    if can_merge(last_region.frags[-1], frag, args.distance):
                        last_region.frags[-1].merge_frags(frag)
                    else:
                        last_region.frags.append(frag)

    for region in regions:
        line = "\t".join(region.from_chr, region.from_start, region.from_end, region.from_strand, region.from_summit, region.from_signal)
        print(line)


if __name__ == "__main__":
    main()


class Region:
    def new(self):
        self.frags = []

    def init(self, line):
        list = line.split()
        self.from_chr = list[0]
        self.from_start = list[1]
        self.from_end = list[2]
        self.from_summit = list[3]
        self.from_signal = list[4]
        self.from_strand = list[5]
        self.frags = []
        frag = Frag.digest_line(line)
        self.frags.append(frag)


class Frag:
    def digest_line(self, line):
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

    def merge_frags(self, frag):
        self.from_start = min(self.from_start, frag.from_start)
        self.from_end = max(self.from_end, frag.from_end)
        self.to_start = min(self.to_start, frag.to_start)
        self.to_end = max(self.to_end, frag.to_end)
