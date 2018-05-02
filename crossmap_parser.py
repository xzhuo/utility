'''
A script to parse crossmap.py results for SV project
'''
import sys
import argparse
import copy


class Region:
    def __init__(self):
        self.frags = []

    def start_line(self, line):
        list = line.split()
        self.from_chr = list[0]
        self.from_start = int(list[1])
        self.from_end = int(list[2])
        self.from_summit = int(list[3])
        self.from_signal = float(list[4])
        self.from_strand = list[5]
        self.frags = []
        frag = Frag()
        frag.digest_line(line)
        self.frags.append(frag)

    def summit_frag_list(self):  # Return a list of all frags with summit (theoratically should only be one)
        return [x for x in self.frags if x.from_chr == self.from_chr and x.from_start <= self.from_summit and x.from_end >= self.from_summit]

    def length(self, assembly):  # assembly has to be "from" or "to"
        length = 0
        for frag in self.frags:
            length += frag.length(assembly)
        return length

    def combine_frags(self, size_limit):
        ''' Combine frags based on chr, strand, and within size_limit.'''
        frag_dict = {}

        for frag in self.frags:
            key = frag.to_chr + frag.to_strand
            if key in frag_dict:
                frag_dict[key].append(frag)
            else:
                frag_dict[key] = [frag]
        combined_frag_dict = {}
        for key in frag_dict:
            split_pos = []
            for i in range(1, len(frag_dict[key])):
                if not Region.can_merge(frag_dict[key][i - 1], frag_dict[key][i], size_limit):
                    split_pos.append(i)
            if not len(split_pos):
                combined_frag_dict[key] = frag_dict[key]
            for i in range(len(split_pos)):
                key_index = key + str(i)
                if i == 0:
                    combined_frag_dict[key_index] = frag_dict[key][:split_pos[i]]
                else:
                    combined_frag_dict[key_index] = frag_dict[key][split_pos[i - 1]:split_pos[i]]
                if i + 1 not in range(len(split_pos)):
                    key_index = key + str(i + 1)
                    combined_frag_dict[key_index] = frag_dict[key][split_pos[i]:]
        self.frags = []
        for key_index in combined_frag_dict:
            merged_frag = Region.merge_all_frags(combined_frag_dict[key_index])
            self.frags.append(merged_frag)

    def keep_primary(self, perc):
        '''If primary fragments found (frag with summit), then delele frags that are less then given perctage of total length'''
        self.from_summit
        length = self.length("from")
        can_remove = []
        found_primary = False
        for i, frag in enumerate(self.frags):
            if frag.length("from") / length < perc:
                can_remove.append(i)
            if self.from_summit > frag.from_start and self.from_summit < frag.from_end:
                found_primary = True
        if found_primary:
            self.frags = [i for j, i in enumerate(self.frags) if j not in can_remove]

    @staticmethod
    def merge_all_frags(frag_list):
        new_frag = Frag()
        new_frag.from_chr = [x.from_chr for x in frag_list][0]
        new_frag.from_strand = [x.from_strand for x in frag_list][0]
        new_frag.from_start = min([x.from_start for x in frag_list])
        new_frag.from_end = max([x.from_end for x in frag_list])
        new_frag.to_chr = [x.to_chr for x in frag_list][0]
        new_frag.to_strand = [x.to_strand for x in frag_list][0]
        new_frag.to_start = min([x.to_start for x in frag_list])
        new_frag.to_end = max([x.to_end for x in frag_list])
        return new_frag

    def can_merge(frag1, frag2, distance):
        return (frag1.from_chr == frag2.from_chr and
                frag1.from_strand == frag2.from_strand and
                abs(max(frag1.from_start, frag2.from_start) - min(frag1.from_end, frag2.from_end) < distance) and
                frag1.to_chr == frag2.to_chr and
                frag1.to_strand == frag2.to_strand and
                abs(max(frag1.to_start, frag2.to_start) - min(frag1.to_end, frag2.to_end)) < distance and
                (frag2.to_start >= frag1.to_start if (frag1.to_strand == "+") else frag2.to_end <= frag1.to_end))


class Frag:
    def __init__(self):
        pass

    def digest_line(self, line):
        list = line.split()
        self.to_chr = list[7]
        self.to_start = int(list[8])
        self.to_end = int(list[9])
        self.to_strand = list[12]
        if list[6] == '->':
            self.from_chr = list[0]
            self.from_start = int(list[1])
            self.from_end = int(list[2])
            self.from_strand = list[5]
        else:
            from_frag_list = list[6].strip('()').split(":")
            self.from_chr = from_frag_list[1]
            self.from_start = int(from_frag_list[2])
            self.from_end = int(from_frag_list[3])
            self.from_strand = from_frag_list[4]

    def merge_frags(self, frag):
        self.from_start = min(self.from_start, frag.from_start)
        self.from_end = max(self.from_end, frag.from_end)
        self.to_start = min(self.to_start, frag.to_start)
        self.to_end = max(self.to_end, frag.to_end)

    def length(self, assembly):  # assembly has to be "from" or "to"
        if assembly == "from":
            length = self.from_end - self.from_start
        elif assembly == "to":
            length = self.to_end - self.to_start
        return length


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
    parser.add_argument(
        '--max',
        '-m',
        action="store",
        dest="max",
        default=50000,
        help='The distance used to define the max of peak size',
    )
    parser.add_argument(
        '--perc',
        '-p',
        action="store",
        dest="perc",
        default=0.2,
        help='Threshold to remove non-primary fragments (fraction of total region length)',
    )
    return parser.parse_args()


def main():
    args = _get_args()
    with open(args.input, 'r') as Fh:
        regions = []
        last_region = Region()
        for line in Fh:
            line = line.rstrip()
            split = line.split()[6]
            if split != 'Fail':
                if split == '->' or split.startswith('(split.1:'):
                    if len(last_region.frags) > 0:
                        regions.append(copy.deepcopy(last_region))
                    last_region.start_line(line)
                else:
                    frag = Frag()
                    frag.digest_line(line)
                    if Region.can_merge(last_region.frags[-1], frag, args.distance):
                        last_region.frags[-1].merge_frags(frag)
                    else:
                        last_region.frags.append(frag)

    for region in regions:
        region.combine_frags(args.max)
        region.keep_primary(args.perc)

        for frag in region.frags:
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                  (frag.from_chr, frag.from_start, frag.from_end, frag.from_strand, region.from_summit, region.from_signal,
                   frag.to_chr, frag.to_start, frag.to_end, frag.to_strand))


if __name__ == "__main__":
    main()
