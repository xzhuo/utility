import sys
import re
import copy

if len(sys.argv) != 3:
    print("python axt.split.py <axt file> <output file>")
    sys.exit()


def sub_align(Align, relative_start, relative_end):
    ref_seq = Align['seqs'][0]
    query_seq = Align['seqs'][1]
    Sub_align = {}
    Sub_align['seqs'] = []
    Sub_align['ref_chr'] = Align['ref_chr']
    Sub_align['query_chr'] = Align['query_chr']
    Sub_align['strand'] = Align['strand']

    Sub_align['seqs'].append(ref_seq[relative_start:relative_end])
    prev_ref_seq = ref_seq[0:relative_start]
    prev_ref_gaps = len(re.findall("-", prev_ref_seq))
    Sub_align['ref_start'] = Align['ref_start'] + relative_start - prev_ref_gaps
    Sub_align['ref_end'] = Sub_align['ref_start'] + len(re.findall(r'[^-]', Sub_align['seqs'][0])) - 1

    Sub_align['seqs'].append(query_seq[relative_start:relative_end])
    prev_query_seq = query_seq[0:relative_start]
    prev_query_gaps = len(re.findall("-", prev_query_seq))
    if Align['strand'] == "+":
        Sub_align['query_start'] = Align['query_start'] + relative_start - prev_query_gaps
        Sub_align['query_end'] = Sub_align['query_start'] + len(re.findall(r'[^-]', Sub_align['seqs'][1])) - 1
    else:
        Sub_align['query_end'] = Align['query_end'] - relative_start + prev_query_gaps
        Sub_align['query_start'] = Sub_align['query_end'] - len(re.findall(r'[^-]', Sub_align['seqs'][1])) + 1
    return Sub_align


def split(Align):
    gaps = []
    for gap in re.finditer(r'-{10,}', Align['seqs'][0]):
        gaps.append(gap.span())
    for gap in re.finditer(r'-{10,}', Align['seqs'][1]):
        gaps.append(gap.span())
    gaps.sort(key=lambda x: x[0])
    relative_start = 0
    aligns = []
    Right_align = copy.deepcopy(Align)
    for gap in gaps:
        Left_align = sub_align(Align, relative_start, gap[0])
        Right_align = sub_align(Align, gap[1], None)
        relative_start = gap[1]
        aligns.append(Left_align)
    aligns.append(Right_align)
    return aligns


all_aligns = []
with open(sys.argv[1]) as fin:
    Align = {}
    for line in fin:
        if re.match("#", line):
            continue

        list = line.rstrip().split()
        if len(list) >= 8:
            if 'ref_chr' in Align:
                aligns = split(Align)
                all_aligns.extend(aligns)
            Align.update({'item': list[0], 'ref_chr': list[1], 'ref_start': int(list[2]), 'ref_end': int(list[3]), 'query_chr': list[4],
                         'query_start': int(list[5]), 'query_end': int(list[6]), 'strand': list[7]})
            Align['seqs'] = []
        elif len(list) == 1:
            Align['seqs'].append(list[0])
    aligns = split(Align)
    all_aligns.extend(aligns)

with open(sys.argv[2], 'w') as fout:
    for index, Align in enumerate(all_aligns):
        fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(index, Align['ref_chr'], Align['ref_start'], Align['ref_end'], 
                   Align['query_chr'], Align['query_start'], Align['query_end'], Align['strand']))
        fout.write('{}\n'.format(Align['seqs'][0]))
        fout.write('{}\n'.format(Align['seqs'][1]))
        fout.write("\n")
