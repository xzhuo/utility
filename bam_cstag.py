import pysam
import re
import argparse

parser = argparse.ArgumentParser(description='simple arguments')
parser.add_argument(
    '--fasta',
    '-f',
    action="store",
    dest="fasta",
    help='The input fasta reference file.',
)
parser.add_argument(
    '--bam',
    '-b',
    action="store",
    dest="bam",
    help='The input bam file',
)
parser.add_argument(
    '--out',
    '-o',
    action="store",
    dest="out",
    help='The output sam file',
)
args = parser.parse_args()

# bam = "HG00741.maternal.CHM13Y_EBV.bam"
# out = "HG00741.maternal.CHM13Y_EBV.cslong.sam"
# fasta = "/taproom/data/xiaoyu/genomes/chm13chrY_EBV.fasta"

bam = pysam.AlignmentFile(args.bam, "rb")
outfile = pysam.AlignmentFile(args.out, "w", template=bam)
fa = pysam.FastaFile(args.fasta)


for read in bam.fetch():
    short_cs = read.get_tag('cs')
    blocks = read.get_blocks()
    chr = read.reference_name
    start = read.reference_start
    end = read.reference_end
    seq = fa.fetch(reference=chr, start=start, end=end)
    cs_array = re.split(r'([\:\-\+\*])', short_cs)
    cs_array.pop(0)
    symbol_list = [cs_array[i] for i in range(0, len(cs_array), 2)]
    value_list = [cs_array[i] for i in range(1, len(cs_array), 2)]
    curr_start = start
    cs_list = [{'symbol': key, 'value': value} for key, value in zip(symbol_list, value_list)]
    # print(short_cs)
    for i in cs_list:
        if i['symbol'] == ':':
            ref_length = int(i['value'])
            item_seq = fa.fetch(reference=chr, start=curr_start, end=curr_start + ref_length)
            i['symbol'] = '='
            i['value'] = item_seq
        elif i['symbol'] == '-':
            ref_length = len(i['value'])
            # test_seq = fa.fetch(reference=chr, start=curr_start, end=curr_start + ref_length)
            # print(test_seq)
        elif i['symbol'] == '+':
            ref_length = 0
        else:
            ref_length = len(i['value']) / 2
            # test_seq = fa.fetch(reference=chr, start=curr_start, end=curr_start + ref_length)
            # print(test_seq)
        curr_start += ref_length
    long_cs = ''
    for i in cs_list:
        long_cs += i['symbol']
        long_cs += i['value']

    read.set_tag('cs', long_cs)
    outfile.write(read)

outfile.close()
bam.close()
fa.close()
