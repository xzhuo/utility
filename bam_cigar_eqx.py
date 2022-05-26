import pysam
import argparse
from paf_cigar_eqx import cigar_remove_eqx


def modify_cigar(bam_file, out_file):
    bam = pysam.AlignmentFile(bam_file, threads = 8)
    out = pysam.AlignmentFile(out_file, "wb", template=bam, threads = 8)
    for read in bam.fetch(until_eof=True):
        try:
            cigar = read.cigarstring
            new_cigar = cigar_remove_eqx(cigar)
            read.cigarstring = new_cigar
        except:
            pass
        out.write(read)

    out.close()
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='Convert cigar string from eqx to regular M in the bam/sam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam/sam file with =|X instead of M in cigar string')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bam file with M in cigar string')

    args = parser.parse_args()
    bam_file = args.bam
    out_file = args.out
    modify_cigar(bam_file, out_file)


if __name__ == '__main__':
    main()