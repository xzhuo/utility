import os
import argparse
import pysam
import re
from Bio.Seq import Seq

def match_end(position, read_name, seq, outfile, failed_outfile, extend):
    """
    From the extended index looking upstream, print the first T found.
    """
    clean_seq = seq
    for i in range(extend, -1, -1):
        if seq[i] == 'T':
            clean_seq = seq[i:]
            break
    match = None
    if len(clean_seq) < 5:
        pattern = re.compile('TGAAA'[:len(clean_seq)])
        match = pattern.match(str(clean_seq))
    else:
        match = re.search('TGAAA.*', str(clean_seq))
    if match:
        outfile.write(f">{read_name}_{position}\n{match.group(0)}\n")
    else:
        # If no match, write the extended sequence
        soft_clip = seq[extend:]
        failed_outfile.write(f">{read_name}_{position}\n{soft_clip}\n")


def extract_softclip(bam_file, out, failed, extend):
    bam = pysam.AlignmentFile(bam_file, "rb")
    with open(out, 'w') as outfile:
        with open(failed, 'w') as failed_outfile:
            # Iterate through each read in the BAM file 
            for read in bam:
                if read.is_supplementary or read.is_secondary or read.is_unmapped:
                    continue
                else:
                    if read.query_alignment_start > 0:
                        soft_left = read.query_sequence[:read.query_alignment_start + extend]
                        revcom_left = Seq(soft_left).reverse_complement()
                        match_end("left", read.query_name, revcom_left, outfile, failed_outfile, extend)
                    if read.query_alignment_end < read.query_length:
                        soft_right = read.query_sequence[read.query_alignment_end - extend:]
                        match_end("right", read.query_name, soft_right, outfile, failed_outfile, extend)
                # breakpoint()


def main():
    parser = argparse.ArgumentParser(description='Extract soft clipped sequences from all reads of a bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output extracted seq in fasta format')
    parser.add_argument('-f', '--failed', type=str, required=True,
                        help='output seq does not match the pattern in fasta format')
    parser.add_argument('-e', '--extend', type=int, default=5,
                        help='extend the region by this many bp before extracting the sequence')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    extract_softclip(bam_file, args.out, args.failed, args.extend)

if __name__ == '__main__':
    main()
