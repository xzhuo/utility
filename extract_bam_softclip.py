import os
import argparse
import pysam
import re
# from Bio.Seq import Seq


def extract_softclip(bam_file, out, failed, extend):
    bam = pysam.AlignmentFile(bam_file, "rb")
    with open(out, 'w') as outfile:
        with open(failed, 'w') as failed_outfile:
            # Iterate through each read in the BAM file 
            for read in bam:
                if read.is_supplementary or read.is_secondary or read.is_unmapped:
                    continue
                else:
                    if (read.is_reverse and read.query_alignment_start > 0) or (read.is_forward and read.query_alignment_end < read.query_length):
                        offset = read.query_length - read.query_alignment_start if read.is_reverse else read.query_alignment_end
                        forward_sequence = read.get_forward_sequence()
                        adj_offset = offset
                        for i in range(offset, offset - extend, -1):
                            if forward_sequence[i] == 'T':
                                adj_offset = i
                                break
                        soft_clip = forward_sequence[adj_offset:]
                        # breakpoint()
                        pattern = re.compile('TGAAA'[:len(soft_clip)]) if len(soft_clip) < 5 else re.compile('TGAAA.*')
                        match = pattern.search(str(soft_clip))
                        if match:
                            outfile.write(f">{read.query_name}_clip\n{match.group(0)}\n")
                        else:
                            # If no match, write the clipped sequence
                            soft_clip = forward_sequence[offset:]
                            failed_outfile.write(f">{read.query_name}_clip\n{soft_clip}\n")

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
