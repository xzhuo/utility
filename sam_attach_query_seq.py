from turtle import right
import pysam
import argparse


def revcom(seq):
    tab = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tab)[::-1]

def attach_query_seq(bam_file, fasta_file, out_file):
    fasta = pysam.FastaFile(fasta_file)
    bam = pysam.AlignmentFile(bam_file)
    out = pysam.AlignmentFile(out_file, "w", template=bam)
    for read in bam.fetch():
        cigar = read.cigartuples
        left_clip_length = cigar[0][1] if cigar[0][0] == 5 else 0
        right_clip_length = cigar[-1][1] if cigar[-1][0] == 5 else 0
        if not read.is_reverse:
            read.query_sequence = fasta.fetch(reference=read.query_name, start=left_clip_length, end=left_clip_length + read.infer_query_length())
        else:
            read.query_sequence = revcom(fasta.fetch(reference=read.query_name, start=right_clip_length, end=right_clip_length + read.infer_query_length()))
        out.write(read)
    out.close()
    bam.close()

def main():
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
        help='The input bam/sam file without query seq',
    )
    parser.add_argument(
        '--out',
        '-o',
        action="store",
        dest="out",
        help='The output sam file with hard clipping seq',
    )
    args = parser.parse_args()

    attach_query_seq(args.bam, args.fasta, args.out)

if __name__ == '__main__':
    main()