import pysam
import argparse


def revcom(seq):
    tab = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tab)[::-1]

def attach_query_seq(bam_file, fasta_file, out_file):
    fasta = pysam.FastaFile(fasta_file)
    bam = pysam.AlignmentFile(bam_file)
    out = open(out_file, "w")
    for idx, read in enumerate(bam.fetch()):
        pairs = read.get_aligned_pairs()
        strand = "-" if read.is_reverse else "+"
        # query_seq = ''
        # ref_seq = ''
        # for pos in pairs:
        #     query_seq += "-" if pos[0] == None else read.query_sequence[pos[0]]
        #     ref_seq += "-" if pos[1] == None else fasta.fetch(reference=read.reference_name, start=pos[1], end=pos[1]+1)
        query_seq = ''.join(["-" if pos[0] == None else read.query_sequence[pos[0]] for pos in pairs])
        ref_seq = ''.join(["-" if pos[1] == None else fasta.fetch(reference=read.reference_name, start=pos[1], end=pos[1]+1) for pos in pairs])
        cols = '\t'.join([str(idx), read.reference_name, read.reference_start, read.reference_end, read.query_name, read.query_alignment_start, read.query_alignment_end, strand, "60"])
        out = cols+"\n"+query_seq+"\n"+ref_seq+"\n"
        out.write(out)
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
        help='The output axt file',
    )
    args = parser.parse_args()

    attach_query_seq(args.bam, args.fasta, args.out)

if __name__ == '__main__':
    main()