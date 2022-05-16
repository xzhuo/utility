import sys
import pysam
import argparse
import json


def revcom(seq):
    tab = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tab)[::-1]

def convertAxt(bam_file, ref_file, out_file, format):
    fasta = pysam.FastaFile(ref_file)
    bam = pysam.AlignmentFile(bam_file)
    out = open(out_file, "w")
    for idx, read in enumerate(bam.fetch()):
        cigar = read.cigartuples
        left_clip_length = cigar[0][1] if cigar[0][0] == 5 else 0
        right_clip_length = cigar[-1][1] if cigar[-1][0] == 5 else 0
        start = right_clip_length if read.is_reverse else left_clip_length
        query_start = start
        query_end = start + read.query_length
        pairs = read.get_aligned_pairs()
        strand = "-" if read.is_reverse else "+"
        query_seq = ''.join(["-" if pos[0] == None else read.query_sequence[pos[0]] for pos in pairs])
        ref_seq = ''.join(["-" if pos[1] == None else fasta.fetch(reference=read.reference_name, start=pos[1], end=pos[1]+1) for pos in pairs])
        if format == "axt":
            cols = ' '.join(map(str, [idx, read.reference_name, read.reference_start + 1, read.reference_end, read.query_name, query_start + 1, query_end, strand, "60"]))
            output = cols+"\n"+query_seq+"\n"+ref_seq+"\n\n"
        elif format == "align":
            genomealign = {"chr": read.query_name, "start": query_start + 1, "stop": query_end, "targetseq": ref_seq, "queryseq": query_seq}
            cols = '\t'.join(map(str, [read.reference_name, read.reference_start + 1, read.reference_end, "id:" + str(idx)]))
            output = cols + ",genomealign:" + json.dump(genomealign) + "\n"
        else:
            sys.exit('outout format has to be axt or align.')
        out.write(output)
    out.close()
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--ref',
        '-r',
        action="store",
        dest="ref",
        help='The input reference fasta file.',
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
        help='The output axt/align file',
    )
    parser.add_argument(
        '--format',
        '-f',
        action="store",
        dest="format",
        default="axt",
        help='The output file format, axt or align',
    )
    args = parser.parse_args()

    convertAxt(args.bam, args.ref, args.out, args.format)

if __name__ == '__main__':
    main()