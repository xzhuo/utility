import pysam
import argparse

def modify_tags(bam_file, out_file):
    bam = pysam.AlignmentFile(bam_file, threads = 8)
    out = pysam.AlignmentFile(out_file, "wb", template=bam, threads = 8)
    for read in bam.fetch(until_eof=True):
        try:
            Mm = read.get_tag("Mm")
            if Mm[3]=="?":
                Mm = Mm[:3]+Mm[4:]
            read.set_tag('Mm', Mm)
        except:
            pass
        out.write(read)

    out.close()
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='Change C+m? to C+m in bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bam file with Ml and Mm tags corrected')

    args = parser.parse_args()
    bam_file = args.bam
    out_file = args.out
    modify_tags(bam_file, out_file)


if __name__ == '__main__':
    main()