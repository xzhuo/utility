import os
import argparse
import pysam

def revcom(seq):
    tab = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tab)[::-1]

def attach_tags(bam_file, tag_file, out_file):
    hash = {}
    if os.path.isfile(tag_file):
        tag_files = [tag_file]
    else:
        tag_files = [os.path.join(tag_file, x) for x in os.listdir(tag_file) if len(x) >= 4 and x[-4:] == ".bam"]

    for f in tag_files:
        tag_pysam = pysam.AlignmentFile(f, check_sq=False, threads = 8)
        for read in tag_pysam.fetch(until_eof=True):
            if read.infer_query_length() == read.infer_read_length():
                seq = read.get_forward_sequence()
                try:
                    Mm = read.get_tag("Mm")
                    Ml = read.get_tag("Ml")
                    if Mm[3]=="?":
                        Mm = Mm[:3]+Mm[4:]
                    hash[read.query_name] = {'Mm': Mm, 'Ml': Ml, 'seq': seq}
                except:
                    MM = read.get_tag("MM")
                    ML = read.get_tag("ML")
                    hash[read.query_name] = {'MM': MM, 'ML': ML, 'seq': seq}
        tag_pysam.close()

    bam = pysam.AlignmentFile(bam_file, threads = 8)
    out = pysam.AlignmentFile(out_file, "wb", template=bam, threads = 8)
    for read in bam.fetch():
        query_name = read.query_name
        if query_name in hash:
            try:
                Mm_string = hash[query_name]['Mm']
            except:
                Mm_string = hash[query_name]['MM']
            try:
                Ml_array = hash[query_name]['Ml']
            except:
                Ml_array = hash[query_name]['ML']
            if not read.infer_query_length() == read.infer_read_length():
                cigar = read.cigartuples()
                if cigar[0][0] == 5:
                    left_clip_length = cigar[0][1]
                    # revcom if necessary:
                    seq = revcom(hash[query_name]['seq']) if read.is_reverse() else hash[query_name]['seq']
                    left_clip_seq = seq[:left_clip_length]
                    numC = left_clip_seq.count('C') + left_clip_seq.count('c')
                    Mm_list = Mm_string.splt(",")
                    Mm_type = Mm_list.pop(0)
                    while numC > 0:
                        if Mm_list[0]>=numC:
                            Mm_list[0]-=numC
                            numC = 0
                        else:
                            first = Mm_list.pop(0)
                            numC -= first+1
                            Ml_array.pop(0)
                    Mm_list.insert(0, Mm_type)
                    Mm_string = ','.join(Mm_list)

            read.set_tag('Mm', Mm_string, 'Z')
            read.set_tag('Ml', Ml_array)
        out.write(read)

    out.close()
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='attach Mm and Ml tags from an unmapped bam file to another bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file without the desired tags')
    parser.add_argument('-t', '--tag', type=str, required=True,
                        help='bam file with tags, or dir with these bam files')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bam file with Ml and Mm tags attached')

    args = parser.parse_args()
    tag_file = os.path.abspath(args.tag)
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(tag_file):
        raise ValueError("--tag file does not exist!")
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    
    outfile = args.out
    attach_tags(bam_file, tag_file, outfile)


if __name__ == '__main__':
    main()
