import os
import argparse
import pysam

def methylation_calculation(bam_file, out_file, len_filter):
    """
    calculate the methylation average of the inserted regions in the bam file
    """

    bam = pysam.AlignmentFile(bam_file, threads = 8, check_sq=False)
    out = pysam.AlignmentFile(out_file, "wb", template=bam, threads = 8)
    for read in bam.fetch():
        query_name = read.query_name
        breakpoint()
        if read.infer_query_length == read.infer_read_length:
            read.set_tag('Mm', Mm_string, 'Z')
            read.set_tag('Ml', Ml_array)
        elif all_reads:
            cigar = read.cigartuples
            clip_process = False
            clip_length = 0
            if (not read.is_reverse) and cigar[0][0] == 5:
                clip_length = cigar[0][1]
                clip_process = True
            if read.is_reverse and cigar[-1][0] == 5:
                clip_length = cigar[-1][1]
                clip_process = True

            if clip_process:
                clip_seq = hash[query_name]['seq'][:clip_length]
                numC = clip_seq.count('C') + clip_seq.count('c')

                Mm_string_tail = ""
                if Mm_string[-1:] == ";":
                    Mm_string = Mm_string[:-1]
                    Mm_string_tail = ";"
                
                Mm_list = Mm_string.split(",")
                Mm_type = Mm_list.pop(0)
                Mm_list = list(map(int, Mm_list))
                while numC > 0 and len(Mm_list) > 0:
                    if Mm_list[0]>=numC:
                        Mm_list[0]-=numC
                        numC = 0
                    else:
                        first = Mm_list.pop(0)
                        numC -= first+1
                        Ml_array.pop(0)
                Mm_list = list(map(str, Mm_list))
                Mm_list.insert(0, Mm_type)
                Mm_string = ','.join(Mm_list) + Mm_string_tail
            read.set_tag('Mm', Mm_string, 'Z')
            read.set_tag('Ml', Ml_array)

            read.set_tag('Mm', Mm_string, 'Z')
            read.set_tag('Ml', Ml_array)
        out.write(read)

    out.close()
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regions in the bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-l', '--len', type=int, default=50,
                        help='lenght filter of the bam file. Only reads with length >= len will be considered')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bed like txt file storing the methylation average of the inserted regions')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    
    outfile = args.out
    methylation_calculation(bam_file, args.out, args.len)


if __name__ == '__main__':
    main()
