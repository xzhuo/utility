import os
import argparse
from xmlrpc.client import boolean
import pysam
from copy import copy

def revcom(seq):
    tab = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tab)[::-1]

def attach_tags(bam_file, tag_file, out_file, all_reads=False):
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
                except KeyError:
                    try:
                        MM = read.get_tag("MM")
                        ML = read.get_tag("ML")
                        hash[read.query_name] = {'MM': MM, 'ML': ML, 'seq': seq}
                    except KeyError:
                        pass
                except:
                    pass

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
                Ml_array = copy(hash[query_name]['Ml'])
            except:
                Ml_array = copy(hash[query_name]['ML'])
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
    parser = argparse.ArgumentParser(description='attach Mm and Ml tags from an unmapped bam file to another bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file without the desired tags')
    parser.add_argument('-t', '--tag', type=str, required=True,
                        help='bam file with tags, or dir with these bam files')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bam file with Ml and Mm tags attached')
    parser.add_argument('-a', '--all', required=False, default=False, action='store_true',
                        help='Attach Ml and Mm tags to all reads, including hard-clipped reads')
    args = parser.parse_args()
    tag_file = os.path.abspath(args.tag)
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(tag_file):
        raise ValueError("--tag file does not exist!")
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    
    outfile = args.out
    attach_tags(bam_file, tag_file, outfile, args.all)


if __name__ == '__main__':
    main()
