import argparse
import re
from collections import OrderedDict
# Convert cigar string from eqx to regular M in the paf file 

def cigar_remove_eqx(cigar):
    # convert "=" and "X" to "M" in a cigar string.
    cigar_list = re.split('([M|I|D|N|S|H|P|X|=])',cigar)
    cigar_tuples = [(cigar_list[i+1], cigar_list[i]) for i in range(0, len(cigar_list)-1, 2)]
    new_cigar_tuples = []
    for i in range(len(cigar_tuples)):
        if not last_list:
            last_list = list(cigar_tuples[i])
        else:
            if last_list[0] == 'M' and (cigar_tuples[i][0] == '=' or cigar_tuples[i][0] == 'X'):
                last_list[1] += cigar_tuples[i][1]
            else:
                new_cigar_tuples.append(tuple(last_list))
                if (cigar_tuples[i][0] == '=' or cigar_tuples[i][0] == 'X'):
                    last_list = ["M", cigar_tuples[i][1]]
                else:
                    last_list = list(cigar_tuples[i])
    new_cigar_tuples.append(tuple(last_list))
    new_cigar = ''.join([str(i[1])+str(i[0]) for i in new_cigar_tuples])
    return new_cigar

def convert_paf(in_file, out_file):
    # write to a new paf file:
    with open(in_file, "r") as input_file:
        with open(out_file, "w") as output_file:
            for line in input_file:
                linelist = line.split(maxsplit=13)
                tags = linelist[12]
                tag_list = re.split('[:|\s+]', tags)
                tag_dict = OrderedDict()
                tag_dict = {tag_list[i]+':'+tag_list[i+1]: tag_list[i+2] for i in range(0, len(tag_list), 3)}
                cigar = tag_dict['cg:Z']
                new_cigar = cigar_remove_eqx(cigar)
                tag_dict['cg:Z'] = new_cigar
                new_tags = '\t'.join([str(k)+':'+str(v) for (k,v) in tag_dict.items()])
                linelist[12] = new_tags
                new_line = '\t'.join(linelist)
                output_file.write(new_line)

def main():
    parser = argparse.ArgumentParser(description='Convert cigar string from eqx to regular M in the paf file')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input paf file with =|X instead of M in cigar string')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output paf file with M in cigar string')

    args = parser.parse_args()
    in_file = args.input
    out_file = args.out
    convert_paf(in_file, out_file)