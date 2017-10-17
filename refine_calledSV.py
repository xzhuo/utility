'''
A script to process sv_calls result bed file from zev's DASVC pipeline. https://github.com/zeeev/DASVC
It will combine "insertion" and "deletion" in the same locus into a "replace" SV,
and get exact INDEL positions for query genome (the query starts and ends are starts and ends of alignment block in the sv_calls file).
'''
import sys
import argparse
import ipdb


def _get_args():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--sv',
        '-s',
        action="store",
        dest="sv",
        help='The input sv calling file.',
    )
    parser.add_argument(
        '--bam',
        '-b',
        action="store",
        dest="bam",
        help='The input bam file.',
    )
    parser.add_argument(
        '--distance',
        '-d',
        action="store",
        dest="distance",
        help='The distance extended to both flanking region used for get_flanking_refinedSV.py.',
    )
    parser.add_argument(
        '--target_size',
        '-c',
        action="store",
        dest="target_size",
        help='The target chromosome size file used for get_flanking_refinedSV.py.',
    )
    parser.add_argument(
        '--query_size',
        '-c',
        action="store",
        dest="query_size",
        help='The query chromosome size file used for get_flanking_refinedSV.py.',
    )
    return parser.parse_args()


def main():
    args = _get_args()
    with open(args.sv, 'r') as Fh:
        last_line = None
        for line in Fh:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            linelist = line.split()
            if linelist[3] == "INS:BETWEEN":
                if last_line:
                    print(last_line)
                last_line = line
            elif linelist[3] == "DEL:BETWEEN":
                if last_line:
                    line = merge_line(last_line, line)
                    print(line)
                    last_line = None
                else:
                    print(line)
            else:
                if last_line:
                    print(last_line)
                    last_line = None
                target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence = splitline(line, "INTERNAL")
                query_start, query_end = get_query(args.bam, target_name, target_start, target_end)  # if it is INTERNAL, get the precise insertion/deletion position in query.
                line = "\t".join([target_name, str(target_start), str(target_end), sv_type, str(sv_length), str(per_id), str(matching_bases), query_name, str(query_start), str(query_end), sequence])
                print(line)


def get_query(bam, target_name, target_start, target_end):
    '''
    Parse bam file using pysam, and get position alignment information.
    '''
    import pysam
    target_start -= 1
    samfile = pysam.AlignmentFile(bam, "rb")
    for num, read in enumerate(samfile.fetch(target_name, target_start, target_end)):
        if num == 0:
            position_pairs = read.get_aligned_pairs()
            start_list = [x[0] for x in position_pairs if x[1] == target_start]
            end_list = [x[0] for x in position_pairs if x[1] == target_end]
            if len(start_list) == 1 and len(end_list) == 1:
                final_start = start_list[0] + read.get_tag("QS") + 1
                final_end = end_list[0] + read.get_tag("QS")
        else:
            print("wrong line in %s" % str(read))
            sys.exit("I don't expect multiple hits!")

    return final_start, final_end


def splitline(line, type):
    '''
    Split the line to a list.

    if type == "INTERNAL", per_id and matching_bases are number.
    if type == "BETWEEN", per_id and matching_bases are further splitted using comma.
    '''
    linelist = line.split()
    if len(linelist) == 11:  # if the last sequence column is missing go to else.
        target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence = linelist
    else:
        target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end = linelist
        sequence = ''
    target_start = int(target_start)
    target_end = int(target_end)
    query_start = int(query_start)
    query_end = int(query_end)
    sv_length = int(sv_length)
    if type == "INTERNAL":
        per_id = float(per_id)
        matching_bases = int(matching_bases)
    elif type == "BETWEEN":
        per_id = per_id.split(',')
        per_id[0] = float(per_id[0])
        per_id[1] = float(per_id[1])
        matching_bases = matching_bases.split(',')
        matching_bases[0] = int(matching_bases[0])
        matching_bases[1] = int(matching_bases[1])
    return target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence


def merge_line(last_line, line):
    '''
    Merge 2 lines to 1 line. the first line must be INS:BETWEEN, the 2nd line must be DEL:BETWEEN.

    The 2 lines merged and replaced with 1 line as a new sv_type: REPLACE.
    For lines cannot merge, return last_line\nline
    '''
    target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence = splitline(line, "BETWEEN")
    last_target_name, last_target_start, last_target_end, last_sv_type, last_sv_length, last_per_id, last_matching_bases, last_query_name, last_query_start, last_query_end, last_sequence = splitline(last_line, "BETWEEN")
    # always! last_sv_length == last_query_end - last_query_start
    if target_name == last_target_name and target_start == last_target_start:
        # print("same target position!")
        if query_name == last_query_name and ((query_start == last_query_start and query_end == last_query_end) or query_end < query_start):
            # print("same query position, merge!")
            sv_type = "REPLACE"
            query_start = last_query_start
            query_end = last_query_end
            sv_length = str(sv_length) + "," + str(last_sv_length)
            sequence = sequence + "," + last_sequence
            per_id = str(per_id[0]) + "," + str(per_id[1])
            matching_bases = str(matching_bases[0]) + "," + str(matching_bases[1])
            line = "\t".join([target_name, str(target_start), str(target_end), sv_type, str(sv_length), str(per_id), str(matching_bases), query_name, str(query_start), str(query_end), sequence])
        else:
            print("I don't expect this!")
            ipdb.set_trace()
    else:
        line = last_line + "\n" + line
    return line


    # linelist[3] = "REPLACE"
    # linelist[8] = last_linelist[8]
    # linelist[9] = last_linelist[9]
    # linelist[4] = linelist[4] + "," + last_linelist[4]
    # linelist[10] = linelist[10] + "," + last_linelist[10]
    return line

if __name__ == '__main__':
    main()