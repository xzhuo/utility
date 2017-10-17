import refine_calledSV

SIZE = 50  # the size differential cut off used to filter if there is a sv in flanking region.


def main():
    # the input refinedSV file should be 0 based bed file
    args = refine_calledSV._get_args()
    distance = args.distance
    chrom_size_dict = args.chrom_size
    with open(args.sv, 'r') as Fh:
        for line in Fh:
            line = line.rstrip()
            linelist = line.split()
            target_chr = linelist[0]
            target_start = linelist[1]
            target_end = linelist[2]
            query_chr = linelist[7]
            query_start = linelist[8]
            query_end = linelist[9]
            chrom_size = chrom_size_dict[target_chr]
            target_start_flanking = target_start - distance if target_start > distance else 0
            target_end_flanking = target_end + distance if target_end + distance < chrom_size else chrom_size
            query_start_flanking = query_start - distance if query_start > distance else 0
            query_end_flanking = query_end + distance if query_end + distance < chrom_size else chrom_size
            aln_start_flanking, aln_end_flanking = refine_calledSV.get_query(args.bam, target_chr, target_start_flanking, target_end_flanking)
            if aln_start_flanking - query_start_flanking < SIZE and aln_end_flanking - query_end_flanking < SIZE:
                print("%s\t%d\t%d\t%d\t%d" % (line, target_start_flanking, target_end_flanking, aln_start_flanking, aln_end_flanking))
            else:
                print("Maybe a SV in flanking region?\n%s\t%d\t%d\t%d\t%d" % (line, target_start_flanking, target_end_flanking, aln_start_flanking, aln_end_flanking))
