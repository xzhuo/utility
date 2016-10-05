# Extract 5UTR, genebody, exon, intron, 3UTR, distal promoter TSS-2kb,
# proximal promoter TSS-250bp and core promoter TSS+-35bp from gtf file.
#
# features in the gtf file includes exon, CDS, start_codon and stop_codon.
import sys
import re

gtf = sys.argv[1]
pre = gtf[:-4]
curr_transcript = ''
curr_strand = ''
curr_start = 0
curr_end = 0
curr_start_codon = 0
curr_stop_codon = 0
exon_list = []


def print_to_bed():
    genebody.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_start, curr_end, curr_strand))
    if curr_strand == '+':
        utr5.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_start, curr_start_codon-1, curr_strand))
        utr3.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_stop_codon+1, curr_end, curr_strand))
        distal.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_start-2000, curr_start-1, curr_strand))
        proximal.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_start-350, curr_start-1, curr_strand))
        core.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_start-35, curr_start+34, curr_strand))
    else:
        utr5.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_start_codon+1, curr_end, curr_strand))
        utr3.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_start, curr_stop_codon-1, curr_strand))
        distal.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_end+1, curr_end+2000, curr_strand))
        proximal.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_end+1, curr_end+350, curr_strand))
        core.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, curr_end-34, curr_end+35, curr_strand))
    intron_start = None
    for num, exon_int in enumerate(exon_list):
        exon.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, exon_int[0], exon_int[1], curr_strand))
        if num:  # if num > 0, means it is at least 2nd exon.
            intron.write("%s\t%d\t%d\t.\t.\t%s\n" % (curr_chrom, intron_start, exon_int[0]-1, curr_strand))
        intron_start = exon_int[1]+1


with open(gtf, "r") as input:
    utr5 = open(pre+".5utr.bed", "w")
    utr3 = open(pre+".3utr.bed", "w")
    genebody = open(pre+".genebody.bed", "w")
    exon = open(pre+".exon.bed", "w")
    intron = open(pre+".intron.bed", "w")
    distal = open(pre+".distal_promoter.bed", "w")
    proximal = open(pre+".proximal_promoter.bed", "w")
    core = open(pre+".core_promoter.bed", "w")
    for line in input:
        linelist = line.split('\t')  # tab delimited file
        chrom = linelist[0]
        features = linelist[2]
        start = linelist[3]
        end = linelist[4]
        strand = linelist[6]
        attributes = linelist[8]
        searchobj = re.search(r'gene_id \"(.+)\"; transcript_id \"(.+)\";', attributes)
        gene = searchobj.group(1)
        transcript = searchobj.group(2)
        if features == 'start_codon':
            curr_start_codon = start if strand == '+' else end
        if features == 'stop_codon':
            curr_stop_codon = end if strand == '+' else start
        if features == 'exon':
            exon_list.append((start, end))
        if curr_transcript == transcript:
            curr_start = start if start < curr_start else curr_start
            curr_end = end if end > curr_end else curr_end

        else:
            if curr_transcript.length > 0:
                print_to_bed()
            curr_transcript = transcript
            curr_chrom = chrom
            curr_strand = strand
            curr_start = start
            curr_end = end
    print_to_bed()
    utr5.close()
    utr3.close()
    genebody.close()
    exon.close()
    intron.close()
    distal.close()
    proximal.close()
    core.close()
