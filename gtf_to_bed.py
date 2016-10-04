# Extract 5UTR, genebody, exon, intron, 3UTR, distal promoter TSS-2kb,
# proximal promoter TSS-250bp and core promoter TSS+-35bp from gtf file.
#
# features in the gtf file includes exon, CDS, start_codon and stop_codon.
import sys
import re

gtf = sys.argv[1]
pre = gtf[:-4]
with open(gtf, "r") as input:
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

