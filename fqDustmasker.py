# take a pair of fastq files and run pydustmasker on them
import os
import sys
import argparse
import pydustmasker
from multiprocessing import Pool

def fq_array(fastq_file):
    """
    Read a fastq file and return an array of tuples (header, sequence, quality).
    """
    fq_ary = []
    with open(fastq_file, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            strand = f.readline().strip()  # This line is just to read the '+' line, not used
            quality = f.readline().strip()
            fq_ary.append((header, sequence, strand, quality))
    return fq_ary

def run_fq_dustmasker(fq_ary_1, fq_ary_2):
    """
    Run pydustmasker on a pair of fastq arrays and write the filtered results to files.
    """
    if len(fq_ary_1) != len(fq_ary_2):
        raise ValueError("The two fastq files must have the same number of reads.")

    filtered_fq_1 = []
    filtered_fq_2 = []

    for (header1, seq1, _, qual1), (header2, seq2, _, qual2) in zip(fq_ary_1, fq_ary_2):
        (id1, _, _) = header1.split()
        (id2, _, _) = header2.split()
        if id1 != id2:
            raise ValueError("Unmatched headers in the fastq files: {} vs {}".format(header1, header2))
    
        masker1 = pydustmasker.DustMasker(seq1)
        masker2 = pydustmasker.DustMasker(seq2)
        if masker1.n_masked_bases / len(seq1) < 0.5 and masker2.n_masked_bases / len(seq2) < 0.5:
            # If both sequences are not too low-complexity (< 50%), keep them
            filtered_fq_1.append((header1, seq1, qual1))
            filtered_fq_2.append((header2, seq2, qual2))
    return filtered_fq_1, filtered_fq_2

def main():
    parser = argparse.ArgumentParser(description='Filter a pair of fastq files using pydustmasker, removing low-complexity sequences by pair.')
    parser.add_argument('-1', '--fq1', type=str, required=True,
                        help='input first fastq file')
    parser.add_argument('-2', '--fq2', type=str, required=True,
                        help='input second fastq file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output prefix for the filtered fastq files')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')

    args = parser.parse_args()
    fastq_1 = os.path.abspath(args.fq1)
    fastq_2 = os.path.abspath(args.fq2)
    
    if not os.path.exists(fastq_1) or not os.path.exists(fastq_2):
        raise ValueError("--fq1 or --fq2 file does not exist!")
    fq_ary_1 = fq_array(fastq_1)
    fq_ary_2 = fq_array(fastq_2)
    if args.threads == 1:
        filtered_fq_1, filtered_fq_2 = run_fq_dustmasker(fq_ary_1[0:0], fq_ary_2[0:0])
    else:
        fq_ary_ary_1 = list(map(lambda i:[i], fq_ary_1))
        fq_ary_ary_2 = list(map(lambda i:[i], fq_ary_2))
        with Pool(args.threads) as pool:
            results = pool.starmap(run_fq_dustmasker, zip(fq_ary_ary_1, fq_ary_ary_2))
            nonempty_results = [(r[0][0],r[1][0]) for r in results if len(r[0]) >0 and len(r[1]) > 0]
            filtered_fq_1, filtered_fq_2 = zip(*nonempty_results)

    with open(f"{args.out}_dust_1.fastq", 'w') as out_fq1:
        for header, seq, qual in filtered_fq_1:
            out_fq1.write(f"{header}\n{seq}\n+\n{qual}\n")

    with open(f"{args.out}_dust_2.fastq", 'w') as out_fq2:
        for header, seq, qual in filtered_fq_2:
            out_fq2.write(f"{header}\n{seq}\n+\n{qual}\n")

if __name__ == '__main__':
    main()
