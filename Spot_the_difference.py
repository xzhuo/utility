# with 2 similar but maybe not identical fasta files, which sequences of them are different?

import sys
from Bio import SeqIO

file1 = sys.argv[1]
file2 = sys.argv[2]
dict1 = {}  # the dictionary to store seq.id and seq from file1
dict2 = {}  # the dictionary to stor seq.id and seq from file2


def import_fasta(file, dict):
    for seq_rec in SeqIO.parse(file, "fasta"):
        if seq_rec.id in dict and seq_rec.seq != dict[seq_rec.id]:
            sys.exit("more than 1 sequence %s is contradict to each other in %s!\n" % (seq_rec.id, file))
        dict[seq_rec.id] = seq_rec.seq

import_fasta(file1, dict1)
import_fasta(file2, dict2)

# now all seqs are in 2 dict now:

same = {}  # the dict to store identical sequences in 2 fasta files
diff = {}  # the dict to store sequences with same name but different seq in 2 fasta files
dist = {}  # the dict to store different sequences (different name and different seq)
for id in dict1:
    if id in dict2:
        if dict1[id] == dict2[id]:
            same[id]['seq'] = dict1[id]
        else:
            diff[id]['seq1'] = dict1[id]
            diff[id]['seq2'] = dict2[id]
        del dict2[id]
    else:
        dist['seq1'][id] = dict1[id]
for id in dict2:
    dist['seq2'][id] = dict2[id]

print("there are %d sequences are identical in 2 files.\n." % len(same))
print("there are %d sequences with same name in 2 files are different.\n." % len(diff))
for id in diff:
    print("%s\n" % id)
print("there are %d sequences file1 not in file2.\n." % len(dist['seq1']))
for id in dist['seq1']:
    print("%s\n" % id)
print("there are %d sequences file2 not in file1.\n." % len(dist['seq2']))
for id in dist['seq2']:
    print("%s\n" % id)
