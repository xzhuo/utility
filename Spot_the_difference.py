# with 2 similar but maybe not identical fasta files, which sequences of them are different?

import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import re

file1 = sys.argv[1]
file2 = sys.argv[2]
dict1 = {}  # the dictionary to store seq.id and seq from file1
dict2 = {}  # the dictionary to stor seq.id and seq from file2


def import_fasta(file, dict, format):  # to manipulate seq name. Otherwise better use SeqIO.to_dict
    for seq_rec in SeqIO.parse(file, format, generic_dna):
        if seq_rec.id.find('#') + 1:
            seq_name = seq_rec.id[:seq_rec.id.find('#')]
        else:
            seq_name = seq_rec.id
        if seq_name in dict and seq_rec.seq != dict[seq_name]:
            sys.exit("more than 1 sequence %s is contradict to each other in %s!\n" % (seq_name, file))
        seq_seq = re.match(r'\D+', str(seq_rec.seq)).group()  # trailing number is not deleted under curtain circumstances, delete all the \d here!
        # above: group() convert match object to string
        print("%d" % len(seq_seq))
        dict[seq_name] = str(seq_seq)

import_fasta(file1, dict1, "fasta")
import_fasta(file2, dict2, "embl")

# Or:
# dict1 = SeqIO.to_dict(SeqIO.parse(file1, "fasta", generic_dna))

# now all seqs are in 2 dict now:

same = {}  # the dict to store identical sequences in 2 fasta files
diff = {}  # the dict to store sequences with same name but different seq in 2 fasta files
dist = {}  # the dict to store different sequences (different name and different seq)
dist['seq1'] = {}  # the dict to store different sequences (different name and different seq)
dist['seq2'] = {}  # the dict to store different sequences (different name and different seq)
for id in dict1:
    if id in dict2:
        if dict1[id] == dict2[id]:
            same[id] = {}
            same[id]['seq'] = dict1[id]
        else:
            diff[id] = {}
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
