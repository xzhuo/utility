# with 2 similar but maybe not identical fasta files, which sequences of them are different?

import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import re
import statistics

RMout = sys.argv[1]
RMlib = sys.argv[2]
RMdict = {}  # the dictionary to store seq.id and seq from repeatmasker library

unwant = ("Low_complexity", "Simple_repeat")  # tuple with TEclass I don't care
# unwant = ("Low_complexity", "Simple_repeat", "scRNA", "snRNA", "tRNA", "srpRNA", "Satellite", "Satellite/centr")  # tuple with TEclass I don't care

TEdict = {}  # a dict for unique TE names
infile = open(RMout, "r")
line = infile.readline()
for line in infile:
    line = line.split()
    if not (line == [] or line[0] == 'SW' or line[0] == 'score'):  # for UCSC RMout
        chr = line[5]
        start = int(line[6])
        end = int(line[7])
        name = line[10]
        TEclass = line[11]
        RepEnd = int(line[14])
        score = float(line[1])
        strand = line[9]
        RepStart = int(line[13])
        RepLeft = int(line[15])
        if strand == '-':
            RepStart = int(line[15])
            RepLeft = int(line[13])
            # RepLeft = int(line[11].translate({ord('('): None, ord(')'): None}))  # if I have to go back
        length = end - start + 1
        RepLength = RepLeft + RepEnd
        if TEclass not in unwant:
            if name in TEdict:
                TEdict[name]['length'].append(RepLength)  # append value to the existing list in the dict of dict
            else:
                TEdict[name] = {'length': [RepLength]}  # create the dict of dict




# the constucted dict structure is
# {repname1:{'length':[length1],'div':[divergence1],...}, repname2:{'length':[length2],'div':[divergence2],...},...}
infile.close()


def import_seq(file, dict, format):  # to manipulate seq name. Otherwise better use SeqIO.to_dict
    for seq_rec in SeqIO.parse(file, format, generic_dna):
        if seq_rec.id.find('#') + 1:
            seq_name = seq_rec.id[:seq_rec.id.find('#')]
        else:
            seq_name = seq_rec.id
        if seq_name in dict and seq_rec.seq != dict[seq_name]:
            sys.exit("more than 1 sequence %s is not consistant to each other in %s!\n" % (seq_name, file))
        seq_seq = re.match(r'\D+', str(seq_rec.seq)).group()  # trailing number is not deleted under curtain circumstances, delete all the \d here!
        # above: group() convert match object to string
        dict[seq_name] = str(seq_seq)

import_seq(RMlib, RMdict, "embl")
# import_seq(RBlib, RBdict, "fasta")

# Or:
# RMdict = SeqIO.to_dict(SeqIO.parse(RMlib, "fasta", generic_dna))

for te in TEdict:
    if te in RMdict:
        TEdict[te]['consensuslen'] = len(RMdict[te])
    # elif te in RBdict:
    #     TEdict[te] = RBdict[te]
    else:
        # print("no!!! Could not find %s anywhere! Time to panic!" % te)
        TEdict[te]['consensuslen'] = statistics.median(TEdict[te]['length'])

for te in TEdict:
    print("%s\t%d" % (te, TEdict[te]['consensuslen']))





# # now all seqs are in 2 dict now:

# same = {}  # the dict to store identical sequences in 2 fasta files
# diff = {}  # the dict to store sequences with same name but different seq in 2 fasta files
# dist = {}  # the dict to store different sequences (different name and different seq)
# dist['seq1'] = {}  # the dict to store different sequences (different name and different seq)
# dist['seq2'] = {}  # the dict to store different sequences (different name and different seq)
# for id in RMdict:
#     if id in RBdict:
#         if RMdict[id] == RBdict[id]:
#             same[id] = {}
#             same[id]['seq'] = RMdict[id]
#         else:
#             diff[id] = {}
#             diff[id]['seq1'] = RMdict[id]
#             diff[id]['seq2'] = RBdict[id]
#         del RBdict[id]
#     else:
#         dist['seq1'][id] = RMdict[id]
# for id in RBdict:
#     dist['seq2'][id] = RBdict[id]

# print("there are %d sequences are identical in 2 files.\n." % len(same))
# print("there are %d sequences with same name in 2 files are different.\n." % len(diff))
# for id in diff:
#     print("%s\n" % id)
# print("there are %d sequences RMlib not in RBlib.\n." % len(dist['seq1']))
# for id in dist['seq1']:
#     print("%s\n" % id)
# print("there are %d sequences RBlib not in RMlib.\n." % len(dist['seq2']))
# for id in dist['seq2']:
#     print("%s\n" % id)
