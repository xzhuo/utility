# merge LINE fragments to full length LINEs. Using EMBOSS water to align overlapping part.

import sys
import os
import json
import re
import tempfile
from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline

struc = sys.argv[1]  # LINE structure json file generated using LINE_compile_ruly.py
# RMlib = sys.argv[2]  # embl format repeatmasker library


# water = "/Users/Xiaoyu/EMBOSS/EMBOSS-6.6.0/emboss/water"  # EMBOSS water
embl = "/Users/Xiaoyu/Downloads/Libraries20140131/RepeatMaskerLib.embl"  # RM library embl file

fh = open(struc)
LINE_str = fh.read()
LINEdict = json.loads(LINE_str)

# create a temp embl file because EMBOSS uniform sequence address method can't handle commend lines before 1st entry.
emblFh = open(embl)
with tempfile.NamedTemporaryFile('w+') as tempembl:
    for line in emblFh:
        if line[:2] == 'ID':
            tempembl.write(line)
            break
        # print(line)

    for line in emblFh:
        tempembl.write(line)
    # print(tempembl.readline())
    # print(tempembl.name)
    tempembl.seek(0)
    emblname = tempembl.name
    RMdict = SeqIO.to_dict(SeqIO.parse(tempembl, "embl"))

    def get_seq(seq_list, RMdict):  # the 2nd option: for single seqence, it is the RMdict, for 2 or more sequences list, it is the tempembl
        if len(seq_list) == 1:
            seqname = seq_list[0]
            seq = str(RMdict[seqname].seq)
            return seq

        if len(seq_list) == 2:
            # call water
            with tempfile.NamedTemporaryFile('w+') as tempout:
                aname = seq_list[0]
                bname = seq_list[1]
                aseq = "%s:%s" % (RMdict, aname)
                bseq = "%s:%s" % (RMdict, bname)
                # bstr = "asis::%s" % bseq
                water_cline = WaterCommandline(r"/Users/Xiaoyu/EMBOSS/EMBOSS-6.6.0/emboss/water",
                                               asequence=aseq,
                                               bsequence=bseq,
                                               gapopen=16, gapextend=4, aformat="pair",
                                               outfile=tempout.name)
                print(water_cline)
                stdout, stderr = water_cline()

                # Shitty biopython does not have handy stuff like locatableseq in bioperl. So I have to parse it to find the start and end.
                for eachline in tempout:
                    if not eachline[:1] == "#":
                        linesplit = eachline.split()
                        if len(linesplit) == 4:
                            if linesplit[0] = aname:
                                



    for name in LINEdict:
        if "other" in LINEdict[name] and len(LINEdict[name]) == 1:
            frag = LINEdict[name]['other']['te']
            seq_list = [frag]
            seq = get_seq(seq_list, RMdict)
            print(">%s\n%s" % (name, seq))

        elif "5end" in LINEdict[name] and "orf2" in LINEdict[name] and "3end" in LINEdict[name]:
            frag1 = LINEdict[name]['5end']['te']
            frag2 = LINEdict[name]['orf2']['te']
            frag3 = LINEdict[name]['3end']['te']
            seq_list = [frag1, frag2, frag3]
            seq = get_seq(seq_list, tempembl.name)
            print(">%s\n%s" % (name, seq))
        elif "other" in LINEdict[name] and "3end" in LINEdict[name] and name[:2] == "L2":
            frag1 = LINEdict[name]['other']['te']
            frag2 = LINEdict[name]['3end']['te']
            seq_list = [frag1, frag2]
            seq = get_seq(seq_list, tempembl.name)
            print(">%s\n%s" % (name, seq))
        else:
            print("can't get sequence for %s" % name)