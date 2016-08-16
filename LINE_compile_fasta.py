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

    def Align2_pos(fh, aname, bname):
        for line in fh:
            if not line[:1] == "#":
                linesplit = line.split()
                if len(linesplit) == 4:
                    if linesplit[0] == aname:
                        start = linesplit[1]
                        end = linesplit[3]

    def water_c(RMfile, aname, bname):
        with tempfile.NamedTemporaryFile('w+') as tempout:
            astr = "%s:%s" % (RMfile, aname)
            bstr = "%s:%s" % (RMfile, bname)
            # bstr = "asis::%s" % bseq
            water_cline = WaterCommandline(r"/Users/Xiaoyu/EMBOSS/EMBOSS-6.6.0/emboss/water",
                                           asequence=astr,
                                           bsequence=bstr,
                                           gapopen=16, gapextend=4, aformat="pair",
                                           outfile=tempout.name)
            # print(water_cline)
            stdout, stderr = water_cline()

            # Shitty biopython does not have handy stuff like locatableseq in bioperl.
            # So I have to parse it to find the start and end.
            astart, aend, bstart, bend = Align2_pos(tempout, aname, bname)

    def get_seq(seq_list, RMfile, RMdict):
        if len(seq_list) == 2:
            # call water
            aname = seq_list[0]
            bname = seq_list[1]
            aseq = str(RMdict[aname].seq)
            bseq = str(RMdict[bname].seq)
            astart, aend, bstart, bend = water_c(RMfile, aname, bname)

            if aend == RMdict[aname].length and bstart == 1 and bend <= 180 and aend - astart <= 180: 
                seq = aseq[:astart-1] + bseq
                return seq
            else:
                print("no! terribly wrong!!")
                return "NA"




    for name in LINEdict:
        if "other" in LINEdict[name] and len(LINEdict[name]) == 1:
            frag = LINEdict[name]['other']['te']
            seq = str(RMdict[frag].seq)
            print(">%s\n%s" % (frag, seq))

        elif "5end" in LINEdict[name] and "orf2" in LINEdict[name] and "3end" in LINEdict[name]:
            frag1 = LINEdict[name]['5end']['te']
            frag2 = LINEdict[name]['orf2']['te']
            frag3 = LINEdict[name]['3end']['te']
            seq_list = [frag1, frag2, frag3]
            seq = get_seq(seq_list, tempembl.name, RMdict)
            print(">%s\n%s" % (name, seq))
        elif "other" in LINEdict[name] and "3end" in LINEdict[name] and name[:2] == "L2":
            frag1 = LINEdict[name]['other']['te']
            frag2 = LINEdict[name]['3end']['te']
            seq_list = [frag1, frag2]
            seq = get_seq(seq_list, tempembl.name, RMdict)
            print(">%s\n%s" % (name, seq))
        else:
            print("can't get sequence for %s" % name)