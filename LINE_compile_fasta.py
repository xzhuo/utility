# merge LINE fragments to full length LINEs. Using EMBOSS water to align overlapping part.
import sys
import json
import tempfile
from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline

struc = sys.argv[1]  # LINE structure json file generated using LINE_compile_ruly.py
embl = sys.argv[2]  # embl format repeatmasker library


# water = "/Users/Xiaoyu/EMBOSS/EMBOSS-6.6.0/emboss/water"  # EMBOSS water
# embl = "/Users/Xiaoyu/Downloads/Libraries20090604/RepeatMaskerLib.utf8.embl"  # RM library embl file
# embl = "/Users/Xiaoyu/Downloads/Libraries20140131/RepeatMaskerLib.embl"  # RM library embl file
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
        # Shitty biopython does not have handy stuff like locatableseq in bioperl.
        # So I have to parse it to find the start and end.
        astart = None
        aend = None
        bstart = None
        bend = None
        for line in fh:
            if not line[:1] == "#":
                linesplit = line.split()
                if len(linesplit) == 4:
                    if linesplit[0] == aname:
                        start = int(linesplit[1])
                        end = int(linesplit[3])
                        if astart is None or astart > start:
                            astart = start
                        if aend is None or aend < end:
                            aend = end
                    if linesplit[0] == bname:
                        start = int(linesplit[1])
                        end = int(linesplit[3])
                        if bstart is None or bstart > start:
                            bstart = start
                        if bend is None or bend < end:
                            bend = end
        return(astart, aend, bstart, bend)

    def water_c(RMfile, aname, bname):
        with tempfile.NamedTemporaryFile('w+') as tempout:
            astr = "%s:%s" % (RMfile, aname)
            bstr = "%s:%s" % (RMfile, bname)
            # bstr = "asis::%s" % bseq
            water_cline = WaterCommandline(r"/bar/xzhuo/EMBOSS/EMBOSS-6.6.0/emboss/water",
                                           asequence=astr,
                                           bsequence=bstr,
                                           gapopen=16, gapextend=4, aformat="pair",
                                           outfile=tempout.name)
            # print(water_cline)
            stdout, stderr = water_cline()
            astart, aend, bstart, bend = Align2_pos(tempout, aname, bname)
            return (astart, aend, bstart, bend)

    def get_seq(seq_list, RMfile, RMdict):
        if len(seq_list) == 2:
            aname = seq_list[0]
            bname = seq_list[1]
            aseq = str(RMdict[aname].seq)
            bseq = str(RMdict[bname].seq)
            # call water
            astart, aend, bstart, bend = water_c(RMfile, aname, bname)
            # if aend == len(aseq) and bstart == 1 and bend <= 300 and aend - astart <= 300:
                # seq = aseq + bseq[bend:]
                # return seq

            seq = aseq + bseq[bend:]
            return seq
            # else:
            #     print("no! something wrong with %s and %s!!" % (aname, bname))
            #     return "NA"

        if len(seq_list) == 3:
            aname = seq_list[0]
            bname = seq_list[1]
            cname = seq_list[2]
            aseq = str(RMdict[aname].seq)
            bseq = str(RMdict[bname].seq)
            cseq = str(RMdict[cname].seq)
            # call water
            astart, aend, bstart, bend = water_c(RMfile, aname, bname)
            cstart, cend, dstart, dend = water_c(RMfile, bname, cname)
            # if (aend == len(aseq) and bstart == 1 and bend <= 300 and aend - astart <= 300) and (
            #    cend == len(bseq) and dstart == 1 and dend <= 300 and cend - cstart <= 300):
                # seq = aseq + bseq[bend:cstart] + cseq
                # return seq

            seq = aseq + bseq[bend:cstart] + cseq
            return seq
            # else:
            #     print("no! somthing wrong with %s, %s and %s!!" % (aname, bname, cname))
            #     return "NA"

    outfa = open(struc+'.fa', 'w')
    for name in LINEdict:
        if "other" in LINEdict[name] and len(LINEdict[name]) == 1:
            frag = LINEdict[name]['other']['te']
            seq = str(RMdict[frag].seq)
            outfa.write(">%s\n%s\n" % (frag, seq))

        elif "5end" in LINEdict[name] and "orf2" in LINEdict[name] and "3end" in LINEdict[name]:
            frag1 = LINEdict[name]['5end']['te']
            frag2 = LINEdict[name]['orf2']['te']
            frag3 = LINEdict[name]['3end']['te']
            seq_list = [frag1, frag2, frag3]
            seq = get_seq(seq_list, tempembl.name, RMdict)
            outfa.write(">%s\n%s\n" % (name, seq))
        elif "other" in LINEdict[name] and "3end" in LINEdict[name] and (name[:2] == "L2" or name[:2] == "L3"):
            frag1 = LINEdict[name]['other']['te']
            frag2 = LINEdict[name]['3end']['te']
            seq_list = [frag1, frag2]
            seq = get_seq(seq_list, tempembl.name, RMdict)
            outfa.write(">%s\n%s\n" % (name, seq))
        elif "other" in LINEdict[name] and "5end" in LINEdict[name] and name == "MusHAL1":
            frag1 = LINEdict[name]['5end']['te']
            frag2 = LINEdict[name]['other']['te']
            seq_list = [frag1, frag2]
            seq = get_seq(seq_list, tempembl.name, RMdict)
            outfa.write(">%s\n%s\n" % (name, seq))
        else:
            print("can't get sequence for %s" % name)
