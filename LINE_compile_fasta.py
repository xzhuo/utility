# merge LINE fragments to full length LINEs. Using EMBOSS matcher to align overlapping part.

import sys
import json

struc = sys.argv[1]  # LINE structure json file generated using LINE_compile_ruly.py
# RMlib = sys.argv[2]  # embl format repeatmasker library


fh = open(struc)
LINE_str = fh.read()
LINEdict = json.loads(LINE_str)


def get_seq(seq_list):
    fdfds


for name in LINEdict:
    if "other" in LINEdict[name] and len(LINEdict[name]) == 1:
        frag = LINEdict[name]['other']['te']
        seq_list = [frag]
        seq = get_seq(seq_list)
        print(">%s\n%s" % (name, seq))
        