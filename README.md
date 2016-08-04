# utility
My small scripts

* RMlib_RMout_TElen.py  
parse repeatmasker library and repeatmasker output, get TE subfamily consensus length from repeatmasker library consensus. If not found (due to merging process in RepeatMasker), use median total length of all TE records found in the repeatmasker output file.  
Note there is a buggish problem with biopython parsing RepeatMasker embl file, I have to tweak it a little bit.  
Used for iteres.
