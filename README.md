# utility
My small scripts

* RMlib_RMout_TElen.py  
Parse repeatmasker library and repeatmasker output to get TE subfamily consensus length. Search repeatmasker library first to find sequence consensus. If not found (due to merging process in RepeatMasker), use median total length of all TE records found in the repeatmasker output file as length.  
Note there is a buggish problem parsing RepeatMasker embl file with biopython, I have to tweak it a little bit.  
Used for iteres.
