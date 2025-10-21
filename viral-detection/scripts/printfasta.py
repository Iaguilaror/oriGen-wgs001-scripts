import pysam
import sys
import numpy as np
bam = pysam.AlignmentFile(sys.argv[1], "rc")
fcs = np.zeros(4, dtype=int)
for a in bam.fetch(until_eof=True):
    if not (a.flag & 0x100):
        fcs[((a.flag & (0x4|0x8)) >> 2)] += 1
        if (a.flag & 0x4) or (a.flag & 0x8):
            print(">" + a.qname + f"-{a.flag}:{a.mapping_quality}:{a.reference_name or '*'}:{a.reference_start+1}:{a.next_reference_name or '*'}:{a.next_reference_start+1}:{len(a.query_sequence)}\n{a.query_sequence}")

if sys.argv[2]:
    with open(sys.argv[2], "w") as f:
        print(".".join(sys.argv[1].split("/")[-1].split(".")[:-1]), "\t".join(map(str, fcs)), sep="\t", file=f)