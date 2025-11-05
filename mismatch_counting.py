import subprocess as sp
import shlex
import pysam
import sys
import re

#bamFile = "test.extract.bam"
bamFile = sys.argv[1]
outFile = sys.argv[2]

baseInfo = {}
for read in pysam.AlignmentFile(bamFile):
    readName = read.query_name
    chrom = read.reference_name
    pos = read.pos + 1
    seq = read.seq
    MDZ = read.get_tag("MD")
    cigar = read.cigarstring
    md_match = re.findall(r'\d+', MDZ)
    md_mismatch = re.findall(r'\D+', MDZ)
    cg_num = re.findall(r'\d+',cigar)
    cg_lt = re.findall(r'\D+',cigar)
    baseRM = []
    tmpos = pos
    for n in range(len(cg_lt)):
        if cg_lt[n] == 'M':
            tmpos += int(cg_num[n])
        elif cg_lt[n] == 'D':
            baseRM += range(tmpos - 3, tmpos + int(cg_num[n]) + 4)
            tmpos += int(cg_num[n])
        elif cg_lt[n] == 'I':
            baseRM += range(tmpos - 3, tmpos+4)
    for n in range(len(md_mismatch)):
        ref = md_mismatch[n]
        i = 0
        tmpos2 = pos
        while i <= n:
            base = tmpos2 + int(md_match[i])
            tmpos2 = base + len(re.sub("\^","",md_mismatch[i]))
            i = i + 1
        tmpkey = "_".join([chrom, str(base)])
        if base not in baseRM:
            try:
                baseInfo[tmpkey][1] += 1
            except KeyError:
                baseInfo[tmpkey] = [ref, 1]

with open(outFile,'a') as fout:
    for key in baseInfo:
        tmpkey = key.split("_")
        fout.write("%s\t%s\t%s\t%d\n" %(tmpkey[0],tmpkey[1],baseInfo[key][0],baseInfo[key][1]))
