#!/usr/bin/python
from tinyfasta import FastaParser
import collections
import sys
import os
import argparse
import operator

#usage: python ~/github/GenomeAssembly/src/GenomeAssembly/ScaffoldRename.py --fasta scaffolds_FINAL.fasta --out scaffolds_FINAL_chr.fasta --target scaffold_chromosome.txt

def rev_compl(st):
    st = st.upper()
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(nn[n] for n in reversed(st))


def rename_scaffold(fa_file, out_file, target_file):
    """
    target_file:
    scaffold_1	Vu03	+
    scaffold_2	Vu05	+
    scaffold_3	Vu04	+
    scaffold_4	Vu07	+
    scaffold_5	Vu09	-
    scaffold_6	Vu06	+
    scaffold_7	Vu10	-
    scaffold_8	Vu01	-
    scaffold_9	Vu11	+
    scaffold_10	Vu08	+
    scaffold_11	Vu02	+
    scaffold_12	Vu02	+
    scaffold_13	Vu_un1	+
    scaffold_14	Vu_un2	+
    scaffold_15	Vu_un3	+
    scaffold_16	Vu_un4	+
    scaffold_17	Vu_un5	+
    scaffold_18	Vu_un6	+
    """
    ScaffoldChr = {}
    in_h = open(target_file, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        scaffold, Chr, strand = lines
        ScaffoldChr[scaffold] = (Chr, strand)
    in_h.close()

    DescSeq = {}
    DescLen = {}
    
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.lstrip(">").split()[0]
        seq = str(record.sequence)
        DescSeq[desc] = seq
        # sl = len(seq)
        # DescLen[desc] = sl
    
    ChrSeq = collections.defaultdict()
    for s in ScaffoldChr:
        if s == "scaffold_12":
            c, st = ScaffoldChr[s]
            seq = DescSeq[s]
            if c in ChrSeq:
                ChrSeq[c] = "%s%s%s" % (seq, "N"*100,  ChrSeq[c])
            else:
                ChrSeq[c] = seq
        elif s == "scaffold_11":
            c, st = ScaffoldChr[s]
            seq = DescSeq[s]
            if c in ChrSeq:
                ChrSeq[c] = "%s%s%s" % (ChrSeq[c], "N"*100,  seq)
            else:
                ChrSeq[c] = seq
        else:
            c, st = ScaffoldChr[s]
            seq = DescSeq[s]
            if st == "+":
                newSeq = seq
            elif st == "-":
                newSeq = rev_compl(seq)
            else:
                print("Please make sure that the strand is '-' or '+'.")
                sys.exit(1)
            ChrSeq[c] = newSeq

    chrs = sorted(list(ChrSeq.keys()))
    out_h = open(out_file, "w")
    for c in chrs:
        seq = ChrSeq[c]
        out_h.write(">%s\n%s\n" % (c, seq))
    out_h.close()




# def rename_scaffold(fa_file, out_file):
#     DescSeq = {}
#     DescLen = {}
#     out_h = open(out_file, "w")
#     for record in FastaParser(fa_file):
#         desc = str(record.description)
#         desc = desc.lstrip(">").split()[0]
#         seq = str(record.sequence)
#         DescSeq[desc] = seq
#         sl = len(seq)
#         DescLen[desc] = sl
    
#     sortLength = sorted(DescLen.items(), key=operator.itemgetter(1), reverse=True)
#     print(sortLength)
#     recordNum = len(sortLength)
#     count = 1
#     for i in range(recordNum):
#         d, l = sortLength[i]
#         seq = DescSeq[d]

#         if i != recordNum - 1:
#             out_h.write(">Vu%d\n%s\n" % (count, seq))
#             count += 1
#         else:
#             string = "N" * 100
#             seqs = seq.split(string)
#             c = 1
#             for s in seqs:
#                 out_h.write(">Vu_un%d\n%s\n" % (c, s))
#                 c += 1
#     out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Rename the scaffold name.")
    parser.add_argument("-f", "--fasta", help="The input fasta file")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-t", "--target", help="The file containing scaffold and corresponding chrmosome.")
    args = parser.parse_args()
    rename_scaffold(args.fasta, args.out, args.target)



if __name__ == "__main__":
    main()
            

    