#!/usr/bin/python
from tinyfasta import FastaParser
import sys
import os
import argparse
import collections
import re




def parser_fasta(fa_file):
    DescSeq = {}
    SeqLen = {}
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.strip(">").split(" ")[0]
        seq = str(record.sequence)
        DescSeq[desc] = seq
        sl = len(seq)
        SeqLen[desc] = sl
    return DescSeq, SeqLen


def dump_list(regions):
    uniqRegion = []
    for i in regions:
        for j in i:
            uniqRegion.append(j)
    return uniqRegion


def splice_length(uniqRegion, totalLength):
    SpliceRegions = []
    for i in range(len(uniqRegion)):
        if i == 0:
            r = (0, uniqRegion[i])
            j = i
        else:
            r = (uniqRegion[j], uniqRegion[i])
            j = i
        SpliceRegions.append(r)
    SpliceRegions.append((uniqRegion[-1], totalLength))
    return SpliceRegions

def rev_compl(st):
    st = st.upper()
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))


def replace_scaffold_sequence(sub_seq, qry_seq, region_file, out_file):
    SubSequence, SubLength = parser_fasta(sub_seq)
    QrySequence, QryLength = parser_fasta(qry_seq)


    SubQryRegions = collections.defaultdict(dict)
    SubRegion = collections.defaultdict(list)

    in_h = open(region_file, "r")
    out_h = open(out_file, "w")

    ### N string in subject
    Nstring = "N" 
    for s in SubSequence:
        seq = SubSequence[s]
        matchN = re.finditer(Nstring, seq)
        for m in matchN:
            print(s, m.start(), m.group())


    for line in in_h:
        lines = line.strip().split()
        subject, slength, sstart, send, strand,  query, qlength, qstart, qend = lines[:9]
        sstart = int(sstart)
        send = int(send)
        qstart = int(qstart)
        qend = int(qend)
        SubQryRegions[subject][(sstart, send)] = (query, qstart, qend, strand)
        SubRegion[subject].append([sstart, send])
    in_h.close()

    replacedContigs = set()
    for s in SubRegion:
        replacedContigs.add(s)
        TotalSequence = []
        totalLength = SubLength[s]
        regions = sorted(SubRegion[s])
        uniqRegion = dump_list(regions)
        SpliceRegions = splice_length(uniqRegion, totalLength)
        for r in SpliceRegions:
            seq = SubSequence[s][r[0]:r[1]]
            if r in SubQryRegions[s]:
                qu = SubQryRegions[s][r]
                query, qs, qe, strand = qu
                qseq = QrySequence[query][qs:qe]
                if strand == "+":
                    seq = qseq
                elif strand == "-":
                    seq = rev_compl(qseq)
                else:
                    print("please make sure that the strand is '-' or '+'.")
                    sys.exit(1)
            TotalSequence.append(seq)
        out_h.write(">%s\n%s\n" % (s, "".join(TotalSequence)))

    allSubs = set(list(SubSequence.keys()))
    residueSubs = allSubs - replacedContigs
    residueSubList = sorted(list(residueSubs))
    for s in residueSubList:
        seq = SubSequence[s]
        out_h.write(">%s\n%s\n" % (s, seq))
    
    out_h.close()




def main():
    parser = argparse.ArgumentParser(description="Replace the sequence.")
    parser.add_argument("-s", "--subject", help="The input subject fasta file.")
    parser.add_argument("-q", "--query", help="The input query fasta file.")
    parser.add_argument("-r", "--region", help="The region file containing query aligned to subject.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    replace_scaffold_sequence(args.subject, args.query, args.region, args.out)

if __name__ == "__main__":
    main()

