#!/usr/bin/python
from tinyfasta import FastaParser
from Bio.Seq import Seq
import argparse
import sys


#usage: python ~/github/GenomeAssembly/src/GenomeAssembly/mergeONTcontig.py --input /home/wuzhikun/Project/Vigna/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm.fasta --out ONT_merged_contigs.fasta --target /home/wuzhikun/Project/Vigna/pipeline/mashmap/HiFi2ONTContig_95_50kb/ont_merge_contig.txt

def merged_contig(contig_list):
    AllList = []
    Contig = []
    in_h = open(contig_list, "r")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        infor = tuple(lines[:2])
        if line.startswith("#"):
            continue
        else:
            if lines == ['']:
                if Contig != []:
                    AllList.append(Contig)
                Contig = []
            else:
                Contig.append(infor)
    AllList.append(Contig)
    in_h.close()
    return AllList

def rev_compl(st):
    st = st.upper()
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))


def merge_contig_sequence(contig_file, out_file, target):
    """
    target:
    ctg000070	+
    ctg000190	-
        
    ctg000180	+
    ctg000320	+
    ctg000050	-
        
    ctg000400	+
    ctg000250	+
    #ctg000200	-
        
    ctg000040	+
    ctg000290	-
    ctg000240	-
    ctg000370	-
        
    ctg000360	-
    ctg000010	-
    ctg000160	-
        
    ctg000060	-
    ctg000310	+
    ctg000340	+
    ctg000110	+
        
    ctg000410	-
    ctg000350	+
        
    ctg000100	-
    ctg000120	+
    ctg000090	+
        
    ctg000150	+
    ctg000380	+
        
    ctg000080	-
    ctg000220	+
    ctg000270	-
    ctg000000	-
    ctg000230	-
    ctg000170	+
    ctg000300	-
    """
    AllList = merged_contig(target)
    print(AllList)
    MergeUnit = set()
    DescSeq = {}
    for record in FastaParser(contig_file):
        desc = str(record.description)
        desc = desc.strip(">").split()[0]
        seq = str(record.sequence)
        DescSeq[desc] = seq
    
    out_h = open(out_file, "w")
    
    for contig in AllList:
        Sequences = []
        Names = []
        for a in contig:
            name, strand = a
            Names.append(name)
            MergeUnit.add(name)
            try:
                seq = DescSeq[name]
            except KeyError:
                print("Please check whether the seq %s is in file %s." % (name, contig_file))
                sys.exit(1)
            if strand == "+":
                newSeq = seq.upper()
            elif strand == "-":
                newSeq = rev_compl(seq)
            else:
                print("Please make sure that the strand is '+' or '-'.")
                sys.exit(1)

            if name == "ctg000050":
                newSeq = "%s%s" % (rev_compl(seq[13564006:]), seq[:13564006])
            
            if name == "ctg000200":
                newSeq = rev_compl(seq[38531686:])
            Sequences.append(newSeq)
        out_h.write(">%s\n%s\n" % ("_".join(Names), "N".join(Sequences)))
    

    print(MergeUnit)
    for c in DescSeq:
        if c in MergeUnit:
            continue
        else:
            s = DescSeq[c]
            out_h.write(">%s\n%s\n" % (c, s))

    for name in DescSeq:
        if name == "ctg000200":
            seq = DescSeq[name]
            newSeq = seq[:38531686]
            out_h.write(">%s\n%s\n" % (name, newSeq))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description='Merge some contigs.')
    parser.add_argument("-t", "--target", help="")
    parser.add_argument('-i', '--input', help='The input assembly contigs.')
    parser.add_argument('-o', '--out', help='The output file.')
    args = parser.parse_args()
    merge_contig_sequence(args.input, args.out, args.target)


if __name__ == '__main__':
    main()

