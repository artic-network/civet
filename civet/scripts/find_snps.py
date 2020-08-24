#!/usr/bin/env python3
import os
import argparse
import collections
from Bio import SeqIO

cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Find sequences relative to Wuhan4 reference.')

    parser.add_argument("--input", action="store", type=str, dest="input")
    parser.add_argument("--output", action="store", type=str, dest="output")
    parser.add_argument("--tree", action="store", type=str, dest="tree")

    return parser.parse_args()

def find_snps():

    args = parse_args()

    input_seqs = collections.defaultdict(list)
    outgroup_seq = ""
    
    for record in SeqIO.parse(args.input, "fasta"):
        if record.id == "outgroup":
            outgroup_seq = record.seq.upper()
        else:
            input_seqs[str(record.seq).upper()].append(record.id)

    non_amb = ["A","T","G","C"]

    with open(args.output, "w") as fw:

        for query_seq in input_seqs:
            snps =[]

            for i in range(len(query_seq)):
                bases = [query_seq[i],outgroup_seq[i]]
                if bases[0] != bases[1]:
                    if bases[0] in non_amb and bases[1] in non_amb:
                        
                        snp = f"{i+1}{bases[1]}{bases[0]}" # position-outgroup-query
                        snps.append(snp)
            
            for record_name in input_seqs[query_seq]:
                print(record_name, len(snps))
                snp_str = ";".join(snps)
                fw.write(f"{record_name}\t{args.tree}\t{len(snps)}\t{snp_str}\n")

if __name__ == '__main__':

    find_snps()