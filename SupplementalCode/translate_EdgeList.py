#!/usr/bin/env python
#usage: python .py pSONIC.EdgeList.txt SequenceIDs.txt

import sys 

def translate(edgelist, seqIDs):
    code_gene, gene_code = load_coded_gene_names(seqIDs)
    output = []
    with open(edgelist, "r") as handle:
        for line in handle:
            line = line.strip().split("\t")
            line[1] = code_gene[line[1]]
            line[0] = code_gene[line[0]]
            output.append(line)
    with open("pSONIC.EdgeList.translated.txt", "w") as handle:
        for l in output:
            j = "\t".join(l)
            handle.write('%s\n' % j)

def load_coded_gene_names(sequences_file):
    code_gene = {}
    gene_code = {}
    with open(sequences_file, "r") as handle:
        for line in handle:
            line = line.strip().split(": ")
            gene_code[line[1]] = line[0]
            code_gene[line[0]] = line[1]
    return code_gene, gene_code

def main():
    edges = sys.argv[1]
    seqs = sys.argv[2]
    translate(edges, seqs)

main()
