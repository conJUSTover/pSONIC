#!/usr/bin/python

"""Usage: .py <prefix> 
Also must have files named Orthogroups.csv, all.gff, all.tandem in currect directory."""



from __future__ import print_function
import sys
import os
import shutil
import subprocess
import time
from igraph import Graph
from multiprocessing import Pool
from functools import partial
from pathlib import Path

def SingleCopy_from_Orthofinder(orthogroups_file, tandem_net, tandem_list, species_ploidy):
    singleton_start = time.time()
    print("Starting Singleton Search")
    singletons = []
    with open(orthogroups_file, "r") as handle: #Orthogroups.csv
        first_line = handle.readline().strip() #Species Identifiers
        species_num = len(first_line.split("\t"))
        for line in handle: #starting with line 2, go through all of the orthogroups
            gene_size = []
            line = line.strip().split("\t") #split orthogroup into species-delimited ortholog sets
            gene_size.append(line[0])
            for i in range(len(line[1:]), 0, -1): #line[1:] because this includes OG, whereas species_num does not
                genes = line[i].strip().split(", ") #creates a list of every gene for a given species
                genes = [k.strip() for k in genes if k]
                if genes:
                    group_size, group_genes = get_group_size(genes, tandem_net, tandem_list)#More than one gene for a species 
                    gene_size.insert(1,str(group_size))
                else: gene_size.insert(1,"0")
            while len(gene_size) < species_num + 1:
                gene_size.insert(1,"0")
            singletons.append(gene_size)
    print("Initial filtering done: " + str((time.time() - singleton_start)/60))
    return species_num, singletons

def get_group_size(genes, tandem_net, tandem_list):
    if not genes: return 0, []
    gene_groups = []
    while genes:
        test_gene = genes.pop()
        if find_group(test_gene, gene_groups):
            continue #skip to next gene in list
        potential_tandem_group = test_group(test_gene, tandem_net, tandem_list)
        if potential_tandem_group:
            gene_groups.append(potential_tandem_group)
        else: gene_groups.append([test_gene])
    return_group = [item for sublist in gene_groups for item in sublist]
    return_group = [i for i in return_group if i]
    return len(gene_groups), return_group

def test_group(a, network, edges):
    if a in edges:
        return find_group(a, network)
    else:
        return []

def find_group(a, lists):
    for i in lists:
        if a in i:
            return i
    return []

def edges_to_groups(edge_list, gene_names):
    add_names_part = partial(add_names, gene_list=gene_names)
    with Pool(processes=int(sys.argv[2])) as pool:
        edge_list = pool.map(add_names_part, edge_list)

    master_graph = Graph(n=len(gene_names), edges = [i[-2:] for i in edge_list])
    master_graph.vs["name"] = gene_names

    singletons = master_graph.decompose(minelements=2)
    singletons = [i.vs["name"] for i in singletons]
    return singletons

def print_list(outfile, tandem_out):
    with open(outfile, "w+") as handle:
        for tandem in tandem_out:
            print_string = '\t'.join(tandem)
            handle.write(print_string + '\n')

def load_coded_gene_names(sequences_file):
    code_gene = {}
    gene_code = {}
    with open(sequences_file, "r") as handle:
        for line in handle:
            line = line.strip().split(": ")
            gene_code[line[1]] = line[0]
            code_gene[line[0]] = line[1]
    return code_gene, gene_code

def add_names(edge, gene_list):
    edge.append(gene_list.index(edge[0]))
    edge.append(gene_list.index(edge[1]))
    return edge

def main():
    prefix = sys.argv[1]
    global gff_genes
    start_time = time.time()

    code_gene, gene_code = load_coded_gene_names("SequenceIDs.txt")

    with open(prefix + ".gff", "r") as handle:
        gffs = [i.strip().split('\t') for i in handle]
    gene_names = [i[1] for i in gffs]

    with open(prefix + ".tandem", "r") as handle:
        tandem_genes = [line.strip().split(',') for line in handle]

    tandem_network = edges_to_groups(tandem_genes, gene_names)
    tandem_network.sort(key=len, reverse=True)

    tandem_pairs = []
    for i in tandem_network:
        tandem_pairs += [j for j in i]
    print("Tandem clustering done: " + str((time.time() - start_time)/60))

    ploidy_file = Path(prefix + ".ploidy")
    if ploidy_file.is_file():
        print("Ploidy file found")
        with open(ploidy_file, "r") as handle:
           ploidies = handle.readline().strip().split('\t')
    else:
        print("No ploidy file found. Now assuming all species have ploidy of one.")
        ploidies = False

    species_num, singletons = SingleCopy_from_Orthofinder("Orthogroups.csv", tandem_network, tandem_pairs, ploidies)
    print("Singletons All Found: " + str((time.time() - start_time)/60))

    print_list("Orthogroups_gene_set_size.txt", singletons)

gff_genes = []

main()
