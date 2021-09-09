#!/usr/bin/python


from __future__ import print_function
import sys
import time
import argparse
from igraph import Graph
from multiprocessing import Pool
from functools import partial
from pathlib import Path

def process_group(block, good_bad, good_groups, gff_list):
    if block: #if string has values in it (not valid for header lines in file)
        a_start = [i for i, colour in enumerate(gff_list) if block[0][0] in colour][0]
        b_start = [i for i, colour in enumerate(gff_list) if block[0][1] in colour][0]
        a_end = [i for i, colour in enumerate(gff_list) if block[-1][0] in colour][0]
        b_end = [i for i, colour in enumerate(gff_list) if block[-1][1] in colour][0]
        good_pairs = good_bad.count("G")
        bad_pairs = good_bad.count("N")
        a_range = abs(a_end - a_start) + 1
        a_density = float(a_range) / len(block)
        b_range = abs(b_end - b_start) + 1
        b_density = float(b_range) / len(block)
        return [str(len(block)), str(a_start), str(a_end), str(a_range), str(a_density), str(b_start), str(b_end), str(b_range), str(b_density), block[0][0], block[0][1], block[-1][0], block[-1][1], str(b_density * a_density), str(good_pairs), str(bad_pairs), good_bad]


def SingleCopy_from_Orthofinder(orthogroups_file, tandem_net, tandem_list, species_ploidy, gene_code, species_IDs):
    singleton_start = time.time()
    print("Starting Singleton Search")
    singletons = []
    with open(orthogroups_file, "r") as handle: #Orthogroups.tsv
        first_line = handle.readline().strip() #Species Identifiers
        first_line = first_line.split("\t")
        if not first_line[0] in list(species_IDs.keys()):
            first_line = first_line[1:]
        species_nums = [species_IDs[i] for i in first_line if species_IDs[i]]
        if species_ploidy == False or species_ploidy == None:
            species_ploidy = [1] * len(species_nums)
        for line in handle: #starting with line 2, go through all of the orthogroups
            good_genes = []
            line = line.strip().split("\t") #split orthogroup into species-delimited ortholog sets
            for i in range(len(line[1:]), 0, -1): #line[1:] because this includes OG, whereas species_num does not
                genes = line[i].strip().split(", ") #creates a list of every gene for a given species
                genes = [k.strip() for k in genes if k]
                if genes:
                    genes = [gene_code[gene] for gene in genes]
                    group_size, group_genes = get_group_size(genes, tandem_net, tandem_list)#More than one gene for a species 
                    if group_size <= int(species_ploidy[i - 1]) and group_genes:
                        good_genes += group_genes
            if len(good_genes) > 1: 
                for j in range(1, len(good_genes)):
                    singletons.append([good_genes[j], good_genes[0]])
    print("Initial filtering done: " + str((time.time() - singleton_start)/60))
    return species_nums, singletons

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

def edges_to_groups(edge_list, gene_names, threads):
    add_names_part = partial(add_names, gene_list=gene_names)
    with Pool(processes=threads) as pool:
        edge_list = pool.map(add_names_part, edge_list)

    master_graph = Graph(n=len(gene_names), edges = [i[-2:] for i in edge_list])
    master_graph.vs["name"] = gene_names
    
    singletons = master_graph.decompose(minelements=2)
    singletons = [i.vs["name"] for i in singletons]
    return singletons

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

def trim_block_ends(block, pos_neg):
    block, pos_neg = trim_forward(block, pos_neg, 0)
    block.reverse()
    pos_neg=pos_neg[::-1]
    block, pos_neg = trim_forward(block, pos_neg, 0)
    block.reverse()
    pos_neg=pos_neg[::-1]
    if len(pos_neg) < 6:
        block = []
        pos_neg = []
    return block, pos_neg

def trim_forward(block, pos_neg, size):
    if not block: return [], ""
    if size == len(pos_neg): return [], ""
    if pos_neg[size] == "N":
        block = block[size+1:]
        pos_neg = pos_neg[size+1:]
        block, pos_neg = trim_forward(block, pos_neg, 0)
        return block, pos_neg
    if pos_neg[size] == "-":
        block, pos_neg = trim_forward(block, pos_neg, size + 1)
        return block, pos_neg
    temp_pos_neg = pos_neg.replace('-', '')
    g_start = temp_pos_neg[0:6].count("G")
    n_start = temp_pos_neg[0:6].count("N")
    if g_start > n_start:
        return block, pos_neg
    else:
        block = block[size+1:]
        pos_neg = pos_neg[size+1:]
        block, pos_neg = trim_forward(block, pos_neg, 0)
        return block, pos_neg



def read_collinearity(collin_file):
    with open(collin_file, "r") as handle:
        infile = handle.read().split("## A")
        return infile[1:]

def score_all_blocks(collin_blocks, singleton_groups, good_groups, single_set, gff_list, spec_nums, gene_notandem, threads):
    process_start = time.time()
    edges_trimmed = []
    edge_list = []
    raw_stats = []
    score_part = partial(score_block, singleton_groups=singleton_groups, good_groups=good_groups, single_set=single_set, gff_list=gff_list, spec_nums=spec_nums, gene_notandem=gene_notandem)
    with Pool(processes=threads) as pool:
        scores = pool.map(score_part, collin_blocks)
    scores = [item for item in scores]
    done_time = time.time()
    print("multithreading done: " + str((done_time - process_start)/60))
    for score in scores:
        if len(score) > 1:
            raw_stats.append(score[0])
            edges_trimmed += score[1]
            edge_list += score[2]
            good_groups += score[3:]
        else: raw_stats.append(score[0])
    print("for loop done: " + str((time.time() - done_time)/60))
    return edges_trimmed, raw_stats, edge_list, good_groups

def score_block(block, singleton_groups, good_groups, single_set, gff_list, spec_nums, gene_notandem):
    block = block.strip().split('\n')
    header = block.pop(0).split(' ')
    orientation = header[-1]
    chr_pairs = header[-2]
    align_num = header[1].strip(':')
    good_block = []
    pos_neg = ""
    result = []
    for line in block:
        line = line.strip().split('\t')
        line.append(gene_notandem.index(line[1]))
        line.append(gene_notandem.index(line[2]))
        singleton_match = test_group(line[1], singleton_groups, single_set)
        if singleton_match:
            #Cluster network here, then see if species has any genes in group
            clustered_singletons = cluster_species(singleton_match, spec_nums, False, False)
            species_pos = line[2].split('_')[0]
            species_id = spec_nums.index(species_pos)
#            species_id = int(line[2][0]) ######
            if line[2] in clustered_singletons[species_id]:
                pos_neg += "G"
            elif not clustered_singletons[species_id]: pos_neg += "-"
            else: pos_neg += "N"
        else:
           singleton_match = test_group(line[2], singleton_groups, single_set)
           if singleton_match:
               pos_neg += "N"
           else: pos_neg += "-"
        good_block.append([line[1], line[2], line[3], orientation, pos_neg[-1], chr_pairs, align_num, line[-2], line[-1]])
    result.append([str(len(pos_neg)), str(pos_neg.count("G")), str(pos_neg.count("N")), str(pos_neg.count("-")), pos_neg, orientation])
    if pos_neg.count("G") > 1 and pos_neg.count("G") > pos_neg.count("N"):
        new_block, pos_neg = trim_block_ends(good_block, pos_neg)
        if new_block:
            split, split_pn = split_block(new_block, pos_neg)
            edges_kept = [item for sublist in split for item in sublist]
            result.append([gene for gene in good_block if gene not in edges_kept]) #trimmed edges
            for i,n in enumerate(split):
                if i > 0:
                    for j in range(len(n)):
                        split[i][j][-3] = str(int(split[i][j][-3]) + (int(i)/100))
            edges_kept = [item for sublist in split for item in sublist]
            result.append(edges_kept) #edge list
            for i,n in enumerate(split):
                result.append(process_group(split[i], split_pn[i], good_groups, gff_list))
    return result

def split_block(group, pn):
    scored_pn = pn.replace('-', '')
    if "NNN" not in scored_pn:
        return [group], [pn]
    count_n = 0
    total_count = 0
    for i in pn:
        if count_n == 3: break
        total_count += 1
        if i == "N": count_n +=1
        elif i =="G": count_n = 0
    group_A = group[:total_count]
    group_Apn = pn[:total_count]
    group_A.reverse()
    group_Apn = group_Apn[::-1]
    group_A, group_Apn = trim_forward(group_A, group_Apn, 0)
    group_A.reverse()
    group_Apn = group_Apn[::-1]
    group_B = group[total_count:]
    group_Bpn = pn[total_count:]
    group_B, group_Bpn = trim_forward(group_B, group_Bpn, 0)
    group_B, group_Bpn = split_block(group_B, group_Bpn)
    group_B.insert(0, group_A)
    group_Bpn.insert(0, group_Apn)
    return group_B, group_Bpn

def print_ortho_process(group, spec_num, code_dict, tandem_net, tandem_list):
    group = sorted(group)
    tandem_score = []
    size = []
    orthogroup = []
    trans_orthogroup = []
    for s in range(len(spec_num)):
        j = [n for n in group if n.split('_')[0] == spec_num[s]]
        size.append(str(len(j)))
        orthogroup.append(','.join(j))
        group = [gene for gene in group if gene not in j]
        t = [translate_genes(k, code_dict) for k in j]
        trans_orthogroup.append(','.join(t))
        group_size, all_groups = get_group_size(j, tandem_net, tandem_list)
        tandem_score.append(str(group_size))
    if group:
        size.append(str(len(group)))
        orthogroup.append(','.join(group))
        t = [translate_genes(k, code_dict) for k in group]
        trans_orthogroup.append(','.join(t))
        group_size, all_groups = get_group_size(group, tandem_net, tandem_list)
        tandem_score.append(str(group_size))
    return [tandem_score, size, orthogroup, trans_orthogroup] 
    
def print_orthogroups_quick(network, spec_num, code_dict, tandem_net, tandem_list, threads):
    print_quick = partial(print_ortho_process, spec_num=spec_num, code_dict=code_dict, tandem_net=tandem_net, tandem_list=tandem_list)
    with Pool(processes=threads) as pool:
        result = pool.map(print_quick, network)
    tandem_out = [["OG" + str('{0:07}'.format(i))] + n[0] for i,n in enumerate(result)]
    size = [["OG" + str('{0:07}'.format(i))] + n[1] for i,n in enumerate(result)]
    orthogroup = [["OG" + str('{0:07}'.format(i))] + n[2] for i,n in enumerate(result)]
    trans_ortho = [["OG" + str('{0:07}'.format(i))] + n[3] for i,n in enumerate(result)]
    return tandem_out, size, orthogroup, trans_ortho


def print_list(outfile, tandem_out):
    with open(outfile, "w+") as handle:
        for tandem in tandem_out:
            print_string = '\t'.join(tandem)
            handle.write(print_string + '\n')

def cluster_species(network, spec_num, code_dict, trans_bool):
    #if trans_bool is False, then code_dict has no purpose here
    final_out = []
    for i in range(len(spec_num)):
#        final_out.append(sorted([n for n in network if n.startswith(str(i))]))
        final_out.append(sorted([n for n in network if n.split('_')[0] == spec_num[i]]))
        if trans_bool:
            #translate edge names into genes
            new_genes = [translate_genes(k, code_dict) for k in final_out[i]]
            final_out[i] = new_genes
    return final_out

def print_network(outfile, all_networks, spec_num, code_uncode, translate_bool):
    final_networks = []
    #group edges into species
    for network in all_networks:
        species_network = []
        ordered_genes = cluster_species(network, spec_num, code_uncode, translate_bool)
        for i in ordered_genes:
            species_network.append(', '.join(i)) 
        final_networks.append(species_network)       
    with open(outfile, "w+") as handle:
        for network in final_networks:
            #translate edge names into genes 
            print_string = '\t'.join(network)
            handle.write(print_string + '\n')

def translate_genes(gene, code_uncode):
    return code_uncode[gene]

def load_coded_gene_names(sequences_file):
    code_gene = {}
    gene_code = {}
    with open(sequences_file, "r") as handle:
        for line in handle:
            line = line.strip().split(": ")
            gene_code[line[1]] = line[0]        
            code_gene[line[0]] = line[1]
    return code_gene, gene_code

def singleton_shortcut(infile):
    final_singletons = []
    with open(infile, "r") as handle:
        for line in handle:
            line = line.strip().split('\t')
            final_singletons.append(line)
    return final_singletons

def add_names(edge, gene_list):
    edge.append(gene_list.index(edge[0]))
    edge.append(gene_list.index(edge[1]))
    return edge

def translate(gff, prefix, seqIDs):
    code_gene, gene_code = load_coded_gene_names(seqIDs)
    output = []
    with open(gff, "r") as handle:
        for line in handle:
            line = line.strip().split("\t")
            line[1] = gene_code[line[1]]
            output.append(line)
    output = sorted(output, key = lambda x: (x[0], int(float(x[2]))))
    with open(prefix + ".gff", "w") as handle:
        for l in output: 
            j = "\t".join(l)
            handle.write('%s\n' % j)

def load_species(specID):
    specID_dict = {}
    with open(specID, "r") as handle:
        for line in handle:
            line = line.strip()
            if line[0] == "#":
                line = line[1:]
            line = line.split(": ")
            spec = line[1].split(".")
            line[1] = '.'.join(spec[:-1])
            specID_dict[line[1]] = line[0]
    return specID_dict

def main(prefix, orthogroups, threads, ploidies, sequenceIDs, speciesIDs):
    print("Starting pSONIC")
    start_time = time.time()

    code_gene, gene_code = load_coded_gene_names(sequenceIDs)
    good_groups = [["len(block)", "a_start", "a_end", "a_range", "a_density", "b_start", "b_end", "b_range", "b_density", "a_startgene", "b_startgene", "a_endgene", "b_endgene", "density_product", "positive_pairs", "negative_pairs", "good_bad"]]

    with open(prefix + ".gff", "r") as handle:
        gffs = [i.strip().split('\t') for i in handle]
    gene_names = [i[1] for i in gffs]

    with open(prefix + ".tandem", "r") as handle:
        tandem_genes = [line.strip().split(',') for line in handle]

    tandem_network = edges_to_groups(tandem_genes, gene_names, threads)
    tandem_network.sort(key=len, reverse=True)

    if ploidies:
        with open(ploidies, "r") as handle:
           ploidies = handle.readline().strip().split('\t') 

    tandem_pairs = []
    for i in tandem_network:
        tandem_pairs += [j for j in i]

    species_dict = load_species(speciesIDs)

    species_nums, singletons = SingleCopy_from_Orthofinder(orthogroups, tandem_network, tandem_pairs, ploidies, gene_code, species_dict)
    print("Tether Sets from OrthoFinder All Found: " + str((time.time() - start_time)/60))

    singletons = edges_to_groups(singletons, gene_names, threads)
    singletons.sort(key=len, reverse=True)

    print_network("pSONIC.TetherSetsFromOrthoFinder.csv", singletons, species_nums, code_gene, False)

    gff_genes = []
    for i in gene_names:
        if i not in tandem_pairs:
            gff_genes.append([i])
        else: gff_genes[-1].append(i)

    singleton_set = []
    for i in singletons:
        singleton_set += [j for j in i]        

    print("Reading Collinearity File") 
    raw_groups = read_collinearity(prefix + ".collinearity")
    trimmed_edges, raw_groups, edges, good_groups = score_all_blocks(raw_groups, singletons, good_groups, singleton_set, gff_genes, species_nums, gene_names, threads)

    print("Filtering of Collinear Groups Completed: " + str((time.time() - start_time)/60))
    print_list("pSONIC.GroupsKept.txt", good_groups)
    print_list("pSONIC.EdgeList.txt", [i[:-2] for i in edges]) #don't need 
    print_list("pSONIC.RawGroups.txt", raw_groups) #don't need species headers
    print_list("pSONIC.EdgesToTrim.txt", [i[:-2] for i in trimmed_edges]) #need 

    print("Assembling Collinear Blocks into Syntenic Orthologs.")
    edges += tandem_genes
    output = edges_to_groups(edges, gene_names, threads)

    output.sort(key=len, reverse=True)

    tandem_size, group_size, raw_out, trans_out = print_orthogroups_quick(output, species_nums, code_gene, tandem_network, tandem_pairs, threads)

    print_list("pSONIC.txt", trans_out)
    print_list("pSONIC.Untranslated.txt", raw_out)
    print_list("pSONIC.TandemCheck.csv", tandem_size) ##need
    print_list("pSONIC.GroupSize.csv", group_size) #need 
    print("pSONIC Completed.")
    print("Total Time to run: " + str((time.time() - start_time)/60))


def parse_args():
    parser = argparse.ArgumentParser(usage="pSONIC prefix [translate_gff] [-h] [-og ORTHOGROUPS] [-t THREADS] [-p PLOIDY] [-sID SEQUENCEIDS] [-specID SPECIESIDS] [-gff GFF]")
    parser.add_argument("prefix", help="PREFIX used to run MCScanX. If used with [translate_gff] argument, resulting file will be [PREFIX].gff. Otherwise, pSONIC expects files with the names [PREFIX].collinearity, [PREFIX].gff, and [PREFIX].tandem to be in current directory.")
    parser.add_argument("translate_gff", help="Only translate gff file to fit gene IDs used by OrthoFinder. Requires [-gff] argument.", nargs='?')
    parser.add_argument("-og", "--orthogroups", help="Orthogroups output file from OrthoFinder. File can either have .tsv extension or .csv (for older versions of OrthoFinder). (Default: %(default)s)", default = "Orthogroups.tsv")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads to use. (Default: %(default)s)")
    parser.add_argument("-p", "--ploidy", help="Tab-delimited file of relative ploidies for each species, listed in same order as in ORTHOGROUPS file. (Default: All species are diploid with no WGDs in the tree)")
    parser.add_argument("-sID", "--sequenceIDs", help="SequenceIDs.txt file from OrthoFinder. (Default: SequenceIDs.txt)", default="SequenceIDs.txt")
    parser.add_argument("-gff", help="GFF file to translate before running MCScanX. Only use with [translate_gff] option.")
    parser.add_argument("-specID", "--speciesIDs", help="SpeciesIDs.txt file from OrthoFinder. Hashes in first character can be present.", default="SpeciesIDs.txt")
    args = parser.parse_args()

    if args.translate_gff:
        if args.gff is None:
            parser.error("To use translate_gff function, you must specify -gff option.")
        elif args.gff == args.prefix + ".gff":
            parser.error("The gff file you want translated has the same file prefix as the MCScanX prefix. Please change the name of this file before proceeding to prevent this file from being overwritten.")
        else: translate(args.gff, args.prefix, args.sequenceIDs)
    else:
        main(args.prefix, args.orthogroups, args.threads, args.ploidy, args.sequenceIDs, args.speciesIDs)

parse_args()
