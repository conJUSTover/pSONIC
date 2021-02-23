#!/usr/bin/python

#usage: .py <out_directory> pSONIC.TetherSetsFromOrthoFinder.csv pSONIC.EdgeList.txt pSONIC.EdgesToTrim.txt pSONIC.Untranslated.txt <list_all_html_files>*.html

import sys
import re


def read_Singletons(singletons_file):
    with open(singletons_file, "r") as handle:
        output = [i.strip().replace(', ', '\t').split("\t") for i in handle]
    return output


def parse_html(out_directory, singleton_list, edges, trimmed_edges, Shit_out, html_file):
#    sys.stdout = open(out_directory + "/" + html_file, "w+")
    p = re.compile('bgcolor=(.*?)</td>')
    final_out = []
    with open(html_file, "r") as handle:
        for line in handle:
            genes = p.findall(line)
            if genes:
                ref_gene = genes.pop(0)
                singleton_group = find_group(ref_gene[10:], singleton_list)
                shit_group = find_group(ref_gene[10:], Shit_out)
                if genes and singleton_group: #ensuring there was originally 2 items in genes and reference gene in singletons_list
                    for i in genes: #for all non-reference genes
                        if i[10:] in singleton_group: #if the gene is same singleton group as reference gene
                            j = i.replace("#ffff99", "#4575b4") #darker blue
                        else: j = i.replace("#ffff99", "#fc8d59") #orange
                        i += "<"
                        j += "<"
                        line = line.replace(i, j)
                else:
                    for i in genes:
#                        if ",".join([ref_gene[10:], i[10:]]) in trimmed_edges:
#                            j = i.replace("#ffff99", "#77c7ef") #light blue
#                            line = line.replace(i, j)
#                        elif ",".join([i[10:], ref_gene[10:]]) in trimmed_edges:
#                            j = i.replace("#ffff99", "#77c7ef") #light blue
#                            line = line.replace(i, j)

                        #need to find all singletons, not just those in edges


                        if ",".join([ref_gene[10:], i[10:]]) in edges:
                            j = i.replace("#ffff99", "#91bfdb") #light blue
                            line = line.replace(i, j)
                        elif ",".join([i[10:], ref_gene[10:]]) in edges:
                            j = i.replace("#ffff99", "#91bfdb") #light blue
                            i += "<"
                            j += "<"
                            line = line.replace(i, j)
                    if shit_group:
                        if i[10:] in shit_group: #if the gene is same singleton group as reference gene
                            j = i.replace("#ffff99", "#ADD8E6") #light blue
                        #else: j = i.replace("#ffff99", "#786cd1") #purple
                            i += "<"
                            j += "<"
                            line = line.replace(i, j)
            final_out.append(line.strip())
    with open(out_directory + "/" + html_file, "w+") as outfile:
        for print_string in final_out:
            outfile.write(print_string + '\n')        
            #print(line)

def parse_edge_list(edges):
    edge_list = []
    with open(edges, "r") as handle:
        for line in handle:
            line = line.strip().split()
            edge_list.append(",".join(line[:2]))
    return edge_list

def find_group(a, lists):
    for i in lists:
        if a in i:
            return i
    return []

def main():
    trimmed_edges = parse_edge_list(sys.argv[4])
    singletons = read_Singletons(sys.argv[2])
    final_output = read_Singletons(sys.argv[5])
    edges = parse_edge_list(sys.argv[3])
    for i in sys.argv[6:]:
        parse_html(sys.argv[1], singletons, edges, trimmed_edges, final_output, i) 
main()

