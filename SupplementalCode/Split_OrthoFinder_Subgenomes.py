#!/usr/bin/py
#This code split OrthoFinder output if one species is a polyploid. 
#It is split into its respective subgenomes. 
#usage: .py <File> <poly_index> <subA_Prefix> <subB_prefix>

import sys 

OrthoFile = sys.argv[1]
poly_index = int(sys.argv[2])
subA_prefix = sys.argv[3]
subB_prefix = sys.argv[4]

with open(OrthoFile) as handle:
#    firstline = handle.readline()
#    firstline = firstline.rstrip().split('\t')
#    poly_string = firstline.pop(poly_index)
#    firstline.append(poly_string)
#    firstline.append(poly_string)
#    print(*firstline, sep="\t")
    for line in handle:
        line = line.rstrip('\n').split('\t')
        poly = line.pop(poly_index)
        if ',' in poly:
            poly = poly.split(',')
        else: poly = [poly]
        subA = [gene for gene in poly if subA_prefix in gene]
        subB = [gene for gene in poly if subB_prefix in gene]
        subA = ','.join(subA)
        subB = ','.join(subB)
        line.append(subA)
        line.append(subB)
        print(*line, sep="\t")
