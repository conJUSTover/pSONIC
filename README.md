# pSONIC: Ploidy-aware Syntenic Orthologous Networks Identified via Collinearity
The repository serves as a public and official hosting of the pSONIC program (Conover et al., in review). If you have any questions regarding the implementation or running of pSONIC, please submit an issue on GitHub. 

---
## pSONIC is a tool for identifying Orthologs between multiple species, even when polyploidy is involved.    
One of the first steps in comparative genomics studies is to create a comprehenive list of orthologs between all species. Several methods have been developed to do this, but they are designed for and tested on distantly-related species across Orders or even Phyla. With the explosion in number of high-quality genome assemblies from closely related species, these methods designed for deep-phylogenetic identification will become less-than-optimal, as conserved gene order (synteny) can be used as powerful evidence of orthology. 

Additionally, one of the most important findings during the genomics era is the ubiquity of whole genome duplications (WGD) across the tree of life, especially in plants. All plant species have a history of WGD in their evolutionary histories, and inference of orthology is complicated among plants with different histories of WGD as it disrupts the 1:1 expectation of orthologs between two species. 

Here, we use synteny across multiple species, as well as _a priori_ knowledge of ploidy differences between species, to infer a genome-wide set of syntenic orthologs. We combine the outputs of [MCScanX](https://github.com/wyp1125/MCScanX) and [OrthoFinder](https://github.com/davidemms/OrthoFinder) to create a genome-wide set of syntenic orthologs. Details on how pSONIC operates can be found in the manuscript or [preprint](). 

### Getting Started with pSONIC 
pSONIC is a python program that is written and tested using Python v3.7. The only non-standard python package required is the iGraph package (see download and installation instructions [here](https://igraph.org/python/). pSONIC was tested using igraph v0.8.3.  

#### Options for running pSONIC 

```
python3 ../../pSONIC.py -h 
usage: pSONIC prefix [translate_gff] [-h] [-og ORTHOGROUPS] [-t THREADS] [-p PLOIDY] [-sID SEQUENCEIDS] [-gff GFF]

positional arguments:
  prefix                PREFIX used to run MCScanX. If used with
                        [translate_gff] argument, resulting file will be
                        [PREFIX].gff. Otherwise, pSONIC expects files with the
                        names [PREFIX].collinearity, [PREFIX].gff, and
                        [PREFIX].tandem to be in current directory.
  translate_gff         Only translate gff file to fit gene IDs used by
                        OrthoFinder. Requires [-gff] argument.

optional arguments:
  -h, --help            show this help message and exit
  -og ORTHOGROUPS, --orthogroups ORTHOGROUPS
                        Orthogroups output file from OrthoFinder. (Default:
                        Orthogroups.csv)
  -t THREADS, --threads THREADS
                        Number of threads to use. (Default: 1)
  -p PLOIDY, --ploidy PLOIDY
                        Tab-delimited file of relative ploidies for each
                        species, listed in same order as in ORTHOGROUPS file.
                        (Default: All species are diploid with no WGDs in the
                        tree)
  -sID SEQUENCEIDS, --sequenceIDs SEQUENCEIDS
                        SequenceIDs.txt file from OrthoFinder. (Default:
                        SequenceIDs.txt)
  -gff GFF              GFF file to translate before running MCScanX. Only use
                        with [translate_gff] option.

```


### Running pSONIC for the first time 
---

__1. Install and run OrthoFinder on your species of interest__    
You should run OrthoFinder using the '-og' flag (this stops OrthoFinder after the stage of inferring orthogroups. It has many more capabilities than what is required for pSONIC, and you should definitely check them out! It's a great program). We will need three files from the OrthoFinder output: 

```
Blast*.txt
Orthogroups.csv
SequenceIDs.txt
```

You should move these files into the same directory, and concatentate all of the Blast files into a single file named <PREFIX>.blast (where <PREFIX> is any prefix, but must be the same <PREFIX> used throughout). 

------
__2. Prepare files for MCScanX.__     
The gff file should have four tab-separated columns:  

`Sp##	Start_POS End_POS 	GeneID`    

__"Sp##"__ is a two-letter code for each species, followed by the chromosome number (e.g. At01)    
__"Start\_POS"__ and __"End\_POS"__ are the beginning of end of the gene sequence on it's chromosome (i.e. columns 4 and 5 in .gff3 formatted files)     
__"GeneID"__ is a unique name for each gene in the fasta files given to OrthoFinder.     
**_You should ensure that every gene in the .gff file is represented in the .fasta files, and every gene in the .fasta files is represented in the .gff file._      

This gff file should be named anything other than "<PREFIX>.gff". For demostrative purposes, we will name this file "your\_gff\_file.gff". We first must translate the gene names into the coded names used by OrthoFinder. We can do this using pSONIC. 

To run pSONIC for this step, use the following command:    
`python pSONIC.py translate_gff <PREFIX> -gff your_gff_file.gff`    

This will create a file called <PREFIX>.gff, which is identical to your gff files except that all of the gene names have been replaced by the coded gene names listed in "SequenceIDs.txt". 

------
__3. Run MCScanX.__ MCScanX can be downloaded from [here](https://github.com/wyp1125/MCScanX). You should now two files in your current directory named <PREFIX>.gff and <PREFIX>.blast. To run MCScanX, be sure that the MCScanX is executable in your $PATH, then run     
```
MCScanX -b 2 <PREFIX>
```
This runs MCScanX on only interspecific comparisons (-b 2), which will greatly improve the speed of pSONIC. (__NOTE__: If there is a polyploidy event anywhere between any of your species in the anaysis, then run MCScanX without the "-b 2" flag to include intra-specific comparisons)


------     
__4. Run pSONIC__     
Now that we have all our needed files in the same directory, let's run pSONIC as such: 
```
python3 pSONIC.py <PREFIX>
```

You can use the [optional flags](#Options-for-running-pSONIC) to change the settings and input files needed for pSONIC. If there is a polyploidy event anywhere betwen any of your species in the analysis, be sure to provide the `-p` flag.     


## Output Files of pSONIC   
There are several output files generated by pSONIC. The most important file is __pSONIC.txt__ which contains a list of orthogroups. Each tab-separated column of pSONIC.txt contains a list of genes from each species, and the genes within each species are in a comma-separated list. The structure of this output is identical to that of OrthoFinder to ease comparisons of results, and to add utility in pipelines when switching between OrthoFinder and pSONIC results.

The full list of output files of pSONIC include:   
__pSONIC.txt__:  Main output, described above.  
__pSONIC.untranslated.txt__:  Same file as pSONIC.txt, but with the coded gene names used by OrthoFinder.  
__pSONIC.GroupSize.csv__:   A tab-delimimted file describing how many genes from each species is included in each orthogroup.  
__pSONIC.RawGroups.txt__:  A file with six columns, each row describes a different collinear groups identified by MCScanX. The six columns include total collinearity group length, number of gene pairs with "PASS" scores, number of gene pairs with "NOT PASS" scores, number of "No Call" scores, the pattern of every gene pair's calls along the collinear group, and the orientation of the collinear group.    
__pSONIC.GroupsKept.txt__: Describes the collinear blocks that were used for orthogroup construction, after trimming and cutting (see manuscript for details on how these groups are trimmed and cut).   
__pSONIC.TetherSetsFromOrthoFinder.csv__:  A list of tether sets used to scores gene pairs. These tether sets were identified from OrthoFinder output. See the manuscript for more information about how these tether sets are identified.   
__pSONIC.EdgeList.txt__:  A list of all of the gene pairs used to create the orthogroups. Each line represents the two genes in the pair, their BLAST score e-value, the orientation of the collinear block they are present on, the score received by pSONIC, the two chromosomes from which the two genes originate, and an index of which collinear group the gene pair is from (i.e. "##Alignment number" in <PREFIX>.collinearity output from MCScanX.). Collinear groups ending in ".01", ".02", etc. are split from the original collinear group (groups are split when at least three consecutive gene pairs receive "NOT PASS" scores without a "PASS" in between them.)    
__pSONIC.EdgesToTrim.txt__: A list of all the gene pairs that were trimmed from the ends of all collinear blocks.   
__pSONIC.TandemCheck.csv__: A tab-delimited file describing how many genes (or tandemly duplicated gene sets) are present in each Orthogroup. Ideally, the value for each species should be 1 (or equal to the relative ploidy specified in the file provide using the '-p' flag).    
 



