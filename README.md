[![Test](https://github.com/milnus/Corekaburra/actions/workflows/Test.yml/badge.svg)](https://github.com/milnus/Corekaburra/actions/workflows/Test.yml)

# Overview 
Corekaburra looks at the gene synteny across genomes used to build a pan-genome. Using syntenic information Corekaburra 
identifies regions between core gene clusters. Regions are described in terms of their content of accessory gene clusters 
and distance between core genes. Information from neighboring core genes is further used to identify stretches of core 
gene clusters throughout the pan-genome that appear in all genomes given as input. Corekaburra is compatible with outputs 
from standard pan-genome pipelines: [Roary](academic.oup.com/bioinformatics/article/31/22/3691/240757) and [Panaroo](genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02090-4).

# Why and When to use Corekaburra
Corekaburra fits into the existing frameworks of bioinformatics pipelines for pan-genomes. It does not reinvent a new pan-genome pipeline, but leverages the existing ones. Because of this, Corekaburra is build to be a natural extension to the analysis of pan-genomes by summarising information and inferring relationships in the pan-genome otherwise not easily accessible via pan-genome graphs. Other tools provide similar outputs or information, but in their own standalone pan-genome analysis framework or pipeline. Such frameworks/pipelines are [PPanGGolin](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732) and [Panakeia](https://www.biorxiv.org/content/biorxiv/early/2021/03/02/2021.03.02.433540.full.pdf). By building on top of existing tools Corekaburra frees users from potentially cross referencing beteween pan-genomes, which in itself is a challenging task. Corekaburra's workflow also allows it to be extended to any pan-genome tool, with an output similar to the gene_presence_absence.csv produced by Roary, making Corekaburra versatile for future implementations.

# Installation
Corekaburra is writen in Python 3.9, and can be installed via pip and conda. A Docker container is also available.
## pip
```pip install corekaburra```

## building a Conda environment from scratch
```conda create -n Corekaburra python==3.9```
```conda activate Corekaburra```
```pip install corekaburra```

## Conda install
```Comming```

## Docker
See the (Wiki for more information)[https://github.com/milnus/Corekaburra/wiki/Docker.md]

# Help
```
usage: Corekabura -ig file.gff [file.gff ...] -ip path/to/pan_genome [-cg complete_genomes.txt] [-a] [-cc 1.0] [-lc 0.05] [-o path/to/output] [-p OUTPUT_PREFIX] [-d] [-c int] [-l | -q] [-h]

Welcome to Corekaburra! An extension to pan-genome analyses that summarise genomic regions between core genes and segments of neighbouring core genes using gene synteny from a set of input genomes and a pan-genome folder.

Required arguments:
  -ig file.gff [file.gff ...], --input_gffs file.gff [file.gff ...]
                        Path to gff files used for pan-genome
  -ip path/to/pan_genome, --input_pangenome path/to/pan_genome
                        Path to the folder produced by Panaroo or Roary

Analysis modifiers:
  -cg complete_genomes.txt, --complete_genomes complete_genomes.txt
                        text file containing names of genomes that are to be handled as complete genomes
  -a, --no_annotate_refound
                        Flag to toggle off the creation of new gff files, with annotation of refound genes. Only done if input pangenome is detected as coming from Panaroo
  -cc 1.0, --core_cutoff 1.0
                        Percentage of isolates in which a core gene must be present [default: 1.0]
  -lc 0.05, --low_cutoff 0.05
                        Percentage of isolates where genes found in less than these are seen as low-frequency genes [default: 0.05]

Output control:
  -o path/to/output, --output path/to/output
                        Path to where output files will be placed [default: current folder]
  -p OUTPUT_PREFIX, --prefix OUTPUT_PREFIX
                        Prefix for output files, if any is desired
  -d, --discard_corrected
                        Discard gff files corrected with refound genes identified by Panaroo - Only compativle if pan-genome comes from Panaroo [Default: Corrected files are kept]

Other arguments:
  -c int, --cpu int     Give max number of CPUs [default: 1]
  -l, --log             Record program progress in for debugging purpose
  -q, --quiet           Only print warnings
  -h, --help            Show help function
```

# Inputs
## Gff files
Input Gff files must be included in the pan-genome gene_presence_absence.csv-style file.  
The Gffs are also required to contain a ```##FASTA``` dividing the file into annotations at the top and the Fasta genome in the bottom of the file.  
All coding sequences (CDS) annotated in the GFF must also carry an ```ID``` and a ```locus_tag```.  
*** Input Gffs can be gzipped ***

## Pan-genome folder
This is the output folder from a Roary or Panaroo run, or a folder that at minimum contains the gene_presence_absence.csv from Roary or the gene_presence_absence_roary.csv from Panaroo.

## Complete genomes
If some input Gff are to processed as complete or closed genomes, a plain text file can be provided with the filename of these.  
example:
```
complete_genome.gff
complete_genome.gff.gz
/paths/are/allowed/complete_genome.gff
complete_genome
```
All files given in the plain text file of complete genomes must be found in a given gene presence/absence file, but are not required to be among the input gffs, meaning that a single plain text file of complete genomes can be used for analysing subsets of genomes in the pan-genome.

## Adjusting cutoffs
To comply with common practice when handling pan-genomes, the cutoff for when a pan-genome cluster (gene) is perceived as core can be changed using the ```-cc``` arguments with a ratio of gene presence required. By default, this is set to a conservative 100% presence of core gene clusters.  
A second argument dividing accessory genes into two groups (Low frequency and Intermediate frequency) can be controlled using the ```-lc``` argument, with the ratio indicating the maximum presence of a gene cluster to be identified as having a low frequency in the pan-genome. This division of low- and intermediate frequency can be disabled by ```-lc 0```, resulting in all genes being considered as intermediate.

# Outputs
Corekaburra outputs multiple files ranging from summaries to more fine grained outputs. This is aimed at giving the user easy access to information, but still allowing for tailored or deep exploration. 

## Core regions
A Core region is defined by two core gene clusters flanking a stretch of the genome in at least one input genome (Gff). A core region can be described by a distance between the flanking core gene clusters, positive if nucleotides can be found between then, and negative if the two clusters overlap). A region can also be described by the number of encoded accessory genes, low- and intermediate frequency. Using core gene clusters as a reference for a region it is possible to compare the same region across genomes, and in the larger framework of the pan-genome. Additionally, with either or both the distance and number of encoded accessory genes in a region it is possible to identify regions of variability, due to horizontal genetic transfer, deletion or other genomic processes.

```core_pair_summary.csv``` is a file that summarises the identified core regions identified across the input genomes (Gffs). Here information about occurrence and co-occurnece of each core gene pair, and individual core gene occurrences can be found. as well as, distance and accessory gene summary statistics (minimum, maximum, mean, and median).  
This file is a good entery point to the results in most analysis, and should give a good indication of which core regions that could be of interest.

```core_core_accessory_gene_content.tsv``` gives the placement of each accessory genes identified in a core region across all genomes (Gff). It is also given if the accessory gene is identified as a low- or intermediate frequency gene.

```low_frequency_gene_placement.tsv``` summarises each core region across all genomes (Gff) with the distance between core gene clusters, and the number of accessory genes found between them.

## Core segments
The two following files are only given if any core gene is found to have more than two different core gene clusters as neighbours across all input genomes (Gff).

The file ```core_segments.csv``` containg all segments of minimum two core genes identified in a pan-genome, where the start and end of a segments is defined by core gene clusters with more than two neighbours, meaning they could be a potential breakpoint of a genomic rearrangements in at least a single input genome (Gff), or be a misassembly.    

```no_accessory_core_segments.csv``` divides the segments identified in ```core_segments.csv``` into potential smaller segments where core gene clusters must form regions with no accessory genes between them across all genomes. These segments could indicate potential operon structures or other stable genomic feature, that could be disturbed by insertion of accessory genes.

## Core-less contigs
```core_segments.csv``` gives all contigs identified in genomes (Gff) that does not contain a core gene cluster, but only accessory genes. Each contig is given by contig name, its Gff file, and number of low- and intermediate frequency genes found on the contig. 

## Corrected Gffs
A folder containing Gff files that have been corrected by annotating the genes refound by Panaroo. This folder is only expected when a pan-genome from Panaroo is provided, and the ```-a``` or ```-d``` arguments are not given as inputs.  
**Notice this will duplicate your Gff files, meaning that ```-a``` or ```-d``` arguments should be used to avoid this, when dealing with memory issues or large datasets**

# For more info
For more into on Corekaburra, its workings, inputs, outputs and more see the (wiki)[https://github.com/milnus/Corekaburra/wiki]


# Bug reporting and feature requests
Please submit bug reports and feature requests to the issue tracker on GitHub: [Corekaburra issue tracker](https://github.com/milnus/Corekaburra/issues)

# Licence
This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/milnus/Corekaburra/master/LICENSE).
