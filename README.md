[![Test](https://github.com/milnus/Corekaburra/actions/workflows/Test.yml/badge.svg)](https://github.com/milnus/Corekaburra/actions/workflows/Test.yml)
[![codecov](https://codecov.io/gh/milnus/Corekaburra/branch/main/graph/badge.svg?token=090xERhDET)](https://codecov.io/gh/milnus/Corekaburra)

# Overview 
Corekaburra looks at the gene synteny across genomes used to build a pan-genome. Using syntenic information Corekaburra identifies regions between core genes. Regions are described in terms of their content of accessory genes and number of nucleotides between core genes. Information from neighboring core genes is further used to identify stretches of core gene clusters that appear in all genomes given as input. Corekaburra is compatible with outputs from 'standard' pan-genome pipelines: [Roary](academic.oup.com/bioinformatics/article/31/22/3691/240757) and [Panaroo](genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02090-4), and can be extended to others if desired.

# Why and When to use Corekaburra
Corekaburra fits into the existing frameworks of bioinformatics pipelines for pan-genomes as a downstream analysis tool. It does not reinvent a new pan-genome pipeline, but leverages the existing ones. Because of this, Corekaburra is built to be a natural extension to the analysis of pan-genomes by summarising information and inferring relationships in the pan-genome otherwise not easily accessible via pan-genome graphs. Other tools provide similar outputs or information, but in their own stand-alone pan-genome analysis framework or pipeline. Examples of such frameworks/pipelines are [PPanGGolin](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732) and [Panakeia](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08303-3). By building on top of existing tools Corekaburra frees users from potentially cross referencing between pan-genomes, which in itself is a challenging task. Corekaburra's workflow also allows it to be extended to any pan-genome tool, with an output similar to the gene_presence_absence.csv produced by Roary, making Corekaburra versatile for future implementations of pan-genome pipelines.

# Installation
Corekaburra is written in Python 3.9, and can be installed via pip and conda. A Docker container is also available.

## Conda install
```conda install -c bioconda -c conda-forge corekaburra```

## pip
```pip install corekaburra```

## Docker
Pull from [DockerHub](https://hub.docker.com/repository/docker/magnusgj/corekaburra) or see the [Wiki for more information](https://github.com/milnus/Corekaburra/wiki/Docker.md)

# Help
```
usage: Corekaburra -ig file.gff [file.gff ...] -ip path/to/pan_genome [-cg complete_genomes.txt] [-cc 1.0] [-lc 0.05] [-o path/to/output] [-p OUTPUT_PREFIX] [-c int] [-l | -q] [-h] [-v]

Welcome to Corekaburra! An extension to pan-genome analyses that summarise genomic regions between core genes and segments of neighbouring core genes using gene synteny from a set of input genomes and a pan-
genome folder.

Required arguments:
  -ig file.gff [file.gff ...], --input_gffs file.gff [file.gff ...]
                        Path to gff files used for pan-genome
  -ip path/to/pan_genome, --input_pangenome path/to/pan_genome
                        Path to the folder produced by Panaroo or Roary

Analysis modifiers:
  -cg complete_genomes.txt, --complete_genomes complete_genomes.txt
                        text file containing names of genomes that are to be handled as complete genomes
  -cc 1.0, --core_cutoff 1.0
                        Percentage of isolates in which a core gene must be present [default: 1.0]
  -lc 0.05, --low_cutoff 0.05
                        Percentage of isolates where genes found in less than these are seen as low-frequency genes [default: 0.05]

Output control:
  -o path/to/output, --output path/to/output
                        Path to where output files will be placed [default: current folder]
  -p OUTPUT_PREFIX, --prefix OUTPUT_PREFIX
                        Prefix for output files, if any is desired

Other arguments:
  -c int, --cpu int     Give max number of CPUs [default: 1]
  -l, --log             Record program progress in for debugging purpose
  -q, --quiet           Only print warnings
  -h, --help            Show help function
  -v, --version         show program's version number and exit
```

# Example workflow
We have made an [example workflow](https://github.com/milnus/Corekaburra/wiki/Example-workflow) using three *Streptococcus pyogenes* genomes, Panaroo and Corekaburra.

# Inputs
## Gff files
Input Gff files must be included in the pan-genome gene_presence_absence.csv-style file.  
The Gffs are also required to contain a ```##FASTA``` line, dividing the file into annotations at the top and the genome in the bottom of the file.  
All coding sequences (CDS) annotated in the GFF must also carry an ```ID``` and a ```locus_tag```.  
Input Gff files can be in gzipped format, if desired.

## Pan-genome folder
This is the output folder from a Roary or Panaroo run. If designing a gene_presence_absence.csv from another pan-genome tool folder must at minimum contain a gene_presence_absence.csv file with each field quoted as the gene_presence_absence.csv from Roary.

## Complete genomes
If some input Gffs are to be processed as complete or closed genomes, a plain text file can be provided with the filename of these.  
example:
```
complete_genome.gff
complete_genome.gff.gz
/paths/are/allowed/complete_genome.gff
complete_genome
```
All files given in the plain text file as complete genomes must be found in a given gene presence/absence file, but are not required to be among the input gffs, meaning that a single plain text file of complete genomes can be used for analysing subsets of genomes in the pan-genome.

## Adjusting cutoffs
To comply with common practice when handling pan-genomes, the cutoff for when a pan-genome cluster (gene) is perceived as core can be changed using the ```-cc``` arguments with a ratio of gene presence required. By default, this is set to a conservative 100% presence of core genes.  
A second argument dividing accessory genes into two groups (Low frequency and Intermediate frequency) can be controlled using the ```-lc``` argument, with the ratio indicating the maximum presence of a gene cluster to be identified as having a low frequency in the pan-genome. This division of low- and intermediate frequency can be disabled by ```-lc 0```, resulting in all genes being considered as intermediate.

# Outputs
Corekaburra outputs multiple files ranging from summaries to more fine grained outputs. This is aimed at giving the user easy access to information, but still allowing for tailored or deep exploration. See [description of outputs in wiki](https://github.com/milnus/Corekaburra/wiki/Inputs-and-outputs#outputs) and [how to query outputs](https://github.com/milnus/Corekaburra/wiki/Down-stream-analyses#Querying-outputs) 

## Core regions
A Core region is defined by two core genes flanking a stretch of the genome in at least one input genome. A core region can be described by a distance between the flanking core genes, positive if nucleotides can be found between them, and negative if the two genes overlap). A region can also be described by the number of encoded accessory genes. Using core gene clusters as a reference for a region of genomes it is possible to compare the same region or their presence across genomes. Additionally, with either or both the distance and number of encoded accessory genes in a region it is possible to identify regions of variability, due to horizontal genetic transfer, deletion or other genomic processes.

```core_pair_summary.csv``` is a file that summarises the core regions identified across the input genomes (Gff files). Here information about occurrence and co-occurence of each core gene pair, and individual core gene occurrences can be found. Distance and accessory gene summary statistics (minimum, maximum, mean, and median) for each core pair is summarised.  
This file is a good entry point to the results in most analyses, and should give a good indication of which core regions that could be of interest.

```core_core_accessory_gene_content.tsv``` gives the placement of each accessory gene found in a core region across all genomes (Gff). Frequency of accessory genes (low- or intermediate frequency) is given.

```low_frequency_gene_placement.tsv``` summarises each core region across all genomes (Gff) with the distance between core gene clusters, and the number of accessory genes found in the region.

## Core segments
The two following files are only given if any core gene is found to have more than two different core genes as neighbours across all input genomes (Gff), meaning there is structural heterogeneity across the genomes.

The file ```core_segments.csv``` contain all segments of minimum two core genes identified in a pan-genome, where the start and end of a segments is defined by core gene clusters with more than two neighbours, meaning they could be a potential breakpoint of a genomic inversion in at least a single input genome (Gff), or be a misassembly.    

```no_accessory_core_segments.csv``` divides the segments identified in ```core_segments.csv``` into potential smaller segments where core genes must form regions with no accessory genes between them across all genomes. These segments could indicate potential operon structures or other stable genomic features that could be disturbed by insertion of accessory genes.

## Core-less contigs
```coreless_contig_accessory_gene_content.tsv``` gives all contigs identified in genomes (Gff) that do not contain a core gene cluster, but only accessory genes. Each contig is given by contig name, its Gff file, and number of low- and intermediate frequency genes found on the contig. 

# For more info
For more into on Corekaburra, its workings, inputs, outputs and more see the (wiki)[https://github.com/milnus/Corekaburra/wiki]


# Bug reporting and feature requests
Please submit bug reports and feature requests to the issue tracker on GitHub: [Corekaburra issue tracker](https://github.com/milnus/Corekaburra/issues)

# Licence
This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/milnus/Corekaburra/master/LICENSE).
