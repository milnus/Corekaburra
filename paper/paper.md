---
title: 'Corekaburra: pan-genome post-processing using core gene synteny'
tags:  
  - Python
  - microbiology
  - genomics
  - genome
  - pan-genome
authors:
  - name: Magnus G. Jespersen 
    orcid: 0000-0001-9751-9877
    affiliation: 1
  - name: Andrew Hayes
    orcid: 0000-0001-8038-1656
    affiliation: 1
  - name: Mark R. Davies
    orcid: 0000-0001-6141-5179
    affiliation: 1
affiliations:
  - name: Department of Microbiology and Immunology, University of Melbourne at the Peter Doherty Institute for Infection and Immunity, Melbourne, VIC, Australia
    index: 1
date: 02 September 2022  
bibliography: paper.bib
---


<!---References like: [@altschul:1990; @mount:2007]--->

# Summary
Pan-genome analysis enables an assessment of the total gene content from a set of genome sequences. Since the 'first' defined pan-genome [@tettelin:2005] many tools have been developed which improve the methodological process used to construct pan-genomes [@page:2014; @tonkin-hill:2020; @thorpe:2018; @inman:2019; @gautreau:2020]. However, current limitations lay in extracting meaningful interpretations from downstream analyses of pan-genomes. Few tools have been made to address this problem, and often focus on a specific analysis, such as pan-genome association studies [@brynildsrud:2016; @whelan:2020; @lees:2018]. Here we present Corekaburra, a tool to define gene regions based on core gene synteny within a pan-genome. Defining regions by flanking core genes provide context to genes and allow for easier systematic comparisons of genomic features, such as genomic inversions and gene insertions and deletions.

# Statement of need
Bacterial genomes from the same species can vary considerably in their genetic content [@welch:2002]. In population genomics and other studies of bacterial genomes an important piece of information is shared genetic information. Due to this, pan-genomes have become a standard method for investigating variation in genetic content of bacteria [@Medini:2020]. The analysis of pan-genomes is critical to basic research, industrial strain development, and public health surveillance systems. Dispite this, methods to systematically dissect pan-genomes are sparse.  

We propose Corekaburra a program designed to reduce the complexity of outputs from pan-genome pipelines, easing the discovery of regions of variation within a pan-genome. The input for Corekaburra is annotated genomes (Gff3 format with appended genome), similar to those used by existing pan-genome pipelines, and the output folder from a pan-genome tool (currently Roary or Panaroo)[@page:2014; @tonkin-hill:2020]. Corekaburra introduces core gene synteny by scanning and summarizing input Gff files based on the genetic distance of nucleotides and any coding sequences between core genes of the pan-genome. 
Using gene synteny is not a novel concept and is used by many pan-genome tools and comparable methods to facilitate pan-genome accuracy and analysis [@page:2014; @bayliss:2019; @tonkin-hill:2020; @bazin:2020; @beier:2022]. The novelty of Corekaburra is its focus on core genes and defining regions using these. Core genes are 'stable' in occurrence across genomes of a pan-genome, making them good reference points for comparisons across genomes. Additionally, Corekaburra is not associated with a single pan-genome pipeline, but takes a general input defined by presence and absence of genes plus the genome annotations in a standard format used by popular pan-genome pipelines. As long as the two above input formats can be supplied Corekaburra is agnostic to the specifics of the pan-genome tools used. This allows for easy adaptation of Corekaburra to current, future, or custom pan-genome pipelines.

# Acknowledgements
MGJ is supported by The Melbourne Research Scholarship from The University of Melbourne. MRD is supported by a University of Melbourne CR Roper Fellowship.

# References