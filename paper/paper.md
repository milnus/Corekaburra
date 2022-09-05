---
title: 'Corekaburra: Explicit use of core gene synteny to summarise and analyse pan-genomes'
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
Since the 'first' pan-genome analysis [@tettelin:2005] this method has been instrumental in understanding the combined genetic content of multiple genomes. Through time new pan-genome tools have been developed cementing the method and improving it in many ways [@page:2014; @tonkin-hill:2020; @thorpe:2018; @inman:2019; @gautreau:2020]. However, few tools have been made to reduce the complexity of a pan-genome and allow users to analyse them in depth [@brynildsrud:2016; @whelan:2020; @lees:2018]. Here we present Corekaburra, a tool to introduce core gene synteny (order of core genes along the chromosome(s)), that defines regions of genomes by core gene pairs. This gives a summarised version of the pan-genome, and allows for systematic comparisons across genomes, easing analysis of an otherwise complex dataset.    

# Statement of need
Since the early days of genome sequencing it became clear that in bacteria variation of encoded proteins existed even within a species [@welch:2002]. The analysis of the encoded protein sequences and other genes later became the study of pan-genomes, which is now a standard method utilised for bacterial genomics [@Medini:2020]. Even though pan-genomes have become a stable in analysing genome collections, methods to systematically dissect these are sparse.  

We propose Corekaburra as a natural extension to pan-genome pipelines as a way to reduce the complexity of pan-genome outputs, ease analysis, and hence facilitate discovery. Using gene synteny (order of genes across the chromosome(s)) is not a novel concept and is used by many pan-genome tools and comparable methods to facilitate pan-genome accuracy and analysis [@page:2014; @bayliss:2019; @tonkin-hill:2020; @bazin:2020; @beier:2022]. The novelty of Corekaburra is its focus on core genes. As these are 'stable' in occurrence across genomes of a pan-genome, lending themselves as good references for comparisons across genomes. Further, Corekaburra is not associated with a single pan-genome pipeline, but takes a general input defined by presence and absence of genes plus the genome annotations in a standard format used by popular pan-genome pipelines. As long as the two above input formats can be supplied Corekaburra is agnostic to the specifics of the pan-genome tools used. This allows for easy adaptation of Corekaburra to future or custom  pan-genome methods.


# Acknowledgements
MGJ is supported by The Melbourne Research Scholarship from The University of Melbourne. MRD is supported by a University of Melbourne CR Roper Fellowship.

# References