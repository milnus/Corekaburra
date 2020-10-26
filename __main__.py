from .parse_gene_presence_absence import read_gene_presence_absence

#def get_arguments():
    # TODO - Set up arguments parser
    # ARGUMENTS:
    #

def main():
    # TODO - Get arguments

    # TODO Check if all gff files are present in input folder

    # TODO - Open gene_presence_absence file and return dict with a key for each core gene cluster and all locus_tags as the value for each key.
    core_dict, low_freq_dict = read_gene_presence_absence(# TODO - Give arguments)

    # TODO for-loop over each gff - Try to multiprocess
    # TODO Parse gff and extract core and low frequency genes from gffs

    # TODO Get all results and mash together into one.

    ### FUNCTION ###
    # TODO Get the synteny of genes if genome is complete with score 1-n_core_genes
    # TODO Record all contig breaks after core genes
    # TODO Record all neighbouring genes with score 1
    # TODO Record distance between all neighbouring core genes
    # TODO Record number of Accessory genes between neighbouring core genes
    ################

    ### DETERMINE CONSENSUS ###
    # TODO Determine start cluster from possible consensus from complete genomes - else determine relative consensus from connections
    # TODO for-loop - Determine following clusters from connections to first cluster
    # TODO Determine all alternative connections between core genes and their frequency.
    ###########################

    ### DO CALCULATIONS ###
    # TODO mean number length between core genes
    # TODO SD of length between core genes
    #######################

    ### WRITE OUTPUTS ###
    # TODO possibly construct pseudo core with core-core distances
    #####################

    ### WRITE PLOTS ###

    ###################