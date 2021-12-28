import os
from numpy import mean, median


def construct_consensus_alignment(consesus_genome_synteny, alignment_folder, core_neighbour_distance):
    con_gen_aln = {}

    # Add in the first gene and the genomes to the alignment dict
    with open(os.path.join(alignment_folder, f'{consesus_genome_synteny[0]}.aln.fas'), 'r') as alignment:
        for line in alignment.readlines():
            if '>' in line:
                genome_name = line.split(';')[0]
                genome_name = genome_name.split('>')[1]
                con_gen_aln[genome_name] = []
            else:
                con_gen_aln[genome_name].append(line.strip())

        # Add in the intergenic sequence to next gene
        # # Get a random key to get an alignments
        # rand_genome = list(con_gen_aln.keys())[0]

        # Get length to match 60 character fasta lines.
        # miss_length, full_lines, remaining_chars = adjust_character(con_gen_aln, rand_genome, core_neighbour_distance, consesus_genome_synteny=consesus_genome_synteny)
    alignment.close()

    # TODO: Construct the output file
    # with open('consensus_core_gene_alignment', 'w') as XXX:

    # TODO: Constrict a dict that keeps track of the gene, its coordinates and its strand.
    #   Output in tsv format.

    # Loop through the consensus core genome using enumerate
    for i, gene in enumerate(consesus_genome_synteny[1:]):
        # add N as intergenetic distance
        core_gene_pair = sorted([consesus_genome_synteny[i-1], consesus_genome_synteny[i]])
        intergen_dist = int(mean(core_neighbour_distance[f'{core_gene_pair[0]}--{core_gene_pair[1]}']))

        # Check if the distance is negative and change it to 0
        if intergen_dist < 0:
            intergen_dist = 0

        for genome in con_gen_aln.keys():
            con_gen_aln[genome].append("N" * intergen_dist)

        # Add next gene
        with open(os.path.join(alignment_folder, f'{gene}.aln.fas'), 'r') as alignment:
            for line in alignment.readlines():
                if '>' in line:
                    genome_name = line.split(';')[0]
                    genome_name = genome_name.split('>')[1]
                else:
                    con_gen_aln[genome_name].append(line.strip())

        # Add last intergenic sequence
        if i == len(consesus_genome_synteny)-2:
            # add N as intergenetic distance
            core_gene_pair = sorted([consesus_genome_synteny[0], consesus_genome_synteny[-1]])
            intergen_dist = int(mean(core_neighbour_distance[f'{core_gene_pair[0]}--{core_gene_pair[1]}']))
            # print(core_neighbour_distance[f'{core_gene_pair[0]}--{core_gene_pair[1]}'])
            # print(intergen_dist)

            # Check if the distance is negative and change it to 0
            if intergen_dist < 0:
                intergen_dist = 0

            for genome in con_gen_aln.keys():
                con_gen_aln[genome].append("N" * intergen_dist)

    # Write the alignment to file
    # Open the file
    with open('consensus_core_gene_alignment.aln.fas', 'w') as con_aln_file:
        # Loop through genomes in dict
        for genome in con_gen_aln:
            cur_genome_seq = con_gen_aln[genome].copy()
            # Writer header line
            con_aln_file.write(f'>{genome}\n')
            # take a line
            cur_line = cur_genome_seq.pop(0)
            # While there are more lines in the current genome, add them in.
            while len(cur_genome_seq) != 0 or len(cur_line) >= 60:
                # if line current line >= 60 characters, write line else add next set of the sequence
                if len(cur_line) >= 60:
                    cur_line, add_string = cur_line[60:], cur_line[0:60]
                    con_aln_file.write(add_string + '\n')
                else:
                    # pop next line in genome and add to current line
                    cur_line = cur_line + cur_genome_seq.pop(0)
    con_aln_file.close()

    #print(sum([len(x) for x in con_gen_aln[list(con_gen_aln.keys())[0]]]))
    # print(con_gen_aln[list(con_gen_aln.keys())[0]])


        # Check if i is equal to number of genes in the consensus genome
            # if then add the last distance (from last to first gene)

        # Following may be a seperate function to keep it DRY with the above of adding the last to first segment.
        # Take the current gene and the next - sort and paste together using -- to get key to distance dict

        # Find the mean/median or other distance between core genes.

        # Append the number of N's to the file that the distance dictates


## TODO - Panaroo seems to use only the fragment with the least gaps "-" if two fragments of a single gene is found.
#   * can we combine some of them based on genome placement?