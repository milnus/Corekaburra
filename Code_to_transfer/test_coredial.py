import unittest
import warnings
from hypothesis import given
import hypothesis.strategies as st
from gff_parser import get_genome_size_from_gff, segment_genome_content, get_contig_lengths, record_core_core_region, connect_first_n_last_gene_on_contig
from parse_gene_presence_absence import read_gene_presence_absence
from merge_dicts import merge_dicts_lists, merge_dicts_counts
from consesus_core_genome import characterise_rearrangements
from check_inputs import define_input_source, check_gene_data, check_gene_alignments, check_gff_in_pan
from correct_gffs import read_gene_data, extract_genome_fasta
from read_complete_genome_file import parse_complete_genome_file
from consesus_core_genome import construct_core_graph, identify_segments, identify_no_accessory_segments
from random import randint, choices
import os
import glob
from numpy import arange, ceil
import networkx as nx


class TestPresenceAbsenceParser(unittest.TestCase):

    def test_parser_core_genes(self):
        core_genes, accessory_genes, _, _ = read_gene_presence_absence("/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
                                                                 1, 0.05, verbose=False)
        keys = [key for key in core_genes.keys()]
        self.assertEqual(len(core_genes[keys[1]]), 10)

    def test_parser_low_freq_genes(self):
        core_genes, accessory_genes, _, _ = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        low_frew_genes = 0
        for key in accessory_genes.keys():
            low_frew_genes += len(accessory_genes[key])
        self.assertEqual(low_frew_genes, 5)


class TestInputGiven(unittest.TestCase):
    # Test pairing of all files in pan genome
    def test_input_gff_pres_abs_full_pairing(self):
        input_pres_abs = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_pan_folder/gene_presence_absence_roary.csv'
        input_file_list = ['Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff', 'Ajwa_the_Shigella.gff',
                           'Ajwa_the_Legionella.gff', 'Cari_the_Listeria.gff', 'Aman_the_Streptococcus.gff',
                           'Zion_the_Streptococcus.gff', 'Dina_the_Shigella.gff', 'Silas_the_Legionella.gff',
                           'Lilly_the_Shigella.gff', 'Chantal_the_Listeria.gff', 'Cari_the_Shigella.gff',
                           'Cari_the_Legionella.gff', 'Aman_the_Shigella.gff', 'Ajwa_the_Streptococcus.gff',
                           'Aman_the_Legionella.gff', 'Zayan_the_Shigella.gff', 'Chantal_the_Salmonella.gff',
                           'Silas_the_Shigella.gff', 'Zayan_the_Legionella.gff']

        return_bool = check_gff_in_pan(input_file_list, input_pres_abs)

        self.assertEqual(return_bool, True)

    # Test pairing of some files in pan genome - Warning
    def test_input_gff_pres_abs_partial_pairing_catch_warning(self):
        input_pres_abs = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_pan_folder/gene_presence_absence_roary.csv'
        input_file_list = ['Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff', 'Ajwa_the_Shigella.gff']

        with self.assertWarns(Warning):
            check_gff_in_pan(input_file_list, input_pres_abs)

    # Test pairing of some files in pan genome
    def test_input_gff_pres_abs_partial_pairing(self):
        input_pres_abs = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_pan_folder/gene_presence_absence_roary.csv'
        input_file_list = ['Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff', 'Ajwa_the_Shigella.gff']

        return_bool = check_gff_in_pan(input_file_list, input_pres_abs)

        self.assertEqual(return_bool, True)

    # Test when given a file not in pan genome among others that are in the pan genome
    def test_input_gff_pres_abs_file_not_in_pan(self):
        input_pres_abs = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_pan_folder/gene_presence_absence_roary.csv'
        input_file_list = ['Cappuccino.gff', 'Flat_white.gff', 'Doubble_espresso.gff']

        with self.assertRaises(FileNotFoundError):
            check_gff_in_pan(input_file_list, input_pres_abs)

    def test_input_gff_pres_abs_some_file_not_in_pan(self):
        input_pres_abs = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_pan_folder/gene_presence_absence_roary.csv'
        input_file_list = ['Cappuccino.gff', 'Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff']

        with self.assertRaises(FileNotFoundError):
            check_gff_in_pan(input_file_list, input_pres_abs)


class TestGffparser(unittest.TestCase):

    def test_get_genome_size_from_gff(self):
        genome_length = get_genome_size_from_gff(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_locating_inputs/GCA_000006785.gff")

        print(genome_length)
        self.assertEqual(1852433, genome_length)

    def test_segmentation_core_gene_number_all_complete(self):
        core_genes, accessory_genes, acc_gene_dict, _ = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        pre_path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/"
        gff_files = glob.glob(pre_path+"*.gff")

        core_result_dict = {}

        for file in gff_files:
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), i=1, acc_gene_dict=acc_gene_dict)
            core_result_dict = merge_dicts_counts(core_result_dict, core_gene_pairs)

        self.assertEqual(len(core_result_dict.keys()), 10)

    def test_segmentation_distance_number(self):
        core_genes, accessory_genes, acc_gene_dict, _ = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        pre_path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/"
        gff_files = glob.glob(pre_path+"*.gff")

        distance_result_dict = {}

        for file in gff_files:
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), acc_gene_dict=acc_gene_dict, i=1)
            distance_result_dict = merge_dicts_lists(distance_result_dict, core_gene_pair_distance)

        self.assertEqual(len(distance_result_dict.keys()), 10)


    def test_segmentation_accessory_gene_number(self):
        core_genes, accessory_genes, acc_gene_dict, _ = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        pre_path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/"
        gff_files = glob.glob(pre_path+"*.gff")

        accessory_result_dict = {}

        for file in gff_files:
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), acc_gene_dict=acc_gene_dict, i=1)
            accessory_result_dict = merge_dicts_lists(accessory_result_dict, accessory_gene_content)
        accessory_gene_amount = []
        for i, key in enumerate(accessory_result_dict.keys()):
            if i == 10:
                self.assertListEqual([7,6,7,8,7,5,7,9,7,7,6,7,7,6,8,3,7,8,4,5], accessory_result_dict[key])


# TODO - constuct test that tests load a single gff file

class TestLocatingInput_files(unittest.TestCase):

    def test_locate_core_gene_alignments(self):
        input_folder = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_locating_inputs/test_pan_alignment_complete'
        core_dict, _, _, _ = read_gene_presence_absence(os.path.join(input_folder, 'gene_presence_absence_roary.csv'),
                                                        1, 0.01, verbose=False)

        return_value = check_gene_alignments(input_folder, core_dict)
        self.assertEqual(os.path.join(input_folder, 'aligned_gene_sequences'), return_value)


        incomplete_input_folder = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_locating_inputs/test_pan_alignment_incomplete_alignments'
        with self.assertRaises(FileNotFoundError):
            check_gene_alignments(incomplete_input_folder, core_dict)


# TODO - constuct test that loads a mock gene_prensence_absence file
    def test_loading_gene_pres_abs(self):
        # TODO - Construct mock gene_presence_absence_file
        # TODO - Load Panaroo file
        # Load mock Panaroo gene presence absence file
        read_core, read_low_freq, read_acc, _ = read_gene_presence_absence("/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_gene_pres_abs_parser/read_mock_gene_presence_absence.csv",
                                                                  1, 0.05, verbose=False)

        # Test number of keys to see if all genomes are read in
        self.assertEqual(len(read_core), 5)

        # Test accessory genes
        # Extract acc genes found
        ace_genes = [gene for x in read_acc.keys() for y in read_acc[x].keys() for gene in read_acc[x][y]]
        self.assertEqual(len(set(ace_genes)), 9)

        ## Test if the expected number of core genes are found
        for i in arange(0.8, 1, 0.1):
            read_core, _, _, _ = read_gene_presence_absence(
                "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_gene_pres_abs_parser/read_mock_gene_presence_absence.csv",
                i, 0.05, verbose=False)

            # Extract core genes found
            core_genes = [gene for x in read_core.keys() for y in read_core[x].keys() for gene in read_core[x][y]]

            # Test depending on percent presence
            if i == 1:
                self.assertEqual(len(set(core_genes)), 3)
            else:
                self.assertEqual(len(set(core_genes)), 5)


        # Test low frequency genes
        for i in arange(0.1, 0.8, 0.1):
            _, read_low_freq, _, _ = read_gene_presence_absence(
                "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_gene_pres_abs_parser/read_mock_gene_presence_absence.csv",
                1, i, verbose=False)

            low_genes = [gene for x in read_low_freq.keys() for y in read_low_freq[x].keys() for gene in read_low_freq[x][y]]

            low_gene_presence = ceil(5*i)

            # Test depending on percent presence
            if low_gene_presence == 1:
                self.assertEqual(len(set(low_genes)), 4)
            elif low_gene_presence == 2:
                self.assertEqual(len(set(low_genes)), 9)
            elif low_gene_presence == 2:
                self.assertEqual(len(set(low_genes)), 11)
            elif low_gene_presence == 2:
                self.assertEqual(len(set(low_genes)), 13)


        # TODO - Load Roary file

    # TODO -

    # def test_rearrangement_predictions(self):
    #     consus_genome = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    #     alternative_core_pairs = {'D--G': 1, 'E--H': 1}
    #
    #     characterise_rearrangements(consus_genome, alternative_core_pairs)


    def test_soruce_identification(self):
        # Identify Roary file
        path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Test_source_identification/Roray_gene_presence_absence"
        source, file_path = define_input_source(path)
        self.assertEqual(source, "Roary")
        self.assertEqual(file_path,
                         "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Test_source_identification/Roray_gene_presence_absence/gene_presence_absence.csv")

        # Identify Panaroo file
        path = '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Test_source_identification/Panaroo_gene_presence_absence'
        source, file_path = define_input_source(path)
        self.assertEqual(source, "Panaroo")
        self.assertNotEqual(file_path,
                         '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Test_source_identification/Panaroo_gene_presence_absence/gene_presence_absence.csv')
        self.assertEqual(file_path,
                            '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Test_source_identification/Panaroo_gene_presence_absence/gene_presence_absence_roary.csv')

    def test_locate_gene_data(self):
        check_return = check_gene_data('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/same_genome_pan_split')

        self.assertEqual(os.path.isfile(check_return), True)

        with self.assertRaises(FileNotFoundError):
            check_gene_data('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome')

    def test_reading_gene_data_csv(self):
        gene_data = read_gene_data('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_gene_data_reading/test_gene_data.csv')

        expected_dict = {
            'GCA_900636485': {'49_refound_1369': 'ATGACGATAGATGAAGCGTTGCAAAATTTACGTGATAACTTTAATAAAATAATGAATGTCCTAAAAAACGATTGGAAAGCACTATTGTTTCTTGCAATCACAATATTTGGGATGATGGTAACCGTGTCGTATTTTAGCTATCGCGACGCACGACAATATTACGAGTCGCAAATCACAGGACTACGTACACAGCTAAGCAGGACACAAAAGCAGCTTAAACGTGCTAGCGAAGATAGAGCTAGACAGACAAAGCGGATTGCGGAACTTACGCACAACGGAGGGTAG'},
            'GCA_001019635': {'8_refound_250': 'ATGGAACCAAAATTACATCGGCAACTGCGTCAAAAATATGACGACGCTGAAAAACAATATCTTGAAAAGTTTGGAGAAGACTCGCTTGATAGAGTATTTTTTTGGGAGCCAGACGTTTACTTTGATGAGTGGAAAAAGGTTCTACCAGATGCAACACTGGAATTAAACAAAGCTATTAATAGCGGGGTGGCGATTGATCCAGATCCAGAAAACGCAATATATTAA',
                              '8_refound_251': 'ATGAAAAGCTTTTTAAATTTAGTCAAACAAAAGTTGTTTAAACCAGGTCTAAAAAAACTCGTAAAGCTTCACAACTCCCAGAACGTTAATATATGCTTATATATCAACGATTGGAACTAATTTATGGTTCGCACCATGGTTTTTGTGGAAGGATCAAAAGTTGTCCTGAAAATTTCTCTTAACCGTGTTTAA'}
        }

        self.assertEqual(expected_dict, gene_data)

    def test_reading_fasta_from_gff(self):
        single_genome_dict, _, _ = extract_genome_fasta('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Sample_gff_files/single_contig_mock_for_reading_fasta_in_gff.gff')
        self.assertEqual(len(single_genome_dict), 1)
        self.assertEqual(len(single_genome_dict[list(single_genome_dict.keys())[0]]), 540)

        multi_genome_dict, _, _ = extract_genome_fasta('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Sample_gff_files/multi_contig_mock_for_reading_fasta_in_gff.gff')
        self.assertEqual(len(multi_genome_dict), 2)
        self.assertEqual(len(multi_genome_dict[list(multi_genome_dict.keys())[0]]), 540)
        self.assertEqual(len(multi_genome_dict[list(multi_genome_dict.keys())[1]]), 300)

    def test_finding_largest_locus_tag(self):
        _, locus_tag, _ = extract_genome_fasta('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Sample_gff_files/single_contig_mock_for_reading_fasta_in_gff.gff')

        self.assertEqual(locus_tag, 'MONDJAPC_01960')


class TestGffParsing(unittest.TestCase):
    def test_contig_length(self):
        contig_dir = get_contig_lengths('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_length_of_contigs/contigs_length_test_pass.txt')

        contig_lengths = [contig_dir[contig] for contig in contig_dir.keys()]

        true_contig_length = [100, 10, 3, 5000]

        self.assertEqual(contig_lengths, true_contig_length)

    def test_contig_name_identification(self):
        contig_dir = get_contig_lengths(
            '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_length_of_contigs/contigs_length_test_pass.txt')

        true_contig_length = ['contig_1', 'test_contig_2', '3', '4']

        self.assertEqual(true_contig_length, list(contig_dir.keys()))

    def test_duplicate_contig_name(self):
        with self.assertRaises(ValueError):
            get_contig_lengths('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_length_of_contigs/test_duplicate_contig_names.txt')


class TestCompleteGenomeParser(unittest.TestCase):
    def test_parsing_file_no_extension(self):
        expected_names = ['complete_genome1', 'comp.genome', '2']

        complete_genomes = parse_complete_genome_file(
            '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_parsing_complete_genomes/complete_genomes.txt',
            expected_names)

        self.assertEqual(complete_genomes, expected_names)

    def test_parsing_file_with_input_extension(self):
        expected_names = ['complete_genome1', 'comp.genome', '2']
        input_gffs = ['complete_genome1.gff', 'comp.genome.gff', '2.gff']

        complete_genomes = parse_complete_genome_file(
            '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_parsing_complete_genomes/complete_genomes_w_extensions.txt',
            input_gffs)

        self.assertEqual(complete_genomes, expected_names)

    def test_parsing_file_with_extension(self):
        expected_names = ['complete_genome1', 'comp.genome', '2']
        input_gffs = ['complete_genome1.gff', 'comp.genome.gff', '2.gff']

        complete_genomes = parse_complete_genome_file(
            '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_parsing_complete_genomes/complete_genomes.txt',
            input_gffs)

        self.assertEqual(complete_genomes, expected_names)

    def test_parsing_file_with_path(self):
        expected_names = ['complete_genome1', 'comp.genome', '2']
        input_gffs = ['/test/path/complete_genome1', '/all/the/paths/comp.genome', '/no/more/paths/2']

        complete_genomes = parse_complete_genome_file(
            '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_parsing_complete_genomes/complete_genomes.txt',
            input_gffs)

        self.assertEqual(complete_genomes, expected_names)

    def test_parsing_file_with_path_and_extension(self):
        expected_names = ['complete_genome1', 'comp.genome', '2']
        input_gffs = ['/test/path/complete_genome1.gff', '/all/the/paths/comp.genome.gff', '/no/more/paths/2.gff']

        complete_genomes = parse_complete_genome_file(
            '/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_parsing_complete_genomes/complete_genomes.txt',
            input_gffs)

        self.assertEqual(complete_genomes, expected_names)

# class TestCoreGeneSyntenyTypes(unittest.TestCase):
#     def test_assign_core_synteny_types(self):
#         thing = []
#
#         self.assertEqual(len(set(thing)), 1)
#
#         thing_90 = []
#         self.assertEqual(len(set(thing_90)), 2)


class TestRecordingCoreGene(unittest.TestCase):
    def test_regular_core_gene_recording(self):
        expected_previous_core_gene_id ='gene_1'
        expected_previous_core_gene_end_coor = 200
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_cluster_1--pan_cluster_2': 99}
        expected_accessory_gene_content = {'pan_cluster_1--pan_cluster_2': ['Acc_1', 'Acc_2']}
        expected_low_freq_gene_content = {'pan_cluster_1--pan_cluster_2': ['low_1']}
        expected_core_gene_pairs = ['pan_cluster_1--pan_cluster_2']
        expected_num_acc_genes_in_region = {'pan_cluster_1--pan_cluster_2': 3}
        expected_master_info = {'pan_cluster_1--pan_cluster_2--genome_1': ['genome_1', 'pan_cluster_1', 'pan_cluster_2', 99, 3, ['Acc_1', 'Acc_2'], ['low_1']]}

        (previous_core_gene_id,
         previous_core_gene_end_coor,
         acc_genes_in_region,
         low_freq_genes_in_region,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content,
         core_gene_pairs,
         num_acc_genes_in_region,
         master_info) = record_core_core_region(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                gff_name='genome_1',
                                                gff_line=['contig_1', '', '', '100', '200', '', '', '', 'gene_1'],
                                                contig_end=None,
                                                previous_core_gene_id='gene_2',
                                                previous_core_gene_end_coor=0,
                                                acc_genes_in_region=['Acc_1', 'Acc_2'],
                                                low_freq_genes_in_region=['low_1'],
                                                core_gene_pair_distance={},
                                                accessory_gene_content={},
                                                low_freq_gene_content={},
                                                core_gene_pairs=[],
                                                num_acc_genes_in_region={},
                                                master_info={})

        self.assertEqual(expected_previous_core_gene_id, previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, low_freq_genes_in_region)
        self.assertEqual(expected_core_gene_pair_distance, core_gene_pair_distance)
        self.assertEqual(expected_accessory_gene_content, accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, low_freq_gene_content)
        self.assertEqual(expected_core_gene_pairs, core_gene_pairs)
        # self.assertEqual(expected_num_acc_genes_in_region, num_acc_genes_in_region)
        self.assertEqual(expected_master_info, master_info)

    def test_adding_last_core_at_contig_change(self):
        expected_previous_core_gene_id = 'Sequence_break'
        expected_previous_core_gene_end_coor = 100
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_cluster_2--Sequence_break': 99}
        expected_accessory_gene_content = {'pan_cluster_2--Sequence_break': ['Acc_1', 'Acc_2']}
        expected_low_freq_gene_content = {'pan_cluster_2--Sequence_break': ['low_1']}
        expected_core_gene_pairs = ['pan_cluster_2--Sequence_break']
        expected_num_acc_genes_in_region = {'pan_cluster_2--Sequence_break': 3}
        expected_master_info = {
            'pan_cluster_2--Sequence_break--genome_1': ['genome_1', 'pan_cluster_2', 'Sequence_break', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']]}
        (previous_core_gene_id,
         previous_core_gene_end_coor,
         acc_genes_in_region,
         low_freq_genes_in_region,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content,
         core_gene_pairs,
         num_acc_genes_in_region,
         master_info) = record_core_core_region(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                gff_name='genome_1',
                                                gff_line=None,
                                                contig_end=200,
                                                previous_core_gene_id='gene_2',
                                                previous_core_gene_end_coor=100,
                                                acc_genes_in_region=['Acc_1', 'Acc_2'],
                                                low_freq_genes_in_region=['low_1'],
                                                core_gene_pair_distance={},
                                                accessory_gene_content={},
                                                low_freq_gene_content={},
                                                core_gene_pairs=[],
                                                num_acc_genes_in_region={},
                                                master_info={})

        self.assertEqual(expected_previous_core_gene_id, previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, low_freq_genes_in_region)
        self.assertEqual(expected_core_gene_pair_distance, core_gene_pair_distance)
        self.assertEqual(expected_accessory_gene_content, accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, low_freq_gene_content)
        self.assertEqual(expected_core_gene_pairs, core_gene_pairs)
        # self.assertEqual(expected_num_acc_genes_in_region, num_acc_genes_in_region)
        self.assertEqual(expected_master_info, master_info)

    def test_add_first_gene_on_new_contig(self):
        expected_previous_core_gene_id = 'gene_1'
        expected_previous_core_gene_end_coor = 200
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'Sequence_break--pan_cluster_1': 99}
        expected_accessory_gene_content = {'Sequence_break--pan_cluster_1': ['Acc_1', 'Acc_2']}
        expected_low_freq_gene_content = {'Sequence_break--pan_cluster_1': ['low_1']}
        expected_core_gene_pairs = ['Sequence_break--pan_cluster_1']
        expected_num_acc_genes_in_region = {'Sequence_break--pan_cluster_1': 3}
        expected_master_info = {
            'Sequence_break--pan_cluster_1--genome_1': ['genome_1', 'Sequence_break', 'pan_cluster_1', 99, 3,
                                                        ['Acc_1', 'Acc_2'], ['low_1']]}

        (previous_core_gene_id,
         previous_core_gene_end_coor,
         acc_genes_in_region,
         low_freq_genes_in_region,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content,
         core_gene_pairs,
         num_acc_genes_in_region,
         master_info) = record_core_core_region(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                gff_name='genome_1',
                                                gff_line=['contig_1', '', '', '100', '200', '', '', '', 'gene_1'],
                                                contig_end=0,
                                                previous_core_gene_id='Sequence_break',
                                                previous_core_gene_end_coor=200,
                                                acc_genes_in_region=['Acc_1', 'Acc_2'],
                                                low_freq_genes_in_region=['low_1'],
                                                core_gene_pair_distance={},
                                                accessory_gene_content={},
                                                low_freq_gene_content={},
                                                core_gene_pairs=[],
                                                num_acc_genes_in_region={},
                                                master_info={})

        self.assertEqual(expected_previous_core_gene_id, previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, low_freq_genes_in_region)
        self.assertEqual(expected_core_gene_pair_distance, core_gene_pair_distance)
        self.assertEqual(expected_accessory_gene_content, accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, low_freq_gene_content)
        self.assertEqual(expected_core_gene_pairs, core_gene_pairs)
        # self.assertEqual(expected_num_acc_genes_in_region, num_acc_genes_in_region)
        self.assertEqual(expected_master_info, master_info)

    def test_adding_first_gene_of_genome_next_to_sequence_break(self):
        expected_core_gene_pair_distance = {'Sequence_break--pan_cluster_1': 99}
        expected_accessory_gene_content = {'Sequence_break--pan_cluster_1': []}
        expected_low_freq_gene_content = {'Sequence_break--pan_cluster_1': []}
        expected_core_gene_pairs = ['Sequence_break--pan_cluster_1']
        expected_num_acc_genes_in_region = {'Sequence_break--pan_cluster_1': 0}
        expected_master_info = {
            'Sequence_break--pan_cluster_1--genome_1': ['genome_1', 'Sequence_break', 'pan_cluster_1', 99, 0,
                                                        [], []]}

        first_core_gene_id = 'gene_1'
        first_core_gene_start_coor = 100


        (_,
         _,
         _,
         _,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content, core_gene_pairs,
         num_acc_genes_in_region, master_info) = record_core_core_region(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                                         gff_name='genome_1',
                                                                         gff_line=['contig_1', '', '', first_core_gene_start_coor, 100000, '', '', '', first_core_gene_id],
                                                                         contig_end=0,
                                                                         previous_core_gene_id='Sequence_break',
                                                                         previous_core_gene_end_coor=10000000,
                                                                         acc_genes_in_region=[],
                                                                         low_freq_genes_in_region=[],
                                                                         core_gene_pair_distance={},
                                                                         accessory_gene_content={},
                                                                         low_freq_gene_content={},
                                                                         core_gene_pairs=[],
                                                                         num_acc_genes_in_region={},
                                                                         master_info={})

        self.assertEqual(expected_core_gene_pair_distance, core_gene_pair_distance)
        self.assertEqual(expected_accessory_gene_content, accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, low_freq_gene_content)
        self.assertEqual(expected_core_gene_pairs, core_gene_pairs)
        # self.assertEqual(expected_num_acc_genes_in_region, num_acc_genes_in_region)
        self.assertEqual(expected_master_info, master_info)

    def test_adding_last_core_gene_next_to_sequence_break_in_incomplete_genome(self):
        expected_previous_core_gene_id = 'Sequence_break'
        expected_previous_core_gene_end_coor = 200
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_cluster_2--Sequence_break': 99}
        expected_accessory_gene_content = {'pan_cluster_2--Sequence_break': ['Acc_1', 'Acc_2']}
        expected_low_freq_gene_content = {'pan_cluster_2--Sequence_break': ['low_1']}
        expected_core_gene_pairs = ['pan_cluster_2--Sequence_break']
        expected_num_acc_genes_in_region = {'pan_cluster_2--Sequence_break': 3}
        expected_master_info = {
            'pan_cluster_2--Sequence_break--genome_1': ['genome_1', 'pan_cluster_2', 'Sequence_break', 99, 3,
                                                        ['Acc_1', 'Acc_2'], ['low_1']]}

        (previous_core_gene_id,
         previous_core_gene_end_coor,
         acc_genes_in_region,
         low_freq_genes_in_region,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content, core_gene_pairs,
         num_acc_genes_in_region, master_info) = record_core_core_region(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                                         gff_name='genome_1',
                                                                         gff_line=None,
                                                                         contig_end=300,
                                                                         previous_core_gene_id='gene_2',
                                                                         previous_core_gene_end_coor=200,
                                                                         acc_genes_in_region=['Acc_1', 'Acc_2'],
                                                                         low_freq_genes_in_region=['low_1'],
                                                                         core_gene_pair_distance={},
                                                                         accessory_gene_content={},
                                                                         low_freq_gene_content={},
                                                                         core_gene_pairs=[],
                                                                         num_acc_genes_in_region={},
                                                                         master_info={})


        self.assertEqual(expected_previous_core_gene_id, previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, low_freq_genes_in_region)
        self.assertEqual(expected_core_gene_pair_distance, core_gene_pair_distance)
        self.assertEqual(expected_accessory_gene_content, accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, low_freq_gene_content)
        self.assertEqual(expected_core_gene_pairs, core_gene_pairs)
        # self.assertEqual(expected_num_acc_genes_in_region, num_acc_genes_in_region)
        self.assertEqual(expected_master_info, master_info)

    def test_adding_last_core_at_contig_and_first_of_next_in_chain(self):
        expected_previous_core_gene_id = 'gene_1'
        expected_previous_core_gene_end_coor = 200
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_cluster_2--Sequence_break': 99, 'Sequence_break--pan_cluster_1': 999}
        expected_accessory_gene_content = {'pan_cluster_2--Sequence_break': ['Acc_1', 'Acc_2'], 'Sequence_break--pan_cluster_1': []}
        expected_low_freq_gene_content = {'pan_cluster_2--Sequence_break': ['low_1'], 'Sequence_break--pan_cluster_1': []}
        expected_core_gene_pairs = ['pan_cluster_2--Sequence_break', 'Sequence_break--pan_cluster_1']
        expected_num_acc_genes_in_region = {'pan_cluster_2--Sequence_break': 3, 'Sequence_break--pan_cluster_1': 0}
        expected_master_info = {
            'pan_cluster_2--Sequence_break--genome_1': ['genome_1', 'pan_cluster_2', 'Sequence_break', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'Sequence_break--pan_cluster_1--genome_1': ['genome_1', 'Sequence_break', 'pan_cluster_1', 999, 0, [], []]}


        # Add last gene before contig break
        (previous_core_gene_id,
         previous_core_gene_end_coor,
         acc_genes_in_region,
         low_freq_genes_in_region,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content,
         core_gene_pairs,
         num_acc_genes_in_region,
         master_info) = record_core_core_region(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                gff_name='genome_1',
                                                gff_line=None,
                                                contig_end=200,
                                                previous_core_gene_id='gene_2',
                                                previous_core_gene_end_coor=100,
                                                acc_genes_in_region=['Acc_1', 'Acc_2'],
                                                low_freq_genes_in_region=['low_1'],
                                                core_gene_pair_distance={},
                                                accessory_gene_content={},
                                                low_freq_gene_content={},
                                                core_gene_pairs=[],
                                                num_acc_genes_in_region={},
                                                master_info={})


        # Use info from previous gene to add the first core gene on the next contig.
        (previous_core_gene_id,
         previous_core_gene_end_coor,
         acc_genes_in_region,
         low_freq_genes_in_region,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content,
         core_gene_pairs,
         num_acc_genes_in_region,
         master_info) = record_core_core_region(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                gff_name='genome_1',
                                                gff_line=['contig_1', '', '', '1000', '200', '', '', '', 'gene_1'],
                                                contig_end=0,
                                                previous_core_gene_id=previous_core_gene_id,
                                                previous_core_gene_end_coor=previous_core_gene_end_coor,
                                                acc_genes_in_region=acc_genes_in_region,
                                                low_freq_genes_in_region=low_freq_genes_in_region,
                                                core_gene_pair_distance=core_gene_pair_distance,
                                                accessory_gene_content=accessory_gene_content,
                                                low_freq_gene_content=low_freq_gene_content,
                                                core_gene_pairs=core_gene_pairs,
                                                num_acc_genes_in_region=num_acc_genes_in_region,
                                                master_info=master_info)

        # Assess how the test went.
        self.assertEqual(expected_previous_core_gene_id, previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, low_freq_genes_in_region)
        self.assertEqual(expected_core_gene_pair_distance, core_gene_pair_distance)
        self.assertEqual(expected_accessory_gene_content, accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, low_freq_gene_content)
        self.assertEqual(expected_core_gene_pairs, core_gene_pairs)
        # self.assertEqual(expected_num_acc_genes_in_region, num_acc_genes_in_region)
        self.assertEqual(expected_master_info, master_info)

    def test_recording_last_n_first_core_on_closed_contig(self):
        expected_previous_core_gene_id = 'Complete_genome_end_fail'
        expected_previous_core_gene_end_coor = 100
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pairs = ['pan_cluster_1--pan_cluster_2']
        expected_core_gene_pair_distance = {'pan_cluster_1--pan_cluster_2': 100}
        expected_accessory_gene_content = {'pan_cluster_1--pan_cluster_2': ['Acc_2', 'Acc_3', 'Acc_1']}
        expected_low_freq_gene_content = {'pan_cluster_1--pan_cluster_2': ['low_1', 'low_2']}
        expected_master_info = {'pan_cluster_1--pan_cluster_2--genome_1': ['genome_1', 'pan_cluster_1', 'pan_cluster_2', 100, 5,
                                                        ['Acc_2', 'Acc_3', 'Acc_1'], ['low_1', 'low_2']]}

        (previous_core_gene_id,
         previous_core_gene_end_coor,
         acc_genes_in_region,
         low_freq_genes_in_region,
         core_gene_pairs,
         core_gene_pair_distance,
         accessory_gene_content,
         low_freq_gene_content,
         master_info) = connect_first_n_last_gene_on_contig(core_genes={'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'}},
                                                            gff_name='genome_1',
                                                            previous_core_gene_id='gene_2',
                                                            previous_core_gene_end_coor=100,
                                                            first_core_gene_gff_line=['contig_1', '', '', '0', '50', '', '', '', 'gene_1'],
                                                            acc_genes_in_region=['Acc_2', 'Acc_3'],
                                                            first_core_accessory_content=['Acc_1'],
                                                            low_freq_genes_in_region=['low_1'],
                                                            first_core_low_freq_genes=['low_2'],
                                                            contig_size=200,
                                                            core_gene_pairs=[],
                                                            core_gene_pair_distance={},
                                                            accessory_gene_content={},
                                                            low_freq_gene_content={},
                                                            master_info={})

        self.assertEqual(expected_previous_core_gene_id, previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, low_freq_genes_in_region)
        self.assertEqual(expected_core_gene_pairs, core_gene_pairs)
        self.assertEqual(expected_core_gene_pair_distance, core_gene_pair_distance)
        self.assertEqual(expected_accessory_gene_content, accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, low_freq_gene_content)
        self.assertEqual(expected_master_info, master_info)


class TestSegmentationIdentification(unittest.TestCase):
    def test_core_gene_graph_construction_circle_case(self):
        expected_edges = [('pan_cluster_1', 'pan_cluster_2'), ('pan_cluster_1', 'pan_cluster_6'), ('pan_cluster_2', 'pan_cluster_3'),
         ('pan_cluster_3', 'pan_cluster_4'), ('pan_cluster_4', 'pan_cluster_5'), ('pan_cluster_5', 'pan_cluster_6')]

        expected_degrees = [('pan_cluster_1', 2), ('pan_cluster_2', 2), ('pan_cluster_3', 2), ('pan_cluster_4', 2), ('pan_cluster_5', 2), ('pan_cluster_6', 2)]

        expected_edge_weights = [10, 10, 10, 10, 10, 10]

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 10,
                                'pan_cluster_2--pan_cluster_3': 10,
                                'pan_cluster_3--pan_cluster_4': 10,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_1--pan_cluster_6': 10}

        core_graph = construct_core_graph(core_neighbour_pairs)

        # Get edge weights:
        edge_weights = [core_graph.get_edge_data(edge[0], edge[1])['weight'] for edge in list(core_graph.edges)]

        # Assert outputs
        self.assertEqual(expected_edges, list(core_graph.edges))
        self.assertEqual(expected_degrees, list(core_graph.degree))
        self.assertEqual(expected_edge_weights, edge_weights)

    def test_core_gene_graph_construction_circle_case_with_single_break(self):
        expected_edges = [('pan_cluster_1', 'pan_cluster_2'), ('pan_cluster_1', 'pan_cluster_6'), ('pan_cluster_2', 'pan_cluster_3'),
         ('pan_cluster_3', 'pan_cluster_4'), ('pan_cluster_4', 'pan_cluster_5'), ('pan_cluster_5', 'pan_cluster_6')]

        expected_degrees = [('pan_cluster_1', 2), ('pan_cluster_2', 2), ('pan_cluster_3', 2), ('pan_cluster_4', 2), ('pan_cluster_5', 2), ('pan_cluster_6', 2)]

        expected_edge_weights = [9, 10, 9, 10, 10, 10]

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 9,
                                'pan_cluster_1--Sequence_break': 1,
                                'Sequence_break--pan_cluster_2': 1,
                                'pan_cluster_2--pan_cluster_3': 9,
                                'pan_cluster_3--pan_cluster_4': 10,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_1--pan_cluster_6': 10}

        core_graph = construct_core_graph(core_neighbour_pairs)

        # Get edge weights:
        edge_weights = [core_graph.get_edge_data(edge[0], edge[1])['weight'] for edge in list(core_graph.edges)]

        # Assert outputs
        self.assertEqual(expected_edges, list(core_graph.edges))
        self.assertEqual(expected_degrees, list(core_graph.degree))
        self.assertEqual(expected_edge_weights, edge_weights)

    def test_core_gene_graph_construction_three_degree_case(self):
        expected_edges = [('pan_cluster_1', 'pan_cluster_2'), ('pan_cluster_1', 'pan_cluster_6'), ('pan_cluster_1', 'pan_cluster_4'), ('pan_cluster_2', 'pan_cluster_3'),
         ('pan_cluster_3', 'pan_cluster_4'), ('pan_cluster_4', 'pan_cluster_5'), ('pan_cluster_5', 'pan_cluster_6')]

        expected_degrees = [('pan_cluster_1', 3), ('pan_cluster_2', 2), ('pan_cluster_3', 2), ('pan_cluster_4', 3), ('pan_cluster_5', 2), ('pan_cluster_6', 2)]

        expected_edge_weights = [10, 8, 2, 10, 8, 10, 10]

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 10,
                                'pan_cluster_2--pan_cluster_3': 10,
                                'pan_cluster_3--pan_cluster_4': 8,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_1--pan_cluster_6': 8,
                                'pan_cluster_1--pan_cluster_4': 2}

        core_graph = construct_core_graph(core_neighbour_pairs)

        # Get edge weights:
        edge_weights = [core_graph.get_edge_data(edge[0], edge[1])['weight'] for edge in list(core_graph.edges)]

        # Assert outputs
        self.assertEqual(expected_edges, list(core_graph.edges))
        self.assertEqual(expected_degrees, list(core_graph.degree))
        self.assertEqual(expected_edge_weights, edge_weights)

    def test_core_gene_graph_construction_three_degree_n_sequence_breaks_case(self):
        expected_edges = [('pan_cluster_1', 'pan_cluster_2'), ('pan_cluster_1', 'pan_cluster_6'), ('pan_cluster_1', 'pan_cluster_4'), ('pan_cluster_2', 'pan_cluster_3'),
         ('pan_cluster_3', 'pan_cluster_4'), ('pan_cluster_4', 'pan_cluster_5'), ('pan_cluster_5', 'pan_cluster_6')]

        expected_degrees = [('pan_cluster_1', 3), ('pan_cluster_2', 2), ('pan_cluster_3', 2), ('pan_cluster_4', 3), ('pan_cluster_5', 2), ('pan_cluster_6', 2)]

        expected_edge_weights = [10, 8, 2, 10, 8, 10, 10]

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 10,
                                'pan_cluster_2--pan_cluster_3': 10,
                                'pan_cluster_3--pan_cluster_4': 8,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_1--pan_cluster_6': 8,
                                'pan_cluster_1--pan_cluster_4': 2,
                                'pan_cluster_3--Sequence_break': 2,
                                'pan_cluster_6--Sequence_break': 2}

        core_graph = construct_core_graph(core_neighbour_pairs)

        # Get edge weights:
        edge_weights = [core_graph.get_edge_data(edge[0], edge[1])['weight'] for edge in list(core_graph.edges)]

        # Assert outputs
        self.assertEqual(expected_edges, list(core_graph.edges))
        self.assertEqual(expected_degrees, list(core_graph.degree))
        self.assertEqual(expected_edge_weights, edge_weights)

    def test_double_edge_segment_identification_all_2_degree_input(self):
        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 10,
                                'pan_cluster_2--pan_cluster_3': 10,
                                'pan_cluster_3--pan_cluster_4': 10,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_6--pan_cluster_1': 10}

        core_graph = construct_core_graph(core_neighbour_pairs)

        return_1, return_2, return_3 = identify_segments(core_graph, 10)

        self.assertEqual(None, return_1)
        self.assertEqual(None, return_2)
        self.assertEqual(None, return_3)


    def test_double_edge_segment_identification_two_segments(self):
        expected_segments = {'pan_cluster_1--pan_cluster_5': ['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5'], 'pan_cluster_2--pan_cluster_4': ['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 9,
                                'pan_cluster_1--pan_cluster_4': 1,
                                'pan_cluster_2--pan_cluster_3': 10,
                                'pan_cluster_3--pan_cluster_4': 10,
                                'pan_cluster_2--pan_cluster_5': 1,
                                'pan_cluster_4--pan_cluster_5': 9,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_6--pan_cluster_1': 10}

        core_graph = construct_core_graph(core_neighbour_pairs)

        double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(core_graph, 10)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_four_segments(self):
        expected_segments = {'pan_cluster_1--pan_cluster_3': ['pan_cluster_1', 'pan_cluster_2', 'pan_cluster_3'],
                             'pan_cluster_1--pan_cluster_9': ['pan_cluster_1', 'pan_cluster_10', 'pan_cluster_9'],
                             'pan_cluster_3--pan_cluster_6': ['pan_cluster_6', 'pan_cluster_5', 'pan_cluster_4', 'pan_cluster_3'],
                             'pan_cluster_6--pan_cluster_9': ['pan_cluster_6', 'pan_cluster_7', 'pan_cluster_8', 'pan_cluster_9']}

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 9,
                                'pan_cluster_1--pan_cluster_6': 1,
                                'pan_cluster_2--pan_cluster_3': 10,
                                'pan_cluster_3--pan_cluster_4': 9,
                                'pan_cluster_3--pan_cluster_9': 1,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_6--pan_cluster_7': 10,
                                'pan_cluster_7--pan_cluster_8': 10,
                                'pan_cluster_8--pan_cluster_9': 9,
                                'pan_cluster_9--pan_cluster_10': 10,
                                'pan_cluster_1--pan_cluster_10': 10}

        core_graph = construct_core_graph(core_neighbour_pairs)

        double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(core_graph, 10)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_four_degrees(self):
        print(4)
        expected_segments = {'pan_cluster_4--pan_cluster_6': ['pan_cluster_4', 'pan_cluster_5', 'pan_cluster_6']}

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 9,
                                'pan_cluster_2--pan_cluster_3': 9,
                                'pan_cluster_2--pan_cluster_4': 1,
                                'pan_cluster_2--pan_cluster_6': 1,
                                'pan_cluster_3--pan_cluster_4': 9,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_6--pan_cluster_1': 9}

        core_graph = construct_core_graph(core_neighbour_pairs)
        double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(core_graph, 10)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_challenging_paths(self):
        print(4)
        expected_segments = {'pan_cluster_A--pan_cluster_B': ['pan_cluster_A', 'pan_cluster_E', 'pan_cluster_F', 'pan_cluster_G', 'pan_cluster_B']}

        core_neighbour_pairs = {'pan_cluster_A--pan_cluster_C': 4,
                                'pan_cluster_A--pan_cluster_D': 4,
                                'pan_cluster_A--pan_cluster_E': 2,
                                'pan_cluster_B--pan_cluster_C': 5,
                                'pan_cluster_B--pan_cluster_D': 3,
                                'pan_cluster_B--pan_cluster_G': 2,
                                'pan_cluster_C--pan_cluster_D': 1,
                                'pan_cluster_E--pan_cluster_F': 2,
                                'pan_cluster_F--pan_cluster_G': 2,
                                }

        core_graph = construct_core_graph(core_neighbour_pairs)
        print(core_graph.degree)
        double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(core_graph, 10)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_challenging_paths_2(self):
        expected_segments = {'pan_cluster_A--pan_cluster_B': ['pan_cluster_A', 'pan_cluster_F', 'pan_cluster_B'],
                             'pan_cluster_B--pan_cluster_C': ['pan_cluster_B', 'pan_cluster_I', 'pan_cluster_C']}

        core_neighbour_pairs = {'pan_cluster_A--pan_cluster_D': 2,
                                'pan_cluster_A--pan_cluster_E': 1,
                                'pan_cluster_A--pan_cluster_F': 7,
                                'pan_cluster_B--pan_cluster_F': 7,
                                'pan_cluster_B--pan_cluster_I': 8,
                                'pan_cluster_B--pan_cluster_D': 1,
                                'pan_cluster_C--pan_cluster_E': 1,
                                'pan_cluster_C--pan_cluster_D': 1,
                                'pan_cluster_C--pan_cluster_I': 8,
                                'pan_cluster_D--pan_cluster_E': 1
                                }

        core_graph = construct_core_graph(core_neighbour_pairs)
        double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(core_graph, 10)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_all_challenging_paths(self):
        expected_segments = {'pan_cluster_A--pan_cluster_D': ['pan_cluster_A', 'pan_cluster_G', 'pan_cluster_F', 'pan_cluster_E', 'pan_cluster_D'],
                             'pan_cluster_B--pan_cluster_C': ['pan_cluster_B', 'pan_cluster_H', 'pan_cluster_I', 'pan_cluster_J', 'pan_cluster_C']}#,}
                             #'pan_cluster_A--pan_cluster_C': ['pan_cluster_A', 'pan_cluster_K', 'pan_cluster_C'],
                             #'pan_cluster_B--pan_cluster_D': ['pan_cluster_B', 'pan_cluster_L', 'pan_cluster_D']}

        core_neighbour_pairs = {'pan_cluster_A--pan_cluster_B': 4,
                                'pan_cluster_A--pan_cluster_K': 4,
                                'pan_cluster_A--pan_cluster_G': 2,
                                'pan_cluster_B--pan_cluster_H': 2,
                                'pan_cluster_B--pan_cluster_L': 4,
                                'pan_cluster_C--pan_cluster_J': 2,
                                'pan_cluster_C--pan_cluster_K': 4,
                                'pan_cluster_D--pan_cluster_C': 4,
                                'pan_cluster_D--pan_cluster_L': 4,
                                'pan_cluster_D--pan_cluster_E': 2,
                                'pan_cluster_E--pan_cluster_F': 2,
                                'pan_cluster_F--pan_cluster_G': 2,
                                'pan_cluster_H--pan_cluster_I': 2,
                                'pan_cluster_I--pan_cluster_J': 2,
                                'pan_cluster_K--pan_cluster_L': 1
                                }


        core_graph = construct_core_graph(core_neighbour_pairs)
        double_edge_segements, connect_dict, multi_edge_nodes = identify_segments(core_graph, 10)

        print(double_edge_segements)

        self.assertEqual(expected_segments, double_edge_segements)

class TestNoAccessorySegmentIdentifcation(unittest.TestCase):
    def test_no_accessory_genes_in_segment(self):
        expected_sub_sgments = {'pan_cluster_1~~pan_cluster_5': [['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5']],
                                 'pan_cluster_2~~pan_cluster_4': [['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']]}

        double_edge_segements = {'pan_cluster_1~~pan_cluster_5': ['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5'],
                                 'pan_cluster_2~~pan_cluster_4': ['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_6': 0, 'pan_cluster_5--pan_cluster_6': 0, 'pan_cluster_2--pan_cluster_3': 0, 'pan_cluster_3--pan_cluster_4': 0}

        sub_segment_dict = identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        print(sub_segment_dict)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

    def test_accessory_genes_in_segment_first_gene_lonely(self):
        expected_sub_sgments = {'pan_cluster_1~~pan_cluster_5': [['pan_cluster_1'], ['pan_cluster_6', 'pan_cluster_5']]}

        double_edge_segements = {'pan_cluster_1~~pan_cluster_5': ['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_6': 1, 'pan_cluster_5--pan_cluster_6': 0}

        sub_segment_dict = identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

    def test_accessory_genes_in_segment_last_gene_lonely(self):
        expected_sub_sgments = {'pan_cluster_1~~pan_cluster_5': [['pan_cluster_1', 'pan_cluster_6'], ['pan_cluster_5']],
                                'pan_cluster_2~~pan_cluster_4': [['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']]}

        double_edge_segements = {'pan_cluster_1~~pan_cluster_5': ['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5'],
                                 'pan_cluster_2~~pan_cluster_4': ['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_6': 0, 'pan_cluster_5--pan_cluster_6': 1, 'pan_cluster_2--pan_cluster_3': 0, 'pan_cluster_3--pan_cluster_4': 0}

        sub_segment_dict = identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

    def test_accessory_genes_in_segment_middle(self):
        expected_sub_sgments = {'pan_cluster_1~~pan_cluster_4': [['pan_cluster_1', 'pan_cluster_2'], ['pan_cluster_3', 'pan_cluster_4']]}

        double_edge_segements = {'pan_cluster_1~~pan_cluster_4': ['pan_cluster_1', 'pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_2': 0, 'pan_cluster_2--pan_cluster_3': 1, 'pan_cluster_3--pan_cluster_4': 0}

        sub_segment_dict = identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

if __name__ == '__main__':
    unittest.main()
