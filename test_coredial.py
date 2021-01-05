import unittest
from gff_parser import get_genome_size_from_gff, segment_genome_content
from parse_gene_presence_absence import read_gene_presence_absence
from merge_dicts import merge_dicts_lists, merge_dicts_counts
from consesus_core_genome import characterise_rearrangements
from check_inputs import define_input_source, check_gene_data
import os
import glob
from numpy import arange, ceil


class TestPresenceAbsenceParser(unittest.TestCase):

    def test_parser_core_genes(self):
        core_genes, accessory_genes = read_gene_presence_absence("/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
                                                                 1, 0.05, verbose=False)
        keys = [key for key in core_genes.keys()]
        self.assertEqual(len(core_genes[keys[1]]), 10)

    def test_parser_low_freq_genes(self):
        core_genes, accessory_genes = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        low_frew_genes = 0
        for key in accessory_genes.keys():
            low_frew_genes += len(accessory_genes[key])
        self.assertEqual(low_frew_genes, 5)



class TestGffparser(unittest.TestCase):

    def test_get_genome_size_from_gff(self):
        genome_length = get_genome_size_from_gff("/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/Ajwa_the_Streptococcus.gff")

        self.assertEqual(genome_length, 1600)

    def test_segmentation_core_gene_number_all_complete(self):
        core_genes, accessory_genes = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        pre_path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/"
        gff_files = glob.glob(pre_path+"*.gff")

        core_result_dict = {}

        for file in gff_files:
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), i=1)
            core_result_dict = merge_dicts_counts(core_result_dict, core_gene_pairs)

        self.assertEqual(len(core_result_dict.keys()), 10)

    def test_segmentation_distance_number(self):
        core_genes, accessory_genes = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        pre_path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/"
        gff_files = glob.glob(pre_path+"*.gff")

        distance_result_dict = {}

        for file in gff_files:
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), i=1)
            distance_result_dict = merge_dicts_lists(distance_result_dict, core_gene_pair_distance)

        self.assertEqual(len(distance_result_dict.keys()), 10)


    def test_segmentation_accessory_gene_number(self):
        core_genes, accessory_genes = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)

        pre_path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/"
        gff_files = glob.glob(pre_path+"*.gff")

        accessory_result_dict = {}

        for file in gff_files:
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), i=1)
            accessory_result_dict = merge_dicts_lists(accessory_result_dict, accessory_gene_content)
        accessory_gene_amount = []
        for i, key in enumerate(accessory_result_dict.keys()):
            if i == 10:
                self.assertListEqual([7,6,7,8,7,5,7,9,7,7,6,7,7,6,8,3,7,8,4,5], accessory_result_dict[key])


# TODO - constuct test that tests load a single gff file

# TODO - constuct test that loads a mock gene_prensence_absence file
    def test_loading_gene_pres_abs(self):
        # TODO - Construct mock gene_presence_absence_file
        # TODO - Load Panaroo file
        # Load mock Panaroo gene presence absence file
        read_core, read_low_freq = read_gene_presence_absence("/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_gene_pres_abs_parser/read_mock_gene_presence_absence.csv",
                                                                  1, 0.05, verbose=False)

        # Test number of keys to see if all genomes are read in
        self.assertEqual(len(read_core), 5)

        ## Test if the expected number of core genes are found
        for i in arange(0.8, 1, 0.1):
            read_core, _ = read_gene_presence_absence(
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
            _, read_low_freq = read_gene_presence_absence(
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
        source = define_input_source("/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Test_source_identification/Roray_gene_presence_absence")
        self.assertEqual(source, "Roary")

        # Identify Panaroo file
        source = define_input_source('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/Test_source_identification/Panaroo_gene_presence_absence')
        self.assertEqual(source, "Panaroo")

    def test_locate_gene_data(self):
        check_return = check_gene_data('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome/same_genome_pan_split')

        self.assertEqual(check_return, True)

        with self.assertRaises(FileNotFoundError):
            check_gene_data('/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/test_pan_genome')

if __name__ == '__main__':
    unittest.main()
