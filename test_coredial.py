import unittest
from gff_parser import get_genome_size_from_gff, segment_genome_content
from parse_gene_presence_absence import read_gene_presence_absence
from merge_dicts import merge_dicts_lists, merge_dicts_counts
import os
import glob

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
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), i=1)
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
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), i=1)
            distance_result_dict = merge_dicts_lists(distance_result_dict, core_gene_pair_distance)

        self.assertEqual(len(distance_result_dict.keys()), 10)


    def test_segmentation_accessory_gene_number(self):
        core_genes, accessory_genes = read_gene_presence_absence(
            "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/mock_gene_presence_absence.csv",
            1, 0.05, verbose=False)
        print(accessory_genes)

        pre_path = "/Users/mjespersen/OneDrive - The University of Melbourne/Phd/Parts/Recombination_hotspots/Code/Between_core_variation/data_for_unit_tests/complete/Single_contig/"
        gff_files = glob.glob(pre_path+"*.gff")

        accessory_result_dict = {}

        for file in gff_files:
            core_gene_pairs, core_gene_pair_distance, accessory_gene_content, low_freq_gene_content, master_info, _ = segment_genome_content(low_freq_genes=accessory_genes, core_genes=core_genes, input_file=os.path.join(pre_path,file), i=1)
            print(file)
            print(low_freq_gene_content)
            accessory_result_dict = merge_dicts_lists(accessory_result_dict, accessory_gene_content)
        accessory_gene_amount = []
        for i, key in enumerate(accessory_result_dict.keys()):
            if i == 10:
                self.assertListEqual([7,6,7,8,7,5,7,9,7,7,6,7,7,6,8,3,7,8,4,5], accessory_result_dict[key])


if __name__ == '__main__':
    unittest.main()
