'''
Unit tests for Corekaburra.

Usage: python -m unittest -v Corekaburra_test
'''

# import
import unittest
import os
import json
from shutil import copyfile
import logging
# pylint: disable=no-name-in-module

# import Corekaburra functions
from Corekaburra import exit_with_error
from Corekaburra import read_complete_genome_file
from Corekaburra import check_inputs
from Corekaburra import parse_gene_presence_absence



# move to folder with mock files. First try Github structure, then try pulled repository structure
try:
    os.chdir('/Corekaburra/unit_tests/unit_test_data/')
except FileNotFoundError:
    os.chdir('unit_test_data/')


class TestExitWithError(unittest.TestCase):
    """ Test for the function carrying out a nice exit """
    def test_exit_w_tmp_folder_deletion(self):
        ''' Test the exit function is able to remove the temporary folder '''

        # copy the placeholder tmp folder to replace it afterwards
        tmp_folder = 'TestExitWithError/tmp_folder'
        tmp_folder_copy = 'TestExitWithError/tmp_folder_copy'
        os.mkdir(tmp_folder_copy)

        tmp_files = os.listdir(tmp_folder)
        for file in tmp_files:
            copyfile(os.path.join(tmp_folder, file), os.path.join(tmp_folder_copy, file))

        with self.assertRaises(SystemExit):
            exit_with_error.exit_with_error(exit_status=2, message='test msg', tmp_folder=tmp_folder)

        os.rename(tmp_folder_copy, tmp_folder)


class TestParsingCompleteGenomes(unittest.TestCase):
    """ Test for the passing of input file containing names of complete genome and checking their presence in the pan-genome """
    def test_all_files_found(self):
        gff_files = ['/path/to/complete_genome_1.gff',
                     '/path/complete_genome_2.gff.gz',
                     'complete_genome_3.gff.gz',
                     'complete_genome_4.gff',
                     'dummy_index_1',
                     'dummy_index_2']

        complete_genome_file = 'TestParsingCompleteGenomes/complete_genomes_file.txt'

        expected_return = ['complete_genome_1',
                           'complete_genome_2',
                           'complete_genome_3',
                           'complete_genome_4']

        return_object = read_complete_genome_file.parse_complete_genome_file(complete_genome_file, gff_files)

        self.assertEqual(return_object, expected_return)

    def test_correct_one_files_not_found(self):
        gff_files = ['/path/complete_genome_2.gff.gz',
                     'complete_genome_3.gff.gz',
                     'complete_genome_4.gff',
                     'dummy_index_1',
                     'dummy_index_2']

        complete_genome_file = 'TestParsingCompleteGenomes/complete_genomes_file.txt'

        with self.assertRaises(SystemExit):
            read_complete_genome_file.parse_complete_genome_file(complete_genome_file,
                                                                 gff_files)


class TestPangenomeSourceProgram(unittest.TestCase):
    """ Test of the function that determines the program from which the pan-genome originated """
    def test_roary_input(self):
        input_folder_path = 'TestPangenomeSourceProgram/Mock_roary'

        return_program, return_path = check_inputs.define_pangenome_program(input_folder_path)

        self.assertEqual("Roary", return_program)
        self.assertEqual(input_folder_path + '/gene_presence_absence.csv', return_path)

    def test_panaroo_input(self):
        input_folder_path = 'TestPangenomeSourceProgram/Mock_panaroo'

        return_program, return_path = check_inputs.define_pangenome_program(input_folder_path)

        self.assertEqual("Panaroo", return_program)
        self.assertEqual(input_folder_path + '/gene_presence_absence_roary.csv', return_path)

    # def test_pirate_input(self): TODO - Make Corekaburra take Pirate input!
    #     pass
    #     input_folder_path = 'TestPangenomeSourceProgram/Mock_pirate'
    #
    #     return_program, return_path = check_inputs.define_pangenome_program(input_folder_path)
    #
    #     self.assertEqual("Pirate", return_program)

    def test_unknown_input(self):
        input_folder_path = 'TestPangenomeSourceProgram/Mock_unknwon'

        with self.assertRaises(SystemExit):
            check_inputs.define_pangenome_program(input_folder_path)


class TestPresenceOfGenedataFile(unittest.TestCase):
    """ Test the function that ensures the presence of the Gene_data.csv file produced by Panaroo """
    def test_Genedata_File_present(self):
        input_folder_path = 'TestPresenceOfGenedataFile/present'
        return_path = check_inputs.check_gene_data(input_folder_path)

        self.assertEqual(return_path, input_folder_path +'/gene_data.csv')

    def test_Genedata_File_absent(self):
        input_folder_path = 'TestPresenceOfGenedataFile/absent'

        with self.assertRaises(SystemExit):
            check_inputs.check_gene_data(input_folder_path)


class TestPresenceOfGffsInPresAbsFile(unittest.TestCase):
    """ Test the function that ensures all gffs given as input are included in the pan-genome provided """
    # Test pairing of all files in pan genome
    def test_input_gff_pres_abs_pairing_all_gffs(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['Silas_the_Salmonella', 'Christina_the_Streptococcus', 'Ajwa_the_Shigella']

        return_bool = check_inputs.check_gff_in_pan(input_file_list, input_pres_abs)

        self.assertEqual(return_bool, True)

    # Test pairing of some files in pan genome - Warning
    def test_input_gff_pres_abs_pairing_some(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff']

        with self.assertWarns(Warning):
            return_bool = check_inputs.check_gff_in_pan(input_file_list, input_pres_abs)

        self.assertEqual(return_bool, True)

    # Test when given a file not in pan genome among others that are in the pan genome
    def test_input_gff_pres_abs_file_not_in_pan(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['not_found.gff', 'Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff']

        with self.assertRaises(SystemExit):
            check_inputs.check_gff_in_pan(input_file_list, input_pres_abs)

    def test_input_gff_pres_abs_some_file_not_in_pan(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['not_found.gff', 'also_not_found.gff', 'definitely_not_found.gff']

        with self.assertRaises(SystemExit):
            check_inputs.check_gff_in_pan(input_file_list, input_pres_abs)


class TestAddingGeneToDict(unittest.TestCase):
    """
    Tests of the function that adds a gene to the dict used to holds that class (core, accessory of low-frequency)
    """
    def test_adding_gene(self):
        main_dict = {'Test_genome': {}}
        gene = 'Test_gene_from_genome'
        pan_gene_name = 'Test_pan_gene'
        genome = 'Test_genome'

        expected_return = {'Test_genome': {'Test_gene_from_genome': 'Test_pan_gene'}}

        return_dict = parse_gene_presence_absence.add_gene_to_dict(main_dict, gene, pan_gene_name, genome)

        self.assertEqual(expected_return, return_dict)

    def test_adding_additional_gene(self):
        main_dict = {'Test_genome': {'Test_gene_from_genome': 'Test_pan_gene'}}
        gene = 'Test_gene_from_genome_2'
        pan_gene_name = 'Test_pan_gene_2'
        genome = 'Test_genome'

        expected_return = {'Test_genome': {'Test_gene_from_genome': 'Test_pan_gene',
                                           'Test_gene_from_genome_2': 'Test_pan_gene_2'}}

        return_dict = parse_gene_presence_absence.add_gene_to_dict(main_dict, gene, pan_gene_name, genome)

        self.assertEqual(expected_return, return_dict)


class TestCheckingFragmentedGenes(unittest.TestCase):
    """
    Test of the function that examines the placement of a potential core gene's placement, if it is fragmented in at least one genome.
    """

    def tearDown(self):
        """ Class to remove created database files of gff files in tmp-folder"""
        for file in os.listdir('test_tmp_folder'):
            if "_db" in file:
                db_path = os.path.join('test_tmp_folder', file)
                os.remove(db_path)

    def test_fragmented_gene_true(self):
        """ Gene is fragmented but found next to each other with nothing in between """
        fragments_in_line = ['Silas_the_Salmonella_tag-1-2.1;Silas_the_Salmonella_tag-1-2.2']
        input_gffs =['TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella.gff',
                      'TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'test_tmp_folder'

        expected_return = [True]

        return_bool = parse_gene_presence_absence.check_fragmented_gene(fragments_in_line, input_gffs, tmp_folder_path)

        self.assertEqual(expected_return, return_bool)

    def test_fragmented_gene_fasle(self):
        """ Gene is fragmented but found next to each other with another gene in between """
        fragments_in_line = ['Silas_the_Salmonella_tag-1-5.1;Silas_the_Salmonella_tag-1-5.2']
        input_gffs = ['TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella.gff',
                      'TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'test_tmp_folder'

        expected_return = [False]

        return_bool = parse_gene_presence_absence.check_fragmented_gene(fragments_in_line, input_gffs, tmp_folder_path)

        self.assertEqual(expected_return, return_bool)

    def test_fragmented_gene_mutiple_genes_fasle(self):
        """ Two genes fragmented with one having nothing and the other having something in between fragments """
        fragments_in_line = ['Silas_the_Salmonella_tag-1-2.1;Silas_the_Salmonella_tag-1-2.2', 'Silas_the_Salmonella_tag-1-5.1;Silas_the_Salmonella_tag-1-5.2']
        input_gffs = ['TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella.gff',
                      'TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'test_tmp_folder'

        expected_return = [True, False]

        return_bool = parse_gene_presence_absence.check_fragmented_gene(fragments_in_line, input_gffs, tmp_folder_path)

        self.assertEqual(expected_return, return_bool)


class TestParsingGenePresenceAbsenceFile(unittest.TestCase):
    """
    Tests for the function that passes the gene presence absence table from pan-genome program
    """
    def test_parsing_w_100_presence(self):
        file_name = 'TestParsingGenePresenceAbsenceFile/gene_presence_absence_roary.csv'
        core_gene_presence = 1
        low_freq_gene = 0.1
        source_program = 'Panaroo'
        input_gffs = ['TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella.gff',
                      'TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'TestParsingGenePresenceAbsenceFile/'

        expected_core_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-1': "A",
                                                            'Silas_the_Salmonella_tag-1-2.1': "B",
                                                            'Silas_the_Salmonella_tag-1-2.2': "B"},
                                   'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-1': "A",
                                                                   'Christina_the_Streptococcus_tag-2-2': "B"},
                                   'Ajwa_the_Shigella': {'Ajwa_the_Shigella_tag-3-1': "A",
                                                         'Ajwa_the_Shigella_tag-3-2': "B"},
                                   'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-1': "A",
                                                           'Ajwa_the_Legionella_tag-4-2': "B"},
                                   'Cari_the_Listeria': {'Cari_the_Listeria_tag-5-1': "A",
                                                         'Cari_the_Listeria_tag-5-2': "B"},
                                   'Aman_the_Streptococcus': {'Aman_the_Streptococcus_tag-6-1': "A",
                                                              'Aman_the_Streptococcus_tag-6-2': "B"},
                                   'Zion_the_Streptococcus': {'Zion_the_Streptococcus_tag-7-1': "A",
                                                              'Zion_the_Streptococcus_tag-7-2': "B"},
                                   'Dina_the_Shigella': {'Dina_the_Shigella_tag-8-1': "A",
                                                         'Dina_the_Shigella_tag-8-2': "B"},
                                   'Silas_the_Legionella': {'Silas_the_Legionella_tag-9-1': "A",
                                                            'Silas_the_Legionella_tag-9-2': "B"},
                                   'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-1': "A",
                                                          'Lilly_the_Shigella_tag-10-2': "B"}}
        expected_low_freq_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-7': "G"},
                                       'Christina_the_Streptococcus': {},
                                       'Ajwa_the_Shigella': {},
                                       'Ajwa_the_Legionella': {},
                                       'Cari_the_Listeria': {},
                                       'Aman_the_Streptococcus': {},
                                       'Zion_the_Streptococcus': {},
                                       'Dina_the_Shigella': {},
                                       'Silas_the_Legionella': {},
                                       'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-6': "F"}}
        expected_acc_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-3': 'C',
                                                           'Silas_the_Salmonella_tag-1-4.1': 'D',
                                                           'Silas_the_Salmonella_tag-1-4.2': 'D',
                                                           'Silas_the_Salmonella_tag-1-5.1': 'E',
                                                           'Silas_the_Salmonella_tag-1-5.2': 'E'},
                                  'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-3': "C",
                                                                  'Christina_the_Streptococcus_tag-2-4': "D",
                                                                  'Christina_the_Streptococcus_tag-2-5': "E"},
                                  'Ajwa_the_Shigella': {"Ajwa_the_Shigella_tag-3-3": "C",
                                                        "Ajwa_the_Shigella_tag-3-4": "D",
                                                        "Ajwa_the_Shigella_tag-3-5": "E"},
                                  'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-3': "C",
                                                          'Ajwa_the_Legionella_tag-4-4': "D",
                                                          'Ajwa_the_Legionella_tag-4-5': "E"},
                                  'Cari_the_Listeria': {"Cari_the_Listeria_tag-5-3": "C",
                                                        "Cari_the_Listeria_tag-5-4": "D",
                                                        "Cari_the_Listeria_tag-5-5": "E"},
                                  'Aman_the_Streptococcus': {"Aman_the_Streptococcus_tag-6-3": "C",
                                                             "Aman_the_Streptococcus_tag-6-4": "D",
                                                             "Aman_the_Streptococcus_tag-6-5": "E"},
                                  'Zion_the_Streptococcus': {"Zion_the_Streptococcus_tag-7-3": "C",
                                                             "Zion_the_Streptococcus_tag-7-4": "D",
                                                             "Zion_the_Streptococcus_tag-7-5": "E"},
                                  'Dina_the_Shigella': {"Dina_the_Shigella_tag-8-3": "C",
                                                        "Dina_the_Shigella_tag-8-4": "D",
                                                        "Dina_the_Shigella_tag-8-5": "E"},
                                  'Silas_the_Legionella': {"Silas_the_Legionella_tag-9-3": "C",
                                                           "Silas_the_Legionella_tag-9-4": "D",
                                                           "Silas_the_Legionella_tag-9-5": "E"},
                                  'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-5': "E"}}
        expected_annotation_dict = {}  # None - Should be done and holds refounds! - TODO Make test for this

        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict, annotation_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
            file_name, core_gene_presence,
            low_freq_gene, source_program,
            input_gffs, tmp_folder_path)

        self.assertEqual(expected_core_gene_dict, core_gene_dict)
        self.assertEqual(expected_low_freq_gene_dict, low_freq_gene_dict)
        self.assertEqual(expected_acc_gene_dict, acc_gene_dict)

    def test_parsing_w_100_presence_roary(self):
        file_name = 'TestParsingGenePresenceAbsenceFile/gene_presence_absence.csv'
        core_gene_presence = 1
        low_freq_gene = 0.1
        source_program = 'Roary'
        input_gffs = ['TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella.gff',
                      'TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'TestParsingGenePresenceAbsenceFile/'



        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict, annotation_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
                file_name, core_gene_presence,
                low_freq_gene, source_program,
                input_gffs, tmp_folder_path)

        expected_core_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-1': "A",
                                                            'Silas_the_Salmonella_tag-1-2.1': "B",
                                                            'Silas_the_Salmonella_tag-1-2.2': "B"},
                                   'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-1': "A",
                                                                   'Christina_the_Streptococcus_tag-2-2': "B"},
                                   'Ajwa_the_Shigella': {'Ajwa_the_Shigella_tag-3-1': "A",
                                                         'Ajwa_the_Shigella_tag-3-2': "B"},
                                   'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-1': "A",
                                                           'Ajwa_the_Legionella_tag-4-2': "B"},
                                   'Cari_the_Listeria': {'Cari_the_Listeria_tag-5-1': "A",
                                                         'Cari_the_Listeria_tag-5-2': "B"},
                                   'Aman_the_Streptococcus': {'Aman_the_Streptococcus_tag-6-1': "A",
                                                              'Aman_the_Streptococcus_tag-6-2': "B"},
                                   'Zion_the_Streptococcus': {'Zion_the_Streptococcus_tag-7-1': "A",
                                                              'Zion_the_Streptococcus_tag-7-2': "B"},
                                   'Dina_the_Shigella': {'Dina_the_Shigella_tag-8-1': "A",
                                                         'Dina_the_Shigella_tag-8-2': "B"},
                                   'Silas_the_Legionella': {'Silas_the_Legionella_tag-9-1': "A",
                                                            'Silas_the_Legionella_tag-9-2': "B"},
                                   'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-1': "A",
                                                          'Lilly_the_Shigella_tag-10-2': "B"}}
        expected_low_freq_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-7': "G"},
                                       'Christina_the_Streptococcus': {},
                                       'Ajwa_the_Shigella': {},
                                       'Ajwa_the_Legionella': {},
                                       'Cari_the_Listeria': {},
                                       'Aman_the_Streptococcus': {},
                                       'Zion_the_Streptococcus': {},
                                       'Dina_the_Shigella': {},
                                       'Silas_the_Legionella': {},
                                       'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-6': "F"}}
        expected_acc_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-3': 'C',
                                                           'Silas_the_Salmonella_tag-1-4.1': 'D',
                                                           'Silas_the_Salmonella_tag-1-4.2': 'D',
                                                           'Silas_the_Salmonella_tag-1-5.1': 'E',
                                                           'Silas_the_Salmonella_tag-1-5.2': 'E'},
                                  'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-3': "C",
                                                                  'Christina_the_Streptococcus_tag-2-4': "D",
                                                                  'Christina_the_Streptococcus_tag-2-5': "E"},
                                  'Ajwa_the_Shigella': {"Ajwa_the_Shigella_tag-3-3": "C",
                                                        "Ajwa_the_Shigella_tag-3-4": "D",
                                                        "Ajwa_the_Shigella_tag-3-5": "E"},
                                  'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-3': "C",
                                                          'Ajwa_the_Legionella_tag-4-4': "D",
                                                          'Ajwa_the_Legionella_tag-4-5': "E"},
                                  'Cari_the_Listeria': {"Cari_the_Listeria_tag-5-3": "C",
                                                        "Cari_the_Listeria_tag-5-4": "D",
                                                        "Cari_the_Listeria_tag-5-5": "E"},
                                  'Aman_the_Streptococcus': {"Aman_the_Streptococcus_tag-6-3": "C",
                                                             "Aman_the_Streptococcus_tag-6-4": "D",
                                                             "Aman_the_Streptococcus_tag-6-5": "E"},
                                  'Zion_the_Streptococcus': {"Zion_the_Streptococcus_tag-7-3": "C",
                                                             "Zion_the_Streptococcus_tag-7-4": "D",
                                                             "Zion_the_Streptococcus_tag-7-5": "E"},
                                  'Dina_the_Shigella': {"Dina_the_Shigella_tag-8-3": "C",
                                                        "Dina_the_Shigella_tag-8-4": "D",
                                                        "Dina_the_Shigella_tag-8-5": "E"},
                                  'Silas_the_Legionella': {"Silas_the_Legionella_tag-9-3": "C",
                                                           "Silas_the_Legionella_tag-9-4": "D",
                                                           "Silas_the_Legionella_tag-9-5": "E"},
                                  'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-5': "E"}}

        self.assertEqual(expected_core_gene_dict, core_gene_dict)
        self.assertEqual(expected_low_freq_gene_dict, low_freq_gene_dict)
        self.assertEqual(expected_acc_gene_dict, acc_gene_dict)

    def test_parsing_w_90_presence(self):
        file_name = 'TestParsingGenePresenceAbsenceFile/gene_presence_absence_roary.csv'
        core_gene_presence = 0.9
        low_freq_gene = 0.1
        source_program = 'Panaroo'
        input_gffs = ['TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella.gff',
                      'TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'TestParsingGenePresenceAbsenceFile/'

        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict, annotation_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
                file_name, core_gene_presence,
                low_freq_gene, source_program,
                input_gffs, tmp_folder_path)

        expected_core_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-1': "A",
                                                            'Silas_the_Salmonella_tag-1-2.1': "B",
                                                            'Silas_the_Salmonella_tag-1-2.2': "B",
                                                            'Silas_the_Salmonella_tag-1-3': 'C',
                                                            'Silas_the_Salmonella_tag-1-4.1': 'D',
                                                            'Silas_the_Salmonella_tag-1-4.2': 'D',},
                                   'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-1': "A",
                                                                   'Christina_the_Streptococcus_tag-2-2': "B",
                                                                   'Christina_the_Streptococcus_tag-2-3': "C",
                                                                   'Christina_the_Streptococcus_tag-2-4': "D"},
                                   'Ajwa_the_Shigella': {'Ajwa_the_Shigella_tag-3-1': "A",
                                                         'Ajwa_the_Shigella_tag-3-2': "B",
                                                         "Ajwa_the_Shigella_tag-3-3": "C",
                                                         "Ajwa_the_Shigella_tag-3-4": "D"},
                                   'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-1': "A",
                                                           'Ajwa_the_Legionella_tag-4-2': "B",
                                                           'Ajwa_the_Legionella_tag-4-3': "C",
                                                           'Ajwa_the_Legionella_tag-4-4': "D"},
                                   'Cari_the_Listeria': {"Cari_the_Listeria_tag-5-3": "C",
                                                         "Cari_the_Listeria_tag-5-4": "D",
                                                         'Cari_the_Listeria_tag-5-1': "A",
                                                         'Cari_the_Listeria_tag-5-2': "B"},
                                   'Aman_the_Streptococcus': {'Aman_the_Streptococcus_tag-6-1': "A",
                                                              'Aman_the_Streptococcus_tag-6-2': "B",
                                                              "Aman_the_Streptococcus_tag-6-3": "C",
                                                              "Aman_the_Streptococcus_tag-6-4": "D"},
                                   'Zion_the_Streptococcus': {"Zion_the_Streptococcus_tag-7-3": "C",
                                                              "Zion_the_Streptococcus_tag-7-4": "D",
                                                              'Zion_the_Streptococcus_tag-7-1': "A",
                                                              'Zion_the_Streptococcus_tag-7-2': "B"},
                                   'Dina_the_Shigella': {"Dina_the_Shigella_tag-8-3": "C",
                                                         "Dina_the_Shigella_tag-8-4": "D",
                                                         'Dina_the_Shigella_tag-8-1': "A",
                                                         'Dina_the_Shigella_tag-8-2': "B"},
                                   'Silas_the_Legionella': {"Silas_the_Legionella_tag-9-3": "C",
                                                            "Silas_the_Legionella_tag-9-4": "D",
                                                            'Silas_the_Legionella_tag-9-1': "A",
                                                            'Silas_the_Legionella_tag-9-2': "B"},
                                   'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-1': "A",
                                                          'Lilly_the_Shigella_tag-10-2': "B"}}
        expected_low_freq_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-7': "G"},
                                       'Christina_the_Streptococcus': {},
                                       'Ajwa_the_Shigella': {},
                                       'Ajwa_the_Legionella': {},
                                       'Cari_the_Listeria': {},
                                       'Aman_the_Streptococcus': {},
                                       'Zion_the_Streptococcus': {},
                                       'Dina_the_Shigella': {},
                                       'Silas_the_Legionella': {},
                                       'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-6': "F"}}
        expected_acc_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-5.1': 'E',
                                                           'Silas_the_Salmonella_tag-1-5.2': 'E'},
                                  'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-5': "E"},
                                  'Ajwa_the_Shigella': {"Ajwa_the_Shigella_tag-3-5": "E"},
                                  'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-5': "E"},
                                  'Cari_the_Listeria': {"Cari_the_Listeria_tag-5-5": "E"},
                                  'Aman_the_Streptococcus': {"Aman_the_Streptococcus_tag-6-5": "E"},
                                  'Zion_the_Streptococcus': {"Zion_the_Streptococcus_tag-7-5": "E"},
                                  'Dina_the_Shigella': {"Dina_the_Shigella_tag-8-5": "E"},
                                  'Silas_the_Legionella': {"Silas_the_Legionella_tag-9-5": "E"},
                                  'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-5': "E"}}

        self.assertEqual(expected_core_gene_dict, core_gene_dict)
        self.assertEqual(expected_low_freq_gene_dict, low_freq_gene_dict)
        self.assertEqual(expected_acc_gene_dict, acc_gene_dict)

    def test_parsing_w_90_presence_roary(self):
        file_name = 'TestParsingGenePresenceAbsenceFile/gene_presence_absence.csv'
        core_gene_presence = 0.90
        low_freq_gene = 0.1
        source_program = 'Roary'
        input_gffs = ['TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella.gff',
                      'TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'TestParsingGenePresenceAbsenceFile/'

        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict, annotation_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
                file_name, core_gene_presence,
                low_freq_gene, source_program,
                input_gffs, tmp_folder_path)

        expected_core_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-1': "A",
                                                            'Silas_the_Salmonella_tag-1-2.1': "B",
                                                            'Silas_the_Salmonella_tag-1-2.2': "B",
                                                            'Silas_the_Salmonella_tag-1-3': 'C',
                                                            'Silas_the_Salmonella_tag-1-4.1': 'D',
                                                            'Silas_the_Salmonella_tag-1-4.2': 'D', },
                                   'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-1': "A",
                                                                   'Christina_the_Streptococcus_tag-2-2': "B",
                                                                   'Christina_the_Streptococcus_tag-2-3': "C",
                                                                   'Christina_the_Streptococcus_tag-2-4': "D"},
                                   'Ajwa_the_Shigella': {'Ajwa_the_Shigella_tag-3-1': "A",
                                                         'Ajwa_the_Shigella_tag-3-2': "B",
                                                         "Ajwa_the_Shigella_tag-3-3": "C",
                                                         "Ajwa_the_Shigella_tag-3-4": "D"},
                                   'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-1': "A",
                                                           'Ajwa_the_Legionella_tag-4-2': "B",
                                                           'Ajwa_the_Legionella_tag-4-3': "C",
                                                           'Ajwa_the_Legionella_tag-4-4': "D"},
                                   'Cari_the_Listeria': {"Cari_the_Listeria_tag-5-3": "C",
                                                         "Cari_the_Listeria_tag-5-4": "D",
                                                         'Cari_the_Listeria_tag-5-1': "A",
                                                         'Cari_the_Listeria_tag-5-2': "B"},
                                   'Aman_the_Streptococcus': {'Aman_the_Streptococcus_tag-6-1': "A",
                                                              'Aman_the_Streptococcus_tag-6-2': "B",
                                                              "Aman_the_Streptococcus_tag-6-3": "C",
                                                              "Aman_the_Streptococcus_tag-6-4": "D"},
                                   'Zion_the_Streptococcus': {"Zion_the_Streptococcus_tag-7-3": "C",
                                                              "Zion_the_Streptococcus_tag-7-4": "D",
                                                              'Zion_the_Streptococcus_tag-7-1': "A",
                                                              'Zion_the_Streptococcus_tag-7-2': "B"},
                                   'Dina_the_Shigella': {"Dina_the_Shigella_tag-8-3": "C",
                                                         "Dina_the_Shigella_tag-8-4": "D",
                                                         'Dina_the_Shigella_tag-8-1': "A",
                                                         'Dina_the_Shigella_tag-8-2': "B"},
                                   'Silas_the_Legionella': {"Silas_the_Legionella_tag-9-3": "C",
                                                            "Silas_the_Legionella_tag-9-4": "D",
                                                            'Silas_the_Legionella_tag-9-1': "A",
                                                            'Silas_the_Legionella_tag-9-2': "B"},
                                   'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-1': "A",
                                                          'Lilly_the_Shigella_tag-10-2': "B"}}
        expected_low_freq_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-7': "G"},
                                       'Christina_the_Streptococcus': {},
                                       'Ajwa_the_Shigella': {},
                                       'Ajwa_the_Legionella': {},
                                       'Cari_the_Listeria': {},
                                       'Aman_the_Streptococcus': {},
                                       'Zion_the_Streptococcus': {},
                                       'Dina_the_Shigella': {},
                                       'Silas_the_Legionella': {},
                                       'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-6': "F"}}
        expected_acc_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-5.1': 'E',
                                                           'Silas_the_Salmonella_tag-1-5.2': 'E'},
                                  'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-5': "E"},
                                  'Ajwa_the_Shigella': {"Ajwa_the_Shigella_tag-3-5": "E"},
                                  'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-5': "E"},
                                  'Cari_the_Listeria': {"Cari_the_Listeria_tag-5-5": "E"},
                                  'Aman_the_Streptococcus': {"Aman_the_Streptococcus_tag-6-5": "E"},
                                  'Zion_the_Streptococcus': {"Zion_the_Streptococcus_tag-7-5": "E"},
                                  'Dina_the_Shigella': {"Dina_the_Shigella_tag-8-5": "E"},
                                  'Silas_the_Legionella': {"Silas_the_Legionella_tag-9-5": "E"},
                                  'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-5': "E"}}

        self.assertEqual(expected_core_gene_dict, core_gene_dict)
        self.assertEqual(expected_low_freq_gene_dict, low_freq_gene_dict)
        self.assertEqual(expected_acc_gene_dict, acc_gene_dict)


if __name__ == '__main__':
    unittest.main()
