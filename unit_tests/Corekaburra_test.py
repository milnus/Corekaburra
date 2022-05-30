'''
Unit tests for Corekaburra.

Usage: python -m unittest -v Corekaburra_test
'''

# import
import unittest
import os
from shutil import copyfile
import logging
from networkx import number_connected_components, connected_components
# pylint: disable=no-name-in-module

# import Corekaburra functions
from Corekaburra import exit_with_error
from Corekaburra import read_complete_genome_file
from Corekaburra import check_inputs
from Corekaburra import parse_gene_presence_absence
from Corekaburra import gff_parser
from Corekaburra import merge_dicts
from Corekaburra import consesus_core_genome
from Corekaburra import summary_table
from Corekaburra import output_writer_functions
from Corekaburra import correct_gffs

# move to folder with mock files. First try Github structure, then try pulled repository structure
try:
    os.chdir('/Corekaburra/unit_tests/unit_test_data/')
except FileNotFoundError:
    os.chdir('unit_test_data/')


class TestExitWithError(unittest.TestCase):
    """ Test for the function carrying out a nice exit """
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

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
            exit_with_error.exit_with_error(exit_status=2, message='test msg', logger=self.logger, tmp_folder=tmp_folder)

        os.rename(tmp_folder_copy, tmp_folder)


class TestCutOffViolations(unittest.TestCase):
    """ Test for the function that examines the cutoffs given for core and low-frequency genes"""
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_low_below_range(self):
        with self.assertRaises(SystemExit):
            check_inputs.check_cutoffs(-0.1, 1, self.logger)

    def test_core_above_range(self):
        with self.assertRaises(SystemExit):
            check_inputs.check_cutoffs(0.05, 1.1, self.logger)

    def test_low_larger_than_core(self):
        with self.assertRaises(SystemExit):
            check_inputs.check_cutoffs(0.6, 0.4, self.logger)


class TestParsingCompleteGenomes(unittest.TestCase):
    """ Test for the passing of input file containing names of complete genome and checking their presence in the pan-genome """
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_all_files_found(self):
        gff_files = ['/path/to/complete_genome_1.gff',
                     '/path/complete_genome_2.gff.gz',
                     'complete_genome_3.gff.gz',
                     'complete_genome_4.gff',
                     'complete_genome_5',
                     'dummy_index_1',
                     'dummy_index_2']

        complete_genome_file = 'TestParsingCompleteGenomes/complete_genomes_file.txt'

        expected_return = ['complete_genome_1',
                           'complete_genome_2',
                           'complete_genome_3',
                           'complete_genome_4',
                           'complete_genome_5']

        return_object = read_complete_genome_file.parse_complete_genome_file(complete_genome_file, gff_files, self.logger)

        self.assertEqual(return_object, expected_return)

    def test_correct_one_files_not_found(self):
        gff_files = ['/path/complete_genome_2.gff.gz',
                     'complete_genome_3.gff.gz',
                     'complete_genome_4.gff',
                     'complete_genome_5',
                     'dummy_index_1',
                     'dummy_index_2']

        complete_genome_file = 'TestParsingCompleteGenomes/complete_genomes_file.txt'

        with self.assertRaises(SystemExit):
            read_complete_genome_file.parse_complete_genome_file(complete_genome_file,
                                                                 gff_files, self.logger)


class TestPangenomeSourceProgram(unittest.TestCase):
    """ Test of the function that determines the program from which the pan-genome originated """
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_roary_input(self):
        input_folder_path = 'TestPangenomeSourceProgram/Mock_roary'

        return_program, return_path = check_inputs.define_pangenome_program(input_folder_path, self.logger)

        self.assertEqual("Roary", return_program)
        self.assertEqual(input_folder_path + '/gene_presence_absence.csv', return_path)

    def test_panaroo_input(self):
        input_folder_path = 'TestPangenomeSourceProgram/Mock_panaroo'

        return_program, return_path = check_inputs.define_pangenome_program(input_folder_path, self.logger)

        self.assertEqual("Panaroo", return_program)
        self.assertEqual(input_folder_path + '/gene_presence_absence_roary.csv', return_path)

    def test_minimal_panaroo_input(self):
        input_folder_path = 'TestPangenomeSourceProgram/Mock_minimal_panaroo'

        return_program, return_path = check_inputs.define_pangenome_program(input_folder_path, self.logger)

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
            check_inputs.define_pangenome_program(input_folder_path, self.logger)


class TestPresenceOfGenedataFile(unittest.TestCase):
    """ Test the function that ensures the presence of the Gene_data.csv file produced by Panaroo """
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_Genedata_File_present(self):
        input_folder_path = 'TestPresenceOfGenedataFile/present'
        return_path = check_inputs.check_gene_data(input_folder_path, self.logger)

        self.assertEqual(return_path, input_folder_path +'/gene_data.csv')

    def test_Genedata_File_absent(self):
        input_folder_path = 'TestPresenceOfGenedataFile/absent'

        with self.assertRaises(SystemExit):
            check_inputs.check_gene_data(input_folder_path, self.logger)


class TestPresenceOfGffsInPresAbsFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    """ Test the function that ensures all gffs given as input are included in the pan-genome provided """
    # Test pairing of all files in pan genome
    def test_input_gff_pres_abs_pairing_all_gffs(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['Silas_the_Salmonella', 'Christina_the_Streptococcus', 'Ajwa_the_Shigella']

        return_bool = check_inputs.check_gff_in_pan(input_file_list, input_pres_abs, self.logger)

        self.assertEqual(return_bool, True)

    # Test pairing of some files in pan genome - Warning
    def test_input_gff_pres_abs_pairing_some(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff']

        return_bool = check_inputs.check_gff_in_pan(input_file_list, input_pres_abs, self.logger)

        self.assertEqual(return_bool, True)

    # Test when given a file not in pan genome among others that are in the pan genome
    def test_input_gff_pres_abs_file_not_in_pan(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['not_found.gff', 'Silas_the_Salmonella.gff', 'Christina_the_Streptococcus.gff']

        with self.assertRaises(SystemExit):
            check_inputs.check_gff_in_pan(input_file_list, input_pres_abs, self.logger)

    def test_input_gff_pres_abs_some_file_not_in_pan(self):
        input_pres_abs = 'TestPresenceOfGffsInPresAbsFile/gene_presence_absence_roary.csv'
        input_file_list = ['not_found.gff', 'also_not_found.gff', 'definitely_not_found.gff']

        with self.assertRaises(SystemExit):
            check_inputs.check_gff_in_pan(input_file_list, input_pres_abs, self.logger)


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
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def tearDown(self):
        """ Class to remove created database files of gff files in tmp-folder"""
        for file in os.listdir('test_tmp_folder'):
            if "_db" in file:
                db_path = os.path.join('test_tmp_folder', file)
                os.remove(db_path)

    def test_fragmented_gene_true(self):
        """ Gene is fragmented but found next to each other with nothing in between """
        fragments_info = [['Silas_the_Salmonella_tag-1-2.1;Silas_the_Salmonella_tag-1-2.2', 'Silas_the_Salmonella']]
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
        gene_data_file = {}
        corrected_dir = ''

        expected_return = [True]

        return_bool = parse_gene_presence_absence.check_fragmented_gene(fragments_info, input_gffs, tmp_folder_path,
                                                                        gene_data_file, corrected_dir, self.logger)

        self.assertEqual(expected_return, return_bool)

    def test_fragmented_gene_fasle(self):
        """ Gene is fragmented but found next to each other with another gene in between """
        fragments_info = [['Silas_the_Salmonella_tag-1-5.1;Silas_the_Salmonella_tag-1-5.2', 'Silas_the_Salmonella']]
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
        gene_data_file = {}
        corrected_dir = ''

        expected_return = [False]

        return_bool = parse_gene_presence_absence.check_fragmented_gene(fragments_info, input_gffs, tmp_folder_path,
                                                                        gene_data_file, corrected_dir, self.logger)

        self.assertEqual(expected_return, return_bool)

    def test_fragmented_gene_mutiple_genes_fasle(self):
        """ Two genes fragmented with one having nothing and the other having something in between fragments """
        fragment_info = [['Silas_the_Salmonella_tag-1-2.1;Silas_the_Salmonella_tag-1-2.2', 'Silas_the_Salmonella'],
                         ['Silas_the_Salmonella_tag-1-5.1;Silas_the_Salmonella_tag-1-5.2', 'Silas_the_Salmonella']]
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
        gene_data_file = {}
        corrected_dir = ''

        expected_return = [True, False]

        return_bool = parse_gene_presence_absence.check_fragmented_gene(fragment_info, input_gffs, tmp_folder_path,
                                                                        gene_data_file, corrected_dir, self.logger)

        self.assertEqual(expected_return, return_bool)

    def test_fragments_on_separate_contigs(self):
        """ One gene fragmented with parts on separate contigs """
        fragments_info = [['Silas_the_Salmonella_tag-1-2.1;Silas_the_Salmonella_tag-1-2.2', 'Silas_the_Salmonella'],
                             ['Silas_the_Salmonella_tag-1-5.1;Silas_the_Salmonella_tag-1-5.2', 'Silas_the_Salmonella']]
        input_gffs = ['TestCheckingFragmentedGenes/Silas_the_Salmonella.gff',
                      'TestCheckingFragmentedGenes/Zion_the_Streptococcus.gff',
                      'TestCheckingFragmentedGenes/Silas_the_Legionella.gff',
                      'TestCheckingFragmentedGenes/Lilly_the_Shigella.gff']
        tmp_folder_path = 'test_tmp_folder'
        gene_data_file = {}
        corrected_dir = ''

        expected_return = [False, False]

        return_bool = parse_gene_presence_absence.check_fragmented_gene(fragments_info, input_gffs, tmp_folder_path,
                                                                        gene_data_file, corrected_dir, self.logger)

        self.assertEqual(expected_return, return_bool)


class TestParsingGenePresenceAbsenceFile(unittest.TestCase):
    """
    Tests for the function that passes the gene presence absence table from pan-genome program
    """
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def tearDown(self):
        try:
            os.remove('TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella_w_refound_db')
        except FileNotFoundError:
            pass

        try:
            for file in os.listdir('TestParsingGenePresenceAbsenceFile/Corrected_gffs/'):
                if '.gff' in file:
                    os.remove(os.path.join('TestParsingGenePresenceAbsenceFile/Corrected_gffs/', file))
        except FileNotFoundError:
            pass

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
        gene_data_file = {}
        corrected_dir = ''

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

        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
            file_name, core_gene_presence,
            low_freq_gene, source_program,
            input_gffs, tmp_folder_path, gene_data_file, corrected_dir, self.logger)

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
        gene_data_file = {}
        corrected_dir = ''


        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
                file_name, core_gene_presence,
                low_freq_gene, source_program,
                input_gffs, tmp_folder_path, gene_data_file, corrected_dir, self.logger)

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
        gene_data_file = {}
        corrected_dir = ''

        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
                file_name, core_gene_presence,
                low_freq_gene, source_program,
                input_gffs, tmp_folder_path, gene_data_file, corrected_dir, self.logger)

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
        gene_data_file = {}
        corrected_dir = ''

        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
                file_name, core_gene_presence,
                low_freq_gene, source_program,
                input_gffs, tmp_folder_path, gene_data_file, corrected_dir, self.logger)

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

    def test_parsign_fragmented_gene_w_refound_component(self):
        file_name = 'TestParsingGenePresenceAbsenceFile/gene_presence_absence_w_refound_fragment.csv'
        core_gene_presence = 0.9
        low_freq_gene = 0.1
        source_program = 'Panaroo'
        input_gffs = ['TestParsingGenePresenceAbsenceFile/Christina_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Ajwa_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Salmonella_w_refound.gff',
                      'TestParsingGenePresenceAbsenceFile/Cari_the_Listeria.gff',
                      'TestParsingGenePresenceAbsenceFile/Aman_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Zion_the_Streptococcus.gff',
                      'TestParsingGenePresenceAbsenceFile/Dina_the_Shigella.gff',
                      'TestParsingGenePresenceAbsenceFile/Silas_the_Legionella.gff',
                      'TestParsingGenePresenceAbsenceFile/Lilly_the_Shigella.gff']
        tmp_folder_path = 'TestParsingGenePresenceAbsenceFile/'
        gene_data_file = {'Silas_the_Salmonella_w_refound': {'0_refound_0': ['CTCTTCCGATCTAATCAAGATTGAGAGGAATTGCTGTTTTTATTGGCAAGACAATTTTACTCTTCCGATCTAATCAAGATTGAGAGGAATT', 'gene_name', 'gene_function']}}
        corrected_dir ='TestParsingGenePresenceAbsenceFile/Corrected_gffs'

        core_gene_dict, low_freq_gene_dict, \
        acc_gene_dict = \
            parse_gene_presence_absence.read_gene_presence_absence(
                file_name, core_gene_presence,
                low_freq_gene, source_program,
                input_gffs, tmp_folder_path, gene_data_file, corrected_dir, self.logger)

        expected_core_gene_dict = {'Silas_the_Salmonella_w_refound': {'Silas_the_Salmonella_tag-1-1': "A",
                                                                      '0_refound_0': "B",
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

        expected_low_freq_gene_dict = {'Silas_the_Salmonella_w_refound': {'Silas_the_Salmonella_tag_2': "G"},
                                       'Christina_the_Streptococcus': {},
                                       'Ajwa_the_Shigella': {},
                                       'Ajwa_the_Legionella': {},
                                       'Cari_the_Listeria': {},
                                       'Aman_the_Streptococcus': {},
                                       'Zion_the_Streptococcus': {},
                                       'Dina_the_Shigella': {},
                                       'Silas_the_Legionella': {},
                                       'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-6': "F"}}

        expected_acc_gene_dict = {'Silas_the_Salmonella_w_refound': {'Silas_the_Salmonella_tag-1-5.1': 'E',
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


class TestReadGeneData(unittest.TestCase):
    """ Function to test the passing of gene_data.csv file from Panaroo """
    def test_read_file(self):
        expected_dict = {'PY_40': {'0_refound_0': ['CTCTTCCGATCTAATCAAGATTGAGAGGAATTGCTGTTTTTATTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAG', 'gene_name', 'gene_function'],
                                   '0_refound_100': ['CTCTTCCGATCTAATCAAGATTGAGAGGAATTGCTTTTTTTTTTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAG', 'gene_name', 'gene_function'],
                                   '0_refound_10': ['CTCTTCCGATCTAATCAAGATTGAGAGGAATTGCGCCTTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAG', 'gene_name', 'gene_function']},
                         'PY_41': {'0_refound_1': ['ATGTTGTAGGAAAATACTTGGAAGAATACGTTGACAGGGGTATTTTTGATAAGGAGCCGTTCCAGACCTTTGATCAGAAAGGGATTGGCCGTCTCTTTAGCCCGTTCGGTTAAGCCTGAGTTGAAACTTGGTATTTGTGGGGAACATGGTGGCGATCCTGCTTCCATTGACTTTTACCACAGCCAAGGCCTGACCTACGTTTCTTGTTCGCCATTTAGAGTGCCGCTTACTCGCTTGGCGGCTGCTCAGGCTGCCATCAAAGCTTCAGGCCACAGTCTTACCCAAGACAAATAG', 'gene_name', 'gene_function']},
                         'PY_42': {'0_refound_2': ['ATGTCACTACTGCATATTCATCACAATAAAAAAAAGACAATAGCCCTAATCGTGCTATTGTCTCAAAATCATTTATTTACTTGAAACTTTATCGTGTTACACCAACAGTTTAA', 'gene_name', 'gene_function']},
                         'PY_43': {'0_refound_4': ['ATGAAACGCTATCAACAAGATGCCCTGCTTTTCAAAAAAAATAGATAAAGAAAAGGCTGCGACAGTATCTGCAAGCAGGGCAAAAGAACTAGAAGATAGGCTCAGTCATCAGCCATTAATTGATGATTATCGAGAAAAGATGCAAGATGCAAGATGCAAGTGATGTGACTCAGTATATCACCAAACGTATAGAAGATCAGTTAAACAAGGAGTTAACAAATGGCAAAAACTAA', 'gene_name', 'gene_function']}}

        return_dict = correct_gffs.read_gene_data('TestReadGeneData/Mock_gene_data.csv')

        self.assertEqual(expected_dict, return_dict)


class TestPrepairForReannotation(unittest.TestCase):
    """ Test for pre-pairing a folder for corrected genomes, and testing if any are present from previous runs """
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def tearDown(self):
        try:
            """ Class to remove created corrected output folder"""
            os.rmdir('TestPrepairForReannotation/Corrected_gff_files')
        except FileNotFoundError:
            pass

    def test_no_files_annotated(self):
        input_gffs = ['Mock_1.gff', 'Mock_2.gff']
        gene_data_dict_return, \
        corrected_gff_out_dir_return, \
        corrected_files_return = correct_gffs.prepair_for_reannotation('TestPrepairForReannotation/Mock_gene_data.csv',
                                                                       'TestPrepairForReannotation/',
                                                                       input_gffs, self.logger)

        self.assertTrue(os.path.isdir('TestPrepairForReannotation/Corrected_gff_files'))
        self.assertEqual(input_gffs, corrected_files_return)

    def test_some_files_annotated(self):
        input_gffs = ['Mock_1.gff', 'Mock_2.gff']
        gene_data_dict_return, corrected_gff_out_dir_return, corrected_files_return = correct_gffs.prepair_for_reannotation(
            'TestPrepairForReannotation/Mock_gene_data.csv',
            'TestPrepairForReannotation/Some_genomes',
            input_gffs, self.logger)

        expected_gffs = ['Mock_2.gff', 'Mock_1_corrected.gff']

        self.assertEqual(expected_gffs, corrected_files_return)

    def test_all_files_annotated(self):
        input_gffs = ['Mock_1.gff', 'Mock_2.gff']
        gene_data_dict_return, corrected_gff_out_dir_return, corrected_files_return = correct_gffs.prepair_for_reannotation(
            'TestPrepairForReannotation/Mock_gene_data.csv',
            'TestPrepairForReannotation/All_genomes',
            input_gffs, self.logger)

        expected_gffs = ['Mock_1_corrected.gff', 'Mock_2_corrected.gff']

        self.assertEqual(expected_gffs, corrected_files_return)


class TestAddGeneToGff(unittest.TestCase):
    """
    Test of the function used to add a gene annotation (line) to a gff file
    """
    # Make a setup and a teardown that copies and renames the mock file
    def setUp(self):
        """ Class to copy the mock gff before modifying"""
        copyfile('TestAddGeneToGff/mocky_test_gff.gff', 'TestAddGeneToGff/mocky_test_gff.gff_copy')

    def tearDown(self):
        """ Class to remove modified gff and rename the original"""
        os.remove('TestAddGeneToGff/mocky_test_gff.gff')
        os.rename('TestAddGeneToGff/mocky_test_gff.gff_copy', 'TestAddGeneToGff/mocky_test_gff.gff')

    def test_adding_a_gene_no_info(self):
        tmp_gff_file = 'TestAddGeneToGff/mocky_test_gff.gff'
        gene_oi = ['TATA', '', '']
        genome_oi = 'CCCCCCCCCCCCTATACCCCCCCC'
        contig = 'test_contig_1'
        strand = '+'
        refound_gene_tag = '0_refound_0'
        largest_locus_tag = 'fer_1432'

        expected_lines = ['##gff-version 3\n', '#test comment line\n', 'test_contig_1\tPanaroo\tCDS\t13\t16\t.\t+\t0\tID=fer_1433;locus_tag=fer_1433;old_locus_tag=0_refound_0\n']

        with open(tmp_gff_file, 'a') as tmp_gff:
            correct_gffs.add_gene_to_gff(tmp_gff, gene_oi[0], genome_oi, contig, strand, refound_gene_tag, gene_oi[1:], largest_locus_tag)

        with open('TestAddGeneToGff/mocky_test_gff.gff', 'r') as added_gff:
            self.assertEqual(expected_lines, added_gff.readlines())

    def test_adding_a_gene_name(self):
        tmp_gff_file = 'TestAddGeneToGff/mocky_test_gff.gff'
        gene_oi = ['TATA', 'Gene_name', '']
        genome_oi = 'CCCCCCCCCCCCTATACCCCCCCC'
        contig = 'test_contig_1'
        strand = '+'
        refound_gene_tag = '0_refound_0'
        largest_locus_tag = 'fer_1432'

        expected_lines = ['##gff-version 3\n', '#test comment line\n', 'test_contig_1\tPanaroo\tCDS\t13\t16\t.\t+\t0\tID=fer_1433;locus_tag=fer_1433;old_locus_tag=0_refound_0;name=Gene_name\n']

        with open(tmp_gff_file, 'a') as tmp_gff:
            correct_gffs.add_gene_to_gff(tmp_gff, gene_oi[0], genome_oi, contig, strand, refound_gene_tag, gene_oi[1:], largest_locus_tag)

        with open('TestAddGeneToGff/mocky_test_gff.gff', 'r') as added_gff:
            self.assertEqual(expected_lines, added_gff.readlines())

    def test_adding_a_gene_annotation(self):
        tmp_gff_file = 'TestAddGeneToGff/mocky_test_gff.gff'
        gene_oi = ['TATA', '', 'Gene_annotation']
        genome_oi = 'CCCCCCCCCCCCTATACCCCCCCC'
        contig = 'test_contig_1'
        strand = '+'
        refound_gene_tag = '0_refound_0'
        largest_locus_tag = 'fer_1432'

        expected_lines = ['##gff-version 3\n', '#test comment line\n', 'test_contig_1\tPanaroo\tCDS\t13\t16\t.\t+\t0\tID=fer_1433;locus_tag=fer_1433;old_locus_tag=0_refound_0;annotation=Gene_annotation\n']

        with open(tmp_gff_file, 'a') as tmp_gff:
            correct_gffs.add_gene_to_gff(tmp_gff, gene_oi[0], genome_oi, contig, strand, refound_gene_tag, gene_oi[1:], largest_locus_tag)

        with open('TestAddGeneToGff/mocky_test_gff.gff', 'r') as added_gff:
            self.assertEqual(expected_lines, added_gff.readlines())

    def test_adding_a_gene_name_and_annotation(self):
        tmp_gff_file = 'TestAddGeneToGff/mocky_test_gff.gff'
        gene_oi = ['TATA', 'Gene_name', 'Gene_annotation']
        genome_oi = 'CCCCCCCCCCCCTATACCCCCCCC'
        contig = 'test_contig_1'
        strand = '+'
        refound_gene_tag = '0_refound_0'
        largest_locus_tag = 'fer_1432'

        expected_lines = ['##gff-version 3\n', '#test comment line\n', 'test_contig_1\tPanaroo\tCDS\t13\t16\t.\t+\t0\tID=fer_1433;locus_tag=fer_1433;old_locus_tag=0_refound_0;name=Gene_name;annotation=Gene_annotation\n']

        with open(tmp_gff_file, 'a') as tmp_gff:
            correct_gffs.add_gene_to_gff(tmp_gff, gene_oi[0], genome_oi, contig, strand, refound_gene_tag, gene_oi[1:], largest_locus_tag)

        with open('TestAddGeneToGff/mocky_test_gff.gff', 'r') as added_gff:
            self.assertEqual(expected_lines, added_gff.readlines())


class TestWriteContig(unittest.TestCase):
    """
    Test of the function used to write a contig in a gff file.
    """
    # Make a setup and a teardown that copies and renames the mock file
    def setUp(self):
        """ Class to copy the mock gff before modifying"""
        copyfile('TestWriteContig/mocky_test_gff.gff', 'TestWriteContig/mocky_test_gff.gff_copy')

    def tearDown(self):
        """ Class to remove modified gff and rename the original"""
        os.remove('TestWriteContig/mocky_test_gff.gff')
        os.rename('TestWriteContig/mocky_test_gff.gff_copy', 'TestWriteContig/mocky_test_gff.gff')

    def test_writing_a_contig(self):
        file_path = 'TestWriteContig/mocky_test_gff.gff'
        contig_name = 'Test_contig_name space'
        sequence = 'AAATAAATGGGCGGGCAAATAAATGGGCGGGCAAATAAATGGGCGGGCAAATAAATGGGCGGGCAAATAAATGGGCGGGCAAATAAATGGGCGGGC'

        expected_lines = ['##gff-version 3\n', '#test comment line\n', '>Test_contig_name space\n',
                          'AAATAAATGGGCGGGCAAATAAATGGGCGGGCAAATAAATGGGCGGGCAAATAAATGGGC\n',
                          'GGGCAAATAAATGGGCGGGCAAATAAATGGGCGGGC\n']

        with open(file_path, 'a') as file:
            correct_gffs.write_contig(file, contig_name, sequence)

        with open('TestWriteContig/mocky_test_gff.gff', 'r') as added_gff:
            self.assertEqual(expected_lines, added_gff.readlines())


class TestAnnotateRefoundGenomes(unittest.TestCase):
    """
    Test of the function used to reannotate refound genes identified by panaroo in a gff file.
    """
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def tearDown(self):
        """ Class to remove modified gff and rename the original"""
        try:
            os.remove('TestAnnotateRefoundGenomes/reannotate_gff_corrected.gff')
        except FileNotFoundError:
            os.remove('TestAnnotateRefoundGenomes/reannotate_gff_tmp.gff')
            os.remove('TestAnnotateRefoundGenomes/reannotate_gff.gff_db')

    def test_annotation_of_pos_stand_gene(self):
        gff_name = 'TestAnnotateRefoundGenomes/reannotate_gff.gff'
        gene_data_dict = {'reannotate_gff': {'0_refound_0': ['CTCTTCCGATCTAATCAAGATTGAGAGGAATTGCTGTTTTTATTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAG', 'gene_name', 'gene_function'],
                                   '0_refound_100': ['CTCTTCCGATCTAATCAAGATTGAGAGGAATTGCTTTTTTTTTTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAG', '', 'gene_function'],
                                   '0_refound_10': ['CTCTTCCGATCTAATCAAGATTGAGAGGAATTGCGCCTTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAG', 'gene_name', '']}}
        tmp_folder_path = 'TestAnnotateRefoundGenomes'
        corrected_gff_out_dir = 'TestAnnotateRefoundGenomes'

        expected_lines = \
            ['##gff-version 3\n',
             '#test comment line\n',
             'test_contig\tProkka\tCDS\t1\t10\t.\t+\t0\tlocus_tag=locus_tag_0097\n',
             'test_contig\tPanaroo\tCDS\t16\t158\t.\t+\t0\tID=locus_tag_0099;locus_tag=locus_tag_0099;old_locus_tag=0_refound_0;name=gene_name;annotation=gene_function\n',
             'test_contig\tPanaroo\tCDS\t174\t316\t.\t+\t0\tID=locus_tag_0100;locus_tag=locus_tag_0100;old_locus_tag=0_refound_100;annotation=gene_function\n',
             'test_contig\tPanaroo\tCDS\t332\t469\t.\t+\t0\tID=locus_tag_0101;locus_tag=locus_tag_0101;old_locus_tag=0_refound_10;name=gene_name\n',
             'test_contig\tProkka\tCDS\t474\t484\t.\t+\t0\tlocus_tag=locus_tag_0098\n',
             '##FASTA\n',
             '>test_contig\n',
             'TTTTTTTTTTTTTTTCTCTTCCGATCTAATCAAGATTGAGAGGAATTGCTGTTTTTATTG\n',
             'GCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAG\n',
             'GACATGGTCAAAGTAACTTTATTTATAATGAATTTTAGTTTTTTTTTTTTTTTCTCTTCC\n',
             'GATCTAATCAAGATTGAGAGGAATTGCTTTTTTTTTTGGCAAGACAATTTTATTTTATCT\n',
             'GATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTAT\n',
             'TTATAATGAATTTTAGTTTTTTTTTTTTTTTCTCTTCCGATCTAATCAAGATTGAGAGGA\n',
             'ATTGCGCCTTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTT\n',
             'TGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAGTTTTTTTTTTT\n',
             'TTTT\n'
             ]

        correct_gffs.annotate_refound_genes(gff_name, gene_data_dict, tmp_folder_path, corrected_gff_out_dir, self.logger)

        with open('TestAnnotateRefoundGenomes/reannotate_gff_corrected.gff', 'r') as added_gff:
            self.assertEqual(expected_lines, added_gff.readlines())

    def test_annotation_of_neg_stand_gene(self):
        gff_name = 'TestAnnotateRefoundGenomes/reannotate_gff.gff'
        gene_data_dict = {'reannotate_gff': {'0_refound_0': ['CTAAAATTCATTATAAATAAAGTTACTTTGACCATGTCCTAGATAACCAAAGAAATAGACCAGCAACATTAAAATCAGATAAAATAAAATTGTCTTGCCAATAAAAACAGCAATTCCTCTCAATCTTGATTAGATCGGAAGAG', 'gene_name', 'gene_function'],
                                             '0_refound_100': ['CTAAAATTCATTATAAATAAAGTTACTTTGACCATGTCCTAGATAACCAAAGAAATAGACCAGCAACATTAAAATCAGATAAAATAAAATTGTCTTGCCAAAAAAAAAAGCAATTCCTCTCAATCTTGATTAGATCGGAAGAG', '', 'gene_function'],
                                             '0_refound_10': ['CTAAAATTCATTATAAATAAAGTTACTTTGACCATGTCCTAGATAACCAAAGAAATAGACCAGCAACATTAAAATCAGATAAAATAAAATTGTCTTGCCAAGGCGCAATTCCTCTCAATCTTGATTAGATCGGAAGAG', 'gene_name', '']}}
        tmp_folder_path = 'TestAnnotateRefoundGenomes'
        corrected_gff_out_dir = 'TestAnnotateRefoundGenomes'
        expected_lines = ['##gff-version 3\n',
                          '#test comment line\n',
                          'test_contig\tProkka\tCDS\t1\t10\t.\t+\t0\tlocus_tag=locus_tag_0097\n',
                          'test_contig\tPanaroo\tCDS\t16\t158\t.\t-\t0\tID=locus_tag_0099;locus_tag=locus_tag_0099;old_locus_tag=0_refound_0;name=gene_name;annotation=gene_function\n',
                          'test_contig\tPanaroo\tCDS\t174\t316\t.\t-\t0\tID=locus_tag_0100;locus_tag=locus_tag_0100;old_locus_tag=0_refound_100;annotation=gene_function\n',
                          'test_contig\tPanaroo\tCDS\t332\t469\t.\t-\t0\tID=locus_tag_0101;locus_tag=locus_tag_0101;old_locus_tag=0_refound_10;name=gene_name\n',
                          'test_contig\tProkka\tCDS\t474\t484\t.\t+\t0\tlocus_tag=locus_tag_0098\n',
                          '##FASTA\n',
                          '>test_contig\n',
                          'TTTTTTTTTTTTTTTCTCTTCCGATCTAATCAAGATTGAGAGGAATTGCTGTTTTTATTG\n',
                          'GCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAG\n',
                          'GACATGGTCAAAGTAACTTTATTTATAATGAATTTTAGTTTTTTTTTTTTTTTCTCTTCC\n',
                          'GATCTAATCAAGATTGAGAGGAATTGCTTTTTTTTTTGGCAAGACAATTTTATTTTATCT\n',
                          'GATTTTAATGTTGCTGGTCTATTTCTTTGGTTATCTAGGACATGGTCAAAGTAACTTTAT\n',
                          'TTATAATGAATTTTAGTTTTTTTTTTTTTTTCTCTTCCGATCTAATCAAGATTGAGAGGA\n',
                          'ATTGCGCCTTGGCAAGACAATTTTATTTTATCTGATTTTAATGTTGCTGGTCTATTTCTT\n',
                          'TGGTTATCTAGGACATGGTCAAAGTAACTTTATTTATAATGAATTTTAGTTTTTTTTTTT\n',
                          'TTTT\n']

        correct_gffs.annotate_refound_genes(gff_name, gene_data_dict, tmp_folder_path, corrected_gff_out_dir, self.logger)

        with open('TestAnnotateRefoundGenomes/reannotate_gff_corrected.gff', 'r') as added_gff:
            self.assertEqual(expected_lines, added_gff.readlines())

    def test_gene_not_found(self):
        gff_name = 'TestAnnotateRefoundGenomes/reannotate_gff.gff'
        gene_data_dict = {'reannotate_gff': {'0_refound_0': [
            'CCCCCCCCCCCCGGGGGGGGGGGGGGGCGGCGCGCGCGCGCGCGGCGCGCGCGGCGCGC',
            'gene_name', 'gene_function']}}

        tmp_folder_path = 'TestAnnotateRefoundGenomes'
        corrected_gff_out_dir = 'TestAnnotateRefoundGenomes'

        with self.assertRaises(SystemExit):
            correct_gffs.annotate_refound_genes(gff_name, gene_data_dict, tmp_folder_path, corrected_gff_out_dir, self.logger)

    # TODO - Add test for annotating of second contig


class TestExtractGenomeFasta(unittest.TestCase):
    def test_extract_genome_fasta(self):
        genome_fasta_dict_expected = {'contig_1': "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}
        largest_locus_tag_expected = 'fer_006'
        header_lines_expected = ['##gff-version3\n', '#test-line\n']

        genome_fasta_dict, largest_locus_tag, header_lines = correct_gffs.extract_genome_fasta('TestExtractGenomeFasta/Mock_gff.gff')

        self.assertEqual(genome_fasta_dict_expected, genome_fasta_dict)
        self.assertEqual(largest_locus_tag_expected, largest_locus_tag)
        self.assertEqual(header_lines_expected, header_lines)


class TestParsingGffFile(unittest.TestCase):
    """ Test of the function that is used to pass a gff file and return a generator object of CDS lines """
    def test_gff_generator_generation_not_corrected(self):
        input_gff_file = 'TestParsingGffFile/Silas_the_Salmonella.gff'

        expected_output = [['contig_1', '.', 'CDS', '1', '90', '.', '.', '.', 'Silas_the_Salmonella_tag-1-1'],
                           ['contig_1', '.', 'CDS', '100', '190', '.', '.', '.', 'Silas_the_Salmonella_tag-1-2.1'],
                           ['contig_1', '.', 'CDS', '200', '290', '.', '.', '.', 'Silas_the_Salmonella_tag-1-2.2'],
                           ['contig_1', '.', 'CDS', '300', '390', '.', '.', '.', 'Silas_the_Salmonella_tag-1-3'],
                           ['contig_1', '.', 'CDS', '400', '490', '.', '.', '.', 'Silas_the_Salmonella_tag-1-4.1'],
                           ['contig_1', '.', 'CDS', '500', '590', '.', '.', '.', 'Silas_the_Salmonella_tag-1-4.2'],
                           ['contig_1', '.', 'CDS', '600', '690', '.', '.', '.', 'Silas_the_Salmonella_tag-1-5.1'],
                           ['contig_1', '.', 'CDS', '700', '790', '.', '.', '.', 'Silas_the_Salmonella_tag-1.7'],
                           ['contig_1', '.', 'CDS', '800', '890', '.', '.', '.', "Silas_the_Salmonella_tag-1-5.2"]]

        return_generator = []
        for line in gff_parser.parse_gff(input_gff_file):
            return_generator += [line]

        for expected, generated in zip(expected_output, return_generator):
            self.assertEqual(expected, generated)

    def test_gff_generator_generation_gzipped_input(self):
        input_gff_file = 'TestParsingGffFile/Silas_the_Salmonella.gff.gz'

        expected_output = [['contig_1', '.', 'CDS', '1', '90', '.', '.', '.', 'Silas_the_Salmonella_tag-1-1'],
                           ['contig_1', '.', 'CDS', '100', '190', '.', '.', '.', 'Silas_the_Salmonella_tag-1-2.1'],
                           ['contig_1', '.', 'CDS', '200', '290', '.', '.', '.', 'Silas_the_Salmonella_tag-1-2.2'],
                           ['contig_1', '.', 'CDS', '300', '390', '.', '.', '.', 'Silas_the_Salmonella_tag-1-3'],
                           ['contig_1', '.', 'CDS', '400', '490', '.', '.', '.', 'Silas_the_Salmonella_tag-1-4.1'],
                           ['contig_1', '.', 'CDS', '500', '590', '.', '.', '.', 'Silas_the_Salmonella_tag-1-4.2'],
                           ['contig_1', '.', 'CDS', '600', '690', '.', '.', '.', 'Silas_the_Salmonella_tag-1-5.1'],
                           ['contig_1', '.', 'CDS', '700', '790', '.', '.', '.', 'Silas_the_Salmonella_tag-1.7'],
                           ['contig_1', '.', 'CDS', '800', '890', '.', '.', '.', "Silas_the_Salmonella_tag-1-5.2"]]

        return_generator = []
        for line in gff_parser.parse_gff(input_gff_file):
            return_generator += [line]

        for expected, generated in zip(expected_output, return_generator):
            self.assertEqual(expected, generated)

    def test_gff_generator_generation_corrected_gff(self):
        input_gff_file = 'TestParsingGffFile/Silas_the_Salmonella_corrected.gff'

        expected_output = [['contig_1', '.', 'CDS', '1', '90', '.', '.', '.', 'Silas_the_Salmonella_tag-1-1'],
                           ['contig_1', '.', 'CDS', '100', '190', '.', '.', '.', 'Silas_the_Salmonella_tag-1-2.1'],
                           ['contig_1', '.', 'CDS', '200', '290', '.', '.', '.', 'Silas_the_Salmonella_tag-1-2.2'],
                           ['contig_1', '.', 'CDS', '300', '390', '.', '.', '.', 'Silas_the_Salmonella_tag-1-3'],
                           ['contig_1', '.', 'CDS', '400', '490', '.', '.', '.', 'Silas_the_Salmonella_tag-1-4.1'],
                           ['contig_1', '.', 'CDS', '500', '590', '.', '.', '.', 'Silas_the_Salmonella_tag-1-4.2'],
                           ['contig_1', '.', 'CDS', '600', '690', '.', '.', '.', 'Silas_the_Salmonella_tag-1-5.1'],
                           ['contig_1', '.', 'CDS', '700', '790', '.', '.', '.', 'Silas_the_Salmonella_tag-1.7'],
                           ['contig_1', '.', 'CDS', '800', '890', '.', '.', '.', "Silas_the_Salmonella_tag-1-5.2"],
                           ['contig_1', 'Panaroo', 'CDS', '900', '1000', '.', '+', '0', 'refound_gene_1']]

        return_generator = []
        for line in gff_parser.parse_gff(input_gff_file):
            return_generator += [line]

        for expected, generated in zip(expected_output, return_generator):
            self.assertEqual(expected, generated)


class TestGetContigLenth(unittest.TestCase):
    """
    Test function that passes a gff file and counts the length of each contig in attached genome
    """
    def test_single_contig(self):
        input_gff_path = 'TestGetContigLenth/single_contig_unwrapped.txt'
        expected_dict = {'contig_1': 1300}

        return_dict = gff_parser.get_contig_lengths(input_gff_path)

        self.assertEqual(expected_dict, return_dict)

    def test_single_wrapped_contig(self):
        input_gff_path = 'TestGetContigLenth/single_contig_wrapped.txt'
        expected_dict = {'contig_1': 1300}

        return_dict = gff_parser.get_contig_lengths(input_gff_path)

        self.assertEqual(expected_dict, return_dict)

    def test_multiple_contigs(self):
        input_gff_path = 'TestGetContigLenth/multi_contig_unwrapped.txt'
        expected_dict = {'contig_1': 1300,
                         'contig_2': 1300}

        return_dict = gff_parser.get_contig_lengths(input_gff_path)

        self.assertEqual(expected_dict, return_dict)

    def test_multiple_wrapped_contigs(self):
        input_gff_path = 'TestGetContigLenth/multi_contig_wrapped.txt'
        expected_dict = {'contig_1': 1300,
                         'contig_2': 1300}

        return_dict = gff_parser.get_contig_lengths(input_gff_path)

        self.assertEqual(expected_dict, return_dict)


class TestRecordCoreCoreRegion(unittest.TestCase):
    """
    Test function that is used to record information of a region identified between two core genes.
    """
    def test_recording_neighbouring_core_genes(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1',
                                   'Core_ID_2': 'pan_gene_2'}}
        gff_name = 'gff_name'
        gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '150', '.', '.', '.', 'Core_ID_2']
        contig_end = 1500
        previous_core_gene_id = 'Core_ID_1'
        previous_core_gene_end_coor = 10
        acc_genes_in_region = []
        low_freq_genes_in_region = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        core_gene_pairs = []
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = 'Core_ID_2'
        expected_previous_core_gene_end_coor = 150
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_gene_1--pan_gene_2': 79}
        expected_accessory_gene_content = {'pan_gene_1--pan_gene_2': []}
        expected_low_freq_gene_content = {'pan_gene_1--pan_gene_2': []}
        expected_core_gene_pairs = ['pan_gene_1--pan_gene_2']
        expected_master_info = {'pan_gene_1--pan_gene_2--gff_name': ['gff_name', 'pan_gene_1', 'pan_gene_2',
                                                                     79, 0, [], []]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pair_distance, return_accessory_gene_content, \
        return_low_freq_gene_content, return_core_gene_pairs, \
        return_master_info = gff_parser.record_core_core_region(core_genes, gff_name, gff_line, contig_end,
                                                                previous_core_gene_id, previous_core_gene_end_coor,
                                                                acc_genes_in_region, low_freq_genes_in_region,
                                                                core_gene_pair_distance, accessory_gene_content,
                                                                low_freq_gene_content, core_gene_pairs,
                                                                master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_recording_neighbouring_core_genes_w_accessory(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1',
                                   'Core_ID_2': 'pan_gene_2'}}
        gff_name = 'gff_name'
        gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '150', '.', '.', '.', 'Core_ID_2']
        contig_end = 1500
        previous_core_gene_id = 'Core_ID_1'
        previous_core_gene_end_coor = 10
        acc_genes_in_region = ['acc_gene_1', 'acc_gene_2']
        low_freq_genes_in_region = ['low_freq_1']
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        core_gene_pairs = []
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = 'Core_ID_2'
        expected_previous_core_gene_end_coor = 150
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_gene_1--pan_gene_2': 79}
        expected_accessory_gene_content = {'pan_gene_1--pan_gene_2': ['acc_gene_1', 'acc_gene_2']}
        expected_low_freq_gene_content = {'pan_gene_1--pan_gene_2': ['low_freq_1']}
        expected_core_gene_pairs = ['pan_gene_1--pan_gene_2']
        expected_master_info = {'pan_gene_1--pan_gene_2--gff_name': ['gff_name', 'pan_gene_1', 'pan_gene_2',
                                                                     79, 3, ['acc_gene_1', 'acc_gene_2'], ['low_freq_1']]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pair_distance, return_accessory_gene_content, \
        return_low_freq_gene_content, return_core_gene_pairs, \
        return_master_info = gff_parser.record_core_core_region(core_genes, gff_name, gff_line, contig_end,
                                                                previous_core_gene_id, previous_core_gene_end_coor,
                                                                acc_genes_in_region, low_freq_genes_in_region,
                                                                core_gene_pair_distance, accessory_gene_content,
                                                                low_freq_gene_content, core_gene_pairs,
                                                                master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_recording_w_fragment_given(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1',
                                   'Core_ID_2': 'pan_gene_1'}}
        gff_name = 'gff_name'
        gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '150', '.', '.', '.', 'Core_ID_2']
        contig_end = 1500
        previous_core_gene_id = 'Core_ID_1'
        previous_core_gene_end_coor = 10
        acc_genes_in_region = []
        low_freq_genes_in_region = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        core_gene_pairs = []
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = 'Core_ID_2'
        expected_previous_core_gene_end_coor = 150
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {}
        expected_accessory_gene_content = {}
        expected_low_freq_gene_content = {}
        expected_core_gene_pairs = []
        expected_master_info = {}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pair_distance, return_accessory_gene_content, \
        return_low_freq_gene_content, return_core_gene_pairs, \
        return_master_info = gff_parser.record_core_core_region(core_genes, gff_name, gff_line, contig_end,
                                                                previous_core_gene_id, previous_core_gene_end_coor,
                                                                acc_genes_in_region, low_freq_genes_in_region,
                                                                core_gene_pair_distance, accessory_gene_content,
                                                                low_freq_gene_content, core_gene_pairs,
                                                                master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_recording_core_gene_before_seqeuncebreak(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1'}}
        gff_name = 'gff_name'
        gff_line = None
        contig_end = 1500
        previous_core_gene_id = 'Core_ID_1'
        previous_core_gene_end_coor = 150
        acc_genes_in_region = []
        low_freq_genes_in_region = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        core_gene_pairs = []
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = 'Sequence_break'
        expected_previous_core_gene_end_coor = 150
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_gene_1--Sequence_break': 1350}
        expected_accessory_gene_content = {'pan_gene_1--Sequence_break': []}
        expected_low_freq_gene_content = {'pan_gene_1--Sequence_break': []}
        expected_core_gene_pairs = ['pan_gene_1--Sequence_break']
        expected_master_info = {'pan_gene_1--Sequence_break--gff_name': ['gff_name', 'pan_gene_1', 'Sequence_break',
                                                                         1350, 0, [], []]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pair_distance, return_accessory_gene_content, \
        return_low_freq_gene_content, return_core_gene_pairs, \
        return_master_info = gff_parser.record_core_core_region(core_genes, gff_name, gff_line, contig_end,
                                                                previous_core_gene_id, previous_core_gene_end_coor,
                                                                acc_genes_in_region, low_freq_genes_in_region,
                                                                core_gene_pair_distance, accessory_gene_content,
                                                                low_freq_gene_content, core_gene_pairs,
                                                                master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_recording_core_gene_before_seqeuncebreak_w_accessory(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1'}}
        gff_name = 'gff_name'
        gff_line = None
        contig_end = 1500
        previous_core_gene_id = 'Core_ID_1'
        previous_core_gene_end_coor = 150
        acc_genes_in_region = ['acc_1']
        low_freq_genes_in_region = ['low_1', "low_2"]
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        core_gene_pairs = []
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = 'Sequence_break'
        expected_previous_core_gene_end_coor = 150
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'pan_gene_1--Sequence_break': 1350}
        expected_accessory_gene_content = {'pan_gene_1--Sequence_break': ['acc_1']}
        expected_low_freq_gene_content = {'pan_gene_1--Sequence_break': ['low_1', "low_2"]}
        expected_core_gene_pairs = ['pan_gene_1--Sequence_break']
        expected_master_info = {'pan_gene_1--Sequence_break--gff_name': ['gff_name', 'pan_gene_1', 'Sequence_break',
                                                                         1350, 3, ['acc_1'], ['low_1', "low_2"]]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pair_distance, return_accessory_gene_content, \
        return_low_freq_gene_content, return_core_gene_pairs, \
        return_master_info = gff_parser.record_core_core_region(core_genes, gff_name, gff_line, contig_end,
                                                                previous_core_gene_id, previous_core_gene_end_coor,
                                                                acc_genes_in_region, low_freq_genes_in_region,
                                                                core_gene_pair_distance, accessory_gene_content,
                                                                low_freq_gene_content, core_gene_pairs,
                                                                master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_recording_first_core_gene_on_contig_as_first_gene(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1'}}
        gff_name = 'gff_name'
        gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'Core_ID_1']
        contig_end = 0
        previous_core_gene_id = 'Sequence_break'
        previous_core_gene_end_coor = 150
        acc_genes_in_region = []
        low_freq_genes_in_region = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        core_gene_pairs = []
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = 'Core_ID_1'
        expected_previous_core_gene_end_coor = 180
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'Sequence_break--pan_gene_1': 89}
        expected_accessory_gene_content = {'Sequence_break--pan_gene_1': []}
        expected_low_freq_gene_content = {'Sequence_break--pan_gene_1': []}
        expected_core_gene_pairs = ['Sequence_break--pan_gene_1']
        expected_master_info = {'Sequence_break--pan_gene_1--gff_name': ['gff_name', 'Sequence_break', 'pan_gene_1',
                                                                         89, 0, [], []]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pair_distance, return_accessory_gene_content, \
        return_low_freq_gene_content, return_core_gene_pairs, \
        return_master_info = gff_parser.record_core_core_region(core_genes, gff_name, gff_line, contig_end,
                                                                previous_core_gene_id, previous_core_gene_end_coor,
                                                                acc_genes_in_region, low_freq_genes_in_region,
                                                                core_gene_pair_distance, accessory_gene_content,
                                                                low_freq_gene_content, core_gene_pairs,
                                                                master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_recording_first_core_gene_on_contig_w_accessory(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1'}}
        gff_name = 'gff_name'
        gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'Core_ID_1']
        contig_end = None
        previous_core_gene_id = 'Sequence_break'
        previous_core_gene_end_coor = 150
        acc_genes_in_region = ['acc_1', 'acc_2', 'acc_3']
        low_freq_genes_in_region = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        core_gene_pairs = []
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = 'Core_ID_1'
        expected_previous_core_gene_end_coor = 180
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pair_distance = {'Sequence_break--pan_gene_1': 89}
        expected_accessory_gene_content = {'Sequence_break--pan_gene_1': ['acc_1', 'acc_2', 'acc_3']}
        expected_low_freq_gene_content = {'Sequence_break--pan_gene_1': []}
        expected_core_gene_pairs = ['Sequence_break--pan_gene_1']
        expected_master_info = {'Sequence_break--pan_gene_1--gff_name': ['gff_name', 'Sequence_break', 'pan_gene_1',
                                                                         89, 3, ['acc_1', 'acc_2', 'acc_3'], []]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pair_distance, return_accessory_gene_content, \
        return_low_freq_gene_content, return_core_gene_pairs, \
        return_master_info = gff_parser.record_core_core_region(core_genes, gff_name, gff_line, contig_end,
                                                                previous_core_gene_id, previous_core_gene_end_coor,
                                                                acc_genes_in_region, low_freq_genes_in_region,
                                                                core_gene_pair_distance, accessory_gene_content,
                                                                low_freq_gene_content, core_gene_pairs,
                                                                master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)


class TestConnectFirstNLastGeneOnContig(unittest.TestCase):
    """
    Test for the function recordning connections between the first and the last gene on a contig in a complete genome
    """
    def test_connect_last_n_first_gene_different_genes_no_accessory(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1',
                                   'Core_ID_2': 'pan_gene_2'}}
        gff_name = 'gff_name'
        previous_core_gene_id = "Core_ID_2"
        previous_core_gene_end_coor = 1450
        first_core_gene_gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'Core_ID_1']
        acc_genes_in_region = []
        first_core_accessory_content = []
        low_freq_genes_in_region = []
        first_core_low_freq_genes = []
        contig_size = 1500
        core_gene_pairs = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = ""
        expected_previous_core_gene_end_coor = 180
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pairs = ['pan_gene_1--pan_gene_2']
        expected_core_gene_pair_distance = {'pan_gene_1--pan_gene_2': 139}
        expected_accessory_gene_content = {'pan_gene_1--pan_gene_2': []}
        expected_low_freq_gene_content = {'pan_gene_1--pan_gene_2': []}
        expected_master_info = {'pan_gene_1--pan_gene_2--gff_name': ['gff_name', 'pan_gene_1', 'pan_gene_2', 139, 0, [], []]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info = gff_parser.connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id,
                                                                            previous_core_gene_end_coor,
                                                                            first_core_gene_gff_line,
                                                                            acc_genes_in_region,
                                                                            first_core_accessory_content,
                                                                            low_freq_genes_in_region,
                                                                            first_core_low_freq_genes, contig_size,
                                                                            core_gene_pairs, core_gene_pair_distance,
                                                                            accessory_gene_content,
                                                                            low_freq_gene_content, master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_connect_last_n_first_gene_different_genes_w_accessory(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1',
                                   'Core_ID_2': 'pan_gene_2'}}
        gff_name = 'gff_name'
        previous_core_gene_id = "Core_ID_2"
        previous_core_gene_end_coor = 1450
        first_core_gene_gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'Core_ID_1']
        acc_genes_in_region = ['acc_1']
        first_core_accessory_content = ['first_acc_1']
        low_freq_genes_in_region = ['low_acc_1']
        first_core_low_freq_genes = ['first_low_1']
        contig_size = 1500
        core_gene_pairs = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = ""
        expected_previous_core_gene_end_coor = 180
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pairs = ['pan_gene_1--pan_gene_2']
        expected_core_gene_pair_distance = {'pan_gene_1--pan_gene_2': 139}
        expected_accessory_gene_content = {'pan_gene_1--pan_gene_2': ['acc_1', 'first_acc_1']}
        expected_low_freq_gene_content = {'pan_gene_1--pan_gene_2': ['first_low_1', 'low_acc_1']}
        expected_master_info = {'pan_gene_1--pan_gene_2--gff_name': ['gff_name', 'pan_gene_1', 'pan_gene_2', 139, 4, ['acc_1', 'first_acc_1'], ['first_low_1', 'low_acc_1']]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info = gff_parser.connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id,
                                                                            previous_core_gene_end_coor,
                                                                            first_core_gene_gff_line,
                                                                            acc_genes_in_region,
                                                                            first_core_accessory_content,
                                                                            low_freq_genes_in_region,
                                                                            first_core_low_freq_genes, contig_size,
                                                                            core_gene_pairs, core_gene_pair_distance,
                                                                            accessory_gene_content,
                                                                            low_freq_gene_content, master_info)

        # Assert expected against returned
        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_connect_same_gene_as_last_n_first_gene_no_accessory(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1',
                                   'Core_ID_2': 'pan_gene_2'}}
        gff_name = 'gff_name'
        previous_core_gene_id = "Core_ID_1"
        previous_core_gene_end_coor = 180
        first_core_gene_gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'Core_ID_1']
        acc_genes_in_region = []
        first_core_accessory_content = []
        low_freq_genes_in_region = []
        first_core_low_freq_genes = []
        contig_size = 1500
        core_gene_pairs = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = ""
        expected_previous_core_gene_end_coor = 180
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pairs = ['pan_gene_1--pan_gene_1']
        expected_core_gene_pair_distance = {'pan_gene_1--pan_gene_1': 1409}
        expected_accessory_gene_content = {'pan_gene_1--pan_gene_1': []}
        expected_low_freq_gene_content = {'pan_gene_1--pan_gene_1': []}
        expected_master_info = {'pan_gene_1--pan_gene_1--gff_name': ['gff_name', 'pan_gene_1', 'pan_gene_1', 1409, 0, [], []]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info = gff_parser.connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id,
                                                                            previous_core_gene_end_coor,
                                                                            first_core_gene_gff_line,
                                                                            acc_genes_in_region,
                                                                            first_core_accessory_content,
                                                                            low_freq_genes_in_region,
                                                                            first_core_low_freq_genes, contig_size,
                                                                            core_gene_pairs, core_gene_pair_distance,
                                                                            accessory_gene_content,
                                                                            low_freq_gene_content, master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)

    def test_connect_same_gene_as_last_n_first_gene_w_accessory(self):
        core_genes = {'gff_name': {'Core_ID_1': 'pan_gene_1',
                                   'Core_ID_2': 'pan_gene_2'}}
        gff_name = 'gff_name'
        previous_core_gene_id = "Core_ID_1"
        previous_core_gene_end_coor = 180
        first_core_gene_gff_line = ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'Core_ID_1']
        acc_genes_in_region = ['acc_2', 'acc_3']
        first_core_accessory_content = ['acc_1']
        low_freq_genes_in_region = ['low_1', 'low_2']
        first_core_low_freq_genes = ['low_3']
        contig_size = 1500
        core_gene_pairs = []
        core_gene_pair_distance = {}
        accessory_gene_content = {}
        low_freq_gene_content = {}
        master_info = {}

        # Set up the expected return values
        expected_previous_core_gene_id = ""
        expected_previous_core_gene_end_coor = 180
        expected_acc_genes_in_region = []
        expected_low_freq_genes_in_region = []
        expected_core_gene_pairs = ['pan_gene_1--pan_gene_1']
        expected_core_gene_pair_distance = {'pan_gene_1--pan_gene_1': 1409}
        expected_accessory_gene_content = {'pan_gene_1--pan_gene_1': ['acc_1', 'acc_2', 'acc_3']}
        expected_low_freq_gene_content = {'pan_gene_1--pan_gene_1': ['low_1', 'low_2', 'low_3']}
        expected_master_info = {'pan_gene_1--pan_gene_1--gff_name': ['gff_name', 'pan_gene_1', 'pan_gene_1', 1409, 6, ['acc_1', 'acc_2', 'acc_3'], ['low_1', 'low_2', 'low_3']]}

        return_previous_core_gene_id, return_previous_core_gene_end_coor, return_acc_genes_in_region, \
        return_low_freq_genes_in_region, return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info = gff_parser.connect_first_n_last_gene_on_contig(core_genes, gff_name, previous_core_gene_id,
                                                                            previous_core_gene_end_coor,
                                                                            first_core_gene_gff_line,
                                                                            acc_genes_in_region,
                                                                            first_core_accessory_content,
                                                                            low_freq_genes_in_region,
                                                                            first_core_low_freq_genes, contig_size,
                                                                            core_gene_pairs, core_gene_pair_distance,
                                                                            accessory_gene_content,
                                                                            low_freq_gene_content, master_info)

        self.assertEqual(expected_previous_core_gene_id, return_previous_core_gene_id)
        self.assertEqual(expected_previous_core_gene_end_coor, return_previous_core_gene_end_coor)
        self.assertEqual(expected_acc_genes_in_region, return_acc_genes_in_region)
        self.assertEqual(expected_low_freq_genes_in_region, return_low_freq_genes_in_region)
        self.assertEqual(expected_accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(expected_low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(expected_core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(expected_core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(expected_master_info, return_master_info)


class TestRecordCorelessContig(unittest.TestCase):
    """
    Test for function that records a contig on which there is no core genes.
    """
    def test_adding_coreless_contig(self):
        coreless_contigs = {}
        acc_genes_in_region = ['acc_1']
        low_freq_genes_in_region = ['low_1']
        gff_name = 'gff_name'
        contig_name = 'gff_contig_1'

        expected_return = {'gff_name--gff_contig_1': [['acc_1'], ['low_1']]}

        return_dict = gff_parser.record_coreless_contig(coreless_contigs, acc_genes_in_region, low_freq_genes_in_region, gff_name, contig_name)

        self.assertEqual(expected_return, return_dict)

    def test_not_adding_coreless_contig(self):
        coreless_contigs = {}
        acc_genes_in_region = []
        low_freq_genes_in_region = []
        gff_name = 'gff_name'
        contig_name = 'gff_contig_1'

        expected_return = {}

        return_dict = gff_parser.record_coreless_contig(coreless_contigs, acc_genes_in_region,
                                                         low_freq_genes_in_region, gff_name, contig_name)

        self.assertEqual(expected_return, return_dict)


class TestSegmentingMockGffs(unittest.TestCase):
    """
    Tests for function that takes in a gff file and segments it into core-core regions
    """
    def test_single_chromosome_complete(self):
        # Set up input
        gff_generator = [['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
                         ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
                         ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
                         ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3']]
        core_genes = {'test_single_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_single_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8'}}
        gff_path = 'TestSegmentingMockGffs/test_single_chromosome.gff'
        acc_genes = {'test_single_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4'}}
        complete_genomes = ['test_single_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7', 'pan_gene_2--pan_gene_7']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1']}

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8']}

        master_info = {'pan_gene_2--pan_gene_5--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_5', 359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
                       'pan_gene_5--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_5', 'pan_gene_7', 269, 1, [], ['pan_gene_6']],
                       'pan_gene_2--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_7', 478, 2, ['pan_gene_1'], ['pan_gene_8']]}

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                              low_freq_genes, gff_path,
                                                                                              acc_genes,
                                                                                              complete_genomes)


        # Evaluate
        self.assertEqual(core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_single_chromosome_draft(self):
        # Set up input
        gff_generator = [['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
                         ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
                         ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
                         ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3']]
        core_genes = {'test_single_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_single_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8'}}
        gff_path = 'TestSegmentingMockGffs/test_single_chromosome.gff'
        acc_genes = {'test_single_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4'}}
        complete_genomes = []

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7', 'Sequence_break--pan_gene_2', 'pan_gene_7--Sequence_break']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'Sequence_break--pan_gene_2': 178,
                                   'pan_gene_7--Sequence_break': 300}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'Sequence_break--pan_gene_2': ['pan_gene_1'],
                                  'pan_gene_7--Sequence_break': []}

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'Sequence_break--pan_gene_2': [],
                                 'pan_gene_7--Sequence_break': ['pan_gene_8']}

        master_info = {
            'pan_gene_2--pan_gene_5--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_5', 359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_5', 'pan_gene_7', 269, 1, [], ['pan_gene_6']],
            'Sequence_break--pan_gene_2--test_single_chromosome': ['test_single_chromosome', 'Sequence_break', 'pan_gene_2', 178, 1, ['pan_gene_1'], []],
            'pan_gene_7--Sequence_break--test_single_chromosome': ['test_single_chromosome',  'pan_gene_7', 'Sequence_break', 300, 1, [], ['pan_gene_8']]}

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_two_chromosomes_complete(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_3'],
            ['gff_name_contig_2', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_5'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
            ['gff_name_contig_2', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_6'],
            ['gff_name_contig_2', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_6']
                         ]
        core_genes = {'test_double_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7',
                                                 'Core_ID_4': 'pan_gene_10',
                                                 'Core_ID_5': 'pan_gene_13',
                                                 'Core_ID_6': 'pan_gene_15'
        }}
        low_freq_genes = {'test_double_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14',
                                                     'low_freq_6': 'pan_gene_16'}}
        gff_path = 'TestSegmentingMockGffs/test_double_chromosome.gff'
        acc_genes = {'test_double_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_3': 'pan_gene_9',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = ['test_double_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7', 'pan_gene_2--pan_gene_7',
                           'pan_gene_10--pan_gene_13', 'pan_gene_13--pan_gene_15', 'pan_gene_10--pan_gene_15']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478,
                                   'pan_gene_10--pan_gene_13': 359,
                                   'pan_gene_13--pan_gene_15': 269,
                                   'pan_gene_10--pan_gene_15': 478}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1'],
                                  'pan_gene_10--pan_gene_13': ['pan_gene_12'],
                                  'pan_gene_13--pan_gene_15': [],
                                  'pan_gene_10--pan_gene_15': ['pan_gene_9']
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8'],
                                 'pan_gene_10--pan_gene_13': ['pan_gene_11'],
                                 'pan_gene_13--pan_gene_15': ['pan_gene_14'],
                                 'pan_gene_10--pan_gene_15': ['pan_gene_16'],
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_5', 359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_5', 'pan_gene_7', 269, 1, [], ['pan_gene_6']],
            'pan_gene_2--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_7', 478, 2, ['pan_gene_1'], ['pan_gene_8']],
            'pan_gene_10--pan_gene_13--test_double_chromosome': ['test_double_chromosome', 'pan_gene_10', 'pan_gene_13', 359, 2, ['pan_gene_12'], ['pan_gene_11']],
            'pan_gene_13--pan_gene_15--test_double_chromosome': ['test_double_chromosome', 'pan_gene_13', 'pan_gene_15', 269, 1, [], ['pan_gene_14']],
            'pan_gene_10--pan_gene_15--test_double_chromosome': ['test_double_chromosome', 'pan_gene_10', 'pan_gene_15', 478, 2, ['pan_gene_9'], ['pan_gene_16']]
        }

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_two_daft_contigs(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_5'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
            ['gff_name_contig_2', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_6'],
        ]
        core_genes = {'test_double_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7',
                                                 'Core_ID_4': 'pan_gene_10',
                                                 'Core_ID_5': 'pan_gene_13',
                                                 'Core_ID_6': 'pan_gene_15'
                                                 }}
        low_freq_genes = {'test_double_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_double_chromosome.gff'
        acc_genes = {'test_double_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = []

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'Sequence_break--pan_gene_2', 'pan_gene_7--Sequence_break',
                           'pan_gene_10--pan_gene_13', 'pan_gene_13--pan_gene_15',
                           'Sequence_break--pan_gene_10', 'pan_gene_15--Sequence_break']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'Sequence_break--pan_gene_2': 178,
                                   'pan_gene_7--Sequence_break': 300,
                                   'pan_gene_10--pan_gene_13': 359,
                                   'pan_gene_13--pan_gene_15': 269,
                                   'Sequence_break--pan_gene_10': 178,
                                   'pan_gene_15--Sequence_break': 300}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'Sequence_break--pan_gene_2': ['pan_gene_1'],
                                  'pan_gene_7--Sequence_break': [],
                                  'pan_gene_10--pan_gene_13': ['pan_gene_12'],
                                  'pan_gene_13--pan_gene_15': [],
                                  'Sequence_break--pan_gene_10': [],
                                  'pan_gene_15--Sequence_break': []
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'Sequence_break--pan_gene_2': [],
                                 'pan_gene_7--Sequence_break': ['pan_gene_8'],
                                 'pan_gene_10--pan_gene_13': ['pan_gene_11'],
                                 'pan_gene_13--pan_gene_15': ['pan_gene_14'],
                                 'Sequence_break--pan_gene_10': [],
                                 'pan_gene_15--Sequence_break': []
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'Sequence_break--pan_gene_2--test_double_chromosome': ['test_double_chromosome', 'Sequence_break', 'pan_gene_2', 178, 1, ['pan_gene_1'], []],
            'pan_gene_7--Sequence_break--test_double_chromosome': ['test_double_chromosome', 'pan_gene_7', 'Sequence_break', 300, 1, [], ['pan_gene_8']],
            'pan_gene_10--pan_gene_13--test_double_chromosome': ['test_double_chromosome', 'pan_gene_10', 'pan_gene_13',
                                                                 359, 2, ['pan_gene_12'], ['pan_gene_11']],
            'pan_gene_13--pan_gene_15--test_double_chromosome': ['test_double_chromosome', 'pan_gene_13', 'pan_gene_15',
                                                                 269, 1, [], ['pan_gene_14']],
            'Sequence_break--pan_gene_10--test_double_chromosome': ['test_double_chromosome', 'Sequence_break', 'pan_gene_10', 178, 0, [], []],
            'pan_gene_15--Sequence_break--test_double_chromosome': ['test_double_chromosome', 'pan_gene_15', 'Sequence_break', 300, 0, [], []]

        }

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_with_coreless_contig_draft_last_contig(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
        ]
        core_genes = {'test_double_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_double_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_double_chromosome.gff'
        acc_genes = {'test_double_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = []

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'Sequence_break--pan_gene_2', 'pan_gene_7--Sequence_break']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'Sequence_break--pan_gene_2': 178,
                                   'pan_gene_7--Sequence_break': 300}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'Sequence_break--pan_gene_2': ['pan_gene_1'],
                                  'pan_gene_7--Sequence_break': []
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'Sequence_break--pan_gene_2': [],
                                 'pan_gene_7--Sequence_break': ['pan_gene_8']
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'Sequence_break--pan_gene_2--test_double_chromosome': ['test_double_chromosome', 'Sequence_break',
                                                                   'pan_gene_2', 178, 1, ['pan_gene_1'], []],
            'pan_gene_7--Sequence_break--test_double_chromosome': ['test_double_chromosome', 'pan_gene_7',
                                                                   'Sequence_break', 300, 1, [], ['pan_gene_8']]
        }

        coreless_contigs = {'test_double_chromosome--gff_name_contig_2': [['pan_gene_12'], ['pan_gene_11', 'pan_gene_14']]}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_with_coreless_contig_complete_last_contig(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
        ]
        core_genes = {'test_double_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_double_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_double_chromosome.gff'
        acc_genes = {'test_double_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = ['test_double_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'pan_gene_2--pan_gene_7']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1']
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8']
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'pan_gene_2--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_7', 478, 2, ['pan_gene_1'], ['pan_gene_8']]
        }

        coreless_contigs = {
            'test_double_chromosome--gff_name_contig_2': [['pan_gene_12'], ['pan_gene_11', 'pan_gene_14']]}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_with_coreless_contig_draft_first_contig(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_2', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_2', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_2', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_2', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3']
        ]
        core_genes = {'test_double_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_double_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_double_chromosome.gff'
        acc_genes = {'test_double_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = []

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'Sequence_break--pan_gene_2', 'pan_gene_7--Sequence_break']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'Sequence_break--pan_gene_2': 178,
                                   'pan_gene_7--Sequence_break': 300}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'Sequence_break--pan_gene_2': ['pan_gene_1'],
                                  'pan_gene_7--Sequence_break': []
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'Sequence_break--pan_gene_2': [],
                                 'pan_gene_7--Sequence_break': ['pan_gene_8']
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'Sequence_break--pan_gene_2--test_double_chromosome': ['test_double_chromosome', 'Sequence_break',
                                                                   'pan_gene_2', 178, 1, ['pan_gene_1'], []],
            'pan_gene_7--Sequence_break--test_double_chromosome': ['test_double_chromosome', 'pan_gene_7',
                                                                   'Sequence_break', 300, 1, [], ['pan_gene_8']]
        }

        coreless_contigs = {
            'test_double_chromosome--gff_name_contig_1': [['pan_gene_12'], ['pan_gene_11', 'pan_gene_14']]}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_with_coreless_contig_middle_contig(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
            # Contig 3 annotations
            ['gff_name_contig_3', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_4'],
            ['gff_name_contig_3', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_5']
        ]
        core_genes = {'test_triple_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7',
                                                 'Core_ID_4': 'pan_gene_10',
                                                 'Core_ID_5': 'pan_gene_13'
                                                 }}
        low_freq_genes = {'test_triple_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_triple_chromosome.gff'
        acc_genes = {'test_triple_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = []

        # Construct expected results
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'Sequence_break--pan_gene_2', 'pan_gene_7--Sequence_break',
                           'Sequence_break--pan_gene_10', 'pan_gene_10--pan_gene_13', 'pan_gene_13--Sequence_break']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'Sequence_break--pan_gene_2': 178,
                                   'pan_gene_7--Sequence_break': 300,
                                   'Sequence_break--pan_gene_10': 178,
                                   'pan_gene_10--pan_gene_13': 359,
                                   'pan_gene_13--Sequence_break': 620}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'Sequence_break--pan_gene_2': ['pan_gene_1'],
                                  'pan_gene_7--Sequence_break': [],
                                  'Sequence_break--pan_gene_10': [],
                                  'pan_gene_10--pan_gene_13': [],
                                  'pan_gene_13--Sequence_break': []
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'Sequence_break--pan_gene_2': [],
                                 'pan_gene_7--Sequence_break': ['pan_gene_8'],
                                 'Sequence_break--pan_gene_10': [],
                                 'pan_gene_10--pan_gene_13': [],
                                 'pan_gene_13--Sequence_break': []
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'Sequence_break--pan_gene_2--test_triple_chromosome': ['test_triple_chromosome', 'Sequence_break',
                                                                   'pan_gene_2', 178, 1, ['pan_gene_1'], []],
            'pan_gene_7--Sequence_break--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_7',
                                                                   'Sequence_break', 300, 1, [], ['pan_gene_8']],
            'Sequence_break--pan_gene_10--test_triple_chromosome': ['test_triple_chromosome', 'Sequence_break', 'pan_gene_10', 178, 0, [], []],
            'pan_gene_10--pan_gene_13--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_10', 'pan_gene_13', 359, 0, [], []],
            'pan_gene_13--Sequence_break--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_13', 'Sequence_break', 620, 0, [], []]
        }

        coreless_contigs = {
            'test_triple_chromosome--gff_name_contig_2': [['pan_gene_12'], ['pan_gene_11', 'pan_gene_14']]}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)
        # Sort expected and returned lists in dicts
        low_freq_gene_content = {x: sorted(low_freq_gene_content[x]) for x in
                                 low_freq_gene_content.keys()}
        accessory_gene_content = {x: sorted(accessory_gene_content[x]) for x in
                                  accessory_gene_content.keys()}
        master_info = {
            x: [sorted(element) if element is list else element for element in master_info[x]] for x in
            master_info.keys()}

        return_low_freq_gene_content = {x: sorted(return_low_freq_gene_content[x]) for x in
                                        return_low_freq_gene_content.keys()}
        return_accessory_gene_content = {x: sorted(return_accessory_gene_content[x]) for x in
                                         return_accessory_gene_content.keys()}
        return_master_info = {x: [sorted(element) if element is list else element for element in return_master_info[x]]
                              for x in return_master_info.keys()}

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_single_core_on_contig(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
            ['gff_name_contig_2', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_4']
        ]
        core_genes = {'test_double_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7',
                                                 'Core_ID_4': 'pan_gene_15'}}
        low_freq_genes = {'test_double_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_double_chromosome.gff'
        acc_genes = {'test_double_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = ['test_double_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'pan_gene_2--pan_gene_7', 'pan_gene_15--pan_gene_15']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478,
                                   'pan_gene_15--pan_gene_15': 1249}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1'],
                                  'pan_gene_15--pan_gene_15': ['pan_gene_12']
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8'],
                                 'pan_gene_15--pan_gene_15': ['pan_gene_11', 'pan_gene_14']
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'pan_gene_2--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_7',
                                                               478, 2, ['pan_gene_1'], ['pan_gene_8']],
            'pan_gene_15--pan_gene_15--test_double_chromosome': ['test_double_chromosome', 'pan_gene_15', 'pan_gene_15',
                                                                 1249, 3, ['pan_gene_12'], ['pan_gene_11', 'pan_gene_14']]
        }

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_segmentation_of_fragmented_core_gene(self):
        # Set up input
        gff_generator = [['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1_1'],
                         ['gff_name_contig_1', '.', 'CDS', '251', '270', '.', '.', '.', 'Core_ID_1_2'],
                         ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
                         ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
                         ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
                         ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3']]
        core_genes = {'test_single_chromosome': {'Core_ID_1_1': 'pan_gene_2',
                                                 'Core_ID_1_2': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_single_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8'}}
        gff_path = 'TestSegmentingMockGffs/test_single_chromosome.gff'
        acc_genes = {'test_single_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4'}}
        complete_genomes = ['test_single_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7', 'pan_gene_2--pan_gene_7']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 339,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1']}

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8']}

        master_info = {
            'pan_gene_2--pan_gene_5--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               339, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'pan_gene_2--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_7',
                                                               478, 2, ['pan_gene_1'], ['pan_gene_8']]}

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_segmentation_of_fragmented_core_gene_lone_contig(self):
        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
            ['gff_name_contig_2', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_4_1'],
            ['gff_name_contig_2', '.', 'CDS', '10020', '1030', '.', '.', '.', 'Core_ID_4_2']
        ]
        core_genes = {'test_double_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7',
                                                 'Core_ID_4_1': 'pan_gene_15',
                                                 'Core_ID_4_2': 'pan_gene_15'}}
        low_freq_genes = {'test_double_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_double_chromosome.gff'
        acc_genes = {'test_double_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4': 'pan_gene_12'}}
        complete_genomes = ['test_double_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'pan_gene_2--pan_gene_7', 'pan_gene_15--pan_gene_15']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478,
                                   'pan_gene_15--pan_gene_15': 1219}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1'],
                                  'pan_gene_15--pan_gene_15': ['pan_gene_12']
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8'],
                                 'pan_gene_15--pan_gene_15': ['pan_gene_11', 'pan_gene_14']
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'pan_gene_2--pan_gene_7--test_double_chromosome': ['test_double_chromosome', 'pan_gene_2', 'pan_gene_7',
                                                               478, 2, ['pan_gene_1'], ['pan_gene_8']],
            'pan_gene_15--pan_gene_15--test_double_chromosome': ['test_double_chromosome', 'pan_gene_15', 'pan_gene_15',
                                                                 1219, 3, ['pan_gene_12'],
                                                                 ['pan_gene_11', 'pan_gene_14']]
        }

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_segmentation_of_fragmented_acc_gene(self):
        # Set up input
        gff_generator = [['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '179', '270', '.', '.', '.', 'Core_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
                         ['gff_name_contig_1', '.', 'CDS', '450', '460', '.', '.', '.', 'acc_ID_2_1'],
                         ['gff_name_contig_1', '.', 'CDS', '470', '500', '.', '.', '.', 'acc_ID_2_2'],
                         ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
                         ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
                         ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3']]
        core_genes = {'test_single_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_single_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8'}}
        gff_path = 'TestSegmentingMockGffs/test_single_chromosome.gff'
        acc_genes = {'test_single_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2_1': 'pan_gene_4',
                                                'acc_ID_2_2': 'pan_gene_4'}}
        complete_genomes = ['test_single_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7', 'pan_gene_2--pan_gene_7']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 339,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1']}

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8']}

        master_info = {
            'pan_gene_2--pan_gene_5--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               339, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'pan_gene_2--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_7',
                                                               478, 2, ['pan_gene_1'], ['pan_gene_8']]}

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_segmentation_of_fragmented_acc_gene_between_first_n_last_gene(self):

        # Set up input
        gff_generator = [['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
                         ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
                         ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
                         ['gff_name_contig_1', '.', 'CDS', '1100', '1150', '.', '.', '.', 'low_freq_3_1'],
                         ['gff_name_contig_1', '.', 'CDS', '1170', '1250', '.', '.', '.', 'low_freq_3_2']]
        core_genes = {'test_single_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_single_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3_1': 'pan_gene_8',
                                                     'low_freq_3_2': 'pan_gene_8'}}
        gff_path = 'TestSegmentingMockGffs/test_single_chromosome.gff'
        acc_genes = {'test_single_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4'}}
        complete_genomes = ['test_single_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7', 'pan_gene_2--pan_gene_7']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1']}

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8']}

        master_info = {'pan_gene_2--pan_gene_5--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_5', 359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
                       'pan_gene_5--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_5', 'pan_gene_7', 269, 1, [], ['pan_gene_6']],
                       'pan_gene_2--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_7', 478, 2, ['pan_gene_1'], ['pan_gene_8']]}

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                              low_freq_genes, gff_path,
                                                                                              acc_genes,
                                                                                              complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_segmentation_of_fragmented_acc_gene_on_coreless_contig(self):

        # Set up input
        gff_generator = [
            # Contig 1 annotations
            ['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_1'],
            ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
            ['gff_name_contig_1', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
            ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
            ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
            ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3'],
            # Contig 2 annotations
            ['gff_name_contig_2', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_4'],
            ['gff_name_contig_2', '.', 'CDS', '450', '500', '.', '.', '.', 'acc_ID_4_1'],
            ['gff_name_contig_2', '.', 'CDS', '501', '503', '.', '.', '.', 'acc_ID_4_2'],
            ['gff_name_contig_2', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_5'],
            # Contig 3 annotations
            ['gff_name_contig_3', '.', 'CDS', '179', '250', '.', '.', '.', 'Core_ID_4'],
            ['gff_name_contig_3', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_5']
        ]
        core_genes = {'test_triple_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7',
                                                 'Core_ID_4': 'pan_gene_10',
                                                 'Core_ID_5': 'pan_gene_13'
                                                 }}
        low_freq_genes = {'test_triple_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8',
                                                     'low_freq_4': 'pan_gene_11',
                                                     'low_freq_5': 'pan_gene_14'}}
        gff_path = 'TestSegmentingMockGffs/test_triple_chromosome.gff'
        acc_genes = {'test_triple_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2': 'pan_gene_4',
                                                'acc_ID_4_1': 'pan_gene_12',
                                                'acc_ID_4_2': 'pan_gene_12'}}
        complete_genomes = []

        # Construct expected results
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7',
                           'Sequence_break--pan_gene_2', 'pan_gene_7--Sequence_break',
                           'Sequence_break--pan_gene_10', 'pan_gene_10--pan_gene_13', 'pan_gene_13--Sequence_break']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 359,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'Sequence_break--pan_gene_2': 178,
                                   'pan_gene_7--Sequence_break': 300,
                                   'Sequence_break--pan_gene_10': 178,
                                   'pan_gene_10--pan_gene_13': 359,
                                   'pan_gene_13--Sequence_break': 620}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': [],
                                  'Sequence_break--pan_gene_2': ['pan_gene_1'],
                                  'pan_gene_7--Sequence_break': [],
                                  'Sequence_break--pan_gene_10': [],
                                  'pan_gene_10--pan_gene_13': [],
                                  'pan_gene_13--Sequence_break': []
                                  }

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'Sequence_break--pan_gene_2': [],
                                 'pan_gene_7--Sequence_break': ['pan_gene_8'],
                                 'Sequence_break--pan_gene_10': [],
                                 'pan_gene_10--pan_gene_13': [],
                                 'pan_gene_13--Sequence_break': []
                                 }

        master_info = {
            'pan_gene_2--pan_gene_5--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               359, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 1, [], ['pan_gene_6']],
            'Sequence_break--pan_gene_2--test_triple_chromosome': ['test_triple_chromosome', 'Sequence_break',
                                                                   'pan_gene_2', 178, 1, ['pan_gene_1'], []],
            'pan_gene_7--Sequence_break--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_7',
                                                                   'Sequence_break', 300, 1, [], ['pan_gene_8']],
            'Sequence_break--pan_gene_10--test_triple_chromosome': ['test_triple_chromosome', 'Sequence_break', 'pan_gene_10', 178, 0, [], []],
            'pan_gene_10--pan_gene_13--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_10', 'pan_gene_13', 359, 0, [], []],
            'pan_gene_13--Sequence_break--test_triple_chromosome': ['test_triple_chromosome', 'pan_gene_13', 'Sequence_break', 620, 0, [], []]
        }

        coreless_contigs = {
            'test_triple_chromosome--gff_name_contig_2': [['pan_gene_12'], ['pan_gene_11', 'pan_gene_14']]}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs.sort(), return_core_gene_pairs.sort())
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_single_fragmented_gene_on_either_side_of_core_gene(self):

        # Set up input
        gff_generator = [['gff_name_contig_1', '.', 'CDS', '90', '180', '.', '.', '.', 'acc_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '179', '270', '.', '.', '.', 'Core_ID_1'],
                         ['gff_name_contig_1', '.', 'CDS', '300', '425', '.', '.', '.', 'low_freq_1'],
                         ['gff_name_contig_1', '.', 'CDS', '450', '460', '.', '.', '.', 'acc_ID_2_1'],
                         ['gff_name_contig_1', '.', 'CDS', '610', '680', '.', '.', '.', 'Core_ID_2'],
                         ['gff_name_contig_1', '.', 'CDS', '685', '690', '.', '.', '.', 'acc_ID_2_2'],
                         ['gff_name_contig_1', '.', 'CDS', '700', '850', '.', '.', '.', 'low_freq_2'],
                         ['gff_name_contig_1', '.', 'CDS', '950', '1000', '.', '.', '.', 'Core_ID_3'],
                         ['gff_name_contig_1', '.', 'CDS', '1100', '1250', '.', '.', '.', 'low_freq_3']]
        core_genes = {'test_single_chromosome': {'Core_ID_1': 'pan_gene_2',
                                                 'Core_ID_2': 'pan_gene_5',
                                                 'Core_ID_3': 'pan_gene_7'}}
        low_freq_genes = {'test_single_chromosome': {'low_freq_1': 'pan_gene_3',
                                                     'low_freq_2': 'pan_gene_6',
                                                     'low_freq_3': 'pan_gene_8'}}
        gff_path = 'TestSegmentingMockGffs/test_single_chromosome.gff'
        acc_genes = {'test_single_chromosome': {'acc_ID_1': 'pan_gene_1',
                                                'acc_ID_2_1': 'pan_gene_4',
                                                'acc_ID_2_2': 'pan_gene_4'}}
        complete_genomes = ['test_single_chromosome']

        # Set up expected outputs
        core_gene_pairs = ['pan_gene_2--pan_gene_5', 'pan_gene_5--pan_gene_7', 'pan_gene_2--pan_gene_7']
        core_gene_pair_distance = {'pan_gene_2--pan_gene_5': 339,
                                   'pan_gene_5--pan_gene_7': 269,
                                   'pan_gene_2--pan_gene_7': 478}

        accessory_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_4'],
                                  'pan_gene_5--pan_gene_7': ['pan_gene_4'],
                                  'pan_gene_2--pan_gene_7': ['pan_gene_1']}

        low_freq_gene_content = {'pan_gene_2--pan_gene_5': ['pan_gene_3'],
                                 'pan_gene_5--pan_gene_7': ['pan_gene_6'],
                                 'pan_gene_2--pan_gene_7': ['pan_gene_8']}

        master_info = {
            'pan_gene_2--pan_gene_5--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_5',
                                                               339, 2, ['pan_gene_4'], ['pan_gene_3'], ],
            'pan_gene_5--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_5', 'pan_gene_7',
                                                               269, 2, ['pan_gene_4'], ['pan_gene_6']],
            'pan_gene_2--pan_gene_7--test_single_chromosome': ['test_single_chromosome', 'pan_gene_2', 'pan_gene_7',
                                                               478, 2, ['pan_gene_1'], ['pan_gene_8']]}

        coreless_contigs = {}

        # Run function
        return_core_gene_pairs, return_core_gene_pair_distance, \
        return_accessory_gene_content, return_low_freq_gene_content, \
        return_master_info, return_coreless_contigs = gff_parser.segment_gff_content(gff_generator, core_genes,
                                                                                     low_freq_genes, gff_path,
                                                                                     acc_genes,
                                                                                     complete_genomes)

        # Evaluate
        self.assertEqual(core_gene_pairs, return_core_gene_pairs)
        self.assertEqual(core_gene_pair_distance, return_core_gene_pair_distance)
        self.assertEqual(accessory_gene_content, return_accessory_gene_content)
        self.assertEqual(low_freq_gene_content, return_low_freq_gene_content)
        self.assertEqual(master_info, return_master_info)
        self.assertEqual(coreless_contigs, return_coreless_contigs)

    def test_something(self): # TODO - What other wired and wonderfull examples can we come up with?
        pass


class TestMergingDicts(unittest.TestCase):
    """ Functions to merge dictionaries and lists into dictionaries """
    # Test merge_dicts_counts
    def test_merge_dicts_counts_list_empty(self):
        input_dict = {}
        input_list = ['x', 'y', 'z']

        expected_dict = {'x': 1,
                         'y': 1,
                         'z': 1}

        return_dict = merge_dicts.merge_dicts_counts(input_dict, input_list)

        self.assertEqual(return_dict, expected_dict)

    def test_merge_dicts_counts_dict_empty(self):
        input_dict = {}
        input_list = {'x': 2, 'y': 2, 'z': 2}

        expected_dict = {'x': 1,
                         'y': 1,
                         'z': 1}

        return_dict = merge_dicts.merge_dicts_counts(input_dict, input_list)

        self.assertEqual(return_dict, expected_dict)

    def test_merge_dicts_counts_list_adding(self):
        input_dict = {'x': 1,
                      'y': 1,
                      'z': 1}
        input_list = ['x', 'y', 'z']

        expected_dict = {'x': 2,
                         'y': 2,
                         'z': 2}

        return_dict = merge_dicts.merge_dicts_counts(input_dict, input_list)

        self.assertEqual(return_dict, expected_dict)

    def test_merge_dicts_counts_dict_adding(self):
        input_dict = {'x': 1,
                      'y': 1,
                      'z': 1}
        input_list = {'x': 1, 'y': 1, 'z': 1}

        expected_dict = {'x': 2,
                         'y': 2,
                         'z': 2}

        return_dict = merge_dicts.merge_dicts_counts(input_dict, input_list)

        self.assertEqual(return_dict, expected_dict)

    def test_merge_dicts_counts_dict_mix(self):
        input_dict = {'x': 1,
                      'y': 1,
                      'z': 1}
        input_list = {'x': 1, 'y': 1}

        expected_dict = {'x': 2,
                         'y': 2,
                         'z': 1}

        return_dict = merge_dicts.merge_dicts_counts(input_dict, input_list)

        self.assertEqual(return_dict, expected_dict)

    def test_merge_dicts_counts_list_mix(self):
        input_dict = {'x': 1,
                      'y': 1}
        input_list = ['x', 'y', 'z']

        expected_dict = {'x': 2,
                         'y': 2,
                         'z': 1}

        return_dict = merge_dicts.merge_dicts_counts(input_dict, input_list)

        self.assertEqual(return_dict, expected_dict)

    # Test merge_dicts_lists
    def test_merge_dicts_lists_empty(self):
        input_dict = {}
        merge_dict = {'x': ['test_3'],
                      'y': ['test_2'],
                      'z': ['test_1']}

        expected_dict = {'x': ['test_3'],
                         'y': ['test_2'],
                         'z': ['test_1']}

        return_dict = merge_dicts.merge_dicts_lists(input_dict, merge_dict)

        self.assertEqual(expected_dict, return_dict)

    def test_merge_dicts_lists_adding(self):
        input_dict = {'x': ['init_3'],
                      'y': ['init_2'],
                      'z': ['init_1']}

        merge_dict = {'x': ['test_3'],
                      'y': ['test_2'],
                      'z': ['test_1']}

        expected_dict = {'x': ['init_3', 'test_3'],
                         'y': ['init_2', 'test_2'],
                         'z': ['init_1', 'test_1']}

        return_dict = merge_dicts.merge_dicts_lists(input_dict, merge_dict)

        self.assertEqual(expected_dict, return_dict)

    def test_merge_dicts_lists_mix(self):
        input_dict = {'x': ['init_3'],
                      'y': ['init_2']}

        merge_dict = {'x': ['test_3'],
                      'y': ['test_2'],
                      'z': ['test_1']}

        expected_dict = {'x': ['init_3', 'test_3'],
                         'y': ['init_2', 'test_2'],
                         'z': ['test_1']}

        return_dict = merge_dicts.merge_dicts_lists(input_dict, merge_dict)

        self.assertEqual(expected_dict, return_dict)


class TestCoreGraphConstruction(unittest.TestCase):
    """
    Test the construction of a network made from core gene pairs and their number of connections.
    """
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

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)

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

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)

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

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)

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

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)

        # Get edge weights:
        edge_weights = [core_graph.get_edge_data(edge[0], edge[1])['weight'] for edge in list(core_graph.edges)]

        # Assert outputs
        self.assertEqual(expected_edges, list(core_graph.edges))
        self.assertEqual(expected_degrees, list(core_graph.degree))
        self.assertEqual(expected_edge_weights, edge_weights)


class TestGeneCoOccurrence(unittest.TestCase):
    """
    Test function that identifies the number of genomes in which two core genes co-occur.
    """
    def test_count_gene_co_occurrence(self):
        core_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-1': "A",
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

        segment = ["B", "A"]

        expected_value = 10
        a_occurrence = 10
        b_occurrence = 10

        return_value, individual_occurrences = consesus_core_genome.count_gene_co_occurrence(core_gene_dict, segment)

        self.assertEqual(expected_value, return_value)
        self.assertEqual(a_occurrence, individual_occurrences["A"])
        self.assertEqual(b_occurrence, individual_occurrences["B"])

    def test_count_gene_co_occurrence_no_occurence(self):
        core_gene_dict = {'Silas_the_Salmonella': {'Silas_the_Salmonella_tag-1-2.1': "B",
                                                   'Silas_the_Salmonella_tag-1-2.2': "B"},
                          'Christina_the_Streptococcus': {'Christina_the_Streptococcus_tag-2-1': "A"},
                          'Ajwa_the_Shigella': {'Ajwa_the_Shigella_tag-3-2': "B"},
                          'Ajwa_the_Legionella': {'Ajwa_the_Legionella_tag-4-2': "B"},
                          'Cari_the_Listeria': {'Cari_the_Listeria_tag-5-1': "A"},
                          'Aman_the_Streptococcus': {'Aman_the_Streptococcus_tag-6-2': "B"},
                          'Zion_the_Streptococcus': {'Zion_the_Streptococcus_tag-7-1': "A"},
                          'Dina_the_Shigella': {'Dina_the_Shigella_tag-8-1': "A"},
                          'Silas_the_Legionella': {'Silas_the_Legionella_tag-9-2': "B"},
                          'Lilly_the_Shigella': {'Lilly_the_Shigella_tag-10-1': "A"}}

        segment = ["B", "A"]

        expected_value = 0
        a_occurrence = 5
        b_occurrence = 5

        return_value, individual_occurrences = consesus_core_genome.count_gene_co_occurrence(core_gene_dict, segment)

        self.assertEqual(expected_value, return_value)
        self.assertEqual(a_occurrence, individual_occurrences["A"])
        self.assertEqual(b_occurrence, individual_occurrences["B"])


class TestSegmentationIdentification(unittest.TestCase):
    """
    Test the function that identifies core gene segments from a pan-genome.
    """
    def test_double_edge_segment_identification_all_2_degree_input(self):
        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 10,
                                'pan_cluster_2--pan_cluster_3': 10,
                                'pan_cluster_3--pan_cluster_4': 10,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_6--pan_cluster_1': 10}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        return_1 = consesus_core_genome.identify_segments(core_graph, 10, {}, num_components)

        self.assertEqual(None, return_1)

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

        core_gene_dict = {'genome_1': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_2': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_3': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_4': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_5': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_6': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_7': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_8': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_9': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},
                          'genome_10': {'tag_1': 'pan_cluster_1', 'tag_2': 'pan_cluster_4', 'tag_3': 'pan_cluster_2', 'tag_4': 'pan_cluster_5'},}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 10, core_gene_dict, num_components)

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

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 10, {}, num_components)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_four_degrees(self):
        # expected_segments = {'pan_cluster_4--pan_cluster_6': ['pan_cluster_4', 'pan_cluster_5', 'pan_cluster_6']}
        expected_segments = {'pan_cluster_2--pan_cluster_4': ['pan_cluster_2',
                                                              'pan_cluster_3',
                                                              'pan_cluster_4'],
                             'pan_cluster_2--pan_cluster_6': ['pan_cluster_2',
                                                              'pan_cluster_1',
                                                              'pan_cluster_6'],
                             'pan_cluster_4--pan_cluster_6': ['pan_cluster_4',
                                                              'pan_cluster_5',
                                                              'pan_cluster_6']}

        core_neighbour_pairs = {'pan_cluster_1--pan_cluster_2': 9,
                                'pan_cluster_2--pan_cluster_3': 9,
                                'pan_cluster_2--pan_cluster_4': 1,
                                'pan_cluster_2--pan_cluster_6': 1,
                                'pan_cluster_3--pan_cluster_4': 9,
                                'pan_cluster_4--pan_cluster_5': 10,
                                'pan_cluster_5--pan_cluster_6': 10,
                                'pan_cluster_6--pan_cluster_1': 9}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 10, {}, num_components)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_challenging_paths(self):
        expected_segments = {'pan_cluster_A--pan_cluster_B': ['pan_cluster_A', 'pan_cluster_E', 'pan_cluster_F', 'pan_cluster_G', 'pan_cluster_B'],
                             'pan_cluster_B--pan_cluster_C': ['pan_cluster_C', 'pan_cluster_B']}

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
        core_gene_dict = {'genome_1': {'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_2': {'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_3': {'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_4': {'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_5': {'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', }}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 5, core_gene_dict, num_components)

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
        core_gene_dict = {'genome_1': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_2': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_3': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_4': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_5': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_6': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_7': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', },
                          'genome_8': {'tag_5': 'pan_cluster_E', 'tag_4': 'pan_cluster_D', 'tag_3': 'pan_cluster_C', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', }}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 8, core_gene_dict, num_components)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_all_challenging_paths(self):
        expected_segments = {'pan_cluster_A--pan_cluster_D': ['pan_cluster_A', 'pan_cluster_G', 'pan_cluster_F', 'pan_cluster_E', 'pan_cluster_D'],
                             'pan_cluster_B--pan_cluster_C': ['pan_cluster_B', 'pan_cluster_H', 'pan_cluster_I', 'pan_cluster_J', 'pan_cluster_C']}#,}

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
        core_gene_dict = {'genome_1': {'tag_5': 'pan_cluster_K', 'tag_4': 'pan_cluster_L', 'tag_3': 'pan_cluster_A', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', 'tag_6': 'pan_cluster_C', 'tag_7': 'pan_cluster_D'},
                          'genome_2': {'tag_5': 'pan_cluster_K', 'tag_4': 'pan_cluster_L', 'tag_3': 'pan_cluster_A', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', 'tag_6': 'pan_cluster_C', 'tag_7': 'pan_cluster_D'},
                          'genome_3': {'tag_5': 'pan_cluster_K', 'tag_4': 'pan_cluster_L', 'tag_3': 'pan_cluster_A', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', 'tag_6': 'pan_cluster_C', 'tag_7': 'pan_cluster_D'},
                          'genome_4': {'tag_5': 'pan_cluster_K', 'tag_4': 'pan_cluster_L', 'tag_3': 'pan_cluster_A', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', 'tag_6': 'pan_cluster_C', 'tag_7': 'pan_cluster_D'},
                          'genome_5': {'tag_5': 'pan_cluster_K', 'tag_4': 'pan_cluster_L', 'tag_3': 'pan_cluster_A', 'tag_2': 'pan_cluster_B', 'tag_1': 'pan_cluster_A', 'tag_6': 'pan_cluster_C', 'tag_7': 'pan_cluster_D'}}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 5, core_gene_dict, num_components)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_less_than_all_present(self):
        expected_segments = {'pan_cluster_B--pan_cluster_D': ['pan_cluster_B', 'pan_cluster_C', 'pan_cluster_D'],
                             'pan_cluster_F--pan_cluster_H': ['pan_cluster_H', 'pan_cluster_G', 'pan_cluster_F'],
                             'pan_cluster_B--pan_cluster_H': ['pan_cluster_B', 'pan_cluster_A', 'pan_cluster_H'],
                             'pan_cluster_D--pan_cluster_F': ['pan_cluster_D', 'pan_cluster_E', 'pan_cluster_F']}

        core_neighbour_pairs = {'pan_cluster_A--pan_cluster_B': 9,
                                'pan_cluster_A--pan_cluster_H': 9,
                                'pan_cluster_B--pan_cluster_H': 1,
                                'pan_cluster_B--pan_cluster_C': 10,
                                'pan_cluster_C--pan_cluster_D': 10,
                                'pan_cluster_D--pan_cluster_E': 9,
                                'pan_cluster_D--pan_cluster_F': 1,
                                'pan_cluster_E--pan_cluster_F': 9,
                                'pan_cluster_F--pan_cluster_G': 10,
                                'pan_cluster_G--pan_cluster_H': 10,
                                'pan_cluster_H--pan_cluster_A': 9,
                                }

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 10, {}, num_components)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_double_edge_segment_identification_segments_node_w_two_gene_segment(self):
        expected_segments = {'pan_cluster_A--pan_cluster_B': ['pan_cluster_A', 'pan_cluster_B'],
                             'pan_cluster_A--pan_cluster_G': ['pan_cluster_A', 'pan_cluster_I', 'pan_cluster_H', 'pan_cluster_G'],
                             'pan_cluster_B--pan_cluster_E': ['pan_cluster_B', 'pan_cluster_C', 'pan_cluster_D', 'pan_cluster_E'],
                             'pan_cluster_E--pan_cluster_G': ['pan_cluster_G', 'pan_cluster_F', 'pan_cluster_E']}

        core_neighbour_pairs = {'pan_cluster_A--pan_cluster_B': 3,
                                'pan_cluster_A--pan_cluster_I': 2,
                                'pan_cluster_A--pan_cluster_G': 1,
                                'pan_cluster_B--pan_cluster_C': 2,
                                'pan_cluster_B--pan_cluster_E': 1,
                                'pan_cluster_C--pan_cluster_D': 3,
                                'pan_cluster_D--pan_cluster_E': 3,
                                'pan_cluster_E--pan_cluster_F': 2,
                                'pan_cluster_F--pan_cluster_G': 2,
                                'pan_cluster_G--pan_cluster_H': 3,
                                'pan_cluster_H--pan_cluster_I': 3
                                }
        core_gene_dict = {'genome_1': {'gene_1': 'pan_cluster_A', 'gene_2': 'pan_cluster_B', 'gene_3': 'pan_cluster_E', 'gene_4': 'pan_cluster_G', 'gene_5': 'pan_cluster_D', 'gene_7': 'pan_cluster_H'},
                          'genome_2': {'gene_1': 'pan_cluster_A', 'gene_2': 'pan_cluster_B', 'gene_3': 'pan_cluster_E', 'gene_4': 'pan_cluster_G', 'gene_5': 'pan_cluster_D', 'gene_7': 'pan_cluster_H'},
                          'genome_3': {'gene_1': 'pan_cluster_A', 'gene_2': 'pan_cluster_B', 'gene_3': 'pan_cluster_E', 'gene_4': 'pan_cluster_G', 'gene_5': 'pan_cluster_D', 'gene_7': 'pan_cluster_H'}}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = consesus_core_genome.identify_segments(core_graph, 3, core_gene_dict, num_components)

        self.assertEqual(expected_segments, double_edge_segements)

    def test_multiple_component_core_graph(self):
        expected_segments = {'pan_cluster_A--pan_cluster_I': ['pan_cluster_A', 'pan_cluster_I'],
                             'pan_cluster_B--pan_cluster_C': ['pan_cluster_C', 'pan_cluster_B'],
                             'pan_cluster_D--pan_cluster_J': ['pan_cluster_D', 'pan_cluster_J'],
                             'pan_cluster_E--pan_cluster_K': ['pan_cluster_E', 'pan_cluster_K'],
                             'pan_cluster_F--pan_cluster_G': ['pan_cluster_G', 'pan_cluster_F'],
                             'pan_cluster_H--pan_cluster_M': ['pan_cluster_H', 'pan_cluster_L', 'pan_cluster_M'],
                             'pan_cluster_Q--pan_cluster_O': ['pan_cluster_Q', 'pan_cluster_P', 'pan_cluster_O']}

        core_neighbour_pairs = {'pan_cluster_A--pan_cluster_B': 1,
                                'pan_cluster_A--pan_cluster_C': 1,
                                'pan_cluster_A--pan_cluster_I': 2,
                                'pan_cluster_B--pan_cluster_C': 2,
                                'pan_cluster_B--pan_cluster_D': 1,
                                'pan_cluster_C--pan_cluster_D': 1,
                                'pan_cluster_D--pan_cluster_J': 2,
                                'pan_cluster_E--pan_cluster_F': 1,
                                'pan_cluster_E--pan_cluster_G': 1,
                                'pan_cluster_E--pan_cluster_K': 2,
                                'pan_cluster_F--pan_cluster_G': 2,
                                'pan_cluster_F--pan_cluster_H': 1,
                                'pan_cluster_G--pan_cluster_H': 1,
                                'pan_cluster_H--pan_cluster_L': 2,
                                'pan_cluster_L--pan_cluster_M': 2,
                                'pan_cluster_O--pan_cluster_P': 2,
                                'pan_cluster_P--pan_cluster_Q': 2,
                                }

        core_gene_dict = {'genome_1': {'tag_1': 'pan_cluster_A', 'tag_2': 'pan_cluster_B', 'tag_3': 'pan_cluster_C',
                                       'tag_4': 'pan_cluster_D', 'tag_5': 'pan_cluster_E', 'tag_6': 'pan_cluster_F',
                                       'tag_7': 'pan_cluster_G', 'tag_8': 'pan_cluster_H', 'tag_9': 'pan_cluster_I',
                                       'tag_10': 'pan_cluster_J', 'tag_11': 'pan_cluster_K', 'tag_12': 'pan_cluster_L',
                                       'tag_13': 'pan_cluster_M', 'tag_14': 'pan_cluster_O', 'tag_15': 'pan_cluster_P',
                                       'tag_16': 'pan_cluster_Q'},
                          'genome_2': {'tag_1': 'pan_cluster_A', 'tag_2': 'pan_cluster_B', 'tag_3': 'pan_cluster_C',
                                       'tag_4': 'pan_cluster_D', 'tag_5': 'pan_cluster_E', 'tag_6': 'pan_cluster_F',
                                       'tag_7': 'pan_cluster_G', 'tag_8': 'pan_cluster_H', 'tag_9': 'pan_cluster_I',
                                       'tag_10': 'pan_cluster_J', 'tag_11': 'pan_cluster_K', 'tag_12': 'pan_cluster_L',
                                       'tag_13': 'pan_cluster_M', 'tag_14': 'pan_cluster_O', 'tag_15': 'pan_cluster_P',
                                       'tag_16': 'pan_cluster_Q'}}

        core_graph = consesus_core_genome.construct_core_graph(core_neighbour_pairs)
        num_components = number_connected_components(core_graph)

        double_edge_segements = {}
        for component in connected_components(core_graph):
            component_graph = core_graph.subgraph(component).copy()
            double_edge_segements = double_edge_segements | consesus_core_genome.identify_segments(component_graph, 2,
                                                                                                   core_gene_dict,
                                                                                                   num_components)

        # comparisons = [True for x in double_edge_segements
        #                if
        #                (x in expected_segments and
        #                (expected_segments[x] == double_edge_segements[x] or expected_segments[x][::-1] == double_edge_segements[x]))
        #                or
        #                (f"{x.split('--')[1]}'--'{x.split('--')[0]}" in expected_segments and
        #                (expected_segments[x] == double_edge_segements[f"{x.split('--')[1]}'--'{x.split('--')[0]}"] or expected_segments[x][::-1] == double_edge_segements[f"{x.split('--')[1]}'--'{x.split('--')[0]}"]))
        #                ]
        key_forward = [x for x in double_edge_segements if x in expected_segments]
        key_reverse = [f"{x.split('--')[1]}--{x.split('--')[0]}" for x in double_edge_segements if f"{x.split('--')[1]}--{x.split('--')[0]}" in expected_segments]
        expected_key_match = key_forward+key_reverse

        # Test if the number of expected segments were returned
        self.assertEqual(len(expected_key_match), len(expected_segments))

        comparisons = [True for returned_key, expected_key in zip(double_edge_segements, expected_key_match)
                       if double_edge_segements[returned_key] == expected_segments[expected_key]
                       or
                       double_edge_segements[returned_key] == expected_segments[expected_key][::-1]]

        # Test of all returned segments look as expected
        self.assertTrue(all(comparisons))

    # TODO - Chat to Andrew about this function how it works and how we can test it more - possibly just run some things to see if it breaks


class TestNoAccessorySegmentIdentifcation(unittest.TestCase):
    """
    Test the function that takes in segments of core genes and divide them into sub-segments based on the accessory content between core genes in segment.
    """
    def test_no_accessory_genes_in_segment(self):
        expected_sub_sgments = {'pan_cluster_1--pan_cluster_5': [['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5']],
                                 'pan_cluster_2--pan_cluster_4': [['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']]}

        double_edge_segements = {'pan_cluster_1--pan_cluster_5': ['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5'],
                                 'pan_cluster_2--pan_cluster_4': ['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_6': 0, 'pan_cluster_5--pan_cluster_6': 0, 'pan_cluster_2--pan_cluster_3': 0, 'pan_cluster_3--pan_cluster_4': 0}

        sub_segment_dict = consesus_core_genome.identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

    def test_accessory_genes_in_segment_first_gene_lonely(self):
        expected_sub_sgments = {'pan_cluster_1--pan_cluster_5': [['pan_cluster_1'], ['pan_cluster_6', 'pan_cluster_5']]}

        double_edge_segements = {'pan_cluster_1--pan_cluster_5': ['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_6': 1, 'pan_cluster_5--pan_cluster_6': 0}

        sub_segment_dict = consesus_core_genome.identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

    def test_accessory_genes_in_segment_last_gene_lonely(self):
        expected_sub_sgments = {'pan_cluster_1--pan_cluster_5': [['pan_cluster_1', 'pan_cluster_6'], ['pan_cluster_5']],
                                'pan_cluster_2--pan_cluster_4': [['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']]}

        double_edge_segements = {'pan_cluster_1--pan_cluster_5': ['pan_cluster_1', 'pan_cluster_6', 'pan_cluster_5'],
                                 'pan_cluster_2--pan_cluster_4': ['pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_6': 0, 'pan_cluster_5--pan_cluster_6': 1, 'pan_cluster_2--pan_cluster_3': 0, 'pan_cluster_3--pan_cluster_4': 0}

        sub_segment_dict = consesus_core_genome.identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

    def test_accessory_genes_in_segment_middle(self):
        expected_sub_sgments = {'pan_cluster_1--pan_cluster_4': [['pan_cluster_1', 'pan_cluster_2'], ['pan_cluster_3', 'pan_cluster_4']]}

        double_edge_segements = {'pan_cluster_1--pan_cluster_4': ['pan_cluster_1', 'pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_2': 0, 'pan_cluster_2--pan_cluster_3': 1, 'pan_cluster_3--pan_cluster_4': 0}

        sub_segment_dict = consesus_core_genome.identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)

    def test_accessory_genes_in_multiple_places(self):
        expected_sub_sgments = {'pan_cluster_1--pan_cluster_4': [['pan_cluster_1'], ['pan_cluster_2', 'pan_cluster_3'], ['pan_cluster_4']]}

        double_edge_segements = {'pan_cluster_1--pan_cluster_4': ['pan_cluster_1', 'pan_cluster_2', 'pan_cluster_3', 'pan_cluster_4']}
        combined_acc_gene_count = {'pan_cluster_1--pan_cluster_2': 1, 'pan_cluster_2--pan_cluster_3': 0, 'pan_cluster_3--pan_cluster_4': 1}

        sub_segment_dict = consesus_core_genome.identify_no_accessory_segments(double_edge_segements, combined_acc_gene_count)

        self.assertEqual(sub_segment_dict, expected_sub_sgments)


class TestSummaryTableConstruction(unittest.TestCase):
    def test_summary_table_calculations(self):
        master_info = {
            'pan_cluster_1--pan_cluster_2--genome_1': ['genome_1', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_1--pan_cluster_2--genome_2': ['genome_2', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_1--pan_cluster_2--genome_3': ['genome_3', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_2--pan_cluster_3--genome_1': ['genome_1', 'pan_cluster_2', 'pan_cluster_3', 100, 2,
                                                       ['Acc_1', 'Acc_2'], []],
            'pan_cluster_2--pan_cluster_4--genome_2': ['genome_2', 'pan_cluster_2', 'pan_cluster_3', 150, 1,
                                                       ['Acc_1', ], []],
            'pan_cluster_2--pan_cluster_3--genome_3': ['genome_3', 'pan_cluster_2', 'pan_cluster_3', 200, 0,
                                                       [], []],
            'pan_cluster_3--pan_cluster_4--genome_1': ['genome_1', 'pan_cluster_2', 'pan_cluster_3', -5, 0,
                                                       [], []],
            'pan_cluster_3--pan_cluster_4--genome_3': ['genome_3', 'pan_cluster_2', 'pan_cluster_3', -10, 0,
                                                       [], []]
        }

        core_gene_dict = {'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2', 'gene_3': 'pan_cluster_3', 'gene_4': 'pan_cluster_4'},
                          'genome_2': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2', 'gene_4': 'pan_cluster_4'},
                          'genome_3': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2', 'gene_3': 'pan_cluster_3', 'gene_4': 'pan_cluster_4'}}

        expected_table = {'pan_cluster_1--pan_cluster_2': ['pan_cluster_1-pan_cluster_2', 3, 3, 3, 3, 99, 99, 99.0, 99.0, 3, 3, 3.0, 3.0],
                          'pan_cluster_2--pan_cluster_3': ['pan_cluster_2-pan_cluster_3', 2, 3, 2, 2, 100, 200, 150.0, 150.0, 0, 2, 1.0, 1.0],
                          'pan_cluster_2--pan_cluster_4': ['pan_cluster_2-pan_cluster_4', 1, 3, 3, 3, 150, 150, 150.0, 150.0, 1, 1, 1.0, 1.0],
                          'pan_cluster_3--pan_cluster_4': ['pan_cluster_3-pan_cluster_4', 2, 2, 3, 2, -10, -5, -7.5, -7.5, 0, 0, 0.0, 0.0]}

        return_table = summary_table.calculate_n_create_summaries(master_info, core_gene_dict)

        self.assertEqual(expected_table, return_table)

    def test_summary_table_calculations_w_sequence_breaks(self):
        master_info = {
            'pan_cluster_1--pan_cluster_2--genome_1': ['genome_1', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_1--Sequence_break--genome_2': ['genome_2', 'pan_cluster_1', 'Sequence_break', 100, 0,
                                                       [], []],
            'Sequence_break--pan_cluster_2--genome_2': ['genome_2', 'Sequence_break', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_1--pan_cluster_2--genome_3': ['genome_3', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_2--pan_cluster_3--genome_1': ['genome_1', 'pan_cluster_2', 'pan_cluster_3', 100, 2,
                                                       ['Acc_1', 'Acc_2'], []],
            'pan_cluster_2--pan_cluster_3--genome_3': ['genome_3', 'pan_cluster_2', 'pan_cluster_3', 200, 0,
                                                       [], []],
            'pan_cluster_3--Sequence_break--genome_1': ['genome_1', 'pan_cluster_3', 'Sequence_break', 100, 2,
                                                       ['Acc_1', 'Acc_2'], []],
            'pan_cluster_3--Sequence_break--genome_3': ['genome_3', 'pan_cluster_3', 'Sequence_break', 200, 0,
                                                       [], []]
        }

        core_gene_dict = {'genome_1': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2', 'gene_3': 'pan_cluster_3'},
                          'genome_2': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2'},
                          'genome_3': {'gene_1': 'pan_cluster_1', 'gene_2': 'pan_cluster_2', 'gene_3': 'pan_cluster_3'}}

        expected_table = {'pan_cluster_1--pan_cluster_2': ['pan_cluster_1-pan_cluster_2', 2, 3, 3, 3, 99, 99, 99.0, 99.0, 3, 3, 3.0, 3.0],
                          'pan_cluster_1--Sequence_break': ['pan_cluster_1-Sequence_break', 1, 3, 0, 0, 100, 100, 100.0, 100.0, 0, 0, 0.0, 0.0],
                          'Sequence_break--pan_cluster_2': ['Sequence_break-pan_cluster_2', 1, 0, 3, 0, 99, 99, 99.0, 99.0, 3, 3, 3.0, 3.0],
                          'pan_cluster_2--pan_cluster_3': ['pan_cluster_2-pan_cluster_3', 2, 3, 2, 2, 100, 200, 150.0, 150.0, 0, 2, 1.0, 1.0],
                          'pan_cluster_3--Sequence_break': ['pan_cluster_3-Sequence_break', 2, 2, 0, 0, 100, 200, 150.0, 150.0, 0, 2, 1.0, 1.0]}

        return_table = summary_table.calculate_n_create_summaries(master_info, core_gene_dict)

        self.assertEqual(expected_table, return_table)


class TestWritingOutputFunction(unittest.TestCase):
    """
    Function to test the creation of output files
    """
    def tearDown(self):
        """ Class to remove created database files of gff files in tmp-folder"""
        for file in os.listdir('TestWritingOutputFunction'):
            if "sv" in file:
                db_path = os.path.join('TestWritingOutputFunction', file)
                os.remove(db_path)

    def test_master_info_writer(self):
        master_info = {
            'pan_cluster_1--pan_cluster_2--genome_1': ['genome_1', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_1--pan_cluster_2--genome_2': ['genome_2', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_1--pan_cluster_2--genome_3': ['genome_3', 'pan_cluster_1', 'pan_cluster_2', 99, 3,
                                                       ['Acc_1', 'Acc_2'], ['low_1']],
            'pan_cluster_2--pan_cluster_3--genome_1': ['genome_1', 'pan_cluster_2', 'pan_cluster_3', 100, 2,
                                                       ['Acc_1', 'Acc_2'], []],
            'pan_cluster_2--pan_cluster_4--genome_2': ['genome_2', 'pan_cluster_2', 'pan_cluster_4', 150, 1,
                                                       ['Acc_1', ], []],
            'pan_cluster_2--pan_cluster_3--genome_3': ['genome_3', 'pan_cluster_2', 'pan_cluster_3', 200, 0,
                                                       [], []],
            'pan_cluster_3--pan_cluster_4--genome_1': ['genome_1', 'pan_cluster_3', 'pan_cluster_4', -5, 0,
                                                       [], []],
            'pan_cluster_3--pan_cluster_4--genome_3': ['genome_3', 'pan_cluster_3', 'pan_cluster_4', -10, 0,
                                                       [], []]
        }
        out_path = 'TestWritingOutputFunction'
        prefix = 'test'

        expected_low_freq = 'TestWritingOutputFunction/low_freq.txt'
        expected_gene_content = 'TestWritingOutputFunction/gene_content.txt'

        output_writer_functions.master_info_writer(master_info, out_path, prefix)

        with open(expected_low_freq, 'r') as expected:
            with open('TestWritingOutputFunction/test_low_frequency_gene_placement.tsv', 'r') as result:
                self.assertEqual(expected.readlines(), result.readlines())

        with open(expected_gene_content, 'r') as expected:
            with open('TestWritingOutputFunction/test_core_core_accessory_gene_content.tsv', 'r') as result:
                self.assertEqual(expected.readlines(), result.readlines())

    def test_summary_info_writer(self):

        input_dict = {
            'pan_cluster_1--pan_cluster_2': ['pan_cluster_1-pan_cluster_2', 3, 3, 3, 3, 99, 99, 99.0, 99.0, 3, 3, 3.0,
                                             3.0],
            'pan_cluster_2--pan_cluster_3': ['pan_cluster_2-pan_cluster_3', 2, 3, 2, 2, 100, 200, 150.0, 150.0, 0, 2,
                                             1.0, 1.0],
            'pan_cluster_2--pan_cluster_4': ['pan_cluster_2-pan_cluster_4', 1, 3, 3, 3, 150, 150, 150.0, 150.0, 1, 1,
                                             1.0, 1.0],
            'pan_cluster_3--pan_cluster_4': ['pan_cluster_3-pan_cluster_4', 2, 2, 3, 2, -10, -5, -7.5, -7.5, 0, 0, 0.0,
                                             0.0]}

        out_path = 'TestWritingOutputFunction'
        prefix = 'test'

        expected_summary_table = 'TestWritingOutputFunction/summary_table.txt'

        output_writer_functions.summary_info_writer(input_dict, out_path, prefix)

        with open(expected_summary_table, 'r') as expected:
            with open('TestWritingOutputFunction/test_core_pair_summary.csv', 'r') as result:
                self.assertEqual(expected.readlines(), result.readlines())

    def test_segment_writer(self):
        input_segments = {'pan_cluster_2--pan_cluster_4': ['pan_cluster_2',
                                                              'pan_cluster_3',
                                                              'pan_cluster_4'],
                             'pan_cluster_2--pan_cluster_6': ['pan_cluster_6',
                                                              'pan_cluster_1',
                                                              'pan_cluster_2'],
                             'pan_cluster_6--pan_cluster_4': ['pan_cluster_6',
                                                              'pan_cluster_5',
                                                              'pan_cluster_4']}

        out_path = 'TestWritingOutputFunction'
        prefix = 'test'

        expected_summary_table = 'TestWritingOutputFunction/core_segments.txt'

        output_writer_functions.segment_writer(input_segments, out_path, prefix)

        with open(expected_summary_table, 'r') as expected:
            with open('TestWritingOutputFunction/test_core_segments.csv', 'r') as result:
                self.assertEqual(expected.readlines(), result.readlines())

    def test_no_acc_segment_writer(self):
        input_segments = {'pan_cluster_2--pan_cluster_4': [['pan_cluster_2'],
                                                           ['pan_cluster_3',
                                                           'pan_cluster_4']],
                          'pan_cluster_6--pan_cluster_2': [['pan_cluster_2'],
                                                           ['pan_cluster_1'],
                                                           ['pan_cluster_6']],
                          'pan_cluster_6--pan_cluster_4': [['pan_cluster_6'],
                                                           ['pan_cluster_5',
                                                            'pan_cluster_4']
                                                           ]}

        out_path = 'TestWritingOutputFunction'
        prefix = 'test'

        expected_summary_table = 'TestWritingOutputFunction/no_acc_segments.txt'

        output_writer_functions.no_acc_segment_writer(input_segments, out_path, prefix)

        with open(expected_summary_table, 'r') as expected:
            with open('TestWritingOutputFunction/test_no_accessory_core_segments.csv', 'r') as result:
                self.assertEqual(expected.readlines(), result.readlines())

    def test_coreless_contig_writer(self):
        coreless_contigs = {'gff_1--contig_x': [['pan_cluster_2'], ['pan_cluster_3']],
                            'gff_1--contig_y': [['pan_cluster_2'], []],
                            'gff_1--contig_z': [[], ['pan_cluster_6']]}

        out_path = 'TestWritingOutputFunction'
        prefix = 'test'

        expected_summary_table = 'TestWritingOutputFunction/no_core_contigs.txt'

        output_writer_functions.non_core_contig_writer(coreless_contigs, out_path, prefix)

        with open(expected_summary_table, 'r') as expected:
            with open('TestWritingOutputFunction/test_coreless_contig_accessory_gene_content.tsv', 'r') as result:
                self.assertEqual(expected.readlines(), result.readlines())


if __name__ == '__main__':
    unittest.main()
