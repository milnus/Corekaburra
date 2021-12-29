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



# move to folder with mock files. First try Github structure, then try pulled repository structure
try:
    os.chdir('/Corekaburra/unit_tests/unit_test_data/')
except FileNotFoundError:
    os.chdir('unit_test_data/')


class TestExitWithError(unittest.TestCase):
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
    def test_Genedata_File_present(self):
        input_folder_path = 'TestPresenceOfGenedataFile/present'
        return_path = check_inputs.check_gene_data(input_folder_path)

        self.assertEqual(return_path, input_folder_path +'/gene_data.csv')

    def test_Genedata_File_absent(self):
        input_folder_path = 'TestPresenceOfGenedataFile/absent'

        with self.assertRaises(SystemExit):
            check_inputs.check_gene_data(input_folder_path)


if __name__ == '__main__':
    unittest.main()
