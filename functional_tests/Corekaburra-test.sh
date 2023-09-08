#!/usr/bin/env bash

# 1. Parse command line arguments
# 2. cd to the test directory
# 3. run tests
# 4. Print summary of successes and failures, exit with 0 if
#    all tests pass, else exit with 1

# Uncomment the line below if you want more debugging information
# about this script.
#set -x

# The name of this test script
this_program_name="Corekaburra-test.sh"
# The program we want to test (either a full path to an executable, or the name of an executable in $PATH)
test_program=""
# Directory containing the test data files and expected outputs
test_data_dir=""
# Number of failed test cases
num_errors=0
# Total number of tests run
num_tests=0

function show_help {
cat << UsageMessage

${this_program_name}: run integration/regression tests for Corekaburra 

Usage:
    ${this_program_name} [-h] [-v] -p program -d test_data_dir 

Example:
    ${this_program_name} -p bin/Corekaburra -d data/tests

-h shows this help message

-v verbose output
UsageMessage
}

# echo an error message $1 and exit with status $2
function exit_with_error {
    printf "${this_program_name}: ERROR: $1\n"
    exit $2
}

# if -v is specified on the command line, print a more verbaose message to stdout
function verbose_message {
    if [ "${verbose}" = true ]; then
        echo "${this_program_name} $1"
    fi
}

# Parse the command line arguments and set the global variables program and test_data_dir 
function parse_args {
    local OPTIND opt

    while getopts "hp:d:v" opt; do
        case "${opt}" in
            h)
                show_help
                exit 0
                ;;
            p)  test_program="${OPTARG}"
                ;;
            d)  test_data_dir="${OPTARG}"
                ;;
            v)  verbose=true
                ;;
        esac
    done

    shift $((OPTIND-1))

    [ "$1" = "--" ] && shift

    if [[ -z ${test_program} ]]; then
        exit_with_error "missing command line argument: -p program, use -h for help" 2
    fi

    if [[ -z ${test_data_dir} ]]; then
        exit_with_error "missing command line argument: -d test_data_dir, use -h for help" 2
    fi
}


# Run a command and check that the output is
# exactly equal the contents of a specified file 
# ARG1: command we want to test as a string
# ARG2: a file path containing the expected output
# ARG3: expected exit status
function test_stdout_exit {
    let num_tests+=1
    output=$(eval $1)
    exit_status=$?
    expected_output_file=$2
    expected_exit_status=$3
    verbose_message "Testing stdout and exit status: $1"
    difference=$(diff <(echo "$output") $expected_output_file)
    if [ -n "$difference" ]; then 
        let num_errors+=1
        echo "Test output failed: $1"
        echo "Actual output:"
        echo "$output"
        expected_output=$(cat $2)
        echo "Expected output:"
        echo "$expected_output"
        echo "Difference:"
        echo "$difference"
    elif [ "$exit_status" -ne "$expected_exit_status" ]; then
        let num_errors+=1
        echo "Test exit status failed: $1"
        echo "Actual exit status: $exit_status"
        echo "Expected exit status: $expected_exit_status"
    fi 
}

# Run a command and check that the output file is
# exactly equal the contents of a specified file
# ARG1: A file returned from program after running
# ARG2: a file path containing the expected output
function test_output_file {
    let num_tests+=1
    output=$1
    expected_output_file=$2
    verbose_message "Testing output file: $1"
    verbose_message "Expected file path: $2"
    difference=$(diff $output $expected_output_file) || let num_errors+=1
    if [ -n "$difference" ]; then
        let num_errors+=1
        echo "Test output failed: $1"
        echo "Actual output:"
        cat $output
        expected_output=$(cat $expected_output_file)
        echo "Expected output:"
        echo "$expected_output"
        echo "Difference:"
        echo "$difference"
    fi
}

# Run a command and check that the exit status is 
# equal to an expected value
# exactly equal the contents of a specified file 
# ARG1: command we want to test as a string
# ARG2: expected exit status
# NB: this is mostly for checking erroneous conditions, where the
# exact output message is not crucial, but the exit status is
# important
function test_exit_status {
    let num_tests+=1
    output=$(eval $1)
    exit_status=$?
    expected_exit_status=$2
    verbose_message "Testing exit status: $1"
    if [ "$exit_status" -ne "$expected_exit_status" ]; then
        let num_errors+=1
        echo "Test exit status failed: $1"
        echo "Actual exit status: $exit_status"
        echo "Expected exit status: $expected_exit_status"
    fi 
}

function call_new_test {
  echo ''
  echo $1
}


# 1. Parse command line arguments.
parse_args $@
# 2. Change to test directory
cd $test_data_dir
# 2. Run tests

call_new_test "Test output for no arguments"
test_stdout_exit "$test_program" no_input.expected 2

call_new_test "Test output for -help argument given"
test_stdout_exit "$test_program -help" no_input.expected 0

call_new_test "Test exit status for a bad command line invocation"
test_exit_status "$test_program --this_is_not_a_valid_argument > /dev/null 2>&1" 2

call_new_test "Test exit status for a bad cutoffs provided - core lower than low-frequency"
test_exit_status "$test_program -ig complete_genome_double_chrom.gff -ip Crash_pan_folder -cg complete_genomes_file -cc 0.1 -lc 0.2 > /dev/null 2>&1" 2

call_new_test "Test exit status for a bad cutoffs provided - core above range"
test_exit_status "$test_program -ig complete_genome_double_chrom.gff -ip Crash_pan_folder -cg complete_genomes_file -cc 1.1 -lc 0.2 > /dev/null 2>&1" 2

call_new_test "Test exit status for a bad cutoffs provided - low-frequency below range"
test_exit_status "$test_program -ig complete_genome_double_chrom.gff -ip Crash_pan_folder -cg complete_genomes_file -cc 1 -lc -0.2 > /dev/null 2>&1" 2

call_new_test "Test exit status for a complete genome not given as input gff file"
test_exit_status "$test_program -ig complete_genome_double_chrom.gff -ip Crash_pan_folder -cg complete_genomes_file > /dev/null 2>&1" 1

call_new_test "Test exit upon unsuccessful identification of source program"
test_exit_status "$test_program -ig complete_genome_double_chrom.gff -ip Crash_pan_folder > /dev/null 2>&1" 1

call_new_test "Test exit upon gff not found in pan is provided as input"
test_exit_status "$test_program -ig complete_genome_single_chrom.gff complete_genome_double_chrom.gff -ip Crash_gff_folder > /dev/null 2>&1" 1

call_new_test "Test Roary input"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_single_chrom_2.gff -ip Roray_run -o test_out_folder > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Simple_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Simple_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Simple_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test Panaroo input"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_single_chrom_2.gff -ip Panaroo_run -o test_out_folder > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Simple_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Simple_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Simple_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder


call_new_test "test gzipped input gffs"
Corekaburra -ig complete_genome_single_chrom.gff.gz complete_genome_single_chrom_2.gff.gz -ip Roray_run -o test_out_folder > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Simple_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Simple_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Simple_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder


call_new_test "test mixed gzipped and non-gzipped input gffs"
Corekaburra -ig complete_genome_single_chrom.gff.gz complete_genome_single_chrom_2.gff -ip Roray_run -o test_out_folder > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Simple_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Simple_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Simple_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder


call_new_test "Test with GFF where reannotated in panaroo is called 'candidate_gene' instead of CDS"
Corekaburra -ig candidate_gene_1.gff candidate_gene_2.gff -ip candidate_gene_pan_input -o test_out_folder > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv candidate_gene_pan_input_expected/core_core_accessory_gene_content.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv candidate_gene_pan_input_expected/low_frequency_gene_placement.expected
test_output_file test_out_folder/core_pair_summary.tsv candidate_gene_pan_input_expected/core_pair_summary.expected
rm -r test_out_folder


call_new_test "Test complete genome with single contig and single complete genome among input"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_single_chrom_2.gff -ip Roray_run -o test_out_folder -cg Complete_single_chromosome.txt > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Single_comple_chromosome_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Single_comple_chromosome_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Single_comple_chromosome_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test complete genome with multiple contigs (Simulate plasmids or two chromosomes)"
Corekaburra -ig complete_genome_double_chrom.gff complete_genome_double_chrom_2.gff -ip complete_double_chromoosme_run -o test_out_folder -cg Complete_double_chromosomes.txt > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv double_comple_chromosome_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv double_comple_chromosome_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv double_comple_chromosome_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test with accessory genes"
Corekaburra -ig genome_single_chrom_larger.gff genome_single_chrom_larger_rearrange.gff -ip Accessory_chrom_run -o test_out_folder -cg complete_larger_genome_list.txt > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Accessory_chrom_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Accessory_chrom_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Accessory_chrom_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test with segments and sub-segments"
Corekaburra -ig genome_single_chrom_larger.gff genome_single_chrom_larger_rearrange.gff -ip Rearrangement_run -o test_out_folder -cg complete_larger_genome_list.txt > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Rearrangement_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Rearrangement_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Rearrangement_run_expected/core_pair_summary.tsv.expected
test_output_file test_out_folder/core_segments.tsv Rearrangement_run_expected/core_segments.tsv.expected
test_output_file test_out_folder/no_accessory_core_segments.tsv Rearrangement_run_expected/no_accessory_core_segments.tsv.expected
rm -r test_out_folder

call_new_test "Test that segmnets can be identified with a core-cutoff that is less than all genomes"
Corekaburra -ip Less_than_all_core_simple/ -ig complete_genome_double_chrom_extra_large.gff complete_genome_double_chrom_larger.gff complete_genome_double_chrom_3_extra_large.gff -o test_out_folder/ -cc 0.9  > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Less_than_all_core_simple_expected/core_core_accessory_gene_content.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Less_than_all_core_simple_expected/low_frequency_gene_placement.expected
test_output_file test_out_folder/core_pair_summary.tsv Less_than_all_core_simple_expected/core_pair_summary.expected
test_output_file test_out_folder/core_segments.tsv Less_than_all_core_simple_expected/core_segments.expected
test_output_file test_out_folder/no_accessory_core_segments.tsv Less_than_all_core_simple_expected/no_accessory_core_segments.expected
rm -r test_out_folder


call_new_test "Test that segments are produced for instances where a core gene may be left out"
Corekaburra -ip Less_than_all_core_complex/ -ig complete_genome_double_chrom_extra_large.gff complete_genome_double_chrom_2_larger.gff complete_genome_double_chrom_3_extra_large.gff -o test_out_folder/ -cc 0.9  > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Less_than_all_core_complex_expected/core_core_accessory_gene_content.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Less_than_all_core_complex_expected/low_frequency_gene_placement.expected
test_output_file test_out_folder/core_pair_summary.tsv Less_than_all_core_complex_expected/core_pair_summary.expected
test_output_file test_out_folder/core_segments.tsv Less_than_all_core_complex_expected/core_segments.expected
test_output_file test_out_folder/no_accessory_core_segments.tsv Less_than_all_core_complex_expected/no_accessory_core_segments.expected
rm -r test_out_folder


call_new_test "Test using two complete chromosomes (complete genome) when >1 contig test"
Corekaburra -ip Multiple_component_graph/ -ig complete_genome_double_chrom_larger.gff complete_genome_double_chrom_2_larger.gff complete_genome_double_chrom_3_larger.gff -o test_out_folder/ -cc 0.9 -cg complete_larger_double_chr_genome_list.txt  > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Multiple_component_graph_complete_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Multiple_component_graph_complete_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Multiple_component_graph_complete_expected/core_pair_summary.tsv.expected
test_output_file test_out_folder/core_segments.tsv Multiple_component_graph_complete_expected/core_segments.tsv.expected
test_output_file test_out_folder/no_accessory_core_segments.tsv Multiple_component_graph_complete_expected/no_accessory_core_segments.tsv.expected
rm -r test_out_folder

call_new_test "Test when core graph forms multiple components - not forming a single 'chromosome' - non circular input gffs"
Corekaburra -ip Multiple_component_graph/ -ig complete_genome_double_chrom_larger.gff complete_genome_double_chrom_2_larger.gff complete_genome_double_chrom_3_larger.gff -o test_out_folder/  > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Multi_component_graph_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Multi_component_graph_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Multi_component_graph_expected/core_pair_summary.tsv.expected
test_output_file test_out_folder/core_segments.tsv Multi_component_graph_expected/core_segments.tsv.expected
test_output_file test_out_folder/no_accessory_core_segments.tsv Multi_component_graph_expected/no_accessory_core_segments.tsv.expected
rm -r test_out_folder

call_new_test "Test when core graph forms multiple components - not forming a single 'chromosome' - circular input gffs"
Corekaburra -ip Multiple_component_graph/ -ig complete_genome_double_chrom_larger.gff complete_genome_double_chrom_2_larger.gff complete_genome_double_chrom_3_larger.gff -o test_out_folder/ -cg complete_larger_double_chr_genome_list.txt > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Multiple_component_graph_complete_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Multiple_component_graph_complete_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Multiple_component_graph_complete_expected/core_pair_summary.tsv.expected
test_output_file test_out_folder/core_segments.tsv Multiple_component_graph_complete_expected/core_segments.tsv.expected
test_output_file test_out_folder/no_accessory_core_segments.tsv Multiple_component_graph_complete_expected/no_accessory_core_segments.tsv.expected
rm -r test_out_folder

call_new_test "Test with decreased core-gene cutoff"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_single_chrom_2.gff genome_single_chrom_larger.gff -ip Change_cutoffs -o test_out_folder -cc 0.9 > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv core_90_cutoff_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv core_90_cutoff_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv core_90_cutoff_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test with increase low cutoff"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_single_chrom_2.gff genome_single_chrom_larger.gff -ip Change_cutoffs -o test_out_folder -lc 0.4 > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Increase_low_cutoff_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Increase_low_cutoff_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Increase_low_cutoff_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test with zero low cutoff"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_single_chrom_2.gff genome_single_chrom_larger.gff -ip Change_cutoffs -o test_out_folder -lc 0 > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv low_freq_cutoff_0_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv low_freq_cutoff_0_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv low_freq_cutoff_0_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test with less than all gffs from pan-genome provided"
Corekaburra -ig complete_genome_single_chrom.gff genome_single_chrom_larger.gff genome_single_chrom_larger_rearrange.gff -ip Less_than_all_gffs -o test_out_folder -cc 0.9 > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Less_than_all_gffs_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Less_than_all_gffs_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Less_than_all_gffs_run_expected/core_pair_summary.tsv.expected
test_output_file test_out_folder/core_segments.tsv Less_than_all_gffs_run_expected/core_segments.tsv.expected
test_output_file test_out_folder/no_accessory_core_segments.tsv Less_than_all_gffs_run_expected/no_accessory_core_segments.tsv.expected
rm -r test_out_folder

call_new_test "Test with a single core gene on a contig that is not complete"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_double_chrom.gff -ip Single_core_contig/ -o test_out_folder/ > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv single_core_contig_draft_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv single_core_contig_draft_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv single_core_contig_draft_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test with a single core gene on a contig that is complete"
Corekaburra -ig complete_genome_single_chrom.gff complete_genome_double_chrom.gff -ip Single_core_contig/ -o test_out_folder/ -cg complete_genomes_file > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv single_core_contig_complete_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv single_core_contig_complete_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv single_core_contig_complete_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test for core genes being fragmented"
Corekaburra -ig complete_genome_single_chrom.gff genome_single_chrom_larger_2.gff genome_single_chrom_larger_rearrange.gff complete_genome_single_chrom_2.gff -ip Fragmented_core_gene_run/ -o test_out_folder/  > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Fragmented_core_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Fragmented_core_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Fragmented_core_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test a fragmented core gene not accepted as core"
Corekaburra -ig complete_genome_single_chrom.gff genome_single_chrom_larger.gff genome_single_chrom_larger_rearrange.gff complete_genome_single_chrom_2.gff -ip Fragmented_core_gene_break_run/ -o test_out_folder/ > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Fragmented_core_gene_break_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Fragmented_core_gene_break_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Fragmented_core_gene_break_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test for accessory genes being fragmented"
Corekaburra -ig complete_genome_single_chrom.gff genome_single_chrom_larger.gff genome_single_chrom_larger_rearrange.gff complete_genome_single_chrom_2.gff -ip Fragmented_accessory_gene_run/ -o test_out_folder/ > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Fragmented_accessory_gene_run_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Fragmented_accessory_gene_run_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Fragmented_accessory_gene_run_expected/core_pair_summary.tsv.expected
rm -r test_out_folder

call_new_test "Test with a core-less contig draft"
Corekaburra -ig complete_genome_double_chrom_2.gff complete_genome_double_chrom.gff -ip Coreless_contig_run/ -o test_out_folder/ > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv coreless_contig_draft_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv coreless_contig_draft_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv coreless_contig_draft_expected/core_pair_summary.tsv.expected
test_output_file test_out_folder/coreless_contig_accessory_gene_content.tsv coreless_contig_draft_expected/coreless_contig_accessory_gene_content.tsv.expected
rm -r test_out_folder

call_new_test "Test with a core-less contig complete"
Corekaburra -ig complete_genome_double_chrom_2.gff complete_genome_double_chrom.gff -ip Coreless_contig_run/ -o test_out_folder/ -cg Complete_double_chromosomes.txt > /dev/null 2>&1
test_output_file test_out_folder/core_core_accessory_gene_content.tsv Coreless_contig_complete_expected/core_core_accessory_gene_content.tsv.expected
test_output_file test_out_folder/low_frequency_gene_placement.tsv Coreless_contig_complete_expected/low_frequency_gene_placement.tsv.expected
test_output_file test_out_folder/core_pair_summary.tsv Coreless_contig_complete_expected/core_pair_summary.tsv.expected
test_output_file test_out_folder/coreless_contig_accessory_gene_content.tsv Coreless_contig_complete_expected/coreless_contig_accessory_gene_content.tsv.expected
rm -r test_out_folder


# 3. End of testing - check if any errors occurrred
if [ "$num_errors" -gt 0 ]; then
    echo "$test_program failed $num_errors out of $num_tests tests"
    exit 1
else
    echo "$test_program passed all $num_tests successfully"
    exit 0
fi
