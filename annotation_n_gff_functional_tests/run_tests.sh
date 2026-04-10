#!/usr/bin/env bash
# Confined functional tests for GFF format handling in Corekaburra.
#
# Tests three GFF flavours that must all produce valid output:
#   1. Prokka raw      — .gff      (baseline)
#   2. Bakta raw       — .gff3     (Bug 1 + Bug 3: extension stripping)
#   3. Panaroo-corrected — _panaroo.gff  (Bug 1 + Bug 2: suffix stripping)
#
# Usage:
#   bash annotation_n_gff_functional_tests/run_tests.sh [-v]
#
# Options:
#   -v  Verbose output (show Corekaburra stdout/stderr)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DATA="${SCRIPT_DIR}/test_data"
PAN_DIR="${TEST_DATA}/panaroo_pan"

num_tests=0
num_passed=0
num_failed=0
failed_tests=()
verbose=false

while getopts "v" opt; do
    case "${opt}" in
        v) verbose=true ;;
        *) echo "Usage: $0 [-v]"; exit 1 ;;
    esac
done

EXPECTED_OUTPUTS=(
    "core_core_accessory_gene_content.tsv"
    "core_pair_summary.tsv"
    "low_frequency_gene_placement.tsv"
)

# ── Helpers ──────────────────────────────────────────────────────────────────

run_test() {
    local test_name="$1"
    local gff_files="$2"
    local output_dir="${SCRIPT_DIR}/test_output_${num_tests}"

    num_tests=$((num_tests + 1))

    echo ""
    echo "────────────────────────────────────────────────────────"
    echo "[INFO] Running: ${test_name}"

    rm -rf "${output_dir}"

    local exit_code=0
    if [ "$verbose" = true ]; then
        Corekaburra -ig ${gff_files} -ip "${PAN_DIR}" -o "${output_dir}" || exit_code=$?
    else
        Corekaburra -ig ${gff_files} -ip "${PAN_DIR}" -o "${output_dir}" > /dev/null 2>&1 || exit_code=$?
    fi

    if [ "$exit_code" -ne 0 ]; then
        num_failed=$((num_failed + 1))
        failed_tests+=("${test_name} (exit code: ${exit_code})")
        echo "[FAIL] ${test_name} — Corekaburra exited with code ${exit_code}"
        rm -rf "${output_dir}"
        return
    fi

    local missing=()
    for expected_file in "${EXPECTED_OUTPUTS[@]}"; do
        if [ ! -s "${output_dir}/${expected_file}" ]; then
            missing+=("${expected_file}")
        fi
    done

    if [ ${#missing[@]} -gt 0 ]; then
        num_failed=$((num_failed + 1))
        failed_tests+=("${test_name} (missing/empty: ${missing[*]})")
        echo "[FAIL] ${test_name} — missing or empty output files: ${missing[*]}"
    else
        num_passed=$((num_passed + 1))
        echo "[PASS] ${test_name}"
    fi

    rm -rf "${output_dir}"
}

# ── Tests ─────────────────────────────────────────────────────────────────────

echo "============================================================"
echo " Corekaburra GFF Format Functional Tests"
echo "============================================================"

# Test 1: Prokka raw GFFs (.gff) — baseline
run_test "Prokka raw GFFs (.gff)" \
    "${TEST_DATA}/genome_A.gff ${TEST_DATA}/genome_B.gff"

# Test 2: Bakta raw GFFs (.gff3) — extension stripping bug fix
run_test "Bakta raw GFFs (.gff3)" \
    "${TEST_DATA}/genome_A.gff3 ${TEST_DATA}/genome_B.gff3"

# Test 3: Panaroo-corrected GFFs (_panaroo.gff) — suffix stripping bug fix
run_test "Panaroo-corrected GFFs (_panaroo.gff)" \
    "${TEST_DATA}/genome_A_panaroo.gff ${TEST_DATA}/genome_B_panaroo.gff"

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo " Summary"
echo "============================================================"
echo " Total:  ${num_tests}"
echo " Passed: ${num_passed}"
echo " Failed: ${num_failed}"

if [ ${#failed_tests[@]} -gt 0 ]; then
    echo ""
    echo " Failed tests:"
    for t in "${failed_tests[@]}"; do
        echo "   - ${t}"
    done
fi

echo "============================================================"

if [ "${num_failed}" -gt 0 ]; then
    exit 1
else
    echo ""
    echo "All ${num_tests} tests passed."
    exit 0
fi
