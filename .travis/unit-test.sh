#!/bin/bash

set -e
errors=0

# Run unit tests
python unit_tests/Corekaburra_test.py || {
    echo "'python python/Corekaburra/Corekaburra_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E Corekaburra/*.py || {
    echo 'pylint -E Corekaburra/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
