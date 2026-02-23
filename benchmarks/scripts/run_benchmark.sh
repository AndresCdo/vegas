#!/bin/bash
# Usage: run_benchmark.sh <config_file> <expected_value> <tolerance> <observable> <description>

set -e

CONFIG=$1
EXPECTED=$2
TOLERANCE=$3
OBSERVABLE=$4
DESC=$5

echo "=== $DESC ==="

./build/vegas run $CONFIG 2>&1 | tail -5

# Get the most recent result file from the config
H5_FILE=$(python -c "import json; c=json.load(open('$CONFIG')); print(c['out'])")

RESULT=$(python benchmarks/scripts/extract_observable.py $H5_FILE $OBSERVABLE)

# Calculate relative error
python -c "
import sys
v=$RESULT
e=$EXPECTED
t=$TOLERANCE
rel_err = abs(v - e) / max(abs(e), 1e-10)
passed = rel_err < t
print(f'Result: $RESULT')
print(f'Expected: $EXPECTED')
print(f'Relative error: {rel_err*100:.2f}%')
print(f'Tolerance: {t*100:.1f}%')
if passed:
    print('PASS')
    sys.exit(0)
else:
    print('FAIL')
    sys.exit(1)
"
