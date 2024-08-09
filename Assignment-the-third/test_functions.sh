#!/bin/bash

set -e
python test_functions.py

diff -q "test_write_results_matched_expected.tsv" "test_write_results_matched.tsv"
diff -q "test_write_results_hopped_expected.tsv" "test_write_results_hopped.tsv"
echo "If no message printed, no errors and files are as expected."