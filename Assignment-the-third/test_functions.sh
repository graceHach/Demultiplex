#!/bin/bash

set -e
python test_functions.py

diff -s "test_write_results_matched_expected.tsv" "test_write_results_matched.tsv"
diff -s "test_write_results_hopped_expected.tsv" "test_write_results_hopped.tsv"
echo "Both files should be identical"