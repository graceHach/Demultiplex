#!/bin/bash

set -e
# defaults are test files
python main.py

diff -q "results/GTAGCGTA_R1.fastq" "../TEST-output_FASTQ/GTAGCGTA_r1.fastq"
diff -q "results/GTAGCGTA_R2.fastq" "../TEST-output_FASTQ/GTAGCGTA_r2.fastq"
diff -q "results/hopped_r1.fastq" "../TEST-output_FASTQ/hopped_r1.fastq"
diff -q "results/hopped_r2.fastq" "../TEST-output_FASTQ/hopped_r2.fastq"
diff -q "results/unknown_r1.fastq" "../TEST-output_FASTQ/unknown_r1.fastq"
diff -q "results/unknown_r2.fastq" "../TEST-output_FASTQ/unknown_r2.fastq"
echo "If no message printed, no errors and files are as expected."