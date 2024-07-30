#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --output=LOG/quality_hists_%j.out
#SBATCH --error=LOG/quality_hists_%j.err
#SBATCH --job-name=quality_hists

set -e

python quality_hist.py -f2_in "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" \
                       -f2_length 8 -f2_title "Read 2" \
                       -f4_in "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" \
                       -f4_length 101 -f4_title "Read 4" \
                       -f3_in "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" \
                       -f3_length 8 -f3_title "Read 3" \
                       -f1_in "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" \
                       -f1_length 101 -f1_title "Read 1"
