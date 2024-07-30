#!/usr/bin/env python

import argparse
import numpy as np
import bioinfo
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f1_in', help="input fastq filename",type=str, default="mini.fq",)
    parser.add_argument('-f1_length', help="length pf reads in f1_in",type=int, default=101)
    parser.add_argument('-f1_title', help="subplot title of f1_in",type=str,default="hello!")
    parser.add_argument('-f2_in', help="input fastq filename",type=str, default="mini.fq",)
    parser.add_argument('-f2_length', help="length pf reads in f2_in",type=int, default=101)
    parser.add_argument('-f2_title', help="subplot title of f2_in",type=str)
    parser.add_argument('-f3_in', help="input fastq filename",type=str, default="mini.fq",)
    parser.add_argument('-f3_length', help="length pf reads in f3_in",type=int, default=101)
    parser.add_argument('-f3_title', help="subplot title of f3_in",type=str)
    parser.add_argument('-f4_in', help="input fastq filename",type=str, default="mini.fq",)
    parser.add_argument('-f4_length', help="length pf reads in f4_in",type=int, default=101)
    parser.add_argument('-f4_title', help="subplot title of f4_in",type=str)
    parser.add_argument('-f_out', help="name of histogram to output. .png will be appended",type=str, default='mini.out')

    # ADD ARGUMENT IF NOT USING PHRED + 33
    return parser.parse_args()

args = parse_args()

# create array to hold all means 
# calcualted as a running mean!
f1_means = np.zeros([args.f1_length])

# open the file
with open(args.f1_in) as fh:
    for index, line in enumerate(fh):
        line=line.strip()
        if index%4==3: # q score lines
            for char_index, char in enumerate(line):
                score = bioinfo.convert_phred(char)
                # new_mean = old_mean + (new_point - old_mean)/n (where n includes new_point)
                # Here, our n is (index-3)/4+1 (1 based index of record)
                n = (index-3)/4+1  # 
                old_mean = f1_means[char_index] # mean before current datapoint, score is added
                f1_means[char_index] = old_mean + (score - old_mean)/n

fig, ax = plt.subplots(2,2)
ax[0,0].errorbar(range(1,102,1),f1_means)
ax[0,0].set_title(args.f1_title)
fig.savefig(args.f_out+".png")


