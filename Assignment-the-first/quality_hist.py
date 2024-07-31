#!/usr/bin/env python

print("hello, python script has begun execution",flush=True)

import argparse
import numpy as np
import bioinfo
import matplotlib.pyplot as plt
import gzip

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
    parser.add_argument('-f_out', help="name of histogram to output. .png will be appended",type=str, default='quality_hist')

    # ADD ARGUMENT IF NOT USING PHRED + 33
    return parser.parse_args()

def get_stats(filename, length):
    """
    Gets np arrays of means, max values, and min values 
    Input:
        filename (str) - string filename of gzipped file
        length (int) - read length
    Output
        all np arrays of length length
        means, mins, maxs
        maxs-mins is used for error bars
    """
    means = np.zeros([length])
    mins = np.zeros([length])
    maxs = np.zeros([length])
    # open the file
    with gzip.open(filename, 'rt') as fh:
        for index, line in enumerate(fh):
            line=line.strip()
            if index%4==3: # q score lines
                for char_index, char in enumerate(line):
                    score = bioinfo.convert_phred(char)
                    if index==0:
                        means[char_index] = score
                        maxs[char_index] = score 
                        mins[char_index] = score 
                    else:
                        if score > maxs[char_index]:
                            maxs[char_index]=score
                        if score < mins[char_index]:
                            mins[char_index]=score
                        # new_mean = old_mean + (new_point - old_mean)/n (where n includes new_point)
                        # Here, our n is (index-3)/4+1 (1 based index of record)
                        n = (index-3)/4+1  # 
                        old_mean = means[char_index] # mean before current datapoint, score is added
                        means[char_index] = old_mean + (score - old_mean)/n

    return means, mins, maxs

args = parse_args()

# create array to hold all means 
# calcualted as a running mean!
print("starting r1 calcs",flush=True)
f1_means, f1_mins, f1_maxs = get_stats(args.f1_in, args.f1_length)
print("calcs done",flush=True)

fig, ax = plt.subplots(2,2, sharey=True, figsize=(12,8),layout='constrained')
ax[0,0].plot(range(1,102,1),f1_means,color='navy',label='Mean')
ax[0,0].fill_between(range(1,102,1), f1_mins, f1_maxs,alpha=.2,label="Full range",color='violet')
ax[0,0].set_title(args.f1_title)
ax[0,0].legend(loc='lower right')
ax[0,0].set_ylabel("Phred quality score")

print("starting r2",flush=True)
f2_means, f2_mins, f2_maxs = get_stats(args.f2_in, args.f2_length)
print("calcs done",flush=True)

ax[1,1].plot(range(1,9,1),f2_means,color='teal',label='Means')
ax[1,1].fill_between(range(1,9,1), f2_mins, f2_maxs,alpha=.2,label="Full range",color='cadetblue')
ax[1,1].set_title(args.f2_title)
ax[1,1].legend(loc='lower right')
ax[1,1].set_xlabel("Read position")

print("starting r3",flush=True)
f3_means, f3_mins, f3_maxs = get_stats(args.f3_in, args.f3_length)
print("calcs done",flush=True)

ax[0,1].plot(range(1,9,1),f3_means,color='seagreen',label='Means')
ax[0,1].fill_between(range(1,9,1), f3_mins, f3_maxs,alpha=.2,label="Full range",color='springgreen')
ax[0,1].set_title(args.f3_title)
ax[0,1].legend(loc='lower right')

print("starting r4",flush=True)
f4_means, f4_mins, f4_maxs = get_stats(args.f4_in, args.f4_length)
print("calcs done",flush=True)

ax[1,0].plot(range(1,102,1),f4_means,color='maroon',label='Means')
ax[1,0].fill_between(range(1,102,1), f4_mins, f4_maxs,alpha=.2,label="Full range",color='firebrick')
ax[1,0].set_title(args.f4_title)
ax[1,0].legend(loc='lower right')
ax[1,0].set_xlabel("Read position")
ax[1,0].set_ylabel("Phred quality score")

fig.suptitle("Mean, Minimum, and Maximum Quality Scores")

fig.savefig(args.f_out+".png")


