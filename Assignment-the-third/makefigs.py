#!/usr/bin/env python
import matplotlib
from matplotlib import pyplot as plt
import argparse 
import functions as f
import numpy as np


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hopped_file', help="input r1 filename",type=str, default="../results_summarized/hopped_frequency.tsv")
    parser.add_argument('-matched_file', help="input r2 filename",type=str, default="../results_summarized/matched_frequency.tsv")
    parser.add_argument('-index_in', help="input filename with indexes",type=str, default='/projects/bgmp/shared/2017_sequencing/indexes.txt')
    parser.add_argument('-results_dir', help="folder to write resulting fastq into",type=str, default='../results_summarized/')
    parser.add_argument('-matched_output_filename', help="filename to write to results. .tsv will be appended",type=str, default="../results_summarized/matched_frequencyTEST")
    parser.add_argument('-hopped_output_filename', help="filename to write to results. .tsv will be appended", type=str, default="../results_summarized/hopped_frequencyTEST")
    args = parser.parse_args()
    return args

args = parse()
# Create a set of the 24 barcodes  
barcodes = f.read_barcodes(args.index_in)

# sort them
bc_list = sorted(list(barcodes))

# Make 2D array

# I'm not using itertools for this beacuse I want them sorted
# first element is r2, second is r3
perms = []
for bar1 in bc_list:
    for bar2 in bc_list:
        perms.append((bar1,bar2))

# keys are the tuples, values are the frequencies, initialized to zero
arr = np.zeros([len(bc_list),len(bc_list)])

inds = {}
for index, bc in enumerate(bc_list):
    inds[bc]=index

with open(args.hopped_file) as hf:
    for line in hf:
        line = line.strip()
        line_el = line.split()
        arr[inds[line_el[0]], inds[line_el[1]]] = line_el[2]

with open(args.matched_file) as mf:
    for index, line in enumerate(mf):
        # this file has headers
        if index>2:
            line = line.strip()
            line_el = line.split()
            arr[inds[line_el[0]], inds[line_el[0]]] = line_el[1]

fig, ax = plt.subplots(figsize=(9,9))
im = ax.imshow(arr,norm=matplotlib.colors.LogNorm())

ax.set_xticks(np.arange(len(bc_list)), labels=bc_list)
ax.set_yticks(np.arange(len(bc_list)), labels=bc_list)
plt.setp(ax.get_xticklabels(), rotation=270, ha="right")
ax.set_title("Frequency of Barcode Pairs")
fig.colorbar(im)

fig.savefig(args.results_dir+"fig1.png")

