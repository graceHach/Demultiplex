#!/usr/bin/env python

# NOTES
# Capitalization r1, r2, r3, r4 indicates the original files R1==r1, R2==r4
# If barcodes are used as keys, they will always be a tuple (R2, R3) 
# If file handles are values, they will always be in a tuple (read1_fh, read2_fh)
# Likewise the barcodes used for the fastq headers will by appended as "[space]R2-R3"

import argparse
import bioinfo
import functions as f
import gzip

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1_in', help="input r1 filename",type=str, default="../TEST-input_FASTQ/R1.fastq.gz")
    parser.add_argument('-r2_in', help="input r2 filename",type=str, default="../TEST-input_FASTQ/R2.fastq.gz")
    parser.add_argument('-r3_in', help="input r3 filename",type=str, default="../TEST-input_FASTQ/R3.fastq.gz")
    parser.add_argument('-r4_in', help="input r4 filename",type=str, default="../TEST-input_FASTQ/R4.fastq.gz")
    parser.add_argument('-index_in', help="input filename with indexes",type=str, default='/projects/bgmp/shared/2017_sequencing/indexes.txt')
    parser.add_argument('-results_dir', help="folder to write resulting fastq into",type=str, default='results/')
    parser.add_argument('-matched_output_filename', help="filename to write to results. .tsv will be appended",type=str, default="../results_summarized/matched_frequencyTEST")
    parser.add_argument('-hopped_output_filename', help="filename to write to results. .tsv will be appended", type=str, default="../results_summarized/hopped_frequencyTEST")
    args = parser.parse_args()
    return args

# INITIALIZE VARIABLES

# Create variables for number of reads that are matched, undetermined, and hopped, initialize the values to zero
num_matched = 0
num_unk = 0
num_hopped = 0

# parse args
args = parse()
# Create a set of the 24 barcodes  
barcodes = f.read_barcodes(args.index_in)
# Initialize a dictionary - this will hold all matched barcode pairs (keys) and a tuple containing the file handles (values). When a new, 
# matched, high quality index pair is seen, the files (R1 and R2) for this pair will be opened, so whether or not a pair is present in the dictionary 
# indiciates whether the file is already open
    # To this dictionary, I'll also add "R1_hopped", "R2_hopped", "R1_unk" and "R2_unk" keys and file handle values
file_handles = {} 
r1_h = open(args.results_dir+"hopped_r1.fastq",'wt')
r4_h = open(args.results_dir+"hopped_r2.fastq",'wt')
r1_u = open(args.results_dir+"unknown_r1.fastq",'wt')
r4_u = open(args.results_dir+"unknown_r2.fastq", 'wt')
file_handles['hopped'] = (r1_h, r4_h)
file_handles["unk"] = (r1_u, r4_u)    # These need to be packed up for later

# Initialize another dictionary to count the number of high-quality reads correspoding to each of the 24 barcodes. The keys are the barcodes, and the 
# values are the number of times a pair has been encountered, so they are initialized to zero.
    # We can initialize this dictionary to contain all matched pairs using the set of known barcodes, initialize values to zero 
    # To this dictionary, we will also add high quality, unmatched pairs as they are encountered, and keep track of the number of occurances of each
    # of these hopped pairs 
    # Note that per Leslie's discussion today, we would like to keep track of the permutations (where order matters) so when these are added to the
    # dictionary as keys, the first barcode will be R2 and the second barcode should be R3
pair_freqs = {}
for code in barcodes:
    pair_freqs[code] = 0

# MAIN ITERATION

# Open all four files simultaneously, using a with open()
    # Note that "with" context manager will automatically close the FastQ files that are opened, but won't close any other files that are opened  
    # during the iteration 
with gzip.open(args.r1_in, 'rt') as r1, gzip.open(args.r2_in, 'rt') as r2, gzip.open(args.r3_in, 'rt') as r3, gzip.open(args.r4_in,'rt') as r4:
    while True:
        # Iterate over the lines of each file in groups of four
        # For each group of four lines:
        # Read into memory the fastQ IDs, the, barcodes the sequence lines, and the Q score lines for one record at a time, for all four filest
        r1_ID = r1.readline().strip()
        r2_ID = r2.readline().strip()
        r3_ID = r3.readline().strip()
        r4_ID = r4.readline().strip()
        
        # Check for EOF
        if r1_ID=='' or r2_ID=='' or r3_ID=='' or r4_ID=='':
             break
        
        r1_seq = r1.readline().strip()
        r2_seq = r2.readline().strip()
        r3_seq = r3.readline().strip()
        r4_seq = r4.readline().strip()
        # Discard plus signs 
        _ = r1.readline()
        _ = r2.readline() 
        _ = r3.readline() 
        _ = r4.readline() 

        r1_qs = r1.readline().strip()
        r2_qs = r2.readline().strip()
        r3_qs = r3.readline().strip()
        r4_qs = r4.readline().strip()

        # Generate the reverse compliment of R3, override the original value of R3_seq with its compliment
        r3_seq = f.rev_comp(r3_seq)

        # Create two new ID lines with fastQ IDs for both R1 and R2 by appending both barcodes to the end of the original fastQ IDS 
        r1_ID_new = r1_ID + " " + r2_seq + "-" + r3_seq
        r4_ID_new = r4_ID + " " + r2_seq + "-" + r3_seq

        # CATEGORIZE THE READS
        r2_mean = bioinfo.qual_score(r2_qs)
        r3_mean = bioinfo.qual_score(r3_qs)
        # Check if the indices are below the quality cutoff, or not within the set of 24
        if r2_mean < 20 or r3_mean < 20 or not r2_seq in barcodes or not r3_seq in barcodes:
            # If so, write R1 and R2, with their respective new ID lines to the undetermined file
            f.write_record(r1_seq,r1_ID_new,r1_qs,file_handles["unk"][0],r4_seq,r4_ID_new,r4_qs,file_handles["unk"][1])
            # increment undetermined count
            num_unk += 1

        # If the previous condition isn't met, check if the barcodes match
        elif not r3_seq == r2_seq:
            # If not, write R1 and R2, with their respective new ID lines to the hopped file
            f.write_record(r1_seq,r1_ID_new,r1_qs,file_handles["hopped"][0],r4_seq,r4_ID_new,r4_qs,file_handles["hopped"][1])
            # increment hopped count
            num_hopped += 1
            # Check if the index pair is present in the dictionary counting pair frequency
            if (r2_seq,r3_seq) in pair_freqs:
                # If so, increment its count
                pair_freqs[(r2_seq,r3_seq)] += 1
            else:
                # If not, add it to the dictionary and set the value to one
                pair_freqs[(r2_seq,r3_seq)] = 1
        
        # If neither of the above conditions is met, this is high quality and matched,
        else:     
            # increment the matched count
            num_matched+=1
            # increment the value in the pair-frequencies dictionary
            pair_freqs[r2_seq] += 1
            # Check if the index-pair is stored in the dictionary of matched pairs and file handles 
            if r2_seq in file_handles:
                # If so, write the records and new ID lines to the appropriate files using the file handles in the dictionary
                f.write_record(r1_seq,r1_ID_new,r1_qs,file_handles[r2_seq][0],r4_seq,r4_ID_new,r4_qs,file_handles[r2_seq][1])
            else:
                # If not, open the file for R1 and R2
                # Store both file handles in dictionary in a tuple
                file_handles[r2_seq] = (open(args.results_dir+r2_seq+"_R1.fastq",'w'),open(args.results_dir+r2_seq+"_R2.fastq",'w'))
                # write the records and new ID lines to the appropriate files using the file handles
                f.write_record(r1_seq,r1_ID_new,r1_qs,file_handles[r2_seq][0],r4_seq,r4_ID_new,r4_qs,file_handles[r2_seq][1])

# WRAP UP
# Iterate through dictionary, close all open files
for key in file_handles:
    file_handles[key][0].close()
    file_handles[key][1].close()

# Close hopped and undetermined files
file_handles["hopped"][0].close()
file_handles["hopped"][1].close()
file_handles["unk"][0].close()
file_handles["unk"][1].close()

# Write results:
print("Number of matched reads:",num_matched)
print("Number of hopped reads:",num_hopped)
print("Number of unknown reads:",num_unk)

f.write_results(num_matched,num_hopped,num_unk,pair_freqs,args.matched_output_filename, args.hopped_output_filename)
