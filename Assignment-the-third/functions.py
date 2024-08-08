#!/usr/bin/env python
# HIGH_LEVEL FUNCTIONS

import bioinfo

def read_barcodes(barcodes_filename):
    '''
    Reads the barcodes file, keeps only the 8 character barcodes
    Input:
        filename (str) - (not file handle) of text file containing barcodes, including file extension
    Output:
        set of barcodes
    '''
    bc_set = set()
    with open(barcodes_filename, 'rt') as fh:
        for index, line in enumerate(fh):
            line=line.strip()
            line_elements = line.split()
            if index: # prevents from reading in header line
                bc_set.add(line_elements[4])
    return bc_set



def rev_comp(sequence):
    '''
    Generates the reverse compliment of a sequence
    If Ns are present, Ns will be used as base pairs
    Input:
        sequence (str) - a DNA sequence, 5' to 3'
    Output:
        rev (str) - its reverse compliment, also 5' to 3'
    '''
    compliment = {'A':"T", "T":'A', 'C':'G','G':'C', 'N':'N', 'U':'A'}
    # comment out following line if runtime is an issue
    # assert bioinfo.validate_base_seq(sequence, RNA_flag='U'in sequence.upper()) == True, "Squence passed into rev_comp is invalid"
    rev_rev = ""
    for char in sequence:
        rev_rev += compliment[char]
    # index [::-1] starts at last character and steps backward
    return rev_rev[::-1]



def write_results(matched_count, hopped_count, unk_count, barcodes_pair_dict, matched_filename, hopped_filename):
    '''
    Writes the results of analysis to two tsv files
    Writes the counts (matched_count, hopped_count and unk_count to the top of matched filename)
    Input:
        matched_count (int) - number of matched reads
        hopped_count (int) - number of hopped reads
        unk_count (int) - number of undetermined reads
        barcodes_pair_dict (dict) -  
            keys - single barcodes (indicating matched pairs)
                OR tuples of barcodes (indicating high-quality hopped reads)
            values - integer frequencies of each pair
        matched_filename (str) - filename to write hopped reads to. .tsv will be appened
        hopped_filename (str) - filename to write hopped reads to. .tsv will be appened
    Output: 
        none
    '''
    matched_filename = matched_filename+'.tsv'
    hopped_filename = hopped_filename+'.tsv'
    with open(matched_filename, 'wt') as m_fh, open(hopped_filename, 'wt') as h_fh:
        counts_strings = ['Matched count:'+'\t'+str(matched_count)+'\n','Hopped count:'+'\t'+str(hopped_count)+'\n','Unknown count:'+'\t'+str(unk_count)+'\n']
        for string in counts_strings:
            m_fh.write(string)
        for key in barcodes_pair_dict:
            # barcode is a single, matched pair
            if len(key)==8:
                line = key+'\t'+str(barcodes_pair_dict[key])+'\n'
                m_fh.write(line)
            # barcode is a tuple of unmacthed pairs
            elif len(key)==2:
                line = key[0]+'\t'+key[1]+'\t'+str(barcodes_pair_dict[key])+'\n'
                h_fh.write(line)
            else:
                print("spurious key encountered:",key)
    return


def write_record(r1_seq, r1_new_ID, r1_qs, r1_fh, r4_seq, r4_new_ID, r4_qs, r4_fh):
    '''
    Writes two records properly to two file handles that are already open!
    Which file handle is appropriate is determined in the central conditional ladder of the main code 
    File handles are retreived from the dictionary and passed into the function if present, or created 
    in the main code and added to the dictionary if not
    Note that barcodes have already been added to the IDS
    Input:
        r1_seq (str) - the full sequence (biological read 1)
        r1_new_ID (str) - fastQ ID with barcodes added (corresponding to biological read 1)
        r1_qs (str) - quality score line (corresponding to biological read 1)
        r1_fh (file handle) - appropriate file handle
        r4_seq (str) - the full sequence (biological read 2)
        r4_new_ID (str) - fastQ ID with barcodes added (corresponding to biological read 1)
        r4_qs (str) - quality score line (corresponding to biological read 2)
        r4_fh (str) - appropriate file handle
    Output:
        None
    '''
    # Implemented such that there will be a trailing newline at the end of file
    r1_fh.write(r1_new_ID+'\n')
    r1_fh.write(r1_seq+'\n')
    r1_fh.write('+'+'\n')   
    r1_fh.write(r1_qs+'\n') 
    r4_fh.write(r4_new_ID+'\n')
    r4_fh.write(r4_seq+'\n')
    r4_fh.write('+'+'\n')   
    r4_fh.write(r4_qs+'\n')  
    return