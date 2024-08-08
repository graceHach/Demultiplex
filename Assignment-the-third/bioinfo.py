#!/usr/bin/env python

# Author: <Grace Hach> <ghach@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = 'ATCGNatcgn'
RNA_bases = 'AUCGNaucgn'

def convert_phred(letter: str) -> int:
    '''
    Converts a single character into a phred score
    Assumes Phred + 33 encoding
    '''
    return ord(letter)-33

def qual_score(phred_string: str) -> float:
    ''' 
    Calculates an average of phred score based on a string of Phred scores
    Assumes Phred + 33 encoding
    '''
    scores = [ord(x) - 33 for x in phred_string]
    return sum(scores)/len(scores)

def validate_base_seq(seq: str, RNA_flag=False):
    '''
    This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.
    '''
    seq = seq.upper()
    return len(seq) == seq.count('U' if RNA_flag else 'T') + seq.count('A') + seq.count('C')  + seq.count('G')


def gc_content(seq: str) -> float:
    '''
    Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.
    '''
    seq = seq.upper()
    assert validate_base_seq(seq), "this is not a DNA string"
    return (seq.count("G")+seq.count("C"))/len(seq)

def calc_median(lst: list) -> float:
    """
    returns the median of an ALREADY SORTED list lst
    """
    # even case
    if len(lst)%2==0:
        middle = len(lst)//2 # upper of two middle data
        return (lst[middle]+lst[middle-1])/2
    # odd case
    else:
        return lst[len(lst)//2]

def oneline_fasta(input_filename: str, output_filename: str):
    '''
    converts a fasta to have one sequence per sequence line
    doesn't store each sequence
    Input:
        in and output filenames, string form
    Output:
        none
    '''
    first=True
    with open(input_filename, 'r') as fh_in, open(output_filename, 'w') as fh_out:
        while True:
            line = fh_in.readline()
            line = line.strip()
            if first:
                fh_out.write(line+'\n')
                first=False
            elif not line:
                break # EOF
            elif line[0]=='>': # sequence ID lines
                fh_out.write('\n'+line+'\n')
            else: # sequence lines
                fh_out.write(line)
    return

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    assert validate_base_seq("ACTG") == True, "incorrect validate_base seq"
    assert validate_base_seq("ACUG", True) == True, "incorrect validate_base seq"
    assert validate_base_seq("actg", False) == True, "incorrect validate_base seq"
    assert validate_base_seq("ATCGCTGCGACGTAGTGTAACTACTCAGAG") == True, "incorrect validate_base seq"
    assert validate_base_seq("ATCGCTGCGACGTAGTGTAACTACTCAGAG", True) == False, "incorrect validate_base seq"
    assert validate_base_seq("") == True, "incorrect validate_base seq"
    print("Validate base seq is working fine!")
