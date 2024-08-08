#!/usr/bin/env python
import functions

#   IMPORTANT
# note that this file must be run from test_function.sh to compare output files

#   TESTING read_barcodes()
correct_set = {"GTAGCGTA", "CGATCGAT", "GATCAAGG", "AACAGCGA", "TAGCCATG", "CGGTAATC", "CTCTGGAT",
    "TACCGGAT", "CTAGCTCA", "CACTTCAC", "GCTACTCT", "ACGATCAG", "TATGGCAC", "TGTTCCGT", 
    "GTCCTAAG", "TCGACAAG", "TCTTCGAC", "ATCATGCG", "ATCGTGGT", "TCGAGAGT", "TCGGATTC", 
    "GATCTTGC", "AGAGTCCA", "AGGATAGC"}
returned_set = functions.read_barcodes('/projects/bgmp/shared/2017_sequencing/indexes.txt')

# make both into ordered list for item by items comparison
correct_set_list = sorted(list(correct_set))
returned_set_list = sorted(list(returned_set))

for i, item in enumerate(correct_set_list):
    assert item==returned_set_list[i], "test of read_barcodes failed"


#   TESTING rev_comp()
seqs = ['GAGCCCTAGCTCACGTTTTT', 'AACAGCGA','TCTGATCTCGCCCCGACAACTGCAAACCCCAA']
correct_rev = ['AAAAACGTGAGCTAGGGCTC','TCGCTGTT','TTGGGGTTTGCAGTTGTCGGGGCGAGATCAGA']

for index, seq in enumerate(seqs):
    assert functions.rev_comp(seq) == correct_rev[index], "test of rev_comp failed"


#   TESTING write_results()
matched_count = 9
hopped_count = 2
unk_count = 3
m_filename = 'test_write_results_matched'
h_filename = 'test_write_results_hopped'
barcodes_pair_dict = {"GTAGCGTA":4, "AACAGCGA":5, ("TCGGATTC","ATCGTGGT"):2}
functions.write_results(matched_count, hopped_count, unk_count, barcodes_pair_dict, m_filename, h_filename)


#   TESTING write_record()
r1_fh = open("r1_test.txt",'wt')
r4_fh = open("r2_test.txt",'wt')
r1_seq = "AAAAA"
r4_seq = "TTTTT"
r1_new_ID = "@I'm the r1 ID GTAGCGTA-GTAGCGTA"
r4_new_ID = "@I'm the r4 ID GTAGCGTA-GTAGCGTA"
r1_qs = "FJJJF"
r4_qs = "J<FJJ"
functions.write_record(r1_seq, r1_new_ID, r1_qs, r1_fh, r4_seq, r4_new_ID, r4_qs, r4_fh)