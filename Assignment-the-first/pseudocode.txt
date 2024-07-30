# NOTES

# If barcodes are used as keys, they will always be a tuple (R2, R3) 
# If file handles are values, they will always be in a tuple (read1_fh, read2_fh)
# Likewise the barcodes used for the fastq headers will by appended as "[space]R2-R3"

# INITIALIZE VARIABLES

# Create variables for number of reads that are matched, undetermined, and hopped, initialize the values to zero
# Create a set of the 24 barcodes  
	# Read this from a file that get arg-parsed
# Initialize a dictionary - this will hold all matched barcode pairs (keys) and a tuple containing the file handles (values). When a new, 
# matched, high quality index pair is seen, the files (R1 and R2) for this pair will be opened, so whether or not a pair is present in the dictionary 
# indiciates whether the file is already open
	# To this dictionary, I'll also add "R1_hopped", "R2_hopped", "R1_unk" and "R2_unk" keys and file handle values 
# Initialize another dictionary to count the number of high-quality reads correspoding to each of the 24 barcodes. The keys are the barcodes, and the 
# values are the number of times a pair has been encountered, so they are initialized to zero.
	# We can initialize this dictionary to contain all matched pairs using the set of known barcodes, initialize values to zero 
	# To this dictionary, we will also add high quality, unmatched pairs as they are encountered, and keep track of the number of occurances of each
	# of these hopped pairs 
	# Note that per Leslie's discussion today, we would like to keep track of the permutations (where order matters) so when these are added to the
	# dictionary as keys, the first barcode will be R2 and the second barcode should be R3

# MAIN ITERATION

# Open all four files simultaneously, using a with open() and a while-True loop 
	# Note that "with" context manager will automatically close the FastQ files that are opened, but won't close any other files that are opened  
	# during the iteration 
# Iterate over the lines of each file in groups of four
# For each group of four lines:
	# Read into memory the fastQ IDs, the, barcodes the sequence lines, and the Q score lines for one record at a time, for all four files
	# Generate the reverse compliment of R3, override the original value of R3_seq with its compliment
	# Create two new ID lines with fastQ IDs for both R1 and R2 by appending both barcodes to the end of the original fastQ IDS  

	# CATEGORIZE THE READS

	# Check if the indices are below the quality cutoff, or not within the set of 24
		# If so, write R1 and R2, with their respective new ID lines to the undetermined file
		# increment undetermined count
	# If the previous condition isn't met, check if the barcodes match
		# If not, write R1 and R2, with their respective new ID lines to the hopped file
		# increment hopped count
			# Check if the index pair is present in the dictionary counting pair frequency
				# If so, increment its count
				# If not, add it to the dictionary and set the value to one
	# If neither of the above conditions is met, this is high quality and matched,
		# increment the matched count
		# increment the value in the pair-frequencies dictionary
		# Check if the index-pair is stored in the dictionary of matched pairs and file handles 
			# If so, write the records and new ID lines to the appropriate files using the file handles in the dictionary
			# If not, open the file for R1 and R2
				# Add the barcode pair to the dictionary as a key and the file handles packed together into a tuple as the value
				# write the records and new ID lines to the appropriate files using the file handles
# break when EOF is reached

# WRAP UP

# Iterate through dictionary, close all open files
# Close hopped and undetermined files
# Using with open(), open a file in order to write results. To this file, write:
	# Number of matched, hopped, and undetermined indices
	# frequencies of each barcode pair by iterating over the dictionary 
