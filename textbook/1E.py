# Minimum Skew Problem
#
# Find a position in a genome minimizing the skew.
#
# Given: A DNA string Genome.
# 
# Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i (from 0 to |Genome|).

# Sample Dataset
# 
# CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG

# Sample Output
# 
# 53 97

import inout 	# my module for handling Rosalind's file I/O
sequence = inout.infilelines[0].strip()

all_min_skew_loc = []
skew, skew_loc, min_skew = 0, 1, 1000
for base in sequence:
	if base == 'C':
		skew = skew - 1
	elif base == 'G':
		skew = skew + 1
	
	if skew < min_skew:
		min_skew = skew
		all_min_skew_loc = [skew_loc]
	elif skew == min_skew:
		all_min_skew_loc.append(skew_loc)
		
	skew_loc = skew_loc + 1

inout.output(" ".join(map(str, all_min_skew_loc)))
