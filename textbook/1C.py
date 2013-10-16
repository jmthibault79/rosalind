# Pattern Matching Problem
#
# Find all occurrences of a pattern in a string.
#
# Given: Two strings, Pattern and Text.
#
# Return: All starting positions where Pattern appears as a substring of Text.

# Sample Dataset
#
# ATAT
# GATATATGCATATACTT

# Sample Output
#
# 1 3 9

import inout 	# my module for handling Rosalind's file I/O
pattern = inout.infilelines[0].strip()
sequence = inout.infilelines[1].strip()

patlen = len(pattern)
seqlen = len(sequence)
output = ''
for idx in range(seqlen - patlen + 1):
	if sequence[idx:idx+patlen] == pattern:
		output = output + str(idx) + " "

inout.output(output)
