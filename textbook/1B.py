# Reverse Complement Problem
#
# Reverse complement a nucleotide pattern.
#
# Given: A DNA string Pattern.
#
# Return: Pattern, the reverse complement of Pattern.

# Sample Dataset
#
# AAAACCCGGT

# Sample Output
#
# ACCGGGTTTT

complement = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' }

import inout 	# my module for handling Rosalind's file I/O
sequence = inout.infilelines[0].strip()

output = ''
for base in reversed(sequence):
	output = output + complement[base]

inout.output(output)
