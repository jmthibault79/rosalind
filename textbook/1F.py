# Approximate Pattern Matching Problem
#
# Find all approximate occurrences of a pattern in a string.
#
# Given: Two strings Pattern and Text along with an integer d.
#
# Return: All positions where Pattern appears in Text with at most d mismatches.

# Sample Dataset
#
# ATTCTGGA
# CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
# 3

# Sample Output
#
# 6 7 26 27 78

import inout 	# my module for handling Rosalind's file I/O
pattern = inout.infilelines[0].strip()
sequence = inout.infilelines[1].strip()
d = int(inout.infilelines[2].strip())

patlen = len(pattern)

def mismatches(s1, s2):
	count = 0
	for loc in range(len(s1)):
		if s1[loc] != s2[loc]:
			count = count + 1
	return count

matches = []
for loc in range(len(sequence) - patlen + 1):
	if mismatches(pattern, sequence[loc:loc + patlen]) <= d:
		matches.append(loc)

inout.output(" ".join(map(str, matches)))
