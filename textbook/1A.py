# Frequent Words Problem
#
# Find the most frequent k-mers in a string.
#
# Given: A string Text and an integer k.
#
# Return: All most frequent k-mers in Text.

# Sample Dataset
#
# ACGTTGCATGTCGCATGATGCATGAGAGCT
# 4

# Sample Output
#
# CATG GCAT

import inout 	# my module for handling Rosalind's file I/O

sequence = inout.infilelines[0].strip()
k = int(inout.infilelines[1].strip())

kmer_counts = {}
max_kmer_count = 0
for idx in range(len(sequence) - k + 1):
	kmer = sequence[idx:idx+k]
	
	if kmer in kmer_counts:
		count = kmer_counts[kmer] + 1
	else:
		count = 1		
	kmer_counts[kmer] = count
	
	if count > max_kmer_count:
		max_kmer_count = count
		max_kmers = kmer
	elif count == max_kmer_count:
		max_kmers = max_kmers + " " + kmer

inout.output(max_kmers)
