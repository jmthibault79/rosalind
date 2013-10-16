# Clump Finding Problem
#
# Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger) string Genome if there is an interval of Genome of length L in which Pattern appears at least t times. For example, TGCA forms a (25,3)-clump in the following Genome: gatcagcataagggtcccTGCAaTGCAtgacaagccTGCAgttgttttac
#
# Find patterns forming clumps in a string.
# 
# Given: A string Genome, and integers k, L, and t.
#
# Return: All distinct k-mers forming (L, t)-clumps in Genome.

# Sample Dataset
#
# CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
# 5 75 4

# Sample Output
#
# CGACA GAAGA AATGT

import inout 	# my module for handling Rosalind's file I/O
sequence = inout.infilelines[0].strip()
k, L, t = map(int, inout.infilelines[1].strip().split())

kmers = []
for window_start in range(len(sequence) - L + 1):
	window_kmers = {}
	for offset in range(L - k + 1):
		kmer = sequence[window_start+offset:window_start+offset+k]
		if kmer in window_kmers:
			count = window_kmers[kmer] + 1
		else:
			count = 1		
		window_kmers[kmer] = count
	
		if count >= t:
			kmers.append(kmer)
			
resultset = " ".join(set(kmers))
inout.output(resultset)
