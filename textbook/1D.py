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

kmers_to_check = []
kmer_locations = {}
for idx in range(len(sequence) - k + 1):
	kmer = sequence[idx:idx+k]
	if kmer in kmer_locations:
		kmer_locations[kmer].append(idx)
	else:
		kmer_locations[kmer] = [idx]
		
	if len(kmer_locations[kmer]) >= t:
		kmers_to_check.append(kmer)
			
kmers_to_keep = []
for kmer in set(kmers_to_check):
	locations = kmer_locations[kmer]

	# once we've passed location[-t+1] then we can't have t locations in the window
	for first_loc_index, first_loc in enumerate(locations[:-t+1]):
		last_in_window = first_loc + L - k	
		locs_in_window = filter(lambda loc: loc <= last_in_window, locations[first_loc_index:])	 			
		if len(locs_in_window) >= t:
			kmers_to_keep.append(kmer)
			break
		
resultset = " ".join(set(kmers_to_keep))
inout.output(resultset)
