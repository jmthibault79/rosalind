# Frequent Words with Mismatches Problem
#
# Find the most frequent k-mers with mismatches in a string.
#
# Given: A string Text as well as integers k and d.
#
# Return: All most frequent k-mers with up to d mismatches in Text.

# Sample Dataset
#
# ACGTTGCATGTCGCATGATGCATGAGAGCT
# 4 1

# Sample Output
#
# GATG ATGC ATGT

import inout 	# my module for handling Rosalind's file I/O
line1 = inout.infilelines[0].strip().split()
sequence, k, d = line1[0], int(line1[1]), int(line1[2])

def enumerate_mismatches (kmer, maxdist):
	if maxdist == 0:
		return [kmer]
	else:
		r = []
		for m_kmer in enumerate_mismatches(kmer, maxdist - 1):
			for loc in range(k):
				for base in ['A', 'C', 'G', 'T']:
					new_kmer = '{}{}{}'.format(m_kmer[:loc], base, m_kmer[loc + 1:])
					r.append(new_kmer)
		return set(r)

kmer_counts = {}
max_kmers = []
max_kmer_count = 0
for idx in range(len(sequence) - k + 1):
	kmer = sequence[idx:idx+k]

	for m_kmer in enumerate_mismatches(kmer, d):
		if m_kmer in kmer_counts:
			count = kmer_counts[m_kmer] + 1
		else:
			count = 1
		kmer_counts[m_kmer] = count

		if count > max_kmer_count:
			max_kmer_count = count
			max_kmers = [m_kmer]
		elif count == max_kmer_count:
			max_kmers.append(m_kmer)

inout.output(' '.join(max_kmers))
