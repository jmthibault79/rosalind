# MEDIANSTRING(Dna, k)
#        BestPattern <- AAA...AA
#        for each k-mer Pattern from AAA...AA to TTT...TT
#            if d(Pattern, Dna) < d(BestPattern, Dna)
#                 BestPattern <- Pattern
#        output BestPattern
#
#     Input: An integer k, followed by a collection of strings Dna.
#     Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern. (Return any one if there are multiple medium strings.)

# Sample Input:
#     3
#     AAATTGACGCAT
#     GACGACCACGTT
#     CGTCAGCGCCTG
#     GCTGAGCACCGG
#     AGTACGGGACAG

# Sample Output:
#     GAC
     
import inout

k = int(inout.infilelines[0].strip())
dna_lines = map(str.strip, inout.infilelines[1:])

bases = ['A', 'C', 'G', 'T']

def enumerate_mismatches(kmer, maxdist):
	if maxdist == 0:
		return [kmer]
	else:
		r = []
		for m_kmer in enumerate_mismatches(kmer, maxdist - 1):
			for loc in range(k):
				for base in bases:
					new_kmer = m_kmer[:loc] + base + m_kmer[loc + 1:]
					r.append(new_kmer)
		return set(r)

def mindist(kmer, sequence):
	for dist in range(k+1):
		for new_kmer in enumerate_mismatches (kmer, dist):
			if new_kmer in sequence:
				return dist
				
	# shouldn't reach here
	return k + 1

def enumerate_all_kmers(k):
	if k == 1:
		return bases
	else:
		kmers = []
		for k_minus_one_mer in enumerate_all_kmers(k - 1):
			for base in bases:
				kmers.append(k_minus_one_mer + base)
		return kmers

best_kmer = ''
best_kmer_dist = 100000000
for kmer in enumerate_all_kmers(k):
	total_dist = 0
	for sequence in dna_lines:
		total_dist = total_dist + mindist(kmer, sequence)
	if total_dist < best_kmer_dist:
		best_kmer = kmer
		best_kmer_dist = total_dist
		
inout.output(best_kmer)