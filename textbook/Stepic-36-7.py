# Implement MOTIFENUMERATION (reproduced below).
#     Input: Integers k and d, followed by a collection of strings Dna.
#     Output: All (k, d)-motifs in Dna.
#    MOTIFENUMERATION(Dna, k, d)
#        for each k-mer a in Dna
#            for each k-mer a* differing from a by at most d mutations
#                if a* appears in each string from Dna with at most d mutations
#                    output a*

# Sample Input:
#     3 1
#     ATTTGGC
#     TGCCTTA
#     CGGTATC
#     GAAAATT

# Sample Output:
#     ATA ATT GTT TTT
     
import inout
import collections

k,d = map(int, inout.infilelines[0].strip().split(' '))
dna_lines = map(str.strip, inout.infilelines[1:])

def enumerate_mismatches (kmer, maxdist):
	if maxdist == 0:
		return [kmer]
	else:
		r = []
		for m_kmer in enumerate_mismatches(kmer, maxdist - 1):
			for loc in range(k):
				for base in ['A', 'C', 'G', 'T']:
					new_kmer = m_kmer[:loc] + base + m_kmer[loc + 1:]
					r.append(new_kmer)
		return set(r)

def motifs (sequence, k, d):
	motifs = []
	for idx in range(len(sequence) - k + 1):
		kmer = sequence[idx:idx+k]
		motifs.extend(enumerate_mismatches(kmer, d))
	return motifs
		
common_motifs = motifs(dna_lines[0], k, d)
for line in dna_lines[1:]:
	common_motifs = set(common_motifs) & set(motifs(line, k, d))

inout.output(' '.join(common_motifs))