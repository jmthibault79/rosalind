# GIBBSSAMPLER(Dna, k, t, N)
#         randomly select k-mers Motifs = (Motif1, ..., Motift) in each string from Dna
#         BestMotifs <- Motifs
#         for i from 1 to N
#             i <- Random(t)
#             construct profile matrix Profile from all strings in Motifs except for Motifi
#             Motifi <- Profile-randomly generated k-mer in the i-th sequence
#             if Score(Motifs) < Score(BestMotifs)
#                 BestMotifs <- Motifs
#         output BestMotifs
#      
# Input: Integers k, t, and N, followed by a collection of strings Dna.
# Output: The strings BestMotifs resulting from running GIBBSSAMPLER(Dna, k, t, N) with
# 20 random starts. Remember to use pseudocounts!

# Sample Input:
#      8 5 100
#      CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
#      GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
#      TAGTACCGAGACCGAAAGAAGTATACAGGCGT
#      TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
#      AATCCACCAGCTCCACGTGCAATGTTGGCCTA

# Sample Output:
#      TCTCGGGG
#      CCAAGGTG
#      TACAGGCG
#      TTCAGGTG
#      TCCACGTG

import inout
import common

k,t,n = map(int, inout.infilelines[0].strip().split(' '))
sequences = map(str.strip, inout.infilelines[1:])

best_motifs = common.repeated_gibbs_sampler(sequences, k, t, n, 20)

inout.output('\n'.join(best_motifs))