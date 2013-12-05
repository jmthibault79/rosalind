# RANDOMIZEDMOTIFSEARCH(Dna, k, t)
#         randomly select k-mers Motifs = (Motif1, ..., Motift) in each string from Dna
#         BestMotifs <- Motifs
#         while forever
#             Profile <- Profile(Motifs)
#             Motifs <- Motifs(Profile, Dna)
#             if Score(Motifs) < Score(BestMotifs)
#                 BestMotifs <- Motifs
#             else
#                 output BestMotifs
#                 return
#                 
# Input: Integers k and t, followed by a collection of strings Dna.
# Output: A collection BestMotifs resulting from running RANDOMIZEDMOTIFSEARCH(Dna, k, t) 1000 times.
# Remember to use pseudocounts!

# Sample Input:
#      8 5
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

k,t = map(int, inout.infilelines[0].strip().split(' '))
sequences = map(str.strip, inout.infilelines[1:])

best_motifs = common.repeated_randomized_motif_search(sequences, k, t, 1000)

inout.output('\n'.join(best_motifs))