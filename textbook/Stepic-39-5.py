# GREEDYMOTIFSEARCH(Dna, k,t)
#        form a set of k-mers BestMotifs by selecting 1st k-mers in each string from Dna
#        for each k-mer Motif in the 1st string from Dna
#            Motif1 <- Motif
#            for i = 2 to t
#                form Profile from motifs Motif1, ..., Motifi - 1
#                Motifi <- Profile-most probable k-mer in the i-th string in Dna
#            Motifs <- (Motif1, ..., Motift)
#            if Score(Motifs) < Score(BestMotifs)
#                BestMotifs <- Motifs
#        output BestMotifs
     
# Input: Integers k and t, followed by a collection of strings Dna.

# Output: A collection of strings BestMotifs resulting from applying GREEDYMOTIFSEARCH(Dna,k,t). If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.

# Sample Input:
#     3 5
#     GGCGTTCAGGCA
#     AAGAATCAGTCA
#     CAAGGAGTTCGC
#     CACGTCAATCAC
#     CAATAATATTCG

# Sample Output:
#     CAG
#     CAG
#     CAA
#     CAA
#     CAA

import inout
import common

k,t = map(int, inout.infilelines[0].strip().split(' '))
sequences = map(str.strip, inout.infilelines[1:])

best_motifs = common.greedy_motif_search(sequences, k, t)

inout.output('\n'.join(best_motifs))