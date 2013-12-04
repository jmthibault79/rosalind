# GREEDYMOTIFSEARCH with pseudocounts

# Sample Input:
#     3 5
#     GGCGTTCAGGCA
#     AAGAATCAGTCA
#     CAAGGAGTTCGC
#     CACGTCAATCAC
#     CAATAATATTCG

# Sample Output:
#      TTC
#      ATC
#      TTC
#      ATC
#      TTC

import inout
import common

k,t = map(int, inout.infilelines[0].strip().split(' '))
sequences = map(str.strip, inout.infilelines[1:])

best_motifs = common.greedy_motif_search_with_pseudocounts(sequences, k, t)

inout.output('\n'.join(best_motifs))