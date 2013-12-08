# Input: An integer k and a string Text.
# Output: Compositionk(Text), where the k-mers are written in lexicographic order.

# Sample Input:
#      5
#      CAATCCAAC

# Sample Output:
#      AATCC
#      ATCCA
#      CAATC
#      CCAAC
#      TCCAA

import inout
import common

k = int(inout.infilelines[0].strip())
sequence = inout.infilelines[1].strip()

kmers = sorted(common.all_kmers(sequence, k))

inout.output('\n'.join(kmers))