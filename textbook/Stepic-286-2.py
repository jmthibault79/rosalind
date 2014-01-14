# Input: A permutation P.
# Output: The sequence of permutations corresponding to applying GREEDYSORTING to P, ending with
# the identity permutation.

# Sample Input:
# (-3 +4 +1 +5 -2)

# Sample Output:
# (-1 -4 +3 +5 -2)
# (+1 -4 +3 +5 -2)
# (+1 +2 -5 -3 +4)
# (+1 +2 +3 +5 +4)
# (+1 +2 +3 -4 -5)
# (+1 +2 +3 +4 -5)
# (+1 +2 +3 +4 +5)

import inout
import common

permutation = common.greedysorting_parse(inout.infilelines[0].strip())
sequence = common.greedysorting(permutation)
sequence_out = common.greedysorting_out(sequence)

inout.output(sequence_out)