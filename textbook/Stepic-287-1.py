# Number of Breakpoints Problem: Find the number of breakpoints in a permutation.
# Input: A permutation P.
# Output: The number of breakpoints in P.

# Sample Input:
# (+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14)

# Sample Output:
# 8

import inout
import common

permutation = common.greedysorting_parse(inout.infilelines[0].strip())

inout.output(str(common.count_breakpoints(permutation)))