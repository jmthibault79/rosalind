# Input: Genomes P and Q.
# Output: The 2-break distance d(P, Q).

# Sample Input:
# (+1 +2 +3 +4 +5 +6)
# (+1 -3 -6 -5)(+2 -4)

# Sample Output:
# 3

import inout
import common

edges = common.breakpoints_parse(inout.infilelines[0].strip() + inout.infilelines[1].strip())

blocks = len(edges) / 2
cycles = common.breakpoints_cycles(edges)
distance = blocks - cycles

inout.output(str(distance))