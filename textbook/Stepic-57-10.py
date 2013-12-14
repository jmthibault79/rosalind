# Input: An integer k.
# Output: A k-universal circular string.

# Sample Input:
#      4

# Sample Output:
#      0000110010111101

import inout
import common

k = int(inout.infilelines[0])

graph = common.k_universal_graph(k)
cycle = common.find_eulerian_cycle(graph)

inout.output(common.assemble_path(cycle[:-1])[:2**k])