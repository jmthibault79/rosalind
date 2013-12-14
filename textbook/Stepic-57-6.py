# Input: The adjacency list of a directed graph that has an Eulerian path.
# Output: An Eulerian path in this graph.

# Sample Input:
#      CTT -> TTA
#      ACC -> CCA
#      TAC -> ACC
#      GGC -> GCT
#      GCT -> CTT
#      TTA -> TAC

# Sample Output:
#      GGCTTACCA           

import inout
import common

edge_strs = map(str.strip, inout.infilelines)

graph = common.parse_graph_edges(edge_strs)
path = common.find_eulerian_path(graph)

inout.output(common.assemble_path(path))