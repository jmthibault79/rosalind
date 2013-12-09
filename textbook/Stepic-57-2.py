# EULERIANCYCLE(Graph)
#         form a cycle Cycle by randomly walking in Graph (never visit an edge twice!)
#         while there are unexplored edges
#             select a node newStart in Cycle with still unexplored edges
#             form Cycle* by traversing Cycle (starting at newStart) and randomly walking
#             Cycle <- Cycle*
#         return Cycle

# Input: The adjacency list of an Eulerian directed graph.
# Output: An Eulerian cycle in this graph.

# Sample Input:
#      0 -> 3
#      1 -> 0
#      2 -> 1,6
#      3 -> 2
#      4 -> 2
#      5 -> 4
#      6 -> 5,8
#      7 -> 9
#      8 -> 7
#      9 -> 6

# Sample Output:
#      6->8->7->9->6->5->4->2->1->0->3->2->6
      
import inout
import common

edge_strs = map(str.strip, inout.infilelines)

graph = common.parse_graph_edges(edge_strs)
cycle = common.find_eulerian_cycle(graph)

inout.output('->'.join(cycle))