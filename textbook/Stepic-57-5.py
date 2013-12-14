# Input: The adjacency list of a directed graph that has an Eulerian path.
# Output: An Eulerian path in this graph.

# Sample Input:
#      0 -> 2
#      1 -> 3
#      2 -> 1
#      3 -> 0,4
#      6 -> 3,7
#      7 -> 8
#      8 -> 9
#      9 -> 6

# Sample Output:
#      6->7->8->9->6->3->0->2->1->3->4
           
import inout
import common

edge_strs = map(str.strip, inout.infilelines)

graph = common.parse_graph_edges(edge_strs)
path = common.find_eulerian_path(graph)

inout.output('->'.join(path))