# Input: An integer representing the source node of a graph, followed by an integer representing the sink
# node of the graph, followed by a list of edges in the graph. The edge notation 0->1:7 indicates that
# an edge connects node 0 to node 1 with weight 7.
# Output: The length of a longest path in the graph, followed by a longest path.

# Sample Input:
# 0
# 4
# 0->1:7
# 0->2:4
# 2->3:2
# 1->4:1
# 3->4:3

# Sample Output:
# 9
# 0->2->3->4

import inout
import common

source = inout.infilelines[0].strip()
sink = inout.infilelines[1].strip()
edges = map(str.strip, inout.infilelines[2:])

dag = common.parse_dag_edges(edges)
ordering = common.wikipedia_depth_first_topological_sort(dag, sink)

weight, backtrack = common.longest_dag_weight(dag, ordering, source, sink)
path = common.output_longest_dag_path(backtrack, source, sink)

inout.output('{}\n{}'.format(weight,path))