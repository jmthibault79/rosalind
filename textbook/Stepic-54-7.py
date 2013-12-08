# Input: A collection of k-mers Patterns.
# Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).

# Sample Input:
#      GAGG
#      GGGG
#      GGGA
#      CAGG
#      AGGG
#      GGAG

# Sample Output:
#      AGG -> GGG
#      CAG -> AGG
#      GAG -> AGG
#      GGA -> GAG
#      GGG -> GGA,GGG    
 
import inout
import common

sequences = map(str.strip, inout.infilelines)

graph = common.debruijn_graph(sequences)		

graph_strs = []
for k,v in graph.iteritems():
	graph_strs.append(common.debruijn_to_str(k,v))

inout.output('\n'.join(graph_strs))