# Input: An integer k and a string Text.
# Output: DeBruijnk(Text).

# Sample Input:
#      4
#      AAGATTCTCTAC

# Sample Output:
#      AAG -> AGA
#      AGA -> GAT
#      ATT -> TTC
#      CTA -> TAC
#      CTC -> TCT
#      GAT -> ATT
#      TCT -> CTA,CTC
#      TTC -> TCT
     
import inout
import common

k = int(inout.infilelines[0].strip())
sequence = inout.infilelines[1].strip()

graph = common.ordered_debruijn_graph(common.all_kmers(sequence, k))		

graph_strs = []
for k,v in graph.iteritems():
	graph_strs.append(common.debruijn_to_str(k,v))

inout.output('\n'.join(graph_strs))