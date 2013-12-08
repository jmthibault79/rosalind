# Input: A collection Patterns of k-mers.
# Output: The overlap graph Overlap(Patterns), in the form of an adjacency list.

# Sample Input:
#      ATGCG
#      GCATG
#      CATGC
#      AGGCA
#      GGCAT

# Sample Output:
#      AGGCA -> GGCAT
#      CATGC -> ATGCG
#      GCATG -> CATGC
#      GGCAT -> GCATG
     
import inout
import common

sequences = map(str.strip, inout.infilelines)
			
inout.output('\n'.join(map(common.overlap_to_str, common.overlap_graph(sequences))))