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

sequences = map(str.strip, inout.infilelines)

def overlap(a, b):
	return a[1:] == b[:-1]
	
def overlap_to_str(o):
	return "{} -> {}".format(*o)

overlap_graph = []
for s1 in sequences:
	for s2 in sequences:
		if overlap(s1, s2):
			overlap_graph.append([s1, s2])
			
inout.output('\n'.join(map(overlap_to_str, overlap_graph)))