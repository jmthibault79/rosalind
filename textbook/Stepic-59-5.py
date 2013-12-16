# Contig Generation Problem: Generate the contigs from a collection of reads (with imperfect coverage).
#      Input: A collection of k-mers Patterns.
#      Output: All contigs in DeBruijn(Patterns).

# Sample Input:
#      ATG
#      ATG
#      TGT
#      TGG
#      CAT
#      GGA
#      GAT
#      AGA

# Sample Output:
#      AGA ATG ATG CAT GAT TGGA TGT

import inout
import common

reads = map(str.strip, inout.infilelines)

graph = common.debruijn_graph_with_duplicates(reads)
paths = common.debruijn_to_contigs(graph)		
contigs = map(common.assemble_path, paths)		

inout.output(' '.join(contigs))