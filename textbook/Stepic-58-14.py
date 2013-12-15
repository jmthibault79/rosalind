# Input: An integer d followed by a collection of paired k-mers PairedReads.
# Output: A string Text with (k, d)-mer composition equal to PairedReads.

# Sample Input:
#      2
#      GAGA|TTGA
#      TCGT|GATG
#      CGTG|ATGT
#      TGGT|TGAG
#      GTGA|TGTT
#      GTGG|GTGA
#      TGAG|GTTG
#      GGTC|GAGA
#      GTCG|AGAT

# Sample Output:
#      GTGGTCGTGAGATGTTGA
     
import inout
import common

d = int(inout.infilelines[0].strip())
pairs = map(str.strip, inout.infilelines[1:])

graph = common.parse_graph_from_pairs(pairs)
path = common.find_eulerian_path(graph)

inout.output(common.assemble_path_from_pairs(path, d))