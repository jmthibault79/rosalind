# Input: A collection of strings Patterns.
# Output: The adjacency list corresponding to Trie(Patterns), in the following format. If Trie(Patterns) has
# n nodes, first label the root with 1 and then label the remaining nodes with the integers 2 through n in
# any order you like. Each edge of the adjacency list of Trie(Patterns) will be encoded by a triple: the first
# two members of the triple must be the integers labeling the initial and terminal nodes of the edge,
# respectively; the third member of the triple must be the symbol labeling the edge.

# Sample Input:
# GGTA
# CG
# GGC

# Sample Output:
# 1 2 G
# 2 3 G
# 3 4 T
# 4 5 A

# 3 6 C

# 1 7 C
# 7 8 G

import inout
import common

strings = map(str.strip, inout.infilelines)

trie = common.create_trie(strings)
trie_out = common.output_trie(trie)

inout.output(trie_out)