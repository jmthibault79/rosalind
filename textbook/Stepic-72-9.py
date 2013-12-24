# Input: Integers n and m, followed by an n * (m + 1) matrix down and an (n + 1) * m matrix right.
# The two matrices are separated by the - symbol.
# Output: The length of a longest path from source (0, 0) to sink (n, m) in the n * m rectangular grid whose
# edges are defined by the matrices down and right.

# Sample Input:
# 4
# 4
# 1 0 2 4 3
# 4 6 5 2 1
# 4 4 5 2 1
# 5 6 8 5 3
# -
# 3 2 4 0
# 3 2 4 2
# 0 7 3 3
# 3 3 0 2
# 1 3 2 2

# Sample Output:
# 34

import inout
import common

n = int(inout.infilelines[0].strip())
m = int(inout.infilelines[1].strip())
if len(inout.infilelines) != 4 + 2 * n:
    raise Exception('Expected {} input lines based on n={}, saw {}'.format(4 + 2 * n, n, inout.infilelines))
downmatrix = common.parse_matrix(map(str.strip, inout.infilelines[2:2 + n]), n, m + 1)
rightmatrix = common.parse_matrix(map(str.strip,inout.infilelines[2 + n + 1:4 + 2 * n]), n + 1, m)
if inout.infilelines[2 + n].strip() != '-':
    raise Exception('Expected - ({}) separating downmatrix from rightmatrix, saw {} ({})'.format(ord('-'), inout.infilelines[2 + n], ord(inout.infilelines[2 + n])))

longest = common.longest_path(n, m, downmatrix, rightmatrix)

inout.output(str(longest))