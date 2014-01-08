# Solve the Middle Edge in Linear Space Problem (for protein strings). Use the BLOSUM62 scoring matrix and a linear indel penalty equal to 5.
# Input: Two amino acid strings.
# Output: A middle edge in the alignment graph in the form (i, j) (k, l), where (i, j) connects to (k, l).
# To compute scores, use the BLOSUM62 scoring matrix and a (linear) indel penalty equal to 5.

# Sample Input:
# PLEASANTLY
# MEASNLY

# Sample Output:
# (4, 3) (5, 4)

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

scoring_matrix = common.parse_scoring_matrix(inout.readlines('BLOSUM62.txt'))
indel_penalty = -5

from_row, from_col, to_row, to_col = common.alignment_middle_edge(scoring_matrix, indel_penalty, str1, str2)
inout.output('({}, {}) ({}, {})'.format(from_row, from_col, to_row, to_col))