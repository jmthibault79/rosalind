# Input: Two amino acid strings v and w (each of length at most 100).
# Output: The maximum alignment score between v and w, followed by an alignment of v and w
# achieving this maximum score. Use the BLOSUM62 scoring matrix, a gap opening penalty of 11, and
# a gap extension penalty of 1.

# Sample Input:
# PRTEINS
# PRTWPSEIN

# Sample Output:
# 8
# PRT---EINS
# PRTWPSEIN-

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

scoring_matrix = common.parse_scoring_matrix(inout.readlines('BLOSUM62.txt'))
gap_open = -11
gap_extension = -1

longest, backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, best_row, best_col = common.scored_longest_common_subsequence_affine(scoring_matrix, gap_open, gap_extension, str1, str2)
aligned1, aligned2 = common.output_longest_common_subsequence_affine_middle(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, str1, str2, best_row, best_col)

inout.output('{}\n{}\n{}'.format(longest, aligned1, aligned2))