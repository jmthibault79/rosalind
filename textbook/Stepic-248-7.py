# Input: Two strings v and w, each of length at most 1000.
# Output: The score of an optimal overlap alignment of v and w, followed by an alignment of a suffix v' of
# v and a prefix w' of w achieving this maximum score. Use an alignment score in which matches count
# +1 and both the mismatch and indel penalties are 2.

# Sample Input:
# PAWHEAE
# HEAGAWGHEE

# Sample Output:
# 1
# HEAE
# HEAG

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

import string
scoring_matrix = common.mismatch_scoring_matrix_overlap(string.ascii_uppercase)
indel_penalty = -2

longest, backtrack_matrix, best_row, best_col = common.scored_longest_common_subsequence_overlap(scoring_matrix, indel_penalty, str1, str2)
aligned1, aligned2 = common.output_longest_common_subsequence_local(backtrack_matrix, str1, str2, best_row, best_col)

inout.output('{}\n{}\n{}'.format(longest, aligned1, aligned2))