# Input: Two nucleotide strings v and w, where v has length at most 1000 and w has length at most 100.
# Output: A highest-scoring fitting alignment between v and w. Use the simple scoring method in which
# matches count +1 and both the mismatch and indel penalties are 1.

# Sample Input:
# GTAGGCTTAAGGTTA
# TAGATA

# Sample Output:
# 2
# TAGGCTTA
# TAGA--TA


import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

scoring_matrix = common.mismatch_scoring_matrix_fitted("ACGT")
indel_penalty = -1

longest, backtrack_matrix, best_row, best_col = common.scored_longest_common_subsequence_fitted(scoring_matrix, indel_penalty, str1, str2)
aligned1, aligned2 = common.output_longest_common_subsequence_local(backtrack_matrix, str1, str2, best_row, best_col)

inout.output('{}\n{}\n{}'.format(longest, aligned1, aligned2))