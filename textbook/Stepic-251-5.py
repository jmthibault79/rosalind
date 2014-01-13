# In the Multiple Longest Common Subsequence Problem, the score of a column of the alignment matrix is equal to 1
# if all of the column's symbols are identical, and 0 if even one symbol disagrees.

# Input: Three DNA strings.
# Output: The length of a longest common subsequence of these three strings, followed by a multiple
# alignment of the three strings corresponding to such an alignment.

# Sample Input:
# ATATCCG
# TCCGA
# ATGTACTG

# Sample Output:
# 3
# ATATCC-G-
# ---TCC-GA
# ATGTACTG-

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()
str3 = inout.infilelines[2].strip()

scoring_matrix = common.mismatch_scoring_matrix_fitted("ACGT")
indel_penalty = -1

score, backtrack_matrix, best_x, best_y, best_z = common.align_3(str1, str2 ,str3)
aligned1, aligned2, aligned3 = common.output_align_3(backtrack_matrix, str1, str2 ,str3, best_x, best_y, best_z)

inout.output('{}\n{}\n{}\n{}'.format(score, aligned1, aligned2, aligned3))