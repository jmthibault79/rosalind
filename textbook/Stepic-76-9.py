# Input: Two protein strings written in the single-letter amino acid alphabet.
# Output: The maximum score of a local alignment of the strings, followed by a local alignment of these
# strings achieving the maximum score. Use the PAM250 scoring matrix and indel penalty sigma = 5.

# Sample Input:
# MEANLY
# PENALTY

# Sample Output:
# 15
# EANL-Y
# ENALTY

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

scoring_matrix = common.parse_scoring_matrix(inout.readlines('PAM250_1.txt'))
indel_penalty = -5

longest, backtrack_matrix, best_row, best_col = common.scored_longest_common_subsequence_local(scoring_matrix, indel_penalty, str1, str2)
aligned1, aligned2 = common.output_longest_common_subsequence_local(backtrack_matrix, str1, str2, best_row, best_col)

inout.output('{}\n{}\n{}'.format(longest, aligned1, aligned2))