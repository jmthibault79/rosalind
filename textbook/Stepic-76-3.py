# Input: Two protein strings written in the single-letter amino acid alphabet.
# Output: The maximum alignment score of these strings followed by an alignment achieving this
# maximum score. Use the BLOSUM62 scoring matrix and indel penalty sigma = 5.

# Sample Input:
# PLEASANTLY
# MEANLY

# Sample Output:
# 8
# PLEASANTLY
# -MEA--N-LY

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

scoring_matrix = common.parse_scoring_matrix(inout.readlines('BLOSUM62.txt'))
indel_penalty = -5

longest, backtrack_matrix = common.scored_longest_common_subsequence(scoring_matrix, indel_penalty, str1, str2)
aligned1, aligned2 = common.output_longest_common_subsequence_aligned(backtrack_matrix, str1, str2, len(str1), len(str2))

inout.output('{}\n{}\n{}'.format(longest, aligned1, aligned2))