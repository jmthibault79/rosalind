# Implement LINEARSPACEALIGNMENT to solve the Global Alignment Problem for a large dataset.
# Input: Two long (10000 amino acid) protein strings written in the single-letter amino acid alphabet.
# Output: The maximum alignment score of these strings, followed by an alignment achieving this
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

score, alignment1, alignment2 = common.linear_space_alignment(scoring_matrix, indel_penalty, str1, str2)
inout.output('{}\n{}\n{}'.format(str(score), alignment1, alignment2))