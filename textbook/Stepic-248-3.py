# Edit Distance Problem: Find the edit distance between two strings.
# Input: Two strings.
# Output: The edit distance between these strings.

# Sample Input:
# PLEASANTLY
# MEANLY

# Sample Output:
# 5

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

# Plan: use the maximum alignment score algorithm from 76-3, setting matches as 0 and mismatches and indels at -1
# Then the edit distance is just the inverse of the score

import string
scoring_matrix = common.mismatch_scoring_matrix(string.ascii_uppercase)
indel_penalty = -1

longest, backtrack_matrix = common.scored_longest_common_subsequence(scoring_matrix, indel_penalty, str1, str2)

inout.output(str(-longest))