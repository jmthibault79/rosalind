# Input: Two strings s and t.
# Output: A longest common subsequence of s and t.
#
# Note: If more than one LCS exists, you may return any one.

# Sample Input:
# AACCTTGG
# ACACTGTGA

# Sample Output:
# AACTGG

import inout
import common

str1 = inout.infilelines[0].strip()
str2 = inout.infilelines[1].strip()

# output_longest_common_subsequence was hitting the default limit of 1000 for the test dataset
# https://class.coursera.org/bioinformatics-001/forum/thread?thread_id=742
import sys
sys.setrecursionlimit(2000)

longest, backtrack_matrix = common.longest_common_subsequence(str1, str2)
lcs = common.output_longest_common_subsequence(backtrack_matrix, str1, len(str1), len(str2))

inout.output(lcs)