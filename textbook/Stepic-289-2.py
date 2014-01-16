# Given two strings, find all their shared k-mers.
# Input: An integer k and two strings.
# Output: All k-mers shared by these strings, in the form of ordered pairs (x, y).

# Sample Input:
# 3
# AAACTCATC
# TTTCAAATC

# Sample Output:
# (0, 4)
# (0, 0)
# (4, 2)
# (6, 6)

import inout
import common

k = int(inout.infilelines[0].strip())
str1, str2 = map(str.strip, inout.infilelines[1:3])

result = common.shared_kmers(k, str1, str2)

def output_one_pair(pair):
    return '({}, {})'.format(pair[0], pair[1])

inout.output('\n'.join(map(output_one_pair, result)))