# Constructing Suffix Array Problem: Construct the suffix array of a string.
# Input: A string Text.
# Output: SuffixArray(Text).

# Sample Input:
# AACGATAGCGGTAGA$

# Sample Output:
# 15, 14, 0, 1, 12, 6, 4, 2, 8, 13, 3, 7, 9, 10, 11, 5

import inout
import common

text = inout.infilelines[0].strip()

array = common.create_suffix_array(text)

inout.output(', '.join(map(str, array.values())))