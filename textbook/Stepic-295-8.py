# Longest Repeat Problem: Find the longest repeat in a string.
# Input: A string Text.
# Output: A longest repeat in Text, i.e., a longest substring of Text that appears in Text more than once.

# Sample Input:
# ATATCGTTTTATCGTT

# Sample Output:
# TATCGTT

import inout
import common

text = inout.infilelines[0].strip()

trie = common.create_suffix_trie(text, 100)
substring = common.find_longest_substring_in_suffix_trie(trie, 1, '')

inout.output(substring)