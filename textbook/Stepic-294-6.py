# Input: A string Text and a collection of strings Patterns.
# Output: All starting positions in Text where a string from Patterns appears as a substring.

# Sample Input:
# AATCGGGTTCAATCGGGGT
# ATCG
# GGGT

# Sample Output:
# 1 4 11 15

import inout
import common

text = inout.infilelines[0].strip()
strings = map(str.strip, inout.infilelines[1:])

trie = common.create_trie(strings)
matches = common.match_trie(trie, text)

inout.output(' '.join(map(str, matches)))