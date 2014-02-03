# Input: A string BWT(Text), followed by a collection of Patterns.
# Output: A list of integers, where the i-th integer corresponds to the number of substring matches of the
# i-th member of Patterns in Text.

# Sample Input:
# TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC
# CCT CAC GAG CAG ATC

# Sample Output:
# 2 1 1 0 1

import inout
import common

bwt_text = inout.infilelines[0].strip()
patterns = inout.infilelines[1].strip().split(' ')

counts = ''
for pattern in patterns:
    counts += str(common.bwt_matching(bwt_text, pattern)) + ' '

inout.output(counts.strip())