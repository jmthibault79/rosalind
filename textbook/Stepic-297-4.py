# Burrows-Wheeler Transform Construction Problem: Construct the Burrows-Wheeler transform of a string.
# Input: A string Text.
# Output: BWT(Text).

# Sample Input:
# GCGTGCCTGGTCA$

# Sample Output:
# ACTGGCT$TGCGGC

import inout
import common

text = inout.infilelines[0].strip()

bwt = common.bwt(text)

inout.output(bwt)