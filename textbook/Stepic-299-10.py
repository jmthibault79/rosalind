# Inverse Burrows-Wheeler Transform Problem: Reconstruct a string from its Burrows-Wheeler transform.
# Input: A string Transform (with a single "$" symbol).
# Output: The string Text such that BWT(Text) = Transform.

# Sample Input:
# TTCCTAACG$A

# Sample Output:
# TACATCACGT$

import inout
import common

text = inout.infilelines[0].strip()

original = common.inv_bwt(text)

inout.output(original)