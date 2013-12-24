# Input: An integer money and an array coins = (coin1, ..., coind).
# Output: The minimum number of coins with denominations coins that changes money.

# Sample Input:
# 40
# 50,25,20,10,5,1

# Sample Output:
# 2

import inout
import common

change = int(inout.infilelines[0].strip())
coins = map(int, inout.infilelines[1].strip().split(','))

numcoins = common.make_change(change, coins)

inout.output(str(numcoins))