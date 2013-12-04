# Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
#     Input: A string Text, an integer k, and a k * 4 matrix Profile.
#     Output: A Profile-most probable k-mer in Text.

# Sample Input:
#     ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
#     5
#     A C G T
#     0.2 0.4 0.3 0.1
#     0.2 0.3 0.3 0.2
#     0.3 0.1 0.5 0.1
#     0.2 0.5 0.2 0.1
#     0.3 0.1 0.4 0.2

# Sample Output:
#     CCGAG
     
import inout
import common

sequence = inout.infilelines[0].strip() 
k = int(inout.infilelines[1].strip())

bases = {}
for idx, base in enumerate(inout.infilelines[2].strip().split(' ')):
	bases[base] = idx

matrix = []
for pos in range(3, k+3):
	matrix.append(map(float, inout.infilelines[pos].strip().split(' ')))

best_kmer = common.profile_most_probable_kmer(sequence, k, bases, matrix)
	
inout.output(best_kmer)