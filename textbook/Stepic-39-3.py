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
bases = inout.infilelines[2].strip().split(' ')

profile_matrix = []
for pos in range(3, k+3):
	probabilities = map(float, inout.infilelines[pos].strip().split(' '))
	profile_matrix.append(dict(zip(bases, probabilities)))
	
print(profile_matrix)

best_kmer = common.profile_most_probable_kmer(sequence, k, profile_matrix)
	
inout.output(best_kmer)