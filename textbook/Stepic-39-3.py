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

sequence = inout.infilelines[0].strip() 
k = int(inout.infilelines[1].strip())

bases = {}
for idx, base in enumerate(inout.infilelines[2].strip().split(' ')):
	bases[base] = idx

matrix = []
for pos in range(3, k+3):
	matrix.append(map(float, inout.infilelines[pos].strip().split(' ')))

def probability(kmer):
	total = 0.0
	for idx, base in enumerate(kmer):
		total = total + matrix[idx][bases[base]]
	return total

best_kmer = ''
best_kmer_probability = 0
for idx in range(len(sequence) - k + 1):
	kmer = sequence[idx:idx+k]
	prob = probability(kmer)
	if prob > best_kmer_probability:
 		best_kmer = kmer
 		best_kmer_probability = prob
 		
inout.output(best_kmer)