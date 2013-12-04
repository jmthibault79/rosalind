# Functions for Stepic assignments

# 39-3
def profile_probability(kmer, bases, profile):
	total = 0.0
	for idx, base in enumerate(kmer):
		total = total + profile[idx][bases[base]]
	return total

# 39-3
def profile_most_probable_kmer(sequence, k, bases, profile):
	best_kmer = ''
	best_kmer_probability = 0
	for idx in range(len(sequence) - k + 1):
		kmer = sequence[idx:idx+k]
		prob = profile_probability(kmer, bases, profile)
		if prob > best_kmer_probability:
 			best_kmer = kmer
 			best_kmer_probability = prob
 	return best_kmer
