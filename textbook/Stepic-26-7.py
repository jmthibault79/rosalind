# Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.

# Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).

# Sample Input:
#     20
#     60
#     57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493

# Sample Output:
#     99-71-137-57-72-57

import inout
import collections

M = int(inout.infilelines[0].strip())
N = int(inout.infilelines[1].strip())
spectrum = map(int, inout.infilelines[2].strip().split(' '))

def most_frequent_convolution(m, spectrum):
	convolution = []
	l = len(spectrum)
	for i in range(l):
		for j in range(i+1, l):
			diff = abs(spectrum[i]-spectrum[j])
			if diff >= 57 and diff <= 200:
				convolution.append(diff)
				
	conv_counter = collections.Counter(convolution)
	
	conv_leaderboard = {}
	for conv in conv_counter:
		count = conv_counter[conv]
		if count in conv_leaderboard:
			conv_leaderboard[count].append(conv)
		else:
			conv_leaderboard[count] = [conv]
	
	survivors = []
	survivors_to_choose = m		
	for count in sorted(conv_leaderboard.keys(), reverse=True):
		survivors.extend(conv_leaderboard[count])
		survivors_to_choose = survivors_to_choose - len(conv_leaderboard[count])
		if survivors_to_choose < 0:
			break

	return survivors

amino_acids = most_frequent_convolution(M, spectrum)

def cyclic_spectrum(peptide):
	out_spectrum = [0, sum(peptide)]

	peptide_2 = peptide + peptide	# for easy cyclic access
	for k in range(1, len(peptide)):
		for n in range(len(peptide)):
			subpep = peptide_2[n:n+k]
			out_spectrum.append(sum(subpep))
	return sorted(out_spectrum)

def branch(peptides):
	out_peptides = []
	for p in peptides:
		for amino_acid in amino_acids:
			out_peptides.append(p + [amino_acid])
	return out_peptides

def score(candidate, target):
	import collections
	c_spectrum = cyclic_spectrum(candidate)
	c_counter = collections.Counter(c_spectrum)
	t_counter = collections.Counter(target)
	
	s = 0
	for mass in c_counter:
		if mass in t_counter:
			s = s + min(c_counter[mass],t_counter[mass])
	
	return s

def cut(candidates, target, n):
	if len(candidates) <= n:
		return candidates
		
	leaderboard = {}
	for candidate in candidates:
		s = score(candidate, target)
		if s in leaderboard:
			leaderboard[s].append(candidate)
		else:
			leaderboard[s] = [candidate]
	
	survivors = []
	survivors_to_choose = n		
	for s in sorted(leaderboard.keys(), reverse=True):
		survivors.extend(leaderboard[s])
		survivors_to_choose = survivors_to_choose - len(leaderboard[s])
		if survivors_to_choose < 0:
			break
			
	return survivors
	
# I'm sure there's a better way to do this but I don't know enough Python yet
def mklist(item):
	return [item]
		
candidates = map(mklist,amino_acids)
winner = ''
winner_score = 0
while candidates:
	candidates = branch(candidates)
	new_candidates = []
	for candidate in candidates:
		c_mass = sum(candidate)
		t_mass = max(spectrum)

		# if the mass of the candidate peptide equals the mass of the target peptide
		if c_mass == t_mass:
			new_candidates.append(candidate)
			c_score = score(candidate, spectrum)
			if c_score > winner_score:
				winner = candidate
				winner_score = c_score
		elif c_mass < t_mass:
			new_candidates.append(candidate)
		# else: the candidate mass is too large, so it does not go on to the next round
	candidates = cut(new_candidates, spectrum, N) 
	
inout.output('-'.join(map(str,winner)))