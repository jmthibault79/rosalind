#  LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)
#        Leaderboard <- {0-peptide}
#        LeaderPeptide <- 0-peptide
#        while Leaderboard is non-empty
#            Leaderboard <- Expand(Leaderboard)
#            for each Peptide in Leaderboard
#                if Mass(Peptide) = ParentMass(Spectrum)
#                    if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
#                        LeaderPeptide <- Peptide
#                else if Mass(Peptide) > ParentMass(Spectrum)
#                    remove Peptide from Leaderboard
#            Leaderboard <- Cut(Leaderboard, Spectrum, N)
#        output LeaderPeptide

# Input: Integer N and a collection of integers Spectrum.

# Output: LeaderPeptide after running LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N).

# Sample Input:
#     10
#     0 71 113 129 147 200 218 260 313 331 347 389 460

# Sample Output:
#     113-147-71-129

# Alternate solution without using amino acids as letters.
# Necessary for adaptation to nonstandard amino acids (Stepic 26-7).

# peptides are now represented as lists of ints

import inout

N = int(inout.infilelines[0].strip())
spectrum = map(int, inout.infilelines[1].strip().split(' '))

amino_acids = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]		

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
		t_mass = spectrum[-1]

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