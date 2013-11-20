#  CYCLOPEPTIDESEQUENCING(Spectrum)
#        List <- {0-peptide}
#        while List is nonempty
#            List <- Expand(List)
#            for each peptide Peptide in List
#                if Cyclospectrum(Peptide) = Spectrum
#                    output Peptide
#                    remove Peptide from List
#                else if Peptide is not consistent with Spectrum
#                    remove Peptide from List

# Sample Input:
#     0 113 128 186 241 299 314 427

# Sample Output:
#     186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186

#W IL KQ

import inout
import peptide

spectrum = map(int, inout.infilelines[0].strip().split(' '))

def branch(peptides):
	out_peptides = []
	for p in peptides:
		for amino_acid in peptide.amino_acids:
			out_peptides.append(p + amino_acid)
	return out_peptides

def consistent(candidate, target):
	import collections
	c_counter = collections.Counter(candidate)
	t_counter = collections.Counter(target)
	
	for mass in c_counter:
		if mass not in t_counter:
			return False
		if c_counter[mass] > t_counter[mass]:
			return False	
	return True

def output_format(pep):
	masses = []
	for amino_acid in pep:
		masses.append(peptide.mass_table[amino_acid])
	return '-'.join(map(str,masses))
		
candidates = peptide.amino_acids
winners = []
while candidates:
	candidates = branch(candidates)
	new_candidates = []
	for candidate in candidates:
		c_spectrum = peptide.cyclic_spectrum(candidate)
		l_spectrum = peptide.linear_spectrum(candidate)
		if c_spectrum == spectrum:
			winners.append(candidate)
		elif consistent(l_spectrum, spectrum):
			new_candidates.append(candidate)
	candidates = new_candidates 

inout.output(' '.join(set(map(output_format, winners))))