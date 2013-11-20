# Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
#     Input: An amino acid string Peptide.
#     Output: Cyclospectrum(Peptide).

# Sample Input:
#     LEQN

# Sample Output:
#     0 113 114 128 129 227 242 242 257 355 356 370 371 484

import inout

peptide = inout.infilelines[0].strip()

mass_table = {}
for line in inout.readlines('integer_mass_table.txt'):
	amino_acid, mass = line.strip().split(' ')
	mass_table[amino_acid] = int(mass)

def total_mass(peptide):
	total = 0
	for amino_acid in peptide:
		total = total + mass_table[amino_acid]
	return total
	
spectrum = [0, total_mass(peptide)]

peptide_2 = peptide + peptide	# for easy cyclic access
for k in range(1, len(peptide)):
	for n in range(len(peptide)):
		subpep = peptide_2[n:n+k]
		spectrum.append(total_mass(subpep))
		
inout.output(' '.join(map(str, sorted(spectrum))))