import inout

mass_table = {}
for line in inout.readlines('integer_mass_table.txt'):
	amino_acid, mass = line.strip().split(' ')
	mass_table[amino_acid] = int(mass)

def total_mass(peptide):
	total = 0
	for amino_acid in peptide:
		total = total + mass_table[amino_acid]
	return total

def cyclic_spectrum(peptide):
	out_spectrum = [0, total_mass(peptide)]

	peptide_2 = peptide + peptide	# for easy cyclic access
	for k in range(1, len(peptide)):
		for n in range(len(peptide)):
			subpep = peptide_2[n:n+k]
			out_spectrum.append(total_mass(subpep))
	return sorted(out_spectrum)

def linear_spectrum(peptide):
	out_spectrum = [0]

	for i in range(0, len(peptide)):
		for j in range(i, len(peptide)):
			subpep = peptide[i:j+1]
			out_spectrum.append(total_mass(subpep))
	return sorted(out_spectrum)

amino_acids = mass_table.keys()