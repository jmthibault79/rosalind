# Protein Translation Problem: Translate an RNA string into an amino acid string.

# Sample Input:
#     AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

#Sample Output:
#     MAMAPRTEINSTRING

import inout

stop_codons = []
codon_map = {}
for line in inout.readlines('RNA_codon_table_1.txt'):
	codon, amino_acid = line.split(' ')
	if amino_acid:
		codon_map[codon] = amino_acid.strip()
	else:
		stop_codons.append(codon)

sequence = inout.infilelines[0].strip()
outstr = ''
while sequence:
	codon, sequence = sequence[:3], sequence[3:]
	if codon in stop_codons:
		break
	else:
		outstr = outstr + codon_map[codon]
		
inout.output(outstr)
