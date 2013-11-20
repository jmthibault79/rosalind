import inout

stop_codons = []
codon_map = {}
for line in inout.readlines('RNA_codon_table_1.txt'):
	tokens = line.strip().split(' ')
	codon = tokens[0]
	if len(tokens) == 1:
		stop_codons.append(codon)
	else:
		amino_acid = tokens[1]
		codon_map[codon] = amino_acid

def transcribe(sequence):
	outstr = ''
	while sequence:
		codon, sequence = sequence[:3], sequence[3:]
		if codon in stop_codons:
			break
		else:
			outstr = outstr + codon_map[codon]
	return outstr
