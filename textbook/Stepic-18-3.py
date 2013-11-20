# Protein Translation Problem: Translate an RNA string into an amino acid string.

# Sample Input:
#     AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

#Sample Output:
#     MAMAPRTEINSTRING

import inout
import codon

sequence = inout.infilelines[0].strip()
		
inout.output(codon.transcribe(sequence))
