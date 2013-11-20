# Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
#     Input: A DNA string Text and an amino acid string Peptide.
#     Output: All substrings of Text encoding Peptide (if any such substrings exist).

# Sample Input:
#     ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
#     MA

# Sample Output:
#     ATGGCC
#     GGCCAT
#     ATGGCC

# Note: The solution may contain repeated strings if the same string occurs more than once as a substring of Text and encodes Peptide.

import inout
import codon

sequence = inout.infilelines[0].strip()
pattern = inout.infilelines[1].strip()

complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(kmer):
	r = ''
	for base in kmer:
		r = complement[base] + r
	return r

def DNA_to_RNA(kmer):
	return kmer.replace('T', 'U') 

k = len(pattern) * 3

outstr = ''
for idx in range(len(sequence) - k + 1):
	kmer = sequence[idx:idx+k]
	if codon.transcribe(DNA_to_RNA(kmer)) == pattern or codon.transcribe(DNA_to_RNA(reverse_complement(kmer))) == pattern:
		outstr = outstr + kmer + '\n'

inout.output(outstr)
