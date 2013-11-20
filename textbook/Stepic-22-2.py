
# How many subpeptides does a linear peptide of given length n have? (Include the 0-peptide and the entire peptide.)

# Note that the spectrum of a linear peptide does not contain as many masses as the spectrum of a cyclic peptide. For instance, the theoretical spectrum of the cyclic peptide NQEL contains 14 masses (corresponding to "", N, Q, E, L, LN, NQ, QE, EL, ELN, LNQ, NQE, QEL, and NQEL); however, the theoretical spectrum of the linear peptide NQEL, shown below, does not contain masses corresponding to LN, LNQ, or ELN, since these subpeptides "wrap around" the end of the linear peptide.

def subpep(n):
	return n * (n+1) / 2 + 1

print(subpep(4))
print(subpep(44415))