# How many subpeptides does a cyclic peptide of length n have?

# For example, the cyclic peptide NQEL has 12 subpeptides: N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, and LNQ. Subpeptides may occur more than once if an amino acid occurs multiple times in the peptide (e.g., ELEL also has 12 subpeptides).

def subpep(n):
	return n * (n-1)

print(subpep(4))
print(subpep(31315))
print(subpep(46121))