import inout 	# module for handling Rosalind's file I/O
sequence = inout.infilelines[0].strip()
reversed_seq = sequence[::-1]

complements = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

rc = [complements[x] for x in reversed_seq]
inout.output(''.join(rc))
