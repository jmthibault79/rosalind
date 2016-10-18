import inout 	# module for handling Rosalind's file I/O
sequence = inout.infilelines[0].strip()

d = {}
for char in sequence:
	if char in d:
		d[char] += 1
	else:
		d[char] = 1

counts = (d['A'], d['C'], d['G'], d['T'])
inout.output(' '.join(map(str, counts)))
