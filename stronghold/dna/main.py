def readlines(filename):
    import os
    infile = open(os.path.expanduser(filename), 'r')
    infilelines = infile.readlines()
    infile.close()
    return infilelines

infile = 'rosalind_dna.txt'
line = (readlines(infile)[0]).strip()

d = {}
for char in line:
	if char in d:
		d[char] += 1
	else:
		d[char] = 1

print d['A'], d['C'], d['G'], d['T'],
