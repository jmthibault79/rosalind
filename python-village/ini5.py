def save_even_lines (infilename, outfilename):
	import os
	infile = open(os.path.expanduser(infilename), 'r')
	outfile = open(os.path.expanduser(outfilename), 'w')
	printme = False
	for line in infile:
		if printme:
			outfile.write(line)
		printme = not printme
	outfile.close()
	infile.close()
	
save_even_lines('~/Downloads/rosalind_ini5.txt', '~/Downloads/rosalind_ini5_even.txt')
