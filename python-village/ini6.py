def string_token_counter (s):
	dict = {}
	for token in s.split():
		if token in dict:
			dict[token] = dict[token] + 1
		else:
			dict[token] = 1
	for token in dict:
		print token, dict[token]
		
import os
infile = open(os.path.expanduser('~/Downloads/rosalind_ini6.txt'), 'r')
infilelines = ' '.join(infile.readlines())
infile.close()

string_token_counter(infilelines)
