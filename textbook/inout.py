import sys
import os

infilename = os.path.expanduser(sys.argv[1])
outfilename = os.path.expanduser(sys.argv[2])

with open(infilename, 'r') as infile:
	infilelines = infile.readlines()
	
def output(outstr):
	with open(outfilename, 'w') as outfile:
		outfile.write(outstr)
