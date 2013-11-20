import sys
import os

infilename = os.path.expanduser(sys.argv[1])
outfilename = os.path.expanduser(sys.argv[2])

def readlines(filename):
	with open(os.path.expanduser(filename), 'r') as infile:
		return infile.readlines()

def writelines(filename, outstr):
	with open(filename, 'w') as outfile:
		outfile.write(outstr)

def output(outstr):
	writelines(outfilename, outstr)
	
infilelines = readlines(infilename)
