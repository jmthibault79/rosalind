import inout 	# module for handling Rosalind's file I/O
sequence = inout.infilelines[0].strip()

transcribed = sequence.replace('T', 'U')
inout.output(transcribed)
