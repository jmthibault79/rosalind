# Spectral Convolution Problem: Compute the convolution of a spectrum.
#     Input: A collection of integers Spectrum.
#     Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k, it should
#     appear exactly k times; you may return the elements in any order.

# Sample Input:
#     0 137 186 323

# Sample Output:
#     137 137 186 186 323 49

import inout

spectrum = map(int, inout.infilelines[0].strip().split(' '))

convolution = []
l = len(spectrum)
for i in range(l):
	for j in range(i+1, l):
		diff = spectrum[i]-spectrum[j]
		if diff != 0:
			convolution.append(abs(diff))	
	
inout.output(' '.join(map(str,sorted(convolution))))