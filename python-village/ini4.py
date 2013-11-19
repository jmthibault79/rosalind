def sum_all_odd_inclusive (a, b):
	acc = 0
	for candidate in range(a, b+1):
		if candidate % 2 == 1:
			acc = acc + candidate
	return acc
	
print sum_all_odd_inclusive(4876, 9779)
