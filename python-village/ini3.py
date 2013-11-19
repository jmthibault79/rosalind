def two_inclusive_slices (str, begin1, end1, begin2, end2):
	return str[begin1:end1+1] + ' ' + str[begin2:end2+1]
	
teststr = 'HumptyDumptysatonawallHumptyDumptyhadagreatfallAlltheKingshorsesandalltheKingsmenCouldntputHumptyDumptyinhisplaceagain.'
print two_inclusive_slices(teststr, 22, 27, 97, 102)

teststr2 = 'CfolD0CTeratoscincusIJyQusJq6QDn1MZh8hmvuPwrO45Dsbq7hVxZ1clgsjIZTavUmvX1dfLzDFfYOSy5xZ9n6KWcNktHvxPjJWXMNhV1cpMpLVf5hr8WOAKm4Fx1Dw7KSCGalbatrusOc5nOE1g7pCzoNWruDqAk6M.'
print two_inclusive_slices(teststr2, 7, 19, 135, 142)
