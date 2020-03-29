"""
Package split_mf.py contains two functions.

split_mf(mf) transforms a molecular formula such as XxnYZo+p
into the data structure ([('Xx', n), ('Y', 1), ('Z', o)], p)

match_MS_MplusH(mf1, mf2) return True if mf1 is like mf2 but with one more H atom,
as if mf1 were the [M+H] species derived from M, as frequently observed in Mass Spectrometry.
Molecular charges of species of molecular formulas mf1 et mf2 are not considered.
"""
def split_mf(mf):
	"""
	split_mf() transforms a molecular formula mf such as XxnYZo+p
	into the data structure ([('Xx', n), ('Y', 1), ('Z', o)], p)
	"""
	d = {}
# d is a dictionary with element names as keys
	charge = 0
# electric charge is zero by default
	signpos = mf.find('+')
# look for a + sign for a positively charged molecule
	if signpos > -1:
# the molecule is positively charged
		strcharge = mf[signpos+1:]
# string for the positive electric charge
		charge = int(strcharge) if strcharge else 1
# charge string might be empty, meaning a +1 charge
		mf = mf[0:signpos]
# keep non charged part of mf as mf
	signpos = mf.find('-')
# look for a - sign for a negatively charged molecule
	if signpos > -1:
# the molecule is negatively charged
		strcharge = mf[signpos+1:]
# string for the negative electric charge
		charge = -int(strcharge) if strcharge else -1
# charge string might be empty, meaning a -1 charge
		mf = mf[0:signpos]
# keep non charged part of mf as mf
	numchars = len(mf)
# number of characters in mf for splitting in mf into parts, one part per element
	uppers = [i for i in range(numchars) if mf[i].isupper()] + [numchars]
# an element name always starts with an uppercase letter, get their location,
# append length of mf to indicate stop position for the last element in mf
	numelts = len(uppers)-1
# number of elements in md, -1 due to the stop-value in uppers
	pieces = [mf[uppers[i]:uppers[i+1]] for i in range(numelts)]
# split mf into pieces related to each element
	for piece in pieces:
# for each piece of mf related to an element
		if len(piece) == 1:
# a one-letter element symbol, no occurence number meaning 1
			d[piece] = 1
# one occurence for the current element
		else:
# two or more characteres for the current element
			ipos = 2 if piece[1].islower() else 1
# if the second character in piece is lowercase, it means that the symbol of the element
# is a two-character symbol (like Br). pos is the position of the number of occurences
			key = piece[0:ipos]
# key is the element name
			strvalue = piece[ipos:]
# strvalue, is the number of occurences, as tring
			value = int(strvalue) if strvalue else 1
# if the string for the number of occurence is empty, the number of occurences is 1
			d[key] = value
# associate element symbol with number of occurences
	a = []
# array of (symbol, number of occurences) pairs, start with C, then H,
# then other elements in alphabetic order
	if 'C' in d:
# start with C, if any
		a.append( ('C', d['C']) )
		del d['C']
	if 'H' in d:
# continue with H, if any
		a.append( ('H', d['H']) )
		del d['H']
	if bool(d):
# continue with other elements, if any
		a += [(k, d[k]) for k in sorted(d.keys())]
# append paits to array a with first pair elements sorted in alphabetic order
	return (a, charge)
# return pair : (array of (symbol, number of occurences) pairs, electric charge)
	
def match_MS_MplusH(mf1, mf2):
	"""
	match_MS_MplusH() return True if mf1 is like mf2 but with one more H atom,
	as if mf1 were the [M+H] species derived from M, as frequently observed in Mass Spectrometry.
	Molecular charges of species of molecular formulas mf1 et mf2 are not considered.
	"""
	mf1mf = split_mf(mf1)[0]
# get array of (symbol, number of occurences) pairs for mf1
	mf2mf = split_mf(mf2)[0]
# get array of (symbol, number of occurences) pairs for mf2
	if len(mf1mf) != len(mf2mf):
# different number of elements, cannot be [M+H] and [M]
		return False
	if len(mf1mf) < 2:
# at least two elements present in mf1 and mf2
		return False
	elts1 = [pair[0] for pair in mf1mf]
# elements in mf1
	elts2 = [pair[0] for pair in mf2mf]
# elements in mf2
	if elts1 != elts2:
# the lists of elements much be identical
		return False
	if elts1[1] != 'H':
# the second element in mf1 and mf2 must be H, meaning that mf1 and mf2 contain C!
		return False
	nums1 = [pair[1] for pair in mf1mf if pair[0] != "H"]
# number of occurences in mf1, for H excepted
	nums2 = [pair[1] for pair in mf2mf if pair[0] != "H"]
# number of occurences in mf2, for H excepted
	idents = [p[0] == p[1] for p in zip(nums1, nums2)]
# list of identity as Booleans for number of occurences of elements, H exceepted
	if False in idents:
# for all non-H elements there must be no difference in the number of occurences
		return False
	numH1 = mf1mf[1][1]
# number of H atoms in mf1
	numH2 = mf2mf[1][1]
# number of H atoms in mf2
	if numH1 != numH2 + 1:
# mf1 id [M+H] if mf2 is M
		return False
	return True
 # return True if all tests passed successfully
 
if __name__ == "__main__":
	formulas = ["C2H6O", "CH5N+", "CH6N2+2", "CCl4"]
# testing for split_mf()
	for mf in formulas:
		print("%s: %s" % (mf, str(split_mf(mf))))
	print()
	formulas = [("CH6N", "CH5N"), ("CH7N", "CH5N"), ("CH6NO", "CH5NO2"), ("CH5NBr", "CH5N")]
# testing for match_MS_MplusH()
	for mf1, mf2 in formulas:
		print("%s and %s: %s" % (mf1, mf2, str(match_MS_MplusH(mf1, mf2))))
