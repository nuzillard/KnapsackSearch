"""
fakeACD.py defines a transform() function that is used by sdfrw.sdfTransform(filenameIn, filenameOut, transform)
in order to add CNMR_SHIFTS in "ACD CNMR Predictor" style tags to an input sdf file 
using 13C NMR chemical shift values calculated by nmrshiftdb2, from KnapsackSearch or addnmrsdb.py.
The resulting .sdf file is correctly imported in an ACD Database created from ACD CNMR Predictor (+DB).

example: python fakeACD.py nmrsdb_quercetin2D.sdf
result: creates fake_acd_nmrsdb_quercetin2D.sdf
note: nmrsdb_quercetin2D.sdf comes from the processing of quercetin2D.sdf by addnmrsdb.py
"""
import os
import sys
import sdfrw

molid = 0
# sorry for the global variable

def acd_string_format(s):
	"""
	formats strings that are more than 200 characters long to fit with the way
	ACD formats long strings in sdf values
	"""
	limit = 200
# maximum length of a line
	lens = len(s)
# get length of string sent by the caller
	if lens <= limit:
# is it a small string?
		return s
# nothing to do for strings with less then 200 characters (included)
	lensm1 = lens-1
# lensm1 is lens minus one
	ranges = [(x*limit, (x+1)*limit) for x in range(lensm1 // limit)]
# calulates ranges at pairs of character indexes for string slicing
	ranges.append( (ranges[-1][1], lens) )
# update ranges with the indexes of the last (most probably incomplete) slice
	return '\n'.join(s[r[0]:r[1]] for r in ranges)
# get the list of slices, joins with newlines and return the result to the caller

def transform(mol):
	"""
	transform()
	"""
	global molid
	nmrskey = "NMRSHIFTDB2_ASSIGNMENT"
	acdkey = "CNMR_SHIFTS"
	idkey = "ID_fakeACD"
# sdf tag in use here
	molid += 1
# prepare index for the current molecule	
	if not sdfrw.sdfHasProp(mol, nmrskey):
# is NMRSHIFTDB2_ASSIGNMENT missing in mol?
		return mol
# do not anything with mol if no NMRSHIFTDB2_ASSIGNMENT tag
	nmrs_val = sdfrw.sdfGetProp(mol, nmrskey)
# get nmrshiftdb assignment
	nmrlines = [line[:-2] for line in nmrs_val.split('\n')]
	nmrs_parts = [line.split(', ') for line in nmrlines]
	fake = acd_string_format(';'.join([str(n)+':'+d[0]+'|'+d[1] for (n, d) in enumerate(nmrs_parts)]))
# make fake CNMR_SHIFTS value from NMRSHIFTDB2_ASSIGNMENT value
	idline = str(molid)
# build molecule index in file
	ok, mol = sdfrw.sdfSetChangeProp(mol, acdkey, fake)
# create the tag-value pair in the list of tag-values for ACD-style chemical shifts
	ok, mol = sdfrw.sdfSetNoChangeProp(mol, idkey, idline)
# create the tag-value pair in the list of tag-values for sequential molecule id
	return mol
# return new molecule

def fakeACD(filenameIn):
	head, tail = os.path.split(filenameIn)
# prepare the creation of the output file name
	filenameOut = os.path.join(head, 'fake_acd_' + tail)
# prepend fake_acd_ to input file name to obtain the output file name
	sdfrw.sdfTransform(filenameIn, filenameOut, transform)
# simply call sdfTransform() from the sdfrw library, with the locally defined transform() function.
		
if __name__ == "__main__":
	argc = len(sys.argv)
# get number of command line arguments + 1
	if argc != 2:
# only one argument possible, the name of the .sdf input file
		print("Usage: %s: filename.sdf" % (argv[0],))
		sys.exit(1)
	filenameIn = sys.argv[1]
# get input .sdf file name
	fakeACD(filenameIn)
# just do it!
