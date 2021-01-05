"""
smi2ADC.py converts a .smi file to a .sdf file ready for importation into an ACD database.
Such a database may be validated for chemical shifts (Database->Tools->Check Chemical Shifts)
and exported to a second .sdf file. The latter can be changed to incorporate the calculated shifts
as if they were experimental chemical shifts.

Usage: python smi2ACD.py [input-file-name.smi [output-file-name.sdf]]

If no file name is given, the defaut input file name "small.smi" will be used.
If no output file name is given, its value is deduced from the input file name,
so that input file name aaa\\bbb.smi is transformed into output file name aaa\\fake_acd_bbb.sdf

The input file has no title line. SMILES chains and following compound names are separated by a single space character
"""

import sys
import os

from rdkit import Chem
try:
	from rdkit.Chem import rdCoordGen
	drawingFunction = rdCoordGen.AddCoords
# recent, nice (slow?) 2D coordinates generator
except:
	from rdkit.Chem import AllChem
	drawingFunction = AllChem.Compute2DCoords
# historical 2D coordinates generator

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

def one_mol(mol):
	"""
	one_mol() adds 2D coordinates to the current molecules mol, assigns
	a defaut chemical shift (99.99) to all carbons atoms and stores
	chemical shift information as an sdf tag.
	"""
	drawingFunction(mol)
# 2D coordinate generation
	catoms = [a.GetIdx()+1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
# list of indexes (start at 1) of C atoms
	ctag = 'CNMR_SHIFTS'
# SDF tag name for ACD string of 13C NMR chemical shift data
	default_shift = '99.99'
# defaut value of 13C NMR chemical shifts
	cvalue = acd_string_format(';'.join(['%d:%d|%s' % (i, c, default_shift) for (i, c) in enumerate(catoms)]))
# build ACD-formatted string of 13C NMR chemical shifts for molecule mol
	mol.SetProp(ctag, cvalue)
# store in molecule mol the ACD-formatted string of 13C NMR chemical shifts under tag CNMR_SHIFTS

def err_msg(argv):
	"""
	err_msg() prints the usage message in case of incorrect command line arguments
	"""
	print("Usage: python", argv[0], "[input-file-name.smi [output-file-name.sdf]]")
	sys.exit(1)

def run(filenameIn, filenameOut=None):
	"""
	run() processes the SMARTS-Name pairs in input file to create a .sdf file of the
	corresponding molecular structures with added 2D coordinates and deep-fake chemical
	shift values (all 99.99) formatted for ACD software.
	If no output file name is given, a new one is deduced from the input file name (vide supra).
	"""
	rootIn, extIn = os.path.splitext(filenameIn)
# get root and extension of input file
	if extIn != '.smi' and extIn != '.SMI':
# extension is not the expected one for a file with SMILES chains (.smi)
		err_msg(argv)
# bad extension, print usage message and exit
	head, tail = os.path.split(rootIn)
# separate access path from file name
	if not filenameOut:
# no output file name given
		filenameOut = os.path.join(head, "fake_acd_" + tail + '.sdf')
# deduce output file name from input file name
	else:
# an output file name was provided
		rootOut, extOut = os.path.splitext(filenameOut)
# get extension of proposed file name
		if extOut != '.sdf' and extOut != '.SDF':
# the extension of the proposed output file name is not what expected (.sdf)
			err_msg(argv)
# print usage message if the extension of the output file name is not .sdf
	supplier = Chem.SmilesMolSupplier(filenameIn, titleLine=False)
# get handle to SMI file reader
	writer = Chem.SDWriter(filenameOut)
# get handle to SDF file writer
	converted = 0
# count successfully converted structures
	for imol, mol in enumerate(supplier, start=1):
# iterate through SMILES chains
		if mol is None:
# error while reading/analyzing the current SMILES chain
			print("Skipped:", imol)
# message for invalid SMILES chains
			continue
# skip to next line in the input file
		one_mol(mol)
# supplements current structure with 2D coordinates and fake CNMR data
		writer.write(mol)
# write current molecule
		converted += 1
# increment number of molecules written to output file
	print("Converted from", filenameIn, "to", filenameOut+':', converted)
# summary of conversion with file names and number of converted molecules

if __name__ == '__main__':
# what else?
	argv = sys.argv
# get command line arguments
	argc = len(argv)
# get number of command line arguments
	defaultFilenameIn = 'small.smi'
# default input file name
	filenameOut=None
# default output file name

	ok = True
# validity of command line arguments
	if argc == 1:
# no input file name given
		filenameIn = defaultFilenameIn
# no input file name given, use default name
	elif argc == 2:
# no output file name given
		filenameIn = argv[1]
# no output file name given, get input file name only
	elif argc == 3:
# input and output file names given
		filenameIn = argv[1]
		filenameOut = argv[2]
# get input and output file names
	else:
# too many command line arguments
		ok = False
# remember that command line arguments are not as expected
	if not ok:
# wrong command line arguments
		err_msg(argv)
# print usage message and exit
	else:
# compliant command line arguments
		run(filenameIn, filenameOut)
# process .smi input file to produce .sdf output file
