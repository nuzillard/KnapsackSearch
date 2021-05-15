"""
python rdcharge.py inputfile.sdf --> inputfile.sdf
python rdcharge.py inputfile.sdf outputfile.sdf --> outputfile.sdf

transforms the input file so that non-zero valence data associated by RDKit to electrically charged atoms
is reset to zero. The input file has to be created by RDKit with the V2000 format.
Does not use RDKit for compound reading and writing.

This module has the same function as novalence.py but has not to interprete the content of the file.
"""


import sdfrw
import tempfile
import sys
import os


emptyMol = sdfrw.sdfEmptyMol()

def correct(old_txt):
	"""
	correct() sets to 0 the atom valence indication from the mol block in old_text and
	returns the corrected text
	"""
	new_lines = []
# initialize list of corrected lines in old_txt
	old_lines = old_txt.split('\n')
# get lines from old_txt
	creator_line = old_lines[1].strip()
# line 2, indexed 1 in old_txt, starts with text producer
	creator = creator_line.split()[0]
# get creatot of old_txt
	if creator != 'RDKit':
# creator of old_txt is not RDKit
		return ''
# denies text correction
	count_line = old_lines[3].rstrip()
# line 4, indexed 3 in old_txt, gives the number of atom lines and the file type (V2000 or V3000)
	if 'V2000' not in count_line:
# file type is not V2000
		return ''
# denies text correction 
	numAtom = int(count_line[0:3])
# number of atoms is the first field of the count line
	for i in range(0, 4):
# loop over lines of old_txt preamble
		new_lines.append(old_lines[i])
# no change for preamble lines
	for i in range(4, 4+numAtom):
# loop over atom lines
		line = old_lines[i]
# get current atom line
		pieces = line.split()
# put current atom line into pieces
		valence = pieces[9]
# valence in the 10th field, indexed 9
		if valence != '0':
# a non-zero valence was set, atom is charged
			line = line[0:48]+'  0'+line[51:]
# line surgery to remove valence indication
		new_lines.append(line)
# accounts for the modified atom line
	for i in range(4+numAtom, len(old_lines)):
# loop over the lines that follow the atom lines
		new_lines.append(old_lines[i])
# no change for the lines that follow the atom lines
	new_txt = '\n'.join(new_lines)
# rebuild the text with modified atom lines
	return new_txt
# returns the sdf text without valence idications in the atom lines


def transform(mol):
	"""
	transforms the sdfrw molecule so that the valence data on electrically charged atoms is reset to 0
	"""
	old_txt = sdfrw.sdfGetMolBlock(mol)
# get mol block
	new_txt = correct(old_txt)
# apply correction to old mol block to obtain the new mol block
	if len(new_txt) == 0:
# the mol block is not valid
		return emptyMol
# return an emplty mol if the old mol block is not valid, that will not be written on output file
	mol = sdfrw.sdfSetMolBlock(mol, new_txt)
# change the old mol block to the new one
	return mol
# return chaned compound

def rdcharge(pathIn, pathOut=""):
	"""
	rdcharge() reads an SDF file from pathIn and writes valence-corrected molécules in place if pathOut is empty
	or to pathOut if the latter is not empty.
	"""
	inplace = len(pathOut) == 0
# determine if the transformation looks as if operated in-place
	if inplace:
# in-place-like file transformation
		pathOut = tempfile.mktemp('.sdf')
# get a temporary name for the content of the output SDF file
		
	sdfrw.sdfTransform(pathIn, pathOut, transform)
# apply compound transformation, so that atom valence data is reset to 0	
	if inplace:
# pseudo in-place processing
		os.replace(pathOut, pathIn)
# for in-place processing, copy (possibly overwrite) the temporary output file with the name of the input file

if __name__ == '__main__':
# what else?
	argv = sys.argv
	argc = len(argv)
# get command line arguments and the number of arguments
	if (argc != 3) and (argc != 2):
# unexpected number of arguments
		print(argv[0]+': bad number of command line arguments', file=sys.stderr)
		print('Usage:', argv[0], 'path_to_input_file path_to_output_file', file=sys.stderr)
# error message on standard error
		sys.exit(1)
# exit
	pathIn = argv[1]
	pathOut = argv[2] if argc == 3 else ""
	rdcharge(pathIn, pathOut)
# does the job
