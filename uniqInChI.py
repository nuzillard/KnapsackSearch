"""
Dedicated to sdf files that refer to the LOTUS database, with lotus_id tag defined.

python uniqInChI.py inputfile.sdf --> inputfile.sdf
python uniqInChI.py inputfile.sdf outputfile.sdf --> outputfile.sdf

keeps only from input file the non-duplicated (or unique) compounds.
Prints InChI string and lotus_id about duplicated compounds.
Does not use RDKit for structure reading and writing.
"""
import sys
import os
import tempfile
import sdfrw
from rdkit import Chem

def uniqInChI(pathIn, pathOut=""):
	"""
	uniqInChI() reads an SDF file from pathIn and writes unique molécules in place if pathOut is empty
	or to pathOut if the latter is not empty. Supplements the SDF tags with InChI_rdkit and InChIKey_RDKit.
	The mol blocks are not modified. rdkit is there to calculate the InChI strings.
	"""
	inplace = len(pathOut) == 0
# determine if the transformation looks as if operated in-place
	if inplace:
# in-place-like file transformation
		pathOut = tempfile.mktemp('.sdf')
# get a temporary name for the content of the output SDF file

	molEmpty = sdfrw.sdfEmptyMol()
# prepare an empty mol
	molDict = {}
# molDict is a dict whose keys an InChI string and values are arrays of lotus_id.
# the input file has to originate from LOTUS, http://lotus.naturalproducts.net
	count = 0
# initial count of processed structures
	uniqcount = 0
#initial count of unique molecules
	with open(pathIn) as fhIn, open(pathOut, "w") as fhOut:
# open for reading input file and for writing unique structures
		for molIn in sdfrw.sdfReader(fhIn):
# read the current structure
			mb = sdfrw.sdfGetMolBlock(molIn)
# get the MOL block of the current structure
			rdmol = Chem.MolFromMolBlock(mb)
# get the RDKit molecule from molblock of the current structure
			inchi = Chem.MolToInchi(rdmol)
# calculate the InChI of the current structure
			inchikey = Chem.InchiToInchiKey(inchi)
# calculate the InChIKey of the current structure
			lts = sdfrw.sdfGetProp(molIn, "lotus_id")
# get the lotus_id of the current structure
			count += 1
# one more molecule from input
			if inchi in molDict:
# the current molecule is a duplicate (or even worse)
				molDict[inchi].append(lts)
# append the new lotus_id to the list of lotus_ids associated to the current InChI string
				molOut = molEmpty
# prepare to write an empty molecule (that is, skipping the writing of the current structure in output file)
			else:
# the current molecule is not a duplicate
				molDict[inchi] = [lts]
# associate the current InChI string with its list of lotus_ids, a singleton for the moment
				molOut = molIn
# molecule on output is the current molecule
				ok, molOut = sdfrw.sdfSetChangeProp(molOut, "InChI_rdkit", inchi)
# remember the InChI of the current compound in the output file
				ok, molOut = sdfrw.sdfSetChangeProp(molOut, "InChIKey_rdkit", inchikey)
# remember the InChIKey of the current compound in the output file
				uniqcount += 1
# one more unique compound
			sdfrw.sdfWrite(fhOut, molOut)
# write out compound


	if inplace:
# pseudo in-place processing
		os.replace(pathOut, pathIn)
# for in-place processing, copy (possibly overwrite) the temporary output file with the name of the input file

	"""
	for inchi, idList in molDict.items():
# loop through the inchi-lotus_ids dict
		if len(idList) > 1:
# the current inchi has been found more than once
			print(inchi)
# print InChI of the currently duplicated structure
			for id in idList:
# loop over lotus_id of duplicated compounds
				print('\t' + id)
# print duplicated lotus_id, shifted by a tab
	"""
	print("Read: %d -- Written: %d -- Discarded: %d" % (count, uniqcount, count-uniqcount))
# print statistics

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
# exit on error
	pathIn = argv[1]
	pathOut = argv[2] if argc == 3 else ""
	uniqInChI(pathIn, pathOut)
# does the job
