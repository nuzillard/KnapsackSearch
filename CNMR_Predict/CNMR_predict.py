"""
Yet another way to build an ACD database with predicted chemical shifts
replacing experimental chemical shifts.
The input .sdf to CNMR_predict.py may have been created by an exportation from an ACD database
in which "experimental" 13C NMR chemical shifts were predicted by nmrshiftdb2 (by fakeACD.py
or by KnapsackSearck) and validated by ACD (ACD Database->Tools->Check Chemical Shifts).
CNMR_predict.py replaces the "experimental" chemical shifts by those calculated for validation.
and writes the result in the output .sdf file.
The output file may be imported in a new ACD Database. Chemical shift validation
indicates that "experimental" chemical shifts are identical to those calculated for validation.
CNMR_predict.py also formats chemical shift values in the style of NMReDATA and of MixONat.
"""

import sdfrw
import sys
from rdkit import Chem

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

def cnum1(s):
	"""
	x[y] -> y
	z    -> z
	"""
	p1 = s.find('[')
	if p1 < 0:
		return s
	p2 = s.find(']')
	return s[p1+1:p2]

def cnum2(s):
	"""
	x[y] -> x[y]
	z    -> z[z]
	"""
	return s if '[' in s else s + '[' + s + ']'

def get_nmredata(simple_ids, shifts):
	"""
	combines atom ids from MOL block and chemical shifts to produce an nmredata-formatted list of chemical shifts
	"""
	return '\n'.join(['%s, %s, %s \\' % x for x in zip(simple_ids, shifts, simple_ids)])

def get_onat_shifts(full_ids, shifts):
	"""
	combines explicit atom ids "ACD_id[MOL block_id]" and chemical shifts to produce a mixonat-formatted list of chemical shifts
	"""
	return '\n'.join(["%s \t%s \t" % x for x in zip(full_ids, shifts)])

def get_onat_mults(mol, simple_ids, shifts):
	"""
	creates the tag-value pairs for the identification of carbon atoms with different multiplicities
	and their associated chemical shift values
	"""
	molH = Chem.AddHs(Chem.MolFromMolBlock(sdfrw.sdfGetMolBlock(mol), removeHs=False))
# transform sdfrw Mol block into an RDKit molecule with added H atoms.
	cnmr_dict = dict(zip(map(int, simple_ids), shifts))
# create a dictionary with atom index as key and chemical shift as value
	mults = [[], [], [], []]
# prepare a list indexed by a number of H atoms bore by C atoms, with values as
# 2-element lists made of a C atome index and the associated chemical shift value
	onat_mults = []
# prepare a list indexed by a number of H atoms bore by C atoms, with text 
# to be printed about the index and chemical shift value pairs associated to the numbers of H
	for c, cshift in cnmr_dict.items():
# loop over C atoms indexes and associated chemical shift values
		a_idx = c - 1
# convert MOL block atom index (start at 1) to RDKit index (start at 0)
		atom = molH.GetAtomWithIdx(a_idx)
# retrieve atom with RDKit index a_idx
		neighbors = atom.GetNeighbors()
# get list of neighbors of current C atom
		isH = [1 if x.GetAtomicNum() == 1 else 0 for x in neighbors]
# affects 1 to H atoms among the neighbors of current C atom and 0 to other atoms
		mult = sum(isH)
# calculates the number of H atoms bound to the current C atom (multiplicity)
		if mult < 4:
# current structure is not the one of methane, the only compound with a C and four H atoms
			mults[mult].append([c, cshift])
# update mults list
	tags = ['Quaternaries', 'Tertiaries', 'Secondaries', 'Primaries']
# array of tags associated to the multiplicity values from 0 to 3
	for mult, tag in enumerate(tags):
# loop through number of H and associated textual tags
		items = mults[mult]
# for the current multiplicity value, get pairs of atom indexes and chemical shifts
		value = '\n'.join(["%d\t%s" % tuple(x) for x in items])
# format string about the current multiplicity
		onat_mult = (tag, value)
# tag-value pair for the current multiplicity
		onat_mults.append(onat_mult)
# update onat_mults
	return onat_mults
# return data necessary to create the four sdf tags related to C atom multiplicity, 
# formatted for MixONat

def get_exported_shifts(ids, shifts):
	"""
	combines ACD implicit atom ids (either x[y] or z) and chemical shifts to produce an ACD-formatted list of chemical shifts
	"""
	new_exported_pieces = [ str(i)+':'+p[0]+'|'+p[1] for (i, p) in enumerate(zip(ids, shifts)) ]
# reconstruct exported data from atom numbers and calculated chemical shifts
	return acd_string_format(';'.join(new_exported_pieces))
# build and return a text from the reconstructed data as if they had been exported

def reshape_from_acd(s):
	"""
	s is what ACD makes during database exportation from multi-line data.
	newlines were replaced by semicolon. reshape_from_acd() reverses this transformation
	"""
	one_line = s.replace('\n', '')
# the 200 character long lines are joined together
	all_lines = one_line.replace('\\;', '\\\n')
# replace semicolon by new-line
	return all_lines

def transform(mol):
	"""
	transform() transforms a molecule mol so that the 13C NMR data under the CNMR_CALC_SHIFTS shifts
	replace those in CNMR_SHIFTS. All tags created by a previous exportation from a .NMRUDB
	file to a .sdf file are removed.	
	"""
	calc = sdfrw.sdfGetProp(mol, "CNMR_CALC_SHIFTS")
# get values of ACD-calculated 13C NMR chemical shifts encoded by ACD in a text
	exported = sdfrw.sdfGetProp(mol, "CNMR_SHIFTS")
# get values of ACD-exported 13C NMR chemical shifts encoded by ACD in a text
	if calc and exported:
# the molecule has ACD-calculated and ACD-exported 13C NMR chemical shift data

		calc = calc.replace('\n', '')
# make a single text line from text of ACD-calculated data, possibly cut over several lines
		calc_pieces = calc.split(';')
# make a list of ACD-calculated data
		shifts_raw = [ x.split(':')[1].split(',')[0][8:] for x in calc_pieces ]
# extract the list of chemical shifts from ACD-calculated values
		shifts = [ '%.2f' % (float(shift),) for shift in shifts_raw ]
# list of formatted chemical shifts

		exported = exported.replace('\n', '')
# make a single text line from text of ACD-exported data, possibly cut over several lines
		exported_pieces = exported.split(';')
# make a list of ACD-exported data
		ids = [ x.split(':')[1].split('|')[0] for x in exported_pieces ]
# extract the list of atom indexes from ACD-exported values
# ids may be like x[y] or like z. x(y] means that the atom indexed in the MOL block by y (start at 1)
# has been renumbered by ACD as x (start at 1). z means that either no renumbering occured
# or that the initial index in the MOL block is identical to index of the same atom
# after renumbering (z renumbered as z).

		simple_ids = list(map(cnum1, ids))
# simple C atom identifiers, identical to those of MOL blocks (start at 1)
		full_ids = list(map(cnum2, ids))
# binomial (x[y] or z[z]) C atom identifiers, compatible with atom renumbering by ACD.

		new_nmredata = get_nmredata(simple_ids, shifts)
# get NMReDATA-formatted text relative to newly calculated chemical shifts
		onat_shifts = get_onat_shifts(full_ids, shifts)
# get MixONat-formatted text relative to newly calculated chemical shifts
		onat_mults = get_onat_mults(mol, simple_ids, shifts)
# get MixONat-formatted text relative to newly calculated chemical shifts separated according to C atom multiplicity
		new_exported_shifts = get_exported_shifts(ids, shifts)
# get ACD-formatted text relative to newly calculated chemical shifts, to be possibly imported by ACD/Labs DB

		old_nmrshiftdb = sdfrw.sdfGetProp(mol, "NMRSHIFTDB2_ASSIGNMENT")
# get value associated to NMRSHIFTDB2_ASSIGNMENT, if any.
		if old_nmrshiftdb:
# there were an NMRSHIFTDB2_ASSIGNMENT but it was flattened by the exportation from an ACD DB
			done, mol = sdfrw.sdfSetChangeProp(mol, "NMRSHIFTDB2_ASSIGNMENT", reshape_from_acd(old_nmrshiftdb))
# remove ACD formatting of NMRSHIFTDB2_ASSIGNMENT, if any, that occurs during database exportation

		done, mol = sdfrw.sdfSetChangeProp(mol, "CNMR_SHIFTS", new_exported_shifts)
# change molecule according to the new NMR data, as if they were exported by ACD database

		done, mol = sdfrw.sdfSetChangeProp(mol, "Predicted 13C shifts", onat_shifts)
# change molecule with new chemical shifts formated for mixonat

		done, mol = sdfrw.sdfSetChangeProp(mol, "NMREDATA_ASSIGNMENT", new_nmredata)
# change molecule with new chemical shifts formatted in the nmredata style
		origin_lines = ['Source=Calculation',  'Method=Database->Tools->Check Chemical Shifts', 'Software=ACD/Labs C+H NMR Predictors and DB', 'Version=2020.1.0']
		origin = '\n'.join(line + ' \\' for line in origin_lines)
		done, mol = sdfrw.sdfSetChangeProp(mol, "NMREDATA_ORIGIN", origin)
# change molecule with new data about the origin of chemical shifts

		for tag, value in onat_mults:
# loop over carbone types, according to multiplicity
			done, mol = sdfrw.sdfSetChangeProp(mol, tag, value)
# change molecule for atom ids and chemical shifts for current carbon type

		done, mol = sdfrw.sdfClearProp(mol, "HNMR_SHIFTS")
		done, mol = sdfrw.sdfClearProp(mol, "CNMR_CALC_SHIFTS")
		done, mol = sdfrw.sdfClearProp(mol, "HNMR_CALC_SHIFTS")
		done, mol = sdfrw.sdfClearProp(mol, "NMR_CONSTANTS")
		done, mol = sdfrw.sdfClearProp(mol, "CHEMICAL_SHIFTS.13C")
		done, mol = sdfrw.sdfClearProp(mol, "CHEMICAL_SHIFTS.1H")
# remove sdf tags added by acd DB during database exportation
	return mol
# return the transformed molecule

def run(inpufilename, outputfilename):
	"""
	transform all molecules from input file and writes the result to output file
	"""
	sdfrw.sdfTransform(inpufilename, outputfilename, transform)
# apply molecule transformation transform() on input .sdf file to produce the .sdf output file


if __name__ == "__main__":
# what else?
	argc = len(sys.argv)
# get number of command line arguments
	if argc not in [1, 3]:
# only allowed (not the case here): no arguments to use defaults as input and output file names or 2 file name arguments
		print("Usage:", sys.argv[0]+':', "0 or 2 arguments required")
# print usage message
	else:
# the number of command line arguments is as expected
		if argc == 1:
# use default file names
			inpufilename = "fake_acd_first_v1_exported.sdf"
			outputfilename = "pnmrnp3_first.sdf"
# define default file names
		else:
# explicit input and output file names
			inpufilename = sys.argv[1]
			outputfilename = sys.argv[2]
# get input and output file names from command line
		run(inpufilename, outputfilename)
# transfer ACD-calculated chemical shifts from input file as experimental chemical shifts in output file

