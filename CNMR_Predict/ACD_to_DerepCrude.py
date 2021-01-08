"""
This script transforms a .sdf file (aaa.sdf) with 13C NMR chemical shift values formatted
for being read by the ACD/Labs software (i. e. with a CNMR_SHIFTS tag) into
a .sdf file (aaa_DC.sdf) for DerepCrude, in which
the 13C NMR chemical shifts are given by the MY_CNMR tag.
"""
from rdkit import Chem
import sys
import os.path

def cnum(s):
# x[y] -> y
# z    -> z
	p1 = s.find('[')
	if p1 < 0:
		return s
	p2 = s.find(']')
	return s[p1+1:p2]


def proc_one_mol(mol):

	cnmr_string = mol.GetProp("CNMR_SHIFTS")
	# get ACD 13C chemical shift string
	
	cnmr_lines = cnmr_string.split('\n')
	# get list of lines in cnmr_string
	
	cnmr_data = ''.join(cnmr_lines)
	# build 13C NMR information as a single string
	
	cnmr_pieces = cnmr_data.split(';')
	# list of 13C NMR data, one element per atom
	
	cnmr_pieces = map(lambda x: x.split(':')[1], cnmr_pieces)
	# get rid of useless order number before ':' in each 13C NMR element
	
	cnmr_pairs = [x.split('|') for x in cnmr_pieces]
	# separate each 13C NMR element between atom number info and chemical shift value
	
	cnmr_lists = list(map(list, zip(*cnmr_pairs)))
	# get two lists, one for atom numbers, one for chemical shift values
	
	cnmr_lists[0] = map(cnum, cnmr_lists[0])
	# replace the list of atom numbers by true atom numbers
	
	cnmr_shifts_string = map(list, zip(*cnmr_lists))
	# build a list of two-element lists for all C atoms, with atom number and chemical shift value
	
	cnmr_shifts = map (lambda x: [int(x[0]), float(x[1])], cnmr_shifts_string)
	# 13C NMR data with atom number as int (start at 1) and chemical shift as float, instead of strings
	
	cnmr_shifts = sorted(cnmr_shifts, key=lambda x: x[0], reverse=False)
	
	my_cnmr = '\n'.join(["%3d %6.2f" % tuple(x) for x in cnmr_shifts])
	# formatted 13C NMR data with one line per atom number (start at 1) and chemical shift value
	
	symbols = ''.join([mol.GetAtomWithIdx(s[0]-1).GetSymbol() for s in cnmr_shifts])
	# string of atoms symbols from 1ACD 13C NMR data
	
	validity = [x == 'C' for x in symbols]
	# list of Booleans, True if atom is ACD 13C NMR data is C, False otherwise
	
	valid = False not in validity
	# ACD data valid only if 13C NMR data is about carbon atoms
	
	return valid, symbols, my_cnmr

def run(in_filename, out_filename):

	reader = Chem.SDMolSupplier(in_filename, removeHs = False)
	writer = Chem.SDWriter(out_filename)
	
	for mol in reader:
		valid, symbols, my_cnmr = proc_one_mol(mol)
				
		if valid:
			mol.SetProp("MY_CNMR", my_cnmr)
		else:
			name = mol.GetProp("Name")
			print("Invalid:", name+":", symbols)
		writer.write(mol)
	writer.close()

if __name__ == "__main__":
	argv = sys.argv
	in_filename = argv[1] if len(argv) == 2 else ''
	if not in_filename:
		print("Usage:", argv[0]+':', "input-file-name.sdf")
		sys.exit(1)
	head, tail = os.path.splitext(in_filename)
	out_filename = head + "_DC" + tail

    run(in_filename, out_filename)
