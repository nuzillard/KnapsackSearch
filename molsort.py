"""
molsort.py:
	reads: familyname_2D_nmr.txt
	writes: familyname_2D_nmr_sorted.txt

familyname_2D_nmr.txt: written by predict_sdf.bat by redirection of the standard output
familyname_2D_nmr_sorted.txt: read by familyname_2D_nmr_sorted.txt

called by: python molsort.py familyname

Each line in file familyname_2D_nmr.txt contains a molecule index (start at 0),
a carbon atom index (start at 0) and the corresponding 13C NMR chemical shift, all separated by a space.
Molecules are indexed according to their position in file familyname_2D_tagged.sdf.
Carbon atoms are indexed according to their position in the MolBlock in file familyname_2D_tagged.sdf.
The molecule indexes do not follow the native order and carbon indexes for a given molecule
do not follow the native order.
molsort.py reorders carbon atoms in molecules and molecules in file familyname_2D_nmr.txt
to produce file familyname_2D_nmr_sorted.txt.
File familyname_2D_nmr_sorted.txt contains one line per compound, made of 3 ' | '-separated fields:
field 1 is the molecule index (start at 1), field 2 is a space-separated list of carbon atom indexes (start at 1),
and field 3 is the list of chemical shifts for carbon atoms of indexes ordered as given by field 2
"""

import sys

def molsort(inputfilename):
	"""
	inputfilename is the name of a file that is created by predict_sdf.bat.
	molsort retains only the lines that contain mol_id-atom_id-chemical_shift triplets.
	molsort creates a list of pairs. Each pair contains the molecule sequence index (start at 0)
	and an list of pairs. Each pair contains a carbon atom index (start at 0) and a chemical shift value as a float.
	"""
	curmolid = -1
# no molecule yet
	curmol = []
# current molecule is initially empty
	mols = []
# the list of molecules is initially empty
	with open(inputfilename) as f:
		for line in f:
# open input file and prepares for sequential line reading
			if not (line and line[0].isdigit()):
				continue
# skip to next line if line is empty or does not start with a digit
			fields = line.split()
			molid = int(fields[0])
			atomid = int(fields[1])
			chemshift = float(fields[2])
# get molecule index, atom index and chemical shift as strings and converts to numerics
			if molid != curmolid:
# first line for a new molecule
				if len(curmol) != 0:
# the (old) current molecule has to be stored in the list of molecules only if the molecule is not empty.
					appendmol(curmolid, curmol, mols)
# store the current molecule
					curmol = []
# the current molecule is cleared
				curmolid = molid
# remember the index of the current molecule
			curmol.append((atomid, chemshift))
# update the current molecule with the new atom index and its chemical shift
	appendmol(curmolid, curmol, mols)
# store the very last molecule
	return sorted(mols, key = lambda p: p[0])
# return the list of molecules, sorted according to molecule index

def appendmol(curmolid, curmol, mols):
	"""
	sort the atom index and associated chemical shifts from curmol
	in the order of indexes and appends the curmolid (molecule index)
	and molecule in mols.
	"""
	mols.append((curmolid, sorted(curmol, key = lambda p: p[0])))

def print_formatted_1(output):
	"""
	textual output of molecules in the style of the input file, without the heading lines.
	but with indexes starting at 1.
	not used here but kept for other possible applications
	"""
	for (molid, curmol) in output:
#		for (atomid, chemshift) in curmol:
#			print("%d %d %.2f" % (molid+1, atomid+1, chemshift))
		return ''.join(["%d %d %.2f\n" % (molid+1, atomid+1, chemshift) for (atomid, chemshift) in curmol])

def print_formatted_2(output):
	"""
	textual output of mecules in 3 |-separated field per line (and per molecule).
	molecule id | list of atom indexes | list of chemical shifts in the same order as in the list of atom indexes
	"""
	output_formatted = ""
	for (molid, curmol) in output:
		atomids = ' '.join([str(curatom[0]+1) for curatom in curmol])
		chemshifts = ' '.join(["%.2f" % (curatom[1],) for curatom in curmol])
#		print(" | ".join([str(molid+1), atomids, chemshifts]))
		output_formatted += (" | ".join([str(molid+1), atomids, chemshifts]) + '\n')
	return output_formatted

if __name__ == "__main__":
	family = sys.argv[1]
# get family name from command line
	inputfilename = family + "_2D_nmr.txt"
	outputfilename = family + "_2D_nmr_sorted.txt"
# build names from family name for input and output files
	output = molsort(inputfilename)
# get list of (mol_id, list of (carbon atom index, chemical shift))
	output_formatted = print_formatted_2(output)
# get formatted output
	with open(outputfilename, 'w') as fp:
		fp.write(output_formatted)
# write formatted output to file familyname_2D_nmr_sorted.txt
