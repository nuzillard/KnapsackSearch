"""
The NMRSDBSort.py module provides the molsort(inputfilename, outputfilename) function.

The input file has been created by predictSdf, that runs the java program Demo.class,
that embeds predictorc.jar, the main part of the 13C NMR chemical shift predictor of the nmrshiftdb2 project.

The output of Demo.class is made of lines with triplets that contain a molecule index, a carbon atom index,
and a chemical shift value. Indexes might not appear in the increasing order.

The goal of this module is to reorder data written by Demo.class and to reformat them in the output file
for future needs.
"""

def molsort(inputfilename, outputfilename):
	"""
	inputfilename is the name of a file that is created by predict_sdf.bat.
	molsort retains only the lines that contain mol_id-atom_id-chemical_shift triplets.
	molsort creates a list of pairs. Each pair contains the molecule sequence index (start at 0)
	and an list of pairs. Each pair contains a carbon atom index (start at 0) and a chemical shift value as a float.
	molsorts writes a line per molecule with rordered atom indexes and in the order of molecule indexes.
	Each line contains threee parts separated by ' | '. Part one is the molecule index (start at 1)
	Part two is the sorted list of carbon atoms indexes (start at 1) with values separated by ' '.
	Part three is the list of chemical shifts in the order defined by atom indexes in part two. Values are separated by ' '.
	"""
	curmolid = -1
# no molecule yet
	curmol = []
# current molecule is initially empty
	mols = []
# the list of molecules is initially empty
	with open(inputfilename) as fpIn:
		for line in fpIn:
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
	sortedmols = sorted(mols, key = lambda p: p[0])
# list of molecules, sorted according to molecule index

	with open(outputfilename, 'w') as fpOut:
# open output file for writing
		expected = 0
# expected index for the current molecule, start at 0
		for (molid, curmol) in sortedmols:
# iterate through NMR data of molecules, molid is part one of the line in output
			while expected != molid :
# the expected molecule index is less than the one given by nmrshiftdb2
# because there was a molecule without any carbon atom
				fpOut.write('\n')
# write a blank line instead of a line with atom and chemical shift data
				expected += 1
# increment expected so that it can reach molid, given by nmrshiftdb2
			expected += 1
# get expected ready for the next molecule
			atomids = ' '.join([str(curatom[0]+1) for curatom in curmol])
# collect atom indexes, atomids is part two of the line in output
			chemshifts = ' '.join(["%.2f" % (curatom[1],) for curatom in curmol])
# collect chemical shift values, chemshifts is part three of the line in output
			output_formatted = " | ".join([str(molid+1), atomids, chemshifts]) + '\n'
# build line for output
			fpOut.write(output_formatted)
# write file to file named outputfilename


def appendmol(curmolid, curmol, mols):
	"""
	sort the atom index and associated chemical shifts from curmol
	in the order of indexes and appends the curmolid (molecule index)
	and molecule in mols.
	"""
	mols.append((curmolid, sorted(curmol, key = lambda p: p[0])))

"""
if __name__ == "__main__":
# demo for NMRSDBSort package
	inputfilename = "foobar_nmr.txt"
	outputfilename = "foobar_nmr_sorted.txt"
	molsort(inputfilename, outputfilename)
"""

