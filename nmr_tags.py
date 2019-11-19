"""
nmr_tags.py:
	reads: familyname_2D_tagged.sdf
	reads: familyname_2D_nmr_sorted.txt
	writes: familyname.nmredata.sdf

familyname_2D_tagged.sdf: written by basic_tags.py
familyname_2D_nmr_sorted.txt: wwritten by molsort.py
familyname.nmredata.sdf: read by user, possibly using EdiSdf

called by: python nmr_tags.py familyname

The file familyname_2D_tagged.sdf contains structural data about compounds in KNApSAcK
related to family familyname.
The file familyname_2D_nmr_sorted.txt contains 13C NMR chemical shift data for all molecules
in familyname_2D_tagged.sdf
nmr_tags.py appends NMR data to familyname_2D_tagged.sdf according to the NMReDATA format
and writes the result to file familyname.nmredata.sdf
"""

from rdkit import Chem
import sys

family = sys.argv[1]
# get family name from command line
inputfilename = family + "_2D_nmr_sorted.txt"
sdfnamein = family + "_2D_tagged.sdf"
sdfnameout = family + ".nmredata.sdf"
# build names from family name for input and output files

reader = Chem.SDMolSupplier(sdfnamein)
writer = Chem.SDWriter(sdfnameout)
# open input and output SDF files for reading and for writing

with open(inputfilename) as fpin:
	for molindex, m in enumerate(reader, start=1):
# iterate through compounds in file familyname_2D_tagged.sdf
		carbons = [a for a in m.GetAtoms() if a.GetSymbol() == 'C']
# make list of carbon atoms
		if not carbons:
			continue
# skip to next molecule if the current one has no carbon atom
		nmrindex, atomstring, cshiftstring = fpin.readline().strip().split('|')
# extract carbon atom index and correspond chemical shifts, as strings
		atoms = atomstring.strip().split()
# get list of carbon atom indexes
		cshifts = [float(x) for x in cshiftstring.strip().split()]
# get list of chemical shifts, in the order of carbon atom indexes
		triplets = zip(atoms, cshifts, atoms)
# prepare write according to NMReDATA format
		nmrlines = ["%s, %.1f, %s \\" % triplet for triplet in triplets]
# NMR data formaating for each carbon atom
		assignment = '\n'.join(nmrlines)
# create property value for SDF
		m.SetProp('NMREDATA_ASSIGNMENT', assignment)
# append NMReDATA assignment tag to the current compound
		writer.write(m)
# write the current NMReDATA record to file familyname.nmredata.sdf
