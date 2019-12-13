"""
make_2D_sdf.py:
	reads: familyname_pickled
	writes: familyname_2D.sdf

familyname_pickled: written by compounds.py
familyname_2D.sdf: read by basic_tags.py

called by: python make_2D_sdf.py familyname

File familyname_pickled contains data about all compounds mentioned
in KNApSAcK about the genera of family familyname, according to familyname_genera.txt.
File familyname_2D.sdf is an SDF file that contains all these compounds,
described only with 2D MolBlocks, including 3D geometry hints.
"""

import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

family = sys.argv[1]
# get family name from command line
inputfilename = family + "_pickled"
outputfilename = family + "_2D.sdf"
# build names from family name for input and output files

with open(inputfilename, "rb") as fpin, open(outputfilename, 'w') as sdfile:
	data = pickle.load(fpin)
# get data for the whole set of compounds
	for compound in data:
# each compound is a dictionary whose keys were defined in source code file compounds.py
		cid = compound["C_ID"]
		smiles = compound["SMILES"]
		inchi = compound["InChICode"]
# get compound_id, SMILES, and InChI for the current compound
		if not smiles:
			print(cid, "SMILES is empty")
			continue
# next compound if no SMILES (should not happen?)
		m = Chem.MolFromSmiles(smiles)
# RDKit molecule construction from SMILES decoding
		if m is None:
			print(cid, "Cannot build molecule from SMILES")
			continue
# next compound if molecule construction failed. May happen.
		inchi_again = Chem.MolToInchi(m)
# determine InChI from molecule
		if(inchi_again != inchi):
# check for the identity of the InChI provided by KNApSAcK and the one recalculated by RDKit from SMILES
			print(cid, "Roundtrip failed")
			print(cid, "InChI_1", inchi)
			print(cid, "InChI_2", inchi_again)
			continue
# next compound if the two InChIs are not identical. This is the less tolerant option.
		m.SetProp("_Name", cid)
# give the compound_id as title to the MolBlock
		try:
			AllChem.Compute2DCoords(m)
# calculate 2D coordinates
		except:
			print("Cannot generate structure diagram for molecule")
			continue
# next compound if 2D coordinate generation failed (why should this happen?)
		sdfile.write(Chem.MolToMolBlock(m) + "$$$$\n")
# write MolBlock of the current compound to file familyname_2D.sdf and make an SDF record of it.
