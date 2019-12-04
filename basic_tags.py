"""
basic_tags.py:
	reads: familyname_2D.sdf
	reads: familyname_pickled
	reads: familyname_cid_all_species.txt
	writes: familyname_2D_tagged.sdf

familyname_2D.sdf: written by make_2D_sdf.py
familyname_pickled: written by compounds.py
familyname_cid_all_species.txt: written by family_cids.py
familyname_2D_tagged.sdf: read by familyname_2D_tagged.sdf

called by: python basic_tags.py familyname

File familyname_2D.sdf contains all the MolBlocks of the compounds associated to
a family as reported in file familyname_genera.txt.
basic_tags.py supplements familyname_2D.sdf with binomial_data from file familyname_cid_all_species
and with structural data from file familyname_pickled.
"""

import pickle
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import sys

family = sys.argv[1]
# get family name from command line
inputfilename1 = family + "_pickled"
inputfilename2 = family + "_cid_all_species.txt"
sdfnamein = family + "_2D.sdf"
sdfnameout = family + "_2D_tagged.sdf"
# build names from family name for input and output files

with open(inputfilename1, 'rb') as fpin:
	data = pickle.load(fpin)
# get compound data from pickled file
d = {compound["C_ID"]: compound for compound in data}
# make a directory with compound_ids as keys and the corresponding dictionary as value

with open(inputfilename2) as fpin:
	cids_all_species_dict = dict([line.strip().split(' : ') for line in fpin])
# make a dictionary with compound_ids as keys and binomial_names as values

reader = Chem.SDMolSupplier(sdfnamein)
writer = Chem.SDWriter(sdfnameout)
# SDF file input and output

for m in reader:
# get molecule for input SDF file
	cid = m.GetProp("_Name")
# get compound_id stored as compound title
	tags = d[cid]
# get structural data from familyname_pickled as stored in d
	for key, value in tags.items():
		m.SetProp(key, value)
# transfer key-value pairs from familyname_pickled as SDF tag-value pairs
	nitrogens = [a for a in m.GetAtoms() if a.GetSymbol() == 'N']
	m.SetProp('Contains_Nitrogen', 'True' if nitrogens else 'False')
# append a new SDF tag with True/False value for the presence of at least one nitrogen atom
	formula = CalcMolFormula(m)
	m.SetProp('Formula', formula)
# replace formula fom KNApSAcK by the one from RDKit, because the former might be wrong (4 Dec 2019)
	if cid in cids_all_species_dict:
		m.SetProp('In_Species', cids_all_species_dict[cid])
# append a new SDF tag for the presence of the compound in organism binomial_name
	writer.write(m)
# store current compound in familyname_2D_tagged.sdf
