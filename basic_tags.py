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
from rdkit.Chem.Descriptors import ExactMolWt
import sys
import os
import tempfile

import split_mf

def compare_MF(oldmf, newmf):
	"""
	Returns a Boolean value that indicates whether molecular formulas
	oldmf (KNApCAcK) and newmf (from RDMol from SMILES from KNApCAcK are compatible.
	"""
	if oldmf == newmf:
		return True
# equality is perfect
	signpos = newmf.find('+')
	if signpos > -1:
		shortnewmf = newmf[0:signpos]
		if oldmf == shortnewmf:
			print("Fixed problem: positive charge was missing in formula from KNApSAcK")
			return True
# KNApSAcK MF systematically "forgets" the molecular electric charge. Check only for elements.
	signpos = newmf.find('-')
	if signpos > -1:
		shortnewmf = newmf[0:signpos]
		if oldmf == shortnewmf:
			print("Fixed problem: negative charge was missing in formula from KNApSAcK")
			return True
# KNApSAcK MF systematically "forgets" the molecular electric charge. Check only for elements.
	ok_MS = split_mf.match_MS_MplusH(oldmf, newmf)
	if ok_MS:
		print("Fixed problem: formula from KNApSAcK was M+H")
	return ok_MS

# first stage

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

tempfn = tempfile.mktemp('.sdf')
# tempfn is a name for a temp file

reader = Chem.SDMolSupplier(sdfnamein)
writer = Chem.SDWriter(tempfn)
# SDF file input and output. Temp file contains molecules fixed in a first stage for molecular formula errors
# but not for duplicates. Second stage processes file entries concerning identical structures.

inchis = []
# list of the InChIs of the molecules stored in the order they appear in the temp file
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
	old_formula = m.GetProp('Formula')
# recall formula from KNApSAcK
	new_formula = CalcMolFormula(m)
# calculate formula from molecule
	if old_formula != new_formula:
# process compounds with mismatch between stored and recalculated molecular formula
		print("%s: Formula mismatch. Stored: %s. Recalculated: %s" % (cid, old_formula, new_formula))
# print warning message
	if not compare_MF(old_formula, new_formula):
		continue
# proceed to next compound
	m.SetProp('Formula', new_formula)
# fix molecular formula for charge ou M+H error
	m.SetProp('Mw', str(ExactMolWt(m)))
# recalculate exact molecular weight, according to RDKit
	if cid in cids_all_species_dict:
		m.SetProp('In_Species', cids_all_species_dict[cid])
# append a new SDF tag for the presence of the compound in organism binomial_name
	writer.write(m)
# store current compound in temporary .sdf file
	inchis.append(m.GetProp('InChIKey'))
# apprnd to InChI of the stored molecule to the list of InChIs of stored molecules
writer.close()
# close temporary output sdf file

# second stage

inchidict = {}
# prepare an empty dictioanry with InChiKeys as key an a list of compound index as value.
# If the were no replicated structures, all values should be one-element lists.
for i, inchi in enumerate(inchis, start=1):
# scan through the list of InChIKey prepared during first stage
	inchidict[inchi] = inchidict[inchi] + [i] if inchi in inchidict else [i]
# create or extend index list related to the current InChiKey
irregular_dict = {inchi: list for (inchi, list) in inchidict.items() if len(list) != 1}
# select from inchidict only the items that signal replicated compounds
has_irregulars = bool(irregular_dict)
# Boolean that signals the existence of replicates

writer = Chem.SDWriter(sdfnameout)
with open(tempfn, 'rb') as tmp:
	reader = Chem.ForwardSDMolSupplier(tmp)
# read ftom temp file with possible replicates, write to familyname_2D_tagged.sdf without replicates

	if not has_irregulars:
# no replicates found
		for m in reader:
# scan through compounds in temp file
			writer.write(m)
# copy molecules in temp file to familyname_2D_tagged.sdf, compound by compound

	else:
# replicated InChIKeys were found
		irregular_inchis = list(irregular_dict)
# list of InChIKeys for which replicated structures were found
		irregular_lists = [irregular_dict[key] for key in irregular_inchis]
# list of the list of indexes that correspond to replicated compounds, InChIKey per InChIKey
		irregular_indexes = sum(irregular_lists, [])
# list of all index of compounds that have replicates
		irregular_mols_dict = {inchi: [] for inchi in irregular_inchis}
# prepare a dictionary with InChIKeys as keys and list of correponding replicated compounds
		for (i, m) in enumerate(reader, start=1):
# scan through compounds and compound indexes from temp file
			if i in irregular_indexes:
# the current molecule has replicates
				inchi = m.GetProp('InChIKey')
# get the InChIKey of the current compound, which has replicates
				irregular_mols_dict[inchi].append(m)
# append the current compound to the list of replicates and according to InChIKey
# processing for this compound is delayed (see below)
			else:
# regular compound
				writer.write(m)
# write current regular compound to output file

		print("Replicates:")
# signals the processing of the compounds whose processing was delayed because thet have
# replicated structures within the file
		for inchi in irregular_inchis:
# scan through InChIKey for which there are replicated compounds
			mols = irregular_mols_dict[inchi]
# get list of compounds that share the current InChIKey
			cids = '||'.join([m.GetProp('C_ID') for m in mols])
# join the KNApSAcK ids of the replicated compounds
			print('\t'+inchi+': '+cids.replace('||', ', '))
# signals which InChIKey is currently processed, along with the list of corresponding KNApSAcK ids
			names = '||'.join(sorted(list(set(sum([m.GetProp('Name').split('||') for m in mols], [])))))
# join the names of the replicated compounds
			species = '||'.join(sorted(list(set(sum([m.GetProp('In_Species').split('||') for m in mols], [])))))
# join the binomial names of the replicated compounds
#			cas_rns = '||'.join(sorted(list(set(sum([m.GetProp('CAS RN').split('||') for m in mols], [])))))
			cas_rns = '||'.join([m.GetProp('CAS RN') for m in mols])
# join the CAS Registry Numbers of the replicated compounds
			mol = mols[0]
# new compound attributes will be assigned to the first compound in the list of replicates
			mol.SetProp('_Name', cids)
			mol.SetProp('C_ID', cids)
			mol.SetProp('Name', names)
			mol.SetProp('In_Species', species)
			mol.SetProp('CAS RN', cas_rns)
# assign new property values to the first compound in the list of replicates
			writer.write(mol)
# writes the first compound in the list of replicates to the output file.
# The other compounds in the the list are left over
writer.close()
# close familyname_2D_tagged.sdf
os.unlink(tempfn)
# unlink temp file
