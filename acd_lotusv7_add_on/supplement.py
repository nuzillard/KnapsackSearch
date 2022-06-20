import sdfrw
import csv
from rdkit import Chem

sdfIn = "acd_lotusv7_tmp.sdf"
csvIn = "220525_frozen_metadata.csv"
sdfOut = "acd_lotusv7.sdf"

def get_row(csvreader, Q_target):
	"""
	search in csvfile for a line with a QID identical to Q_target
	"""
	search = True
	while search:
		row = next(csvreader)
		Q = row[0][31:]
		search = (Q != Q_target)
	return row
# return a line splitted by the CSV reader as a list of fields
		

with open(sdfIn) as fhIn, open(csvIn, encoding='utf-8') as fhCsv, open(sdfOut, "w", encoding='utf-8') as fhOut:
# open all files
	csvreader = csv.reader(fhCsv)
# get CSV reader
	firstrow = next(csvreader)
# ignore first line with headers inside
	for mol in sdfrw.sdfReader(fhIn):
# iterate through mols in input SDF file
		title = sdfrw.sdfGetTitle(mol)
		Q_ID, lineno = title.split('_')
# get QID of current compound
		row = get_row(csvreader, Q_ID)
# find the CSV data for the current compound
		inchikey = row[1]
		inchi = row[2]
		smiles = row[3]
		exact_mass = row[5]
		xlogp = row[6]
		npclassifier_01pathway = row[12]
		npclassifier_02superclass = row[13]
		npclassifier_03class = row[14]
		classyfire_chemontid = row[15]
		classyfire_01kingdom = row[16]
		classyfire_02superclass = row[17]
		classyfire_03class = row[18]
		classyfire_04directparent = row[19]
# extract molecular data from CSV LOTUS file
		ok, mol = sdfrw.sdfSetChangeProp(mol, "inchikey", inchikey)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "inchi", inchi)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "smiles", smiles)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "exact_mass", exact_mass)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "xlogp", xlogp)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "npclassifier_01pathway", npclassifier_01pathway)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "npclassifier_02superclass", npclassifier_02superclass)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "npclassifier_03class", npclassifier_03class)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "classyfire_chemontid", classyfire_chemontid)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "classyfire_01kingdom", classyfire_01kingdom)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "classyfire_02superclass", classyfire_02superclass)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "classyfire_03class", classyfire_03class)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "classyfire_04directparent", classyfire_04directparent)
# copy current molecular data to current molecule
		mb = sdfrw.sdfGetMolBlock(mol)
		rdmol = Chem.MolFromMolBlock(mb)
# get molblock of current compound from input SDF file
		smiles_recalculated = Chem.MolToSmiles(rdmol)
		inchi_recalculated = Chem.MolToInchi(rdmol)
		inchikey_recalculated = Chem.InchiToInchiKey(inchi_recalculated)
# recalculate molecular identifier
		inchikey_match = str(inchikey == inchikey_recalculated)
# check for molecular identity
		ok, mol = sdfrw.sdfSetChangeProp(mol, "smiles_recalculated", smiles_recalculated)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "inchi_recalculated", inchi_recalculated)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "inchikey_recalculated", inchikey_recalculated)
		ok, mol = sdfrw.sdfSetChangeProp(mol, "inchikey_match", inchikey_match)
# copy recalculated molecular identifiers to current molecule
		sdfrw.sdfWrite(fhOut, mol)
# write current molecule to SDF output file
