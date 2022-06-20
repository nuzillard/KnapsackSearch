import csv

fileIn = "220525_frozen_metadata.csv"
d = dict()
# a dictionnary with QIDs as keys. '' as values. for QID unicity enforcement.
with open(fileIn, encoding='utf-8') as csvfile:
	csvreader = csv.reader(csvfile)
	for irow, row in enumerate(csvreader, start=1):
		if irow == 1:
# row 1 contains CSV headers
			continue
# skip line with headers
		entity = row[0][31:]
# field 0 (index start at 0) contains web link to Wikidata.
# the characters after 31 (index start at 0) constitute the wikidata QID (or entity)
		SMILES_orig = row[3]
# the SMILES of the compound
		if entity in d:
# QID already there
			continue
# nothing to do with the current line in fileIn
		d[entity] = ''
# remember the current QID was already processed
		print (SMILES_orig + ' \t' + entity + '_' + str(irow))
# format for the reading of a .smi file
