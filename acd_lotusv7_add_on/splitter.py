"""
Splits a .sdf file in ../ into pieces in the current directory.
The input file name and the root name of the pieces are hard-coded hereafter.
The size of the pieces is also hard-coded.
"""

import sdfrw

filenameIn = '..\\fake_acd_lotusv7.sdf'
rootOut = 'fake_acd_lotusv7_'
sizePiece = 10000

with open(filenameIn) as fhIn:
# open input file for reading
	for imol, mol in enumerate(sdfrw.sdfReader(fhIn)):
# iterate through molecules and associated index (start at 0)
		ilocal = imol % sizePiece
# index in current output file (possibly not created yet)
		ipiece = imol // sizePiece
# index of current output file (possibly not created yet)
		if ilocal == 0:
# a new output file will be created 
			if ipiece != 0:
# a new ouput file will be created and an old one still exist
				fhOut.flush()
				fhOut.close()
# flush and close the completed existing output file
			fhOut = open(rootOut + ("%02d" % (ipiece,)) + '.sdf', 'w')
# create a new output file with a new name deduced from its index (start at 0)
			print(ipiece)
# user information about the index of the new output file (start at 0)
		sdfrw.sdfWrite(fhOut, mol)
# write current molecule in the current output file
	fhOut.flush()
	fhOut.close()
# flush and close the last output file.

