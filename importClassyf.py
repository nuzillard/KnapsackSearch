"""
importClassyf.py inserts classification data from ClassyFire into the .sdf file it originates from.

3 arguments: python importClassyf filenameIn filenameClass filenameOut: filenameIn and filenameClass -> filenameOut
2 arguments: python importClassyf filenameIn filenameClass: filenameIn and filenameClass -> filenameIn
"""
import sdfrw
import sys
import os
import tempfile

def importClassyf(filenameIn, filenameClass, filenameOut=""):
	"""
	importClassyf() imports the classification data from filenameClass into filenameOut
	If filenameOut is empty or points to the same file as filenameIn, the result is placed in file filenameIn
	
	"""
	if not filenameOut or (os.path.isfile(filenameOut) and os.path.samefile(filenameIn, filenameOut)):
# 2 arguments case
		filenameInterm = tempfile.NamedTemporaryFile().name
# get a name for a temp file
	else:
# 3 arguments case
		filenameInterm = filenameOut
# use thirs argument as output file name
	with open(filenameIn) as fhIn, open(filenameClass) as fhCl, open(filenameInterm, 'w') as fhOut:
# open all files
		for (mol, molCl) in zip(sdfrw.sdfReader(fhIn), sdfrw.sdfReaderVariant(fhCl)):
# read molecules simulanously from filenameIn and filenameClass
			props = sdfrw.sdfGetPropList(molCl)
# get list of molecule properties from file with classification data
			for prop in props:
# loop through classification property names
				value = sdfrw.sdfGetProp(molCl, prop)
# get current property value
				if value:
# the current property value is not empty
					value = ' | '.join(value.split('\t'))
# replace tabs by vertical bars in property value
				prop = "ClassyFire_" + prop
# change the property name to make its origin apparent
				done, mol = sdfrw.sdfSetChangeProp(mol, prop, value)
# set a new property in the current molcule
			sdfrw.sdfWrite(fhOut, mol)
# write the current molecule to output (3 arguments cas) or to the temp file (2 arguments case)
	if filenameInterm != filenameOut:
# this is the 2 arguments case, the temp file must be copied to the input file
		with open(filenameInterm, 'r') as fhIn, open(filenameIn, 'w') as fhOut:
# open files for copy
			for line in fhIn:
# get line from input (temp file)
				fhOut.write(line)
# write current line to input file
		if os.path.isfile(filenameInterm):
# is the temp file there?
			os.unlink(filenameInterm)
# delete temp file

if __name__ == "__main__":
#	importClassyf('tazettine_uniq.sdf', 'tazettine_uniq_classified.sdf', 'tazettine_uniq_merged.sdf')
#	importClassyf('tazettine_uniq.sdf', 'tazettine_uniq_classified.sdf')
#	importClassyf('tazettine_uniq.sdf', 'tazettine_uniq_classified.sdf', r"..\LSDos-3.4.11\tazettine_uniq.sdf")
	argc = len(sys.argv)
	if argc != 3 and argc != 4:
# bad number of command line arguments
		print("Usage: python %s file-to-be-completed-by-classification file-from-ClassyFire [file-with-classification]" % (sys.argv[0],))
# print help message
		sys.exit(1)
# exit on error
	importClassyf(*sys.argv[1:])
# import ClassiFire classyfication data into input file and save result.
