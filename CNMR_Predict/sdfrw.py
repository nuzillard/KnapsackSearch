"""
The module sdfrw.py provides:
- sdfReader(fh) yields compounds stored in an opened sdf file (this is an iterator, not a function)
- sdfWrite(fh, mol) that writes a single compound, mol, to a previously opened text file accessed to by file handle fh
- sdfTransform(filenameIn, filenameOut, transformation=defaultTransformation) that reads molecules,
transforms them, and writes them in an other file. The defaultTransformation leaves compounds unchanged.
- sdfHasProp(mol, sdfprop) --> Boolean
- sdfGetPropList(mol) --> List of strings
- sdfGetProp(mol, sdfprop) --> String
- sdfClearProp(mol, sdfprop) --> molecule
- sdfClearAllProps(mol) --> molecule
- sdfSetChangeProp(mol, sdfprop, sdfvalue) --> (Boolean, molecule)
- sdfSetNoChangeProp(mol, sdfprop, sdfvalue) --> (Boolean, molecule)
- sdfEmptyMol() --> molecule
- sdfGetMolBlock(mol) --> String
- sdfSetMolBlock(mol, molblock) --> molecule
- sdfGetTitle(mol) --> String
- sdfChangeTitle(mol, newtitle) --> molecule
sdfrw.py does not depend on any toolkit of cheminformatic functions

Reading and writing .sdf files with molecule storage as a dictionary with two keys:
"molblock", giving access to the molblock without the final '\n' and
"keyvals", giving access to a list of pairs, the first value is a sdf property name
and the second value is a whole property and value text without the final '\n\n'.
The choice of a list to represent sdf tag and value pairs was made to ensure
that their order is preserved during reading and writing operations.

example: python sdfrw.py flavonoids.sdf
result: creates copied_flavonoids.sdf, a copy of flavonoids.sdf
"""

import os
import sys

def getMolFromText(text):
	"""
	getMolFromText() transforms a text which is a section of an .sdf file relative to a molecule
	into a dictionary that represents a molecule.
	"""
	lines = text.rstrip().split('\n')
# lines from text without trailing '\n'
	for i, line in enumerate(lines):
		limit = i
		if line == "M  END":
			break
# find limit, which is the endex of the line that terminates the molblock.
# The list.index() fonction could have been used here
	molblock = '\n'.join(lines[0:limit+1])
# reconstruct the molblock, with no trailing '\n'
	tagsblock = '\n'.join(lines[limit+1:])
# reconstruct the sdf tag block
	items = tagsblock.split('\n\n')[:-1]
# build the list of the key and value text blocks, without the '$$$$' line and without trailing space
	keyvals = []
# extract the tags and build the list of pairs (tag, tag-and-value) in keyvals
	for item in items:
# loop through tag-and-value texts
		less_sign_pos = item.find('<')
# a tag identifier starts with <
		second_more_sign_pos = item.find('>', less_sign_pos)
# a tag identifier stops with the > that follows the <
		keyvals.append((item[less_sign_pos+1:second_more_sign_pos], item))
# extract tag identifier from tag-and-value and updates keyvals with tag identifier and tag-and-value
	return {"molblock": molblock, "keyvals": keyvals}
# return the dictionary that represents the molecule described by text

def sdfReader(fh):
	"""
	sdfReader is an iterator for dictionary representation of molecules
	from an opened .sdf file accessible by filehandle fh
	"""
# open input file from its name
	text = ""
# text contains the current text block relative to the current molecule, initially empty
	for line in fh:
# iterate through the lines of the input file
		if not line:
# no line
			break
# exit text file reading when no more line to be read
		text += line
# append current line to current text
		if line.strip() == "$$$$":
# test if the most recently line is the end-of-molecule line
			yield getMolFromText(text)
# send the molecule representation of the current molecule
			text = ""
# reset text for next molecule, if any

def sdfWrite(fh, mol):
	"""
	sdfWrite() writes the text transcription of molecule dictionary representation mol
	in file a previously opened for writing and which operates through filehandle fh
	"""
	fh.write(mol["molblock"])
# write the molblock
	fh.write('\n')
# write the trailing '\n' of the molblock
	for key, val in mol["keyvals"]:
# loop through sdf data in the order of keys
		fh.write(val)
# write the key-and-value text	
		fh.write('\n\n')
# write one '\n' after the key-and-value text and one more '\n' the the blank line that follows it
	fh.write('$$$$\n')
# finally write the end of molecule sign and its trailing '\n'

def sdfHasProp(mol, sdfprop):
	"""
	sdfHasProp() returns a boolean that indicates the presence of property sdfprop in molecule mol
	"""
	sdfkeyvals = mol["keyvals"]
	return sdfprop in [pair[0] for pair in sdfkeyvals] if sdfkeyvals else False

def sdfGetPropList(mol):
	"""
	sdfGetPropList() returns the list of all property names in molecule mol
	"""
	sdfkeyvals = mol["keyvals"]
	return [pair[0] for pair in sdfkeyvals] if sdfkeyvals else []

def sdfGetProp(mol, sdfprop):
	"""
	returns the value of property sdfprop in molecule mol if this property is defined and returns None otherwise
	"""
	sdfkeyvals = mol["keyvals"]
	if not sdfkeyvals:
		return None
	sdfkeys = [pair[0] for pair in sdfkeyvals]
	if sdfprop not in sdfkeys:
		return None
	where = sdfkeys.index(sdfprop)
	value = sdfkeyvals[where][1]
	eol_pos = value.find('\n')
	value = value[eol_pos+1:]
	return value
	
def sdfClearProp(mol, sdfprop):
	"""
	sdfClearProp() tries to remove property sdfprop from molecule mol.
	returns a pair made of a Boolean that indicates success and of the possibly modified molecule
	"""
	sdfkeyvals = mol["keyvals"]
	if not sdfkeyvals:
		return (False, mol)
	sdfkeys = [pair[0] for pair in sdfkeyvals]
	if sdfprop not in sdfkeys:
		return (False, mol)
	where = sdfkeys.index(sdfprop)
	cleared = sdfkeyvals.pop(where)
	mol["keyvals"] = sdfkeyvals
	return (True, mol)

def sdfClearAllProps(mol):
	"""
	sdfClearAllProps() returns molecule mol without any property inside
	"""
	mol["keyvals"] == []
	return mol

def sdfSetChangeProp(mol, sdfprop, sdfvalue):
	"""
	sdfSetChangeProp() sets property sdfprop to value sdfvalue for molecule mol.
	If the sdfprop is already defined in mol, the value is changed.
	returns a pair made of Boolean True and of the new molecule
	"""
	sdfkeyvals = mol["keyvals"]
	sdfkeys = [pair[0] for pair in sdfkeyvals] if sdfkeyvals else []
	newpair = (sdfprop, ">  <"+sdfprop+'>\n'+sdfvalue)
	if sdfprop not in sdfkeys:
		sdfkeyvals.append(newpair)
	else:
		where = sdfkeys.index(sdfprop)
		sdfkeyvals[where] = newpair
	mol["keyvals"] = sdfkeyvals
	return (True, mol)

def sdfSetNoChangeProp(mol, sdfprop, sdfvalue):
	"""
	sdfSetNoChangeProp() set property sdfprop to value sdfvalue for molecule mol.
	If the sdfprop is already defined in mol, the value is not changed.
	returns a pair made of Boolean and of the new molecule.
	The Boolean is False only if the property sdfprop was already defined in mol 
	"""
	sdfkeyvals = mol["keyvals"]
	sdfkeys = [pair[0] for pair in sdfkeyvals] if sdfkeyvals else []
	if sdfprop not in sdfkeys:
		newpair = (sdfprop, ">  <"+sdfprop+'>\n'+sdfvalue)
		sdfkeyvals.append(newpair)
		mol["keyvals"] = sdfkeyvals
		return (True, mol)
	else:
		return (False, mol)

def sdfEmptyMol():
	"""
	sdfEmptyMol() returns a molecule without any MOL block and  any property
	"""
	return {"molblock": "", "keyvals": []}
	

def sdfGetMolBlock(mol):
	"""
	sdfGetMolBlock() returns the MOL block of the molecule
	"""
	return mol["molblock"]
	
def sdfSetMolBlock(mol, molblock):
	"""
	sdfSetMolBlock() sets the MOL block of molecule mol to molblock
	"""
	mol["molblock"] = molblock
	return mol

def sdfGetTitle(mol):
	"""
	sdfGetTitle() returns the title of molecule mol
	"""
	molblock = mol["molblock"]
	if not molblock:
		return None
	first_eol_pos = molblock.find('\n')
	return molblock[0:first_eol_pos]

def sdfChangeTitle(mol, newtitle):
	"""
	sdfChangeTitle() returns mol with title newtitle
	"""
	molblock = mol["molblock"]
	if not molblock:
		return None
	first_eol_pos = molblock.find('\n')
	molblock = newtitle + molblock[first_eol_pos:]
	mol["molblock"] = molblock
	return mol

def defaultTransformation(mol):
	"""
	defaultTransformation() returns the mol it receives as parameter.
	For testing, corrsponds to a .sdf copy action without any molecule transformation
	"""
	return mol

def sdfTransform(filenameIn, filenameOut, transformation=defaultTransformation):
	"""
	sdfTransform() copies .sdf file named filenameIn into a .sdf file named filenameOut.
	The molecules in input are transformed according to function transformation()
	Not really useful but was used as part of the testing of sdfReader() and sdfWrite()
	"""
	with open(filenameIn) as fhIn, open(filenameOut, "w") as fhOut:
# open for writing transformed molecules
		for molIn in sdfReader(fhIn):
# read input molecule
			molOut = transformation(molIn)
# transform input molecule to output molecule
			sdfWrite(fhOut, molOut)
# write out molecule

if __name__ == "__main__":
# testing module functions for sdf file reading and writing
	argc = len(sys.argv)
# get number of command line arguments + 1
	if argc != 2:
# only one argument possible, the name of the .sdf input file
		print("Usage: %s: filename.sdf" % (argv[0],))
		sys.exit(1)
	filenameIn = sys.argv[1]
# get input .sdf file name
	head, tail = os.path.split(filenameIn)
# prepare the creation of the output file name
	filenameOut = os.path.join(head, 'copied_' + tail)
# prepend copied_ to input file name to obtain the output file name
	sdfTransform(filenameIn, filenameOut)
# make copy
