"""
python tautomer.py inputfile.sdf --> inputfile.sdf
python tautomer.py inputfile.sdf outputfile.sdf --> outputfile.sdf

transforms the input file so that "strange tautomers" are converted to "usual tautomeric form".
Mainly alihatic iminols are converted to amides.

RDKit is used to read and write SDF files, so that molecules with electrically charged 
atoms have their mol block containing non-zero valence information ACD software does not understand.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import os
import sys
import tempfile

try:
	from rdkit.Chem import rdCoordGen
	drawingFunction = rdCoordGen.AddCoords
# recent, nice (slow?) 2D coordinates generator
except:
	from rdkit.Chem import AllChem
	drawingFunction = AllChem.Compute2DCoords
# historical 2D coordinates generator

#  Courtesy of Richy Leroy

# iminol
smarts1 = "[NH0]=C([OH])([!O])" # secondary iminol
target1 = Chem.MolFromSmarts(smarts1)
smarts11 = "[NH1]=C([OH])([!O])" # primary iminol
target11 = Chem.MolFromSmarts(smarts11)
smarts12 = "[A&!H][NH0]=C([OH])O" # secondary carbamate
target12 = Chem.MolFromSmarts(smarts12)
# smarts13 = "[NH]=C([OH])O" # primary carbamate
smarts13 = "[NH1]=C([OH])O" # primary carbamate
target13 = Chem.MolFromSmarts(smarts13)
smarts14 = "[CH]([OH])=N" # N-formyl
target14 = Chem.MolFromSmarts(smarts14)

# enol
smarts2 = "[C&!c]([!CX3&!O])([!CX3&!O])=[C&!c]([OH])([!N])"
target2 = Chem.MolFromSmarts(smarts2)

# enethiol
smarts3 = "[C&!c]([!CX3&!SHX2])([!CX3&!SHX2])=[C&!c]([SHX2])"
target3 = Chem.MolFromSmarts(smarts3)
smarts31 = "[N&!n]([!CX3&!SHX2])([!CX3&!SHX2])=[C&!c]([SHX2])"
target31 = Chem.MolFromSmarts(smarts31)


def tautomer(pathIn, pathOut=""):
	"""
	tautomer() reads an SDF file from pathIn and writes tautomerized molecules in place if pathOut is empty
	or to pathOut if the latter is not empty.
	"""
	inplace = len(pathOut) == 0
# determine if transformation looks as if operated in-place
	if inplace:
# in-place-like file transformation
		pathOut = tempfile.mktemp('.sdf')
# get a temporary name for the content of the output SDF file

	fhIn = open(pathIn, 'rb')
# handle to input file
	r = Chem.ForwardSDMolSupplier(fhIn)
# get reader from input file
# I do that in order to be able to close the input, so that it may be overwritten later if tautomer.py is called with 1 argument
	w = Chem.SDWriter(pathOut)
# get writer to output file
	count = 0
# number of compounds read
	new_count = 0
# number of compounds written

	for m in r:
# loop over compounds in input file
		count += 1
# one more compound read
		real = False
# stay False if the current compound do not require tautomer conversion
		if m is None :
# bad molecule
			continue
# jump to next molecule if the current one is bad
		Chem.AddHs(m)
# add H 
		ps = m
# is ps really a copy or simply a synonym of m
		Chem.SanitizeMol(ps)
# cannot harm

# N-formyl
		if ps.HasSubstructMatch(target14):
			real = True
			nb14 = len(ps.GetSubstructMatches(target14))
			#print(nb14, "N-Formyl")
			rxn11 = AllChem.ReactionFromSmarts('[CH1:1]([OH:2])=[N:3]>>[CH1:1](=[OH0D1:2])[NH:3]')
			ps = rxn11.RunReactants((ps,))
			ps = ps[0][0]
			if nb14 != 1 :
				for i in range(1,nb14) :
					ps = rxn11.RunReactants((ps,))
					ps = ps[0][0]
			Chem.SanitizeMol(ps)
			Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
		
# Primary carbamate
		if ps.HasSubstructMatch(target13):
			real = True
			nb13 = len(ps.GetSubstructMatches(target13))
			#print(nb13, "primary carbamate")
			rxn11 = AllChem.ReactionFromSmarts('[#6:1][O:2][CD3:3]([OH1:4])=[NH:5]>>[#6:1][O:2][CD3:3](=[OH0D1:4])[NH2:5]')
			ps = rxn11.RunReactants((ps,))
			ps = ps[0][0]
			if nb13 != 1 :
				for i in range(1,nb13) :
					ps = rxn11.RunReactants((ps,))
					ps = ps[0][0]
			Chem.SanitizeMol(ps)
			Chem.AssignStereochemistry(ps,force=True,cleanIt=True)

# Secondary carbamate
		if ps.HasSubstructMatch(target12):
			real = True
			nb12 = len(ps.GetSubstructMatches(target12))
			#print(nb12, "secondary carbamate")
			rxn11 = AllChem.ReactionFromSmarts('[O:4][CD3:1]([OH:2])=[NH0:3]>>[O:4][CD3:1](=[OH0D1:2])[NH:3]')
			ps = rxn11.RunReactants((ps,))
			ps = ps[0][0]
			if nb12 != 1 :
				for i in range(1,nb12) :
					ps = rxn11.RunReactants((ps,))
					ps = ps[0][0]
			Chem.SanitizeMol(ps)
			Chem.AssignStereochemistry(ps,force=True,cleanIt=True)

# Primary iminol
		if ps.HasSubstructMatch(target11):
			real = True
			nb11 = len(ps.GetSubstructMatches(target11))
			#print(nb11, "primary iminol")
#			rxn11 = AllChem.ReactionFromSmarts('[CD3:1]([OH:2])=[NH:3]>>[CD3:1](=[OH0D1:2])[NH2:3]')
			rxn11 = AllChem.ReactionFromSmarts('[O:2][CD3:3]([OH1:4])=[NH1:5]>>[O:2][CD3:3](=[OH0D1:4])[NH2:5]')
			ps = rxn11.RunReactants((ps,))
			ps = ps[0][0]
			if nb11 != 1 :
				for i in range(1,nb11) :
					ps = rxn11.RunReactants((ps,))
					ps = ps[0][0]
			Chem.SanitizeMol(ps)
			Chem.AssignStereochemistry(ps,force=True,cleanIt=True)			
			if ps.HasSubstructMatch(target11)==True :
				ps = rxn11.RunReactants((ps,))
				ps = ps[0][0]
				Chem.SanitizeMol(ps)
				Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
		
# Scondary iminol
		if ps.HasSubstructMatch(target1)==True :
			real = True
			nb1 = len(ps.GetSubstructMatches(target1))
			#print(nb1, "secondary iminol")
			rxn1 = AllChem.ReactionFromSmarts('[C:1]([OH:2])=[NH0:3]>>[C:1](=[OH0:2])[NH:3]')
			ps = rxn1.RunReactants((ps,))
			ps = ps[0][0]
			if nb1 != 1 :
				for i in range(1,nb1) :
					ps = rxn1.RunReactants((ps,))
					ps = ps[0][0]
			Chem.SanitizeMol(ps)
			Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
			if ps.HasSubstructMatch(target1)==True :
				ps = rxn1.RunReactants((ps,))
				ps = ps[0][0]
				Chem.SanitizeMol(ps)
				Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
			
# enol
		if ps.HasSubstructMatch(target2):
			real = True
			nb2 = len(ps.GetSubstructMatches(target2))
			#print(nb2, "enol")
			rxn2 = AllChem.ReactionFromSmarts('[!c&C:1]([OH:2])=[!c&C:3]>>[C:1](=[OH0:2])[CH:3]')
			ps = rxn2.RunReactants((ps,))
			ps = ps[0][0]
			if nb2 != 1 :
				for i in range(1,nb2) :
					ps = rxn2.RunReactants((ps,))
					ps = ps[0][0]
			Chem.SanitizeMol(ps)
			Chem.AssignStereochemistry(ps,force=True,cleanIt=True)

# enethiol	
		if ps.HasSubstructMatch(target3):
			real = True
			nb3 = len(ps.GetSubstructMatches(target3))
			#print(nb3, "enethiol")
			rxn3 = AllChem.ReactionFromSmarts('[!c&C:1]([SH:2])=[!c&C:3]>>[C:1](=[SH0:2])[CH:3]')
			ps = rxn3.RunReactants((ps,))
			ps = ps[0][0]
			if nb3 != 1 :
				for i in range(1,nb3) :
					ps = rxn3.RunReactants((ps,))
					ps = ps[0][0]
			Chem.SanitizeMol(ps)
			Chem.AssignStereochemistry(ps,force=True,cleanIt=True)

		Chem.AddHs(ps)
# cannot harm
		new_count += 1
# one more compound for writing
		if real :
# a tautomerization took place
#			print("Changed: was", count, "now:", new_count)
			for prop in [(x, m.GetProp(x)) for x in list(m.GetPropNames())]:
# loop over SD tags and values in m 
				ps.SetProp(*prop)
# copy properties from m in ps. Why did they disappear in ps if ps is simply a copy of m?
			realinchi = Chem.MolToInchi(ps,options="-FixedH")
# generate InChI with fixed position of H atoms
			#print(realinchi)
			ps.SetProp("InChI_HFIXED", realinchi)
# save InChI with fixed H
			drawingFunction(ps)
# redraw the tautmeric form
		w.write(ps)
# draw the current molecule, thus introducing the "valence bug" for molecules with electrically charged atoms

	fhIn.close()
	w.close()
# close input and outpul file
	if inplace:
# pseudo in-place processing
		os.replace(pathOut, pathIn)
# for in-place processing, copy (possibly overwrite) the temporary output file with the name of the input file

if __name__ == '__main__':
# what else?
	argv = sys.argv
	argc = len(argv)
# get command line arguments and the number of arguments
	if (argc != 3) and (argc != 2):
# unexpected number of arguments
		print(argv[0]+': bad number of command line arguments', file=sys.stderr)
		print('Usage:', argv[0], 'path_to_input_file path_to_output_file', file=sys.stderr)
# error message on standard error
		sys.exit(1)
# exit on error
	pathIn = argv[1]
	pathOut = argv[2] if argc == 3 else ""
	tautomer(pathIn, pathOut)
# does the job
