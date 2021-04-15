"""
The mySDWriter module defines the mySDWriter class
"""
from rdkit import Chem
from rdkit.six import StringIO
import sys

class mySDWriter:
	"""
	mySDWriter is a class that uses RDKit to write molecules to SDF files
	so that no valence is reported in the atom block for electrically charged atoms
	"""
	def __init__(self, fhOut):
		"""
		The only mySDWriter object member is the handle to the output file
		"""
		self.fhOut = fhOut
	
	def correct(self, old_txt):
		"""
		correct() sets to 0 the atom valence indication from the mol block in old_text and
		returns the corrected text
		"""
		new_lines = []
# initialize list of corrected lines in old_txt
		old_lines = old_txt.split('\n')
# get lines from old_txt
		creator_line = old_lines[1].strip()
# line 2, indexed 1 in old_txt, starts with text producer
		creator = creator_line.split()[0]
# get creatot of old_txt
		if creator != 'RDKit':
# creator of old_txt is not RDKit
			return ''
# denies text correction
		count_line = old_lines[3].rstrip()
# line 4, indexed 3 in old_txt, gives the number of atom lines and the file type (V2000 or V3000)
		if 'V2000' not in count_line:
# file type is not V2000
			return ''
# denies text correction 
		numAtom = int(count_line[0:3])
# number of atoms is the first field of the count line
		for i in range(0, 4):
# loop over lines of old_txt preamble
			new_lines.append(old_lines[i])
# no change for preamble lines
		for i in range(4, 4+numAtom):
# loop over atom lines
			line = old_lines[i]
# get current atom line
			pieces = line.split()
# put current atom line into pieces
			valence = pieces[9]
# valence in the 10th field, indexed 9
			if valence != '0':
# a non-zero valence was set, atom is charged
				line = line[0:48]+'  0'+line[51:]
# line surgery to remove valence indication
			new_lines.append(line)
# accounts for the modified atom line
		for i in range(4+numAtom, len(old_lines)):
# loop over the lines that follow the atom lines
			new_lines.append(old_lines[i])
# no change for the lines that follow the atom lines
		new_txt = '\n'.join(new_lines)
# rebuild the text with modified atom lines
		return new_txt
# returns the sdf text without valence idications in the atom lines

	def write(self, m):
		"""
		writes RDit molecule m
		returns 1 if molecule is invalid
		returns 0 otherwise
		"""
		if m is None:
# invalid molecule
			return 1
# return 1 if molecule is invalid
		sio = StringIO()
		w = Chem.SDWriter(sio)
# get a brand new writer to string
		w.write(m)
# write to string molecule with valence indication for the electrically charged atoms
		w.flush()
# make string complete
		old_txt = sio.getvalue()
# get string
		new_txt = self.correct(old_txt)
# remove valence indication
		self.fhOut.write(new_txt)
# write molecule without valence indication for the electrically charged atoms
		return 0
# nothing went wrong

if __name__ == '__main__':
# what else?
	fhOut = sys.stdout
# write on standard output
	wr = mySDWriter(fhOut)
# get a writer for molecules ACD software can read
	smi = '[O-]C(=O)CC[NH3+]'
# SMILES for testing
	print()
	print('Testing with', smi)
	print()
# show SMILES string used for testing
	m = Chem.MolFromSmiles(smi)
# a simple aminoacid, with two electrically charged atoms
	print('*** Before ***')
	print(Chem.MolToMolBlock(m))
# print RDKit mol block, as usual
	print('**** After ***')
	status = wr.write(m)
# write an sdf file on standard output for this aminoacid
	print()
	print ('*** Failed ***' if status else '*** Done ***')
# print status message
