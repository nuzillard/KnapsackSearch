"""
novalence.py reads .sdf files in V2000 format produced by RDKit and resets to 0
the valence information that is introduced for electrically charged atoms.
The resulting file can be read by ACD software without error.

TO BE REWRITTEN WITH SDFRW
"""
import sys
from rdkit import Chem
import mySDWriter

def novalenceFromHandle(fhIn, fhOut):
	"""
	novalenceFromHandle() reads through file handle fhIn a .sdf file in V2000 format produced by RDKit
	and writes through file handle fhOut the same file but with atom valence set to 0
	"""
	with fhIn, fhOut:
		reader = Chem.ForwardSDMolSupplier(fhIn)
		writer = mySDWriter.mySDWriter(fhOut)
		for m in reader:
			if m:
				writer.write(m)

def novalenceFromPath(pathIn, pathOut):
	"""
	novalenceFromHandle() reads through file path pathIn a .sdf file in V2000 format produced by RDKit
	and writes through file path pathOut the same file but with atom valence set to 0
	"""
	fhIn = open(pathIn, 'rb')
	fhOut = open(pathOut, 'w')
	novalenceFromHandle(fhIn, fhOut)

if __name__ == '__main__':
# what else?
	argv = sys.argv
	argc = len(argv)
# get command line arguments and the number of arguments
	if argc != 3:
# unexpected number of arguments
		print(argv[0]+': bad number of command line arguments', file=sys.stderr)
		print('Usage:', argv[0], 'path_to_input_file path_to_output_file', file=sys.stderr)
# error message on standard error
		sys.exit(1)
# exit
	novalenceFromPath(argv[1], argv[2])
# does the job
