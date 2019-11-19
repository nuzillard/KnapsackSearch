"""
process.py;
	reads: familyname_genera.txt
	writes: familyname.nmredata.sdf

familyname_genera.txt: from wikipedia, interrogation by organism family. See example files and source code of family_from_web.py.
familyname.nmredata.txt: for EdiSdf or any other SDF file reader. Includes MolBlock, InChI, SMILES, 
	formula, Names, molecular mass, presence of nitrogen, organism binomial_name, and 13C NMR chemical shift prediction
	in NMReDATA assignment format.

called by: python process.py familyname

The file familyname_genera.txt is transformed into the file familyname.nmredata.txt
using through-web requests to KNApSAcK, molecular data conversion and drawing by RDKit, and 13C NMR chemical shift calculation by nmrshiftdb2
"""

import sys
import os

def run_command(command):
	"""
	run_command() forwards a command to the OS through a call to os.system()
	"""
	print("Running:",  command)
# print command to run
	code = os.system(command)
# run command
	print("\"%s\" returned with code: %d\n" % (command, code))
# print code returned by OS, 0 for Windows, something else for Linux/MacOS (none tested)

family = sys.argv[1]
# get family name from command line

command = "python family_from_web.py " + family
run_command(command)
"""
family_from_web.py:
	reads: familyname_genera.txt
	writes: familyname_cid_species.txt
"""

command = "python family_cids.py " + family
run_command(command)
"""
family_cids.py:
	reads: familyname_cid_species.txt
	writes: familyname_cids.txt
	writes: familyname_cid_all_species.txt
"""

command = "python compounds.py " + family
run_command(command)
"""
compounds.py:
	reads: familyname_cids.txt
	writes: familyname_pickled
"""

command = "python make_2D_sdf.py " + family
run_command(command)
"""
make_2D_sdf.py:
	reads: familyname_pickled
	writes: familyname_2D.sdf
"""

command = "python basic_tags.py " + family
run_command(command)
"""
basic_tags.py:
	reads: familyname_2D.sdf
	reads: familyname_pickled
	reads: familyname_cid_all_species.txt
	writes: familyname_2D_tagged.sdf
"""

command = "predictSdf " + family + "_2D_tagged.sdf 4 3d 1>" + family + "_2D_nmr.txt 2>errorlog.txt"
run_command(command)
"""
predict_sdf.bat
	reads: familyname_2D_tagged.sdf
	writes: standard output, redirected to familyname_2D_nmr.txt
	writes: standard error, redirected to errorlog.txt
"""

command = "python molsort.py " + family
run_command(command)
"""
molsort.py:
	reads: familyname_2D_nmr.txt
	writes: familyname_2D_nmr_sorted.txt
"""

command = "python nmr_tags.py " + family
run_command(command)
"""
nmr_tags.py:
	reads: familyname_2D_tagged.sdf
	reads: familyname_2D_nmr_sorted.txt
	writes: familyname.nmredata.sdf
"""
