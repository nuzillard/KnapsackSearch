"""
classyfy.py submits a .sdf file to ClassyFire through the pyclassyfire python module
and receives a .sdfile

2 arguments: python classyfy.py filenameIn.sdf filenameOut.sdf: filenameIn.sdf -> filenameOut.sdf
1 argument : python classyfy.py filenameIn.sdf: filenameIn.sdf -> filenameIn_classified.sdf
"""
from pyclassyfire import client
from datetime import datetime
import os
import sys

def get_now():
	"""
	returns formatted localtime with millisecond accuracy
	"""
	return datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

def classyfy(filenameIn, filenameOut=""):
	"""
	classyfy() organizes the call to client.sdf_query(filenameIn, filenameOut)
	"""
	if not filenameOut:
# no filenameOut given
		(root, ext) = os.path.splitext(filenameIn)
# separates filenameOut into root and extension
		root += "_classified"
# appends _classified to filenameIn root if no filenameOut is given
		filenameOut = root + ext
# make filenameOut if not given
	print("Started: " + get_now())
	client.sdf_query(filenameIn, filenameOut)
# submit request to ClassyFire through the pyclassyfire python module
	print("Finished: " + get_now())

if __name__ == '__main__':
	argc = len(sys.argv)
# count command ine arguments (+1)
	if argc != 2 and argc != 3:
# bad number of command line arguments
		print("Usage: python %s file-sent-to-ClassyFire [file-from-ClassyFire]" % (sys.argv[0],))
# print help message
		sys.exit(1)
# exit on error
	classyfy(*sys.argv[1:])
# classyfy with ClassyFire and save result.
