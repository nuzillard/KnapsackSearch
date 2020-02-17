"""
compounds.py:
	reads: familyname_cids.txt
	writes: familyname_pickled

familyname_cids.txt: written by family_cids.py
familyname_pickled: read by make_2D_sdf.py and basic_tags.py

called by: python compounds.py familyname

File familyname_cids.txt contains a list of KNApSAcK compound_ids, one per line.
compounds.py collects data about each compound from the KNApSAcK web site and
stores them pickled in file familyname_pickled.

Modified on February 17, 2020, as a consequence of a change in the
URL for KNApSAK request submission. See also family_from_web.py.
"""

import requests
import pickle
import sys

def getvalue(s):
	"""
	getvalue() takes a string like <aaa>bbbbb<cc> and returns bbbbb
	"""
	p = s.find('>')
	s = s[p+1:]
	s = s[::-1]
# string reverse
	p = s.find('<')
	s = s[p+1:]
	s = s[::-1]
# string reverse, again
	return s

def get_data_from_cid(cid, site, kw, count):
	"""
	get_data_from_cid() takes a compound_id, the beginning of the url from which
	data about the compound is obtained, a list of keys for which values are searched
	and the order number of the current search
	"""
	print(str(count) + ': ' + cid)
# print search order number and compound_id
	url = site + cid
# build the url for the current compound
	reply = requests.get(url)
	ans = reply.text
	lines = ans.split('\n')
# get the list of the text lines returned by the KNApSAcK server
	
	nkw = len(kw)
# number of keys
	ikw = 0
# current key index in the list of keys
	key = kw[ikw]
# current key

	ikeylines = []
# the list of line indexes at which values are present is initially empty
	for iline, line in enumerate(lines):
# scan HTML text lines
		if key in line:
			ikeylines.append(iline+1)
# if the current key is present in the current line, the value that corresponds to the key
# is in the line that follows the current line
			ikw += 1
# ready for next key.
# the order of the keys in kw is the order in which the keys appear in the HTML text.
			if ikw == nkw:
				break
# if all key-value pairs have been determined, text scanning is finished
			key = kw[ikw]
# otherwise, prepare next key
	values = []
# the list of values is initially empty.
	for iline in ikeylines:
# scan through text line indexes for lines that contain values
		line = lines[iline].strip()
# remove surrounding space characters from text line, result should be like <aaa>bbbbb<cc>
		line = getvalue(line).strip()
# getvalue(line) is bbbbb but can be made of space character(s)
		if not line:
# value is empty
			line = "N/A"
# N/A stands for an empty value
		values.append(line)
# the current value is appended to the list of values
	names = values[0].split('<br>')
# values[0] can be either a single compound name or a row of names separated by <br>.
# names is the list of compound names, possibly with a single element
	values[0] = '||'.join(names)
# values[0] is reformatted with '||' instead of '<br>' as name separator
	d = dict(zip(kw, values))
# make a dictionary from keys and values
	d['C_ID'] = cid
# include the compound_id in the dictionary
	return d

def get_data_from_cid_file(cidfilename, site, kw):
	"""
	get_data_from_cid_file() collects dictionaries about compounds whose coumpound_ids
	are in the lines of file cidfilename. kw is the list of dictionary keys.
	data has to be found using site as url after being completed with a compound_id
	"""
	with open(cidfilename) as fpin:
		return [get_data_from_cid(line.strip(), site, kw, count) for count, line in enumerate(fpin, start=1)]

if __name__ == '__main__':
	site = "http://www.knapsackfamily.com/knapsack_core/information.php?sname=C_ID&word="
# invariable part of the request to KNApSAcK for data access by compound_id
	kw = ["Name", "Formula", "Mw", "CAS RN", "InChIKey", "InChICode", "SMILES"]
# list of keys to data items about compounds, in the order they appear in the HTML text
# returned by KNApSAcK for a search by compound_id. "Name" must be the first element of the list
	
	family = sys.argv[1]
# get family name from command line
	inputfilename = family + "_cids.txt"
	ouputfilename = family + "_pickled"
# build names from family name for input and output files
	data = get_data_from_cid_file(inputfilename, site, kw)
# data is a list of dictionaries, one per compound
	print("Number of compounds: %d" % (len(data),))
	with open(ouputfilename, 'wb') as fpout:
		pickle.dump(data, fpout, pickle.HIGHEST_PROTOCOL)
# write pickled data to file familyname_pickled

