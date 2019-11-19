"""
family_from_web.py:
	reads: familyname_genera.txt
	writes: familyname_cid_species.txt

familyname_genera.txt: written by user
familyname_cid_species.txt: read by family_cids.py

called by: python family_from_web.py familyname

The file familyname_genera.txt contains
blank lines, comments with the # sign in column 1, and genus names.
For each genus name, a call to the KNApSAcK search engine returns a web page.
This web page contains a table with one compound_id-binomial_name association per line.
The HTML code of the page is captured and the text is analyzed, retaining only the lines
of interest, in which the string 'class="d1"' specifically appear.
The table has 6 columns, column 1 (index 0) contains the compound_id and 
column 6 (index 5) contains the genus+species string of the binomial_name.
The number of table lines is therefore the number of relevant lines inthe HTML text, divided by 6.

The file familyname_cid_species.txt contains one line per compound_id-binomial_name association
for all the genus names in file familyname_genera.txt.
The compound_id and the binomial_name are separated by a ' | ' string.
The number of associations is printed for each genus.
The first word of the binomial_name is alway identical to the requested genus.
"""

import re
import requests
import sys

def get_cid_spec_pairs(strings, genus):
	"""
	get_cid_spec_pairs() returns a list of (compound_id, binomial_name) pairs for a given genus
	according to strings, the list of the non-HTML parts (the tree leaves) of relevant HTML lines related to the requested genus.
	"""
	pairlist = []
	if len(strings) == 0:
# if no association found for the genus
		return pairlist
# returns an empty list
	n = len(strings) // 6
# the number of associations is the number of table cells divided by 6 (6 coluns in the table)
	pairlist = [(strings[i*6], strings[5+i*6]) for i in range(n)]
# compound_id is in column 1 (index 0) and binomial_name in column 6 (index 5)
	pairlist = [pair for pair in pairlist if pair[1].split()[0] == genus]
# keep only the pairs for which the first word of the binomial_name is identical the species
	return pairlist

def proc_line(line):
	"""
	proc_line() removes the HTML tags from a relevant line from the HTML text bound to a genus.
	Such a line looks like <aaa>bbbbb<cc>. Only bbbbb is returned.
	"""
	return re.sub('<.*?>', '', line)
# use non-greedy (? mark) regular expression to isolate HTML tag content

def proc_text(lines):
	"""
	proc_text() removes the HMTL tags in each line of the lines list and returns the list of relevant data strings.
	"""
	return [proc_line(line) for line in lines]

if __name__ == "__main__":
	family = sys.argv[1]
# get family name from command line
	inputfilename = family + "_genera.txt"
	outputfilename = family + "_cid_species.txt"
# build names from family name for input and output files.
	with open(inputfilename) as infp, open(outputfilename, 'w') as outfp:
		for line in infp:
# get line from genus file
			genus = line.strip()
			if (not genus):
# blank line, skip to next one.
				continue
			if (genus[0] == '#'):
# comment line, skip to next one.
				continue
			url = 'http://www.knapsackfamily.com/knapsack_jsp/result.jsp?sname=organism&word=' + genus
			reply = requests.get(url)
			ans = reply.text
# get HTML text for the current genus.
			lines = [l.strip() for l in ans.split('\n') if (l and ("class=\"d1\"" in l))]
# keep relevant lines from HTML lines, those that are formatted as table cells, of class "d1"
			processed_text = proc_text(lines)
# remove HTML tags from relevant lines
			cid_species_pairs = get_cid_spec_pairs(processed_text, genus)
# get (compound_id, binomial_name) pairs for the current genus
			print("%s: %d" % (genus, len(cid_species_pairs)))
# print the number of associations for the current genus
			for pair in cid_species_pairs:
				outfp.write("%s | %s\n" % pair)
# write compound_id | binomial_name to familyname_cid_species.txt
