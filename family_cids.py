"""
family_cids.py:
	reads: familyname_cid_species.txt
	writes: familyname_cids.txt
	writes: familyname_cid_all_species.txt

familyname_cid_species.txt: written by family_from_web.py
familyname_cids.txt: read by compounds.py
familyname_cid_all_species.txt: read by basic_tags.sdf

called by: python family_cids.py familyname

The file familyname_cid_species.txt contains lines with the compound_id | binomial_name pattern.
The file familyname_cids.txt contains only the compound_ids from file familyname_cid_species.txt
without repetition and listed in increasing alphabetical order
The file familyname_cid_all_species.txt contains lines with the compound_id : binomial_name_list pattern
that associates each compound with a '||'-separated list of binomial_names. It contains one line
per compound_id and the compound_ids are sorted in alphabetical order. The binomial_names
are sorted in alphabetical order in each line.
"""

import sys

family = sys.argv[1]
# get family name from command line
inputfilename = family + "_cid_species.txt"
outputfilename1 = family + "_cids.txt"
outputfilename2 = family + "_cid_all_species.txt"
# build names from family name for input and output files.

cid_dict = {}
# cid_dict associate each compound_id with the list of binomial_names of organisms that contain the compound
with open(inputfilename) as infp:
	for line in infp:
# get line from familyname_cid_species.txt
		cid, species = [x.strip() for x in line.strip().split('|')]
# separate compound_id from binomial_name
		species = ' '.join(species.split())
# normalize binomial_name by removal of internal spaces
		if cid in cid_dict:
			cid_dict[cid].append(species)
# append binomial_name to list if compounds_id has already been encountered
		else:
			cid_dict[cid] = [species]
# new compound_id associated to a single binomial_name
cids = sorted(cid_dict.keys())
# get sorted compound_ids
cids_formatted = ''.join([cid + '\n' for cid in cids])
# build formatted list of compound_ids, one per line.

cids_all_species_dict = {cid: '||'.join(sorted(all_species)) for cid, all_species in cid_dict.items()}
# build the dictionnary that associate each compound_id with the sorted list of binomial_names, with values separated by '||'
cids_all_species_dict_formatted = ''.join(["%s : %s\n" % (cid, cids_all_species_dict[cid]) for cid in cids])
# formatting of familyname_cid_all_species.txt content with compound_id : binomial_name_list pattern

with open(outputfilename1, 'w') as outfp1:
	outfp1.write(cids_formatted)
# write familyname_cids.txt

with open(outputfilename2, 'w') as outfp2:
	outfp2.write(cids_all_species_dict_formatted)
# write familyname_cid_all_species.txt
