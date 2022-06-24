acd_lotus_v7.sdf
================

The acd_lotusv7 database is produced as an acd_lotus_v7.sdf file

The production process involves the action of python scripts
that include the ones in the current directory:

- csv2smi.py
- splitter.py
- supplement.py

and others located in the parent directory:

- ..\CNMR_Predict\smi2ACD.py
- ..\tautomer.py
- ..\rdcharge.py
- ..\CNMR_Predict\CNMR_predict.py
- ..\sdfrw.py

Recipe
======

The LOTUS database version 7 in file 220525_frozen_metadata.csv.gz, was downloaded from the zenodo data repository,
https://zenodo.org/record/6582124 .

The unpacking of 220525_frozen_metadata.csv.gz provided the text file 220525_frozen_metadata.csv from which each line was scanned
to constitute triplets made of the line number, of a wikidata identifier (QID, see wikidata.org) and of a SMILES chain.
There may be more than one line related to the same QID and SMILES; in this case only the first line was considered for future use.
The text produced by the python script csv2smi.py from the triplets was redirected to file lotusv7.smi;
it contained two columns, one for the SMILES chain and the other one for a compound character string formed
by the QID and the line number joined by an underscore character, such as in Q43656_2, and used as compound name.
Wikidata compound Q43656 is cholesterol and it appears in file 220525_frozen_metadata.csv at line 2.
The file lotusv7.smi contains a minimal description for 201022 compounds.

The script smi2ACD.py applied to lotusv7.smi resulted in file fake_acd_lotusv7.sdf in which a fake chemical shift value
was assigned to each carbon atom (99.99).
The resulting chemical shift lists, one per compound, were formatted to be understood by the ACD/Labs software
as if they were experimental chemical shifts. File fake_acd_lotusv7.sdf contains 201022 compounds.
Action of scripts tautomer.py and rdcharge.py achieved in-place tautomer correction and charge/valence correction.
Operations on structure files make use of the RDKit cheminformatics function library for SMILES chain interpretation,
2D structure diagram generation, tautomer correction and SDF structure file handling.

The splitter.py script applied to file fake_acd_lotusv7.sdf created 21 .sdf files stored in a dedicated sub-directory.
These pieces were named fake_acd_lotusv7_xx.sdf with xx ranging from 00 to 20, the 20 first ones contained 10,000 compounds each
and the last one contained the remaining structures.
Chemical shift prediction by ACD/Labs software using the validation algorithm was carried out on these sub-files.
Each sub-file was imported in an ACD/Labs database, validated for chemical shifts and exported as fake_acd_lotusv7_xx_exported.sdf.
Action of the script CNMR_predict on the latter produced the files true_acd_lotusv7_xx.sdf in which
initial fake chemical shift values were replaced by predicted ones.
The content of these files was assembled into a single file, acd_lotusv7_tmp.sdf, containing 201,010 compound descriptions.
The supplement.py script collected chemical taxonomy data present in file 220525_frozen_metadata.csv
and appended them to compound metadata in file acd_lotusv7_tmp.sdf in order to produce file acd_lotusv7.sdf.
Chemical taxonomy data were obtained by the authors of LOTUS from NPclassifier and Classyfire.
The chemical identifiers SMILES, InChI and InChIKey were recalculated from the structures in acd_lotusv7.sdf
and the resulting InChIKeys compared to those provided by 220525_frozen_metadata.csv.
File acd_lotusv7.sdf was finally imported in the ACD database file acd_lotusv7.NMRUDB for future structural dereplication works.

The acd_lotusv7.sdf file can be read by any cheminformatics tool and its compressed version, of size 219 MB,
was stored in a zenodo repository at https://zenodo.org/record/6621129 .
